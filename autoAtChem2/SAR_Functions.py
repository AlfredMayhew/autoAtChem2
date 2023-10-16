"""Contains functions to calculate rates of reactions based on Structure 
Activity Relationships (SARs) used in the MCM."""
#imports
from rdkit import Chem
import numpy as np
import pandas as pd
from collections import defaultdict

def mol_composition(molecule):
    """Get the composition of an RDKit molecule:
    Atomic counts, including hydrogen atoms, and any charge.
    For example, fluoride ion (chemical formula F-, SMILES string [F-])
    returns {9: 1, 0: -1}. (Adapted from https://github.com/rdkit/rdkit/discussions/5339)

    :param molecule: The molecule to analyze
    :type some_input: An RDKit molecule
    :rtype: A dictionary.
    """
    # Check that there is a valid molecule
    if molecule:

        # Add hydrogen atoms--RDKit excludes them by default
        molecule_with_Hs = Chem.AddHs(molecule)
        comp = defaultdict(lambda: 0)

        # Get atom counts
        for atom in molecule_with_Hs.GetAtoms():
            comp[atom.GetSymbol()] += 1

        # If charged, add charge as "atomic number" 0
        charge = Chem.GetFormalCharge(molecule_with_Hs)
        if charge != 0:
            comp[0] = charge
        return comp
    else:
        raise Exception("Invalid molecule.")

def RO2_Class(ro2_mol):
    """Identifies whether a given RO2 is primary, secondary, or tertiary."""
    prim_ro2_smarts = Chem.MolFromSmarts("[#6]([OX2][OX1v1])[#6,$(O-C)]")
    sec_ro2_smarts = Chem.MolFromSmarts("[#6]([OX2][OX1v1])([#6,$(O-C)])[#6,$(O-C)]")
    tert_ro2_smarts = Chem.MolFromSmarts("[#6]([OX2][OX1v1])([#6,$(O-C)])([#6,$(O-C)])[#6,$(O-C)]")
    
    #go through each pattern from tert->sec->prim as a tert RO2 will match the
    #prim pattern
    #first detect methyl-peroxy which will be a special case
    if Chem.MolToSmiles(ro2_mol) == "CO[O]":
        return 1
    elif ro2_mol.HasSubstructMatch(tert_ro2_smarts):
        return 3
    elif ro2_mol.HasSubstructMatch(sec_ro2_smarts):
        return 2
    elif ro2_mol.HasSubstructMatch(prim_ro2_smarts):
        return 1
    else:
        raise Exception(f"The provided molecule does not match to any of tertiary, secondary, or primary RO2. Provided mol = {Chem.MolToSmiles(ro2_mol)}")

def _Jenkin2019_RO2_ref_self_rate(ncon, ro2_class):
    """Returns the reference self reaction rate of a given RO2 as calculated 
    in Jenkin et al 2019."""
    if ro2_class == 1:
        return 10**(-11.7-(3.2*np.exp(-0.55*(ncon-0.52))))
    elif ro2_class == 2:
        return 10**(-12.9-(3.2*np.exp(-0.64*(ncon-2.3))))
    elif ro2_class == 3:
        return 2.1E-17
    else:
        raise Exception(f"Invalid RO2 class entered. Provided RO2 class does not match 1, 2, or 3 (corresponding to primary, secondary, or tertiary): {ro2_class}")
    
def _Jenkin2019_RO2_self_rate(ref_rate, ro2_mol):
    """Calculates the self-reaction rate for a given RO2 according to 
    Equation 17 in Jenkin et al. 2019"""
    #activating factors alpha and beta for different functionalities
    activ_facts = pd.DataFrame([(1, 0, "[#6][OX2v2][OX1v1]"),
                                (8E-5, 0.4, "[OH1][#6][#6][OX2v2][OX1v1]"),
                                (4E-2, 0.15, "[#6]=[#6][#6][OX2v2][OX1v1]"),
                                (5.8E-2, 0.15, "[c][#6][OX2v2][OX1v1]"),
                                (7E-5, 0.4, "[#6][OX2][#6][OX2v2][OX1v1]"),
                                (1.6E-4, 0.4, "[#6](=O)[#6][OX2v2][OX1v1]"),
                                (5.3E-5, 0.4, "[#6](=O)[#6][#6][OX2v2][OX1v1]")],
                               index = ["alkyl", "b-hydroxy", "allylic", 
                                        "b-aryl", "a-alkoxy", "b-oxo", "g-oxo"],
                               columns = ["a", "b", "SMARTS"])
    activ_facts["molSMARTS"] = activ_facts["SMARTS"].apply(Chem.MolFromSmarts)

    #check for matches to any of the activating factors
    matches = {i : ro2_mol.HasSubstructMatch(m) for i,m in activ_facts["molSMARTS"].items()}
    
    #check to see if the RO2 contains either an allylic or aryl RO2 as well as 
    #any other functional groups. If this is the case, then the alpha and beta
    #values must be combined
    matches_add = []
    if matches["b-aryl"]:
        for k,v in matches.items():
            if (v and (k not in ["b-aryl", "allylic", "alkyl"])): #Go through each match and add a combination of the allyl/aryl to the activation factors
                new_i_name =  f"aryl-{k}"   
                activ_facts.loc[new_i_name,:] = [activ_facts.loc["b-aryl", "a"]*activ_facts.loc[k, "a"],
                                                 activ_facts.loc["b-aryl", "b"]+activ_facts.loc[k, "b"],
                                                 "", ""]
                matches_add.append(new_i_name)
    if matches["allylic"]:
        for k,v in matches.items():
            if (v and (k not in ["b-aryl", "allylic", "alkyl"])): #Go through each match and add a combination of the allyl/aryl to the activation factors
                new_i_name =  f"aryl-{k}"   
                activ_facts.loc[new_i_name,:] = [activ_facts.loc["allylic", "a"]*activ_facts.loc[k, "a"],
                                                 activ_facts.loc["allylic", "b"]+activ_facts.loc[k, "b"],
                                                 "", ""]
                matches_add.append(new_i_name)
    for m in matches_add:
        matches[m] = True
        
    #calculate the rate for each functionality
    calc_rate_activ = lambda a, b : a*(ref_rate**(1-b))
    all_rates = {k : calc_rate_activ(a,b) for k,a,b in zip(activ_facts.index, 
                                                           activ_facts["a"].values,
                                                           activ_facts["b"].values)}
    
    #select the fastest rate of the ones corresponding to the present functionality
    max_rate = max([v for k,v in all_rates.items() if matches[k]])

    return max_rate

def Jenkin2019_RO2_Rate(smiles_a, smiles_b):
    """Returns the rate of reaction for two RO2 (input as SMILES strings). 
    Uses the rules from Jenkin et al. 2019 to calculate rates.
    First a reference self-reaction rate is calculated for each RO2, then each 
    reference rate is modified with substituent activation factors, then if the
    two RO2 are different, the geometric mean of the two self-reaction rates is
    taken as the rate of reaction.
    """
    mol_a = Chem.MolFromSmiles(smiles_a)
    mol_b = Chem.MolFromSmiles(smiles_b)
    
    #sum the C, O, and N for use in calculating the reference rate
    composition_a = mol_composition(mol_a)
    composition_b = mol_composition(mol_b)
    ncon_a = composition_a["C"] + composition_a["O"] + composition_a["N"]
    ncon_b = composition_b["C"] + composition_b["O"] + composition_b["N"]
    
    #identify whether each ro2 is primary, secondary, or tertiary
    class_a = RO2_Class(mol_a)
    class_b = RO2_Class(mol_b)
    
    #calculate the self reaction rates
    ref_self_rate_a = _Jenkin2019_RO2_ref_self_rate(ncon_a, class_a)
    ref_self_rate_b = _Jenkin2019_RO2_ref_self_rate(ncon_b, class_b)
    
    self_rate_a = _Jenkin2019_RO2_self_rate(ref_self_rate_a, mol_a)
    self_rate_b = _Jenkin2019_RO2_self_rate(ref_self_rate_b, mol_b)
    
    #see if ro2a and ro2b are the same. If so then this is a self-reaction and
    #the calculated rate can be used. If not, then take the geometric mean of 
    #the two rates to get the cross-reaction rate
    if Chem.CanonSmiles(smiles_a) == Chem.CanonSmiles(smiles_b):
        return_rate = self_rate_a
    else:
        return_rate = 2*(self_rate_a*self_rate_b)**0.5

    return f"{return_rate : .3E}"


