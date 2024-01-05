#imports
import re

def get_reaction_lines(mechanism_path):
    """Returns a list of reactions given in a FACSIMILE format mechanism"""
    with open(mechanism_path) as file:
        mech_txt = file.read()
    
    rxns = mech_txt.split("Reaction definitions")[1].split("\n")
    
    rxn_pattern = re.compile("%.*:.*=.*;")
    
    return [x for x in rxns if rxn_pattern.match(x)]

def get_species_from_lines(mechanism_path):
    """returns all of the species included in a set of FACSIMILE reactions"""
    lines = get_reaction_lines(mechanism_path)
    
    rxn_pattern = re.compile("%.*:(.*)=(.*);")
    
    #pattern to match a species with or without a stoichiometric coefficient
    spec_pattern = re.compile(r"^ *(\d*\.?\d*) *([a-zA-Z_].*) *$")
    
    comps = set()
    for l in lines:
        match = rxn_pattern.match(l)
        reacts = [spec_pattern.match(x)[2] for x in match[1].split("+") if x.strip()]
        prods = [spec_pattern.match(x)[2] for x in match[2].split("+") if x.strip()]
        
        for r in reacts:
            if r:
                comps.add(r)
        for p in prods:
            if p:
                comps.add(p)
    
    return comps

def return_inorganic_species(mechanism_path, 
                             speclist=["N2O5", "H2O2", "NO", "H2", "NA", "HONO", 
                                       "OH", "SO2", "O", "HNO3", "SO3", "O1D", 
                                       "HO2", "HO2NO2", "CO", "SA", "O3",
                                       "HSO3", "NO2", "NO3"]):
    """Function returns list of all species in the mcm inorganic mechanism"""
    
    mech_species = get_species_from_lines(mechanism_path)

    return [x for x in mech_species if x in speclist]

def return_all_species(mechanism_path):
    """Function returns list of all species in the mcm FACSIMILE mechanism"""
    mech_species = get_species_from_lines(mechanism_path)

    return list(mech_species)

