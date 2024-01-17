"""Functions to read the output files from AtChem2"""
#imports
import pandas as pd

def species_concentrations_df(file_path, species="ALL", 
                                  error_for_non_species=False):
    """Reads speciesConcentrations.output files into a pandas dataframe"""
    #Create dataframe from file at the given path
    data = pd.read_csv(file_path, index_col=0, delim_whitespace=True)
    
    #select only the species specified in the "species" variable
    if type(species) == list:
        if (error_for_non_species and 
            (not all(e in data.columns for e in species))):
            raise Exception(f"""Provided species not present in model output: {[", ".join([x for x in species if x not in data.columns])]}""")
        data = data[species]
    elif type(species) == str:
        if species.casefold() != "ALL".casefold():
            if (error_for_non_species and (species not in data.columns)):
                raise Exception(f"""Provided species not present in model output: {species}""")
            data= data[[species]]
    else:
        raise TypeError(f"""Invalid input of species. Species argument must be 
                        a list of species names, a string of the name of a 
                        species, or the string "ALL".
                        Provided input = {species}.""")

    return data

def return_net_0_species(rxn_str):
    """Returns species in a reaction which are present in both the reactants 
    and products"""
    reacts = rxn_str.split("=")[0].split("+")
    prods = rxn_str.split("=")[1].split("+")
    
    intersection = set(reacts) & set(prods)
    
    return [x for x in intersection]

def rate_df(file_path, species="ALL", drop_0=True, drop_net_0=True, 
                 drop_rev=False, error_for_non_species = True):
    """Reads lossRates.output or productionRates.output files into a pandas dataframe"""
    #Create dataframe from file at the given path
    data = pd.read_csv(file_path, index_col=[0,2,3], delim_whitespace=True,
                       keep_default_na=False)
    
    #get df of all reactions for dropping reaction later if needed
    rxns = data.groupby(level=2).first()["reaction"]
    
    #get list of all species present in output to raise error for non-species 
    #if requested
    unique_specs = data.groupby(level=[1]).first().index.unique()   
    
    #select only the species specified in the "species" variable
    if type(species) == list:
        if (error_for_non_species and 
            (not all(e in unique_specs for e in species))):
            raise Exception(f"""Provided species not present in model output: {[", ".join([x for x in species if x not in unique_specs])]}""")
            
        data = data.loc[:,species,:,:]
    elif type(species) == str:
        if species.casefold() != "ALL".casefold():
            if (error_for_non_species and (species not in unique_specs)):
                raise Exception(f"""Provided species not present in model output: {species}""")
            data = data.loc[:,[species],:,:]
    else:
        raise TypeError(f"""Invalid input of species. Species argument must be 
                        a list of species names, a string of the name of a 
                        species, or the string "ALL".
                        Provided input = {species}.""")
                        
    #remove reactions with a rate of 0 throughout the whole model if specified
    if drop_0:
        #get df of the reaction rates summed at every time step
        summed_data = data.groupby(level=[2]).sum()
        #get a list of reaction numbers where the total rate is 0 for the whole model
        zero_rxns = summed_data[summed_data["rate"] == 0].index.tolist()
        
        #remove the 0 reactions from the dataframe
        data = data.drop(zero_rxns, axis=0, level=2)
    
    #remove reactions where a compound is lost and produced simultaniously where specified
    if drop_net_0:
        #get series of net_0 compounds in each reaction
        data["repeated"] = data["reaction"].apply(return_net_0_species)
        #split the list of repeated compounds into a multi-index
        repeated = data.apply(lambda x: x.name[1] in x["repeated"], axis=1)
        remove_idxs = repeated[repeated == True].index
        
        #remove the net 0 reactions from the data
        data = data.drop(remove_idxs)
        #remove the "repeated" column
        data = data.drop("repeated", axis=1)
        
    #remove reactions where the analogous reverse reaction is also present where specified
    if drop_rev:
        #split reactions into columns of reatcants and products
        split_rxns = rxns.str.split('=', n=1, expand=True)
        split_rxns[0] = split_rxns[0].str.split("+")
        split_rxns[1] = split_rxns[1].str.split("+")
        
        #function to see if there are any reverse reactions for a given set of 
        #reactants (r) and products (p) in a dataframe (df)
        has_reversed_rxns = lambda df, r, p : not df[(df[0].apply(lambda x : sorted(x) == sorted(p)) & 
                                                      df[1].apply(lambda x : sorted(x) == sorted(r)))].index.empty
        #apply this function to all rows in the split reactions dataframe
        reversible_reactions = split_rxns.apply(lambda x : has_reversed_rxns(split_rxns, x[0], x[1]), axis=1)
        reversible_reactions = reversible_reactions[reversible_reactions==True].index
        
        #remove the reversible reactions from the data
        data = data.drop(reversible_reactions, level=2)
        
    return data

