#imports
import os
import pandas as pd
from datetime import datetime
import tempfile
import numpy as np
from .species_from_mechanism import return_all_species
import warnings

def wipe_file(file_path):
    open(file_path, 'w').close()
    
def dict_to_config_file(in_dict, filepath):
    if in_dict: #check that the dict isn't empty
        lines = []
        for s in in_dict.keys():
            lines.append(f"{s}   {in_dict[s]:.5e}\n")
        with open(filepath,"w") as file:
            file.writelines(lines)
    else:
        wipe_file(filepath)
    
def list_to_config_file(in_list, filepath):
    if in_list: #check that the dict isn't empty
        lines = [f"{x}\n" for x in in_list]
        with open(filepath,"w") as file:
            file.writelines(lines)
    else:
        wipe_file(filepath)
    
def tuple_list_to_config_file(in_tuple_list, filepath):
    if in_tuple_list: #check that the dict isn't empty
        lines = [f"{x} {y}\n" for x,y in in_tuple_list]
        with open(filepath,"w") as file:
            file.writelines(lines)
    else:
        wipe_file(filepath)
        
def write_config(atchem2_path, initial_concs={}, spec_constrain={}, 
                 spec_constant={}, env_constrain={}, photo_constant = {},
                 photo_constrain = {}, 
                 env_vals = {"TEMP":"298", "PRESS":"1013", "RH":"50",
                             "H2O":"CALC", "DEC":"CALC", "BLHEIGHT":"NOTUSED",
                             "DILUTE":"NOTUSED", "JFAC":"CONSTRAINED",
                             "ROOF":"OPEN", "ASA":"NOTUSED"},
                 spec_output=[]):
    """Prepares model files in specified AtChem2 directory for building and running"""
    #initialConcentrations.config
    dict_to_config_file(initial_concs, 
                        f"{atchem2_path}/model/configuration/initialConcentrations.config")
        
    #speciesConstrained.config
    list_to_config_file(spec_constrain.keys(), 
                        f"{atchem2_path}/model/configuration/speciesConstrained.config")
    #species constraints
    for k,v in spec_constrain.items():
        spec_path = f"{atchem2_path}/model/constraints/species/{k}"
        tuple_list_to_config_file(v,spec_path)
        
    #photolysisConstrained.config
    list_to_config_file(photo_constrain.keys(), 
                        f"{atchem2_path}/model/configuration/photolysisConstrained.config")
    #photolysis constraints
    for k,v in photo_constrain.items():
        spec_path = f"{atchem2_path}/model/constraints/photolysis/{k}"
        tuple_list_to_config_file(v,spec_path)

    #speciesConstant.config
    dict_to_config_file(spec_constant, 
                        f"{atchem2_path}/model/configuration/speciesConstant.config")
    
    #photolysisConstant.config 
    #make lines for config file which has to be in the format "1 1E-4 J1" etc.
    pconst_lines = [f"{k.strip('J')} {v} {k}" for k,v in photo_constant.items()]
    
    list_to_config_file(pconst_lines,
                        f"{atchem2_path}/model/configuration/photolysisConstant.config")
    
    
    #environmentVariables.config
    default_env = {"TEMP":"298",
                   "PRESS":"1013",
                   "RH":"50",
                   "H2O":"CALC",
                   "DEC":"CALC",
                   "BLHEIGHT":"NOTUSED",
                   "DILUTE":"NOTUSED",
                   "JFAC":"CONSTRAINED",
                   "ROOF":"OPEN",
                   "ASA":"NOTUSED"}
    #fill in any missing environment variable values with defaults (if they 
    # aren't supposed to be constrained)
    for k,v in default_env.items():
        if k not in env_vals.keys():
            if k in env_constrain.keys():
                env_vals[k] = "CONSTRAINED"
            else:
                env_vals[k] = v
    
    env_var_lines=f"""1 TEMP			{env_vals['TEMP']}
    2 PRESS			{env_vals['PRESS']}
    3 RH			{env_vals['RH']}
    4 H2O			{env_vals['H2O']}
    5 DEC			{env_vals['DEC']}
    6 BLHEIGHT		{env_vals['BLHEIGHT']}
    7 DILUTE		{env_vals['DILUTE']}
    8 JFAC			{env_vals['JFAC']}
    9 ROOF			{env_vals['ROOF']}
    10 ASA          {env_vals['ASA']}"""  
    
    #append any custom environment variables
    for i,k in enumerate([x for x in env_vals.keys() if x not in default_env.keys()]):
        env_var_lines += f"\n{11+i} {k} {env_vals[k]}"
    
    with open(f"{atchem2_path}/model/configuration/environmentVariables.config","w") as file:
        file.write(env_var_lines)

    #environment constraints
    for k,v in env_constrain.items():
        if env_vals[k] == "CONSTRAINED":
            if k != "JFAC":
                env_path = f"{atchem2_path}/model/constraints/environment/{k}"
            else:
                env_path = f"{atchem2_path}/model/constraints/photolysis/{k}"
            tuple_list_to_config_file(v,env_path)
        else:
            raise Exception(f"Constraint provided for {k}, but value is set as {env_vals[k]} not 'CONSTRAINED'")   
        
    #outputSpecies.config and outputRates.config
    list_to_config_file(spec_output,
                        f"{atchem2_path}/model/configuration/outputSpecies.config")
    list_to_config_file(spec_output,
                        f"{atchem2_path}/model/configuration/outputRates.config")    
    
def write_model_params(atchem2_path, nsteps, model_tstep, tstart, day, month,
                       year, lat=51.51, lon=0.31):
    """Write to model parameters file"""
    model_params_lines=f"""{nsteps}			number of steps
    {model_tstep}			step size (seconds)
    2			species interpolation method (pw constant = 1, pw linear = 2)
    2			conditions interpolation method (pw constant = 1, pw linear = 2)
    {model_tstep}			rates output step size (seconds)
    {tstart}			model start time (seconds)
    0			jacobian output step size (seconds)
    {lat}			latitude (degrees)
    {lon}			longitude (degrees)
    {day:02d}			day
    {month:02d}			month
    {year:04d}			year
    {model_tstep}			reaction rates output step size (seconds)"""

    with open(atchem2_path+"/model/configuration/model.parameters","w") as file:
        file.write(model_params_lines)
    

def build_model(atchem2_path, mechanism_path):
    """Builds specified AtChem2 model, ready for running"""
    script_dir = os.getcwd()
    os.chdir(atchem2_path)
    os.system(f"{atchem2_path}/build/build_atchem2.sh {mechanism_path}")
    os.chdir(script_dir)

def run_model(atchem2_path):
    """Runs specified (pre-built) AtChem2 model"""
    script_dir = os.getcwd()
    os.chdir(atchem2_path)
    os.system(f"{atchem2_path}/atchem2")
    os.chdir(script_dir)
    
                        
def _write_build_run_injections(injection_dict, atchem2_path, mech_path, day, 
                                month, year, t_start, t_end, step_size, 
                                initial_concs, spec_constrain, 
                                spec_constant, env_constrain, photo_constant, 
                                photo_constrain, env_vals, spec_output, 
                                lat, lon):
    """Called by the 'write_build_run' function to configures, build and run
    a specified AtChem2 model including instantaneous increases in 
    concentrations of certain species. 
    
    """
    #dataframe to store model outputs
    stitched_output = pd.DataFrame(dtype=float)
    #dataframes to store rate outputs
    stitched_loss_rates = pd.DataFrame(dtype=float)
    stitched_prod_rates = pd.DataFrame(dtype=float)
    #dataframe to store environment variables
    stitched_env = pd.DataFrame(dtype=float)
    
   
    #make a list of ordered injection times to iterate through
    ordered_times = [list(v.keys()) for v in injection_dict.values()]
    ordered_times = [item for sublist in ordered_times for item in sublist]
    ordered_times = list(set(ordered_times)) #remove duplicates
    ordered_times.sort()
    
    #remove any injections outside of the start and end times and add the start
    #time to the injection list
    ordered_times = [x for x in ordered_times if (x > t_start) and (x < t_end)]
    ordered_times.insert(0, t_start)
    
    all_specs = return_all_species(mech_path)
    
    for i,inj_time in enumerate(ordered_times):
        #copy atchem2 directory
        new_atchem_path = f"{tempfile.gettempdir()}/AtChem2_{year:04d}-{month:02d}-{day:02d}_{i}_{datetime.now()}".replace(' ','_')
        os.system(f"cp -r {atchem2_path} {new_atchem_path}")
        #copy the mechanism to the AtChem directory
        new_mech_path = f"{new_atchem_path}/model/{mech_path.split('/')[-1]}"
        os.system(f"cp {mech_path} {new_mech_path}")
        
        #write config files using data passed
        write_config(new_atchem_path, initial_concs=initial_concs, 
                     spec_constrain=spec_constrain, spec_constant=spec_constant,
                     env_constrain=env_constrain, env_vals=env_vals, 
                     photo_constant = photo_constant, photo_constrain = photo_constrain,
                     spec_output=all_specs) #return all species for now (all needed to set the new start concs)

        
        if (i != (len(ordered_times)-1)): #if this isn't the last iteration then 
        #calculate the next injection time, otherwise the next injection time
        #is just the model end time
            next_injtime = ordered_times[i+1]
        else:
            next_injtime = t_end
            
        if i != 0: #if it isn't the first run, then adjust concentrations based on required injections
        #rewrite initial concentrations file to match the model output from
        #the previous model run
            new_start_concs = stitched_output.loc[min(stitched_output.index, key=lambda x:abs(x-inj_time))] #species concentrations at the closest time to the previous injection time
            
            #change the start concentrations for species injected this time
            specs = [k for k,v in injection_dict.items() if inj_time in v.keys()]
            for s in specs:
                if s != "NOx":
                    new_start_concs.loc[s] = injection_dict[s][inj_time]
                else:
                    #for the NOx constraint, calculate the NO/NO2 ratio 
                    #and change NOx such that the ratio is preserved.
                    old_no = new_start_concs.loc["NO"]
                    old_no2 = new_start_concs.loc["NO2"]
                    
                    old_total_nox = old_no + old_no2
                    nox_deficit = injection_dict[s][inj_time] - old_total_nox                             
                    
                    new_no = (old_no + (nox_deficit*(old_no/old_total_nox)))
                    new_no2 = (old_no2 + (nox_deficit*(old_no2/old_total_nox)))
                    
                    new_start_concs.loc["NO"] = new_no
                    new_start_concs.loc["NO2"] = new_no2
            
            init_lines = [f"{k} {v}\n" for k,v in new_start_concs.to_dict().items()]
            
            with open(new_atchem_path+"/model/configuration/initialConcentrations.config",
                      "w") as file:
                file.writelines(init_lines)
        
    
        #rewrite model parameters file to only run for the length of the 
        #injection of interest
        model_length=(next_injtime+step_size) - inj_time
        nsteps=int(model_length/step_size)
               
        write_model_params(new_atchem_path, nsteps, step_size, inj_time, day, 
                           month, year, lat=lat, lon=lon)
        
        #build and run the model
        build_model(new_atchem_path, new_mech_path)
        run_model(new_atchem_path)
        
        #read the model output and append it to the stitched df
        output = pd.read_csv(f"{new_atchem_path}/model/output/speciesConcentrations.output", 
                             index_col=0, delim_whitespace=True)
        loss_output = pd.read_csv(f"{new_atchem_path}/model/output/lossRates.output", 
                                  delim_whitespace=True,
                                  keep_default_na=False)
        prod_output = pd.read_csv(f"{new_atchem_path}/model/output/productionRates.output", 
                                  delim_whitespace=True,
                                  keep_default_na=False)
        env_output = pd.read_csv(f"{new_atchem_path}/model/output/environmentVariables.output", 
                                 index_col=0, delim_whitespace=True)
        
        #trim off the values that are accounted for by subsequent iterations
        if i != 0:
            output = output.iloc[1:-1,:]
            loss_output = loss_output.loc[loss_output["time"]!=(next_injtime+step_size),:]
            prod_output = prod_output.loc[prod_output["time"]!=(next_injtime+step_size),:]
            env_output = env_output.iloc[1:-1,:]
                
        stitched_output = pd.concat([stitched_output, output])
        stitched_loss_rates = pd.concat([stitched_loss_rates, loss_output])
        stitched_prod_rates = pd.concat([stitched_prod_rates, prod_output])
        stitched_env = pd.concat([stitched_env, env_output])
    
        #remove temporary AtChem2 copy
        os.system(f"rm -r {new_atchem_path}")
        
    #select only the output speices
    stitched_output = stitched_output[spec_output]
    stitched_loss_rates = stitched_loss_rates[stitched_loss_rates["speciesName"].isin(spec_output)]
    stitched_prod_rates = stitched_prod_rates[stitched_prod_rates["speciesName"].isin(spec_output)]

    return (stitched_output,stitched_loss_rates,stitched_prod_rates,stitched_env)
    
    
def _write_build_run_nox_constraint(nox_dict, atchem2_path, mech_path, day, 
                                    month, year, t_start, t_end, step_size, 
                                    initial_concs, spec_constrain, 
                                    spec_constant, env_constrain, photo_constant, 
                                    photo_constrain, env_vals, spec_output, 
                                    lat, lon):
    """Called by the 'write_build_run' function to configures, build and run
    a specified AtChem2 model including a constraint on total NOx, while NO 
    and NO2 are allowed to vary freely.
    """
    warnings.warn("""WARNING. THIS NOX CONSTRAINT FEATURE IS EXPERIMENTAL,
CHECK ANY MODEL OUTPUT THOROUGHLY TO ENSURE THE RESULTS ARE AS EXPECTED.
THE NOX CONSTRAINT FEATURE IS ALSO VERY SLOW AS IT REQUIRES THE REPEATED
BUILDING OF MANY INDIVIDUAL MODELS.""")
    
    #dataframe to store model outputs
    stitched_output = pd.DataFrame(dtype=float)
    #dataframes to store rate outputs
    stitched_loss_rates = pd.DataFrame(dtype=float)
    stitched_prod_rates = pd.DataFrame(dtype=float)
    #dataframe to store environment variables
    stitched_env = pd.DataFrame(dtype=float)
    
   
    #calculate the number of timesteps the model must run for
    model_length = t_end - t_start
    nsteps=int(model_length/step_size)
    
    #make a dictionary of the NOx concentrations at each timestep
    nox_series = pd.Series(index = np.arange(t_start, t_end, step_size))
    for k,v in nox_dict.items():
        nox_series[k] = v

    #fill in empty values
    nox_series = nox_series.sort_index()
    nox_series = nox_series.interpolate()
    nox_series = nox_series.bfill()
    
    all_specs = return_all_species(mech_path)
    
    for istep in range(nsteps):
        step_time = t_start + (istep*step_size)
        
        #copy atchem2 directory
        new_atchem_path = f"{tempfile.gettempdir()}/AtChem2_{year:04d}-{month:02d}-{day:02d}_{istep}_{datetime.now()}".replace(' ','_')
        os.system(f"cp -r {atchem2_path} {new_atchem_path}")
        #copy the mechanism to the AtChem directory
        new_mech_path = f"{new_atchem_path}/model/{mech_path.split('/')[-1]}"
        os.system(f"cp {mech_path} {new_mech_path}")
        
        #write config files using data passed
        write_config(new_atchem_path, initial_concs=initial_concs, 
                     spec_constrain=spec_constrain, spec_constant=spec_constant,
                     env_constrain=env_constrain, env_vals=env_vals, 
                     photo_constant = photo_constant, photo_constrain = photo_constrain,
                     spec_output=all_specs) #return all species for now (all needed to set the new start concs)

        
        if istep != 0: #if it isn't the first run, then adjust concentrations based on required injections
        #rewrite initial concentrations file to match the model output from
        #the previous model run
            new_start_concs = stitched_output.iloc[-1] #species concentrations at the last time step
            
            #for the NOx constraint, calculate the NO/NO2 ratio 
            #and change NOx such that the ratio is preserved.
            old_no = new_start_concs.loc["NO"]
            old_no2 = new_start_concs.loc["NO2"]
            
            old_total_nox = old_no + old_no2
            nox_deficit = nox_series[step_time] - old_total_nox                             
            
            new_no = (old_no + (nox_deficit*(old_no/old_total_nox)))
            new_no2 = (old_no2 + (nox_deficit*(old_no2/old_total_nox)))
            
            new_start_concs.loc["NO"] = new_no
            new_start_concs.loc["NO2"] = new_no2
            
            init_lines = [f"{k} {v}\n" for k,v in new_start_concs.to_dict().items()]
            
            with open(new_atchem_path+"/model/configuration/initialConcentrations.config",
                      "w") as file:
                file.writelines(init_lines)
        else: #just check that we have some NOx in the model if this is the first step
            if not any([x in initial_concs.keys() for x in ["NO","NO2"]]):
                #if there is not initial NO or NO2 specified, then split the
                #given NOx value 50:50 between NO and NO2
                half_val = nox_series[step_time]/2
                
                with open(new_atchem_path+"/model/configuration/initialConcentrations.config",
                          "a") as file:
                    file.write(f"NO2 {half_val}\nNO {half_val}")
            
    
        #rewrite model parameters file to only run for the length of the 
        #injection of interest
               
        write_model_params(new_atchem_path, 1, step_size, step_time, day, 
                           month, year, lat=lat, lon=lon)
        
        #build and run the model
        build_model(new_atchem_path, new_mech_path)
        run_model(new_atchem_path)
        
        #read the model output and append it to the stitched df
        output = pd.read_csv(f"{new_atchem_path}/model/output/speciesConcentrations.output", 
                             index_col=0, delim_whitespace=True)
        loss_output = pd.read_csv(f"{new_atchem_path}/model/output/lossRates.output", 
                                  delim_whitespace=True,
                                  keep_default_na=False)
        prod_output = pd.read_csv(f"{new_atchem_path}/model/output/productionRates.output", 
                                  delim_whitespace=True,
                                  keep_default_na=False)
        env_output = pd.read_csv(f"{new_atchem_path}/model/output/environmentVariables.output", 
                                 index_col=0, delim_whitespace=True)
        
        #trim off the first value (if this isn't the first step)
        if istep != 0:
            output = output.iloc[1:,:]
            env_output = env_output.iloc[1:,:]
                
        stitched_output = pd.concat([stitched_output, output])
        stitched_loss_rates = pd.concat([stitched_loss_rates, loss_output])
        stitched_prod_rates = pd.concat([stitched_prod_rates, prod_output])
        stitched_env = pd.concat([stitched_env, env_output])
    
        #remove temporary AtChem2 copy
        os.system(f"rm -r {new_atchem_path}")
        
    #select only the output speices
    stitched_output = stitched_output[spec_output]
    stitched_loss_rates = stitched_loss_rates[stitched_loss_rates["speciesName"].isin(spec_output)]
    stitched_prod_rates = stitched_prod_rates[stitched_prod_rates["speciesName"].isin(spec_output)]
    
    return (stitched_output,stitched_loss_rates,stitched_prod_rates,stitched_env)
    
def write_build_run(atchem2_path, mech_path, day, month, year, t_start, t_end, 
                    step_size,initial_concs={}, spec_constrain={}, 
                    spec_constant={}, env_constrain={}, photo_constant={}, 
                    photo_constrain={}, 
                    env_vals = {"TEMP":"298", "PRESS":"1013", "RH":"50",
                                "H2O":"CALC", "DEC":"CALC", "BLHEIGHT":"NOTUSED",
                                "DILUTE":"NOTUSED", "JFAC":"CONSTRAINED",
                                "ROOF":"OPEN", "ASA":"NOTUSED"},
                    spec_output=[], lat=51.51, lon=0.31, injection_dict = {},
                    nox_dict = {}):
    """Configures, builds and runs a specified AtChem2 model. 
    
    If spec_inject is specified, then a series of models 
    will be run to simulate a chamber experiments with the introduction of 
    species into the chamber mid-experiment. 
    spec_inject should be a dictionary where keys = species to be injected, 
    and values = second dictionary of times at which injections occur and the
    concentration at the corresponding time e.g. 
    {"C5H8":{36000 : 3E10, 43200 : 2E10}, "NOx":{32400 : 1E11, 43200 : 1E11}} 


    If nox_dict is specified, then a series of one step models 
    will be run with NO and NO2 concentrations adjusted after each, to produce 
    a continuous model output with concentrations of NOx constrained (while NO
    and NO2 can vary freely). The nox_dict variable should be a dictionary with
    keys of model time and values of the desired NOx concentration at each time.
    e.g. {36000 : 1E10, 40000 : 2E10, 45000 : 3E10}
    The NOx concentrations will be linearly interpolated along all of the model
    timesteps.
    WARNING. THIS NOX CONSTRAINT FEATURE IS EXPERIMENTAL AND ALSO VERY SLOW. 
    CHECK ANY MODEL OUTPUT THOROUGHLY TO ENSURE THE RESULTS ARE AS EXPECTED.
    """
    if injection_dict and nox_dict:
        raise Exception("""Cannot run models using both species injections and 
                        NOx constraints. Select either injection_dict or 
                        nox_dict arguments, not both.""")
    elif injection_dict:
        return _write_build_run_injections(injection_dict = injection_dict, 
                                        atchem2_path = atchem2_path, 
                                        mech_path = mech_path, 
                                        day = day, 
                                        month = month, 
                                        year = year, 
                                        t_start = t_start, 
                                        t_end = t_end, 
                                        step_size = step_size,
                                        initial_concs = initial_concs,
                                        spec_constrain = spec_constrain,
                                        spec_constant = spec_constant,
                                        env_constrain = env_constrain,
                                        photo_constant = photo_constant, 
                                        photo_constrain = photo_constrain,
                                        env_vals = env_vals,
                                        spec_output = spec_output,
                                        lat = lat,
                                        lon = lon)
    elif nox_dict:
        return _write_build_run_nox_constraint(nox_dict = nox_dict, 
                                            atchem2_path = atchem2_path, 
                                            mech_path = mech_path, 
                                            day = day, 
                                            month = month, 
                                            year = year, 
                                            t_start = t_start, 
                                            t_end = t_end, 
                                            step_size = step_size,
                                            initial_concs = initial_concs,
                                            spec_constrain = spec_constrain,
                                            spec_constant = spec_constant,
                                            env_constrain = env_constrain,
                                            photo_constant = photo_constant, 
                                            photo_constrain = photo_constrain,
                                            env_vals = env_vals,
                                            spec_output = spec_output,
                                            lat = lat,
                                            lon = lon)
    else:
    
        #copy atchem2 directory
        new_atchem_path = f"{tempfile.gettempdir()}/AtChem2_{year:04d}-{month:02d}-{day:02d}_{datetime.now()}".replace(' ','_')
        os.system(f"cp -r {atchem2_path} {new_atchem_path}")
        #copy the mechanism to the AtChem directory
        new_mech_path = f"{new_atchem_path}/model/{mech_path.split('/')[-1]}"
        os.system(f"cp {mech_path} {new_mech_path}")
        
        #write config files using data passed
        write_config(new_atchem_path, initial_concs=initial_concs, 
                     spec_constrain=spec_constrain, spec_constant=spec_constant,
                     env_constrain=env_constrain, env_vals=env_vals, 
                     photo_constant = photo_constant, photo_constrain = photo_constrain,
                     spec_output=spec_output)
        
        #change model parameters
        model_length=t_end-t_start
        nsteps=int(model_length/step_size)
    
        write_model_params(new_atchem_path, nsteps, step_size, t_start, day, 
                           month, year, lat=lat, lon=lon)
        
        #build and run the model
        build_model(new_atchem_path, new_mech_path)
        run_model(new_atchem_path)
        
        #read the model output and append it to the stitched df
        output = pd.read_csv(f"{new_atchem_path}/model/output/speciesConcentrations.output", 
                             index_col=0, delim_whitespace=True)
        
        
        loss_output = pd.read_csv(f"{new_atchem_path}/model/output/lossRates.output", 
                                  delim_whitespace=True,
                                  keep_default_na=False)
        prod_output = pd.read_csv(f"{new_atchem_path}/model/output/productionRates.output", 
                                  delim_whitespace=True,
                                  keep_default_na=False)
        
    
        env_output = pd.read_csv(f"{new_atchem_path}/model/output/environmentVariables.output", 
                                 index_col=0, delim_whitespace=True)
        
        #remove temporary AtChem2 copy
        os.system(f"rm -r {new_atchem_path}")
        return (output, loss_output, prod_output, env_output)

