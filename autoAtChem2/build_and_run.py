#imports
import os
import pandas as pd
from datetime import datetime

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
                 spec_constant={}, env_constrain={}, jfac_constrain = [(0,1)], 
                 env_vals = {"TEMP":"298","PRESS":"1013","RH":"50","DILUTE":"NOTUSED"},
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
    
    #speciesConstant.config
    dict_to_config_file(spec_constant, 
                        f"{atchem2_path}/model/configuration/speciesConstant.config")
    
    #environmentVariables.config
    env_var_lines=f"""1 TEMP			{env_vals['TEMP']}
    2 PRESS			{env_vals['PRESS']}
    3 RH			{env_vals['RH']}
    4 H2O			CALC
    5 DEC			CALC
    6 BLHEIGHT		NOTUSED
    7 DILUTE		{env_vals['DILUTE']}
    8 JFAC			CONSTRAINED
    9 ROOF			OPEN
    10 ASA          NOTUSED"""  
    with open(f"{atchem2_path}/model/configuration/environmentVariables.config","w") as file:
        file.write(env_var_lines)

    #environament constraints
    for k,v in env_constrain.items():
        if env_vals[k] == "CONSTRAINED":
            env_path = f"{atchem2_path}/model/constraints/environment/{k}"
            tuple_list_to_config_file(v,env_path)
        else:
            raise Exception(f"Constraint provided for {k}, but value is set as {env_vals[k]} not 'CONSTRAINED'")   
    
    if jfac_constrain:
        tuple_list_to_config_file(jfac_constrain,
                                  f"{atchem2_path}/model/constraints/photolysis/JFAC")
    else:
        raise Exception("No JFAC constraint provided")
    
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
    
def write_build_run(atchem2_path, mech_path, day, month, year, t_start, t_end, 
                    step_size,initial_concs={}, spec_constrain={}, 
                    spec_constant={}, env_constrain={}, jfac_constrain = [(0,1)], 
                    env_vals = {"TEMP":"298","PRESS":"1013",
                                "RH":"50","DILUTE":"NOTUSED"},
                    spec_output=[], spec_inject = [], inject_targets = {}, 
                    lat=51.51, lon=0.31):
    """Configures, builds and runs a specified AtChem2 model. 
    
    If spec_inject and inject_targets are specified, then a series of models 
    will be run to simulate a chamber experiments with the introduction of 
    species into the chamber mid-experiment. 
        spec_inject should be a dictionary where keys = species to be injected, 
        and values = a list of times at which injectionsshould occur. e.g. 
        {"C5H8":[36000, 43200], "NOx":[32400, 43200]} 
        
        inject_targets should be a dictionary where keys = species to be injected
        and values = a list of target concentrations at each corresponding 
        time point in spec_inject. e.g. {"C5H8":[3E10, 2E10], 
                                         "NOx":[1E11, 1E11]} """
    #dataframe to store model outputs
    stitched_output = pd.DataFrame(dtype=float)
    #dataframes to store rate outputs
    stitched_loss_rates = pd.DataFrame(dtype=float)
    stitched_prod_rates = pd.DataFrame(dtype=float)
    #dataframe to store environment variables
    stitched_env = pd.DataFrame(dtype=float)
    
    if (not inject_targets) and (not spec_inject): #if species injections ARE NOT needed
        #copy atchem2 directory
        new_atchem_path = f"{atchem2_path}_{year:04d}-{month:02d}-{day:02d}_{datetime.now()}".replace(' ','_')
        os.system(f"cp -r {atchem2_path} {new_atchem_path}")
        #copy the mechanism to the AtChem directory
        new_mech_path = f"{new_atchem_path}/model/{mech_path.split('/')[-1]}"
        os.system(f"cp {mech_path} {new_mech_path}")
        
        #write config files using data passed
        write_config(new_atchem_path, initial_concs=initial_concs, 
                     spec_constrain=spec_constrain, spec_constant=spec_constant,
                     env_constrain=env_constrain, env_vals=env_vals, jfac_constrain = jfac_constrain,
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
        
        stitched_output = pd.concat([stitched_output, output])
        
        loss_output = pd.read_csv(f"{new_atchem_path}/model/output/lossRates.output", 
                                  delim_whitespace=True)
        prod_output = pd.read_csv(f"{new_atchem_path}/model/output/productionRates.output", 
                                  delim_whitespace=True)
        
        stitched_loss_rates = pd.concat([stitched_loss_rates, loss_output])
        stitched_prod_rates = pd.concat([stitched_prod_rates, prod_output])

        env_output = pd.read_csv(f"{new_atchem_path}/model/output/environmentVariables.output", 
                                 index_col=0, delim_whitespace=True)
        stitched_env = pd.concat([stitched_env, env_output])
        
        #remove temporary AtChem2 copy
        os.system(f"rm -r {new_atchem_path}")
        return (stitched_output,stitched_loss_rates,stitched_prod_rates,stitched_env)

    elif (inject_targets) and (spec_inject): #if species injections ARE needed
        #go through each injection in time order and run that portion of the model
       
        for i,(inj_time,specs) in enumerate(spec_inject):
            #copy atchem2 directory
            new_atchem_path = f"{atchem2_path}_{year:04d}-{month:02d}-{day:02d}_{i}_{datetime.now()}".replace(' ','_')
            os.system(f"cp -r {atchem2_path} {new_atchem_path}")
            #copy the mechanism to the AtChem directory
            new_mech_path = f"{new_atchem_path}/model/{mech_path.split('/')[-1]}"
            os.system(f"cp {mech_path} {new_mech_path}")
            
            #write config files using data passed
            write_config(new_atchem_path, initial_concs=initial_concs, 
                         spec_constrain=spec_constrain, spec_constant=spec_constant,
                         env_constrain=env_constrain, env_vals=env_vals, jfac_constrain = jfac_constrain,
                         spec_output=spec_output)

            
            if (i != (len(spec_inject)-1)): #if this isn't the first or last iteration then adjust
                         #the starting concentrations to match the end of the last model run
                next_injtime = spec_inject[i+1][0]
                
                if i != 0: #if it isn't the first run, then adjust concentrations based on required injections
                #rewrite initial concentrations file to match the model output from
                #the previous model run
                    new_start_concs = stitched_output.loc[min(stitched_output.index, key=lambda x:abs(x-inj_time))] #species concentrations at the closest time to the previous injection time
                    
                    #change the start concentrations for species injected this time
                    for s in specs:
                        if s != "NOx":
                            new_start_concs.loc[s] = inject_targets[s][inj_time]
                        else:
                            #for the NOx constraint, calculate the NO/NO2 ratio 
                            #and change NOx such that the ratio is preserved.
                            old_no = new_start_concs.loc["NO"]
                            old_no2 = new_start_concs.loc["NO2"]
                            
                            old_total_nox = old_no + old_no2
                            nox_deficit = inject_targets[s][inj_time] - old_total_nox                             
                            
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
                                          delim_whitespace=True)
                prod_output = pd.read_csv(f"{new_atchem_path}/model/output/productionRates.output", 
                                          delim_whitespace=True)  
                env_output = pd.read_csv(f"{new_atchem_path}/model/output/environmentVariables.output", 
                                         index_col=0, delim_whitespace=True)
                
                #trim off the last value (this will be replaced by the next model run)
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
        return (stitched_output,stitched_loss_rates,stitched_prod_rates,stitched_env)
    
    else:
        raise Exception("""Both 'inject_targets' and 'spec_inject' arguments 
                        must be provided or empty, not just one""")
