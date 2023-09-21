"""Script to plot stacked plot of production and loss rates of species from
AtChem2 output files"""
#imports
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colormaps
from autoAtChem2.read_output import rate_df
from autoAtChem2.utilities import string_to_bool
import sys
import os
import numpy as np
from datetime import date
import pandas as pd

#read in arguments from command line
args=sys.argv
kwarg_dict={"drop_rev" : True, "drop_0" : True, "drop_net_0" : True, "title_page_text":""}

if len(args) < 6:
    raise Exception("""Must provide at least 5 arguments (in this order):
                     - Model Output Path (string, path to a model output directory containing a lossRates.output file and productionRates.output file)
                     - Species of Interest (Comma Separated List e.g. NO2,O3,NO3)
                     - Number of Reactions to List (int, e.g. 10 for Top 10 Reactions, with the rest being lumped into "Other")
                     - Start Time (In Model Time or "START")
                     - End Time (In Model Time or "END")
                     
                     Additional key word arguments (e.g. drop_rev=True) are:
                     - title_page_text (string, Text to add to the title page of the pdf output.)
                     - drop_rev (bool, Ignore reversible reactions? Default: True)
                     - drop_0 (bool, Ignore reactions where the rate is 0 throughout the model. Default: True)
                     - drop_net_0 (bool, Ignore reactions where the species of interest is both a reactant and product. Default: True)
                     """)
else: #if 5 or more arguments passed then get the initial 5
    print(args)    

    out_path=args[1]
    species=args[2].strip("[]").split(",") #strip [] in case it was entered as python syntax list
    top_n=int(args[3])
    
    if args[4].casefold() == "start".casefold():
        start = "START"
    else:
        try:
            start = float(args[4])
        except ValueError:
            raise Exception("Start time must be numeric or 'START'")
    
    if args[5].casefold() == "end".casefold():
        end = "END"
    else:
        try:
            end = float(args[5])
        except ValueError:
            raise Exception("End time must be numeric or 'END'")    

if len(args) > 6: #if there are more than 5 arguments then process the kwargs
    for kwarg in args[6:]:
        kw=kwarg.split("=")[0]
        arg=kwarg.split("=",1)[1] #only split at the first =, may be subsequent = from reaction definitions
        
        if (kw in ["drop_rev", "drop_0", "drop_net_0"]): #bools
            kwarg_dict[kw] = string_to_bool(arg)
        elif kw == "title_page_text": #strings
            kwarg_dict[kw] = arg
            
###############################################################################
#read in loss and production rate files
l_rates = rate_df(f"{out_path}{os.sep}lossRates.output", species = species,
                  drop_0 = kwarg_dict["drop_0"], drop_net_0 = kwarg_dict["drop_net_0"], 
                  drop_rev = kwarg_dict["drop_rev"], 
                  error_for_non_species = False)
p_rates = rate_df(f"{out_path}{os.sep}productionRates.output", species = species,
                  drop_0 = kwarg_dict["drop_0"], drop_net_0 = kwarg_dict["drop_net_0"], 
                  drop_rev = kwarg_dict["drop_rev"], 
                  error_for_non_species = False)


#trim data down based on start and end times provided
#get the index corresponding to the start and end times provided
time_list = p_rates.groupby(level=0).first().index
if start == "START":
    start_num = time_list[0]
else:
    start_num = start
if end == "END":
    end_num = time_list[-1]
else:
    end_num = end
l_rates = l_rates.sort_index(level=0).loc[start_num:end_num,:,:,:]
p_rates = p_rates.sort_index(level=0).loc[start_num:end_num,:,:,:]

#split the data into separate dataframes for each species for further processing
spec_l_dfs = {}
for s in l_rates.index.get_level_values(1).unique():
    spec_l_dfs[s] = l_rates.loc[:,s,:,:]
spec_p_dfs = {}
for s in p_rates.index.get_level_values(1).unique():
    spec_p_dfs[s] = p_rates.loc[:,s,:,:]

#sort dataframe to get the top n reactions for each species, lumping the rest 
#together into "other"
def top_specs(df_dict):
    for s in df_dict.keys():
        #get df of the average rate of each reaction over the whole model
        avg_series = df_dict[s].groupby(level=1)["rate"].median().sort_values(ascending=False)
        #extract top n reactions for this species and then all of the other reactions
        topn_nos = avg_series.iloc[:top_n].index
        other_nos = avg_series.iloc[top_n:].index
        
        #get the dataframe subsets for the top n compounds and a separate subset of
        #the other compounds
        topn_df = df_dict[s].loc[:,topn_nos,:].sort_index(level=0)
        other_df = df_dict[s].loc[:,other_nos,:].sort_index(level=0)
        
        #add the rates of all of the other subsets at each timestep
        other_df = other_df.groupby(level=0).sum()
        #give the other rates df a dummy reaction number and reaction name to 
        #allow for appending (np.inf used so that "other" will be set as the final reaction)
        other_df = other_df.assign(reactionNumber = np.inf, 
                                   reaction = "Other").set_index('reactionNumber', 
                                                                 append=True)
        
        df_dict[s] = pd.concat([topn_df,other_df]).sort_index(level=0)

top_specs(spec_l_dfs)
top_specs(spec_p_dfs)

#calculate the production and loss rates as a percentage of the total loss
#for the percentage plots and add the percentage data to the relavent dataframe
for s, df in spec_l_dfs.items():
    df["%rate"] = df["rate"].divide(df.groupby(level=0).sum()["rate"])*100
for s, df in spec_p_dfs.items():
    df["%rate"] = df["rate"].divide(df.groupby(level=0).sum()["rate"])*100

###############################################################################       
#Plotting
pp=PdfPages("temp_rates_plot.pdf")
cols=colormaps.get_cmap("tab20")


#create title page for pdf output
confirstPage = plt.Figure()
plottxt=f"""Production rates of : {", ".join(spec_p_dfs.keys())}
Loss rates of : {", ".join(spec_l_dfs.keys())}

Reversible Reactions Removed : {kwarg_dict['drop_rev']}
Zero-rate Reactions Removed: {kwarg_dict['drop_0']}
Net-zero Reactions Removed: {kwarg_dict['drop_net_0']}

{kwarg_dict['title_page_text']}

Plotted on {date.today()}"""
confirstPage.text(0.5,0.5,plottxt, size=10, ha="center", wrap=True)
pp.savefig(confirstPage)

#function for plotting rates
def rate_plot(df_dict, fig, col_name, title_end, rate_units):
    n_species = len(df_dict.keys())
    for i,(s,df) in enumerate(df_dict.items()):
        times = df.index.get_level_values(0).unique()
        labels = df["reaction"].unique()
     
        #plot actual plot
        ax=fig.add_subplot(n_species,1,i+1)
        ax.stackplot(times, df[col_name].unstack().values.transpose(),
                     labels = labels,
                     colors = cols([i for i in range(len(labels))]))

        
        ax.set_ylabel(f"Rate ({rate_units})")
        ax.set_xlabel("Time (s)")
        ax.set_title(f"{s} {title_end}")
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    fig.suptitle(f"{title_end} Reactions")
    pp.savefig(fig,bbox_inches="tight")


lfig=plt.Figure(figsize=(10,5*len(spec_l_dfs.keys())))
pfig=plt.Figure(figsize=(10,5*len(spec_p_dfs.keys())))
perc_lfig=plt.Figure(figsize=(10,5*len(spec_l_dfs.keys())))
perc_pfig=plt.Figure(figsize=(10,5*len(spec_p_dfs.keys())))

rate_plot(spec_l_dfs, lfig, "rate", "Loss", "molecule $cm^{-3}\ s^{-1}$")
rate_plot(spec_p_dfs, pfig, "rate", "Production", "molecule $cm^{-3}\ s^{-1}$")
rate_plot(spec_l_dfs, perc_lfig, "%rate", "% Loss", "%")
rate_plot(spec_p_dfs, perc_pfig, "%rate", "% Production", "%")

pp.close()
