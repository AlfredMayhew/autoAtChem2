#imports
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import pandas as pd
from datetime import timedelta
from .reading_concentrations import conc_to_units
import math

def make_axis(model_data, measured_data, ax, name, units, env_df, cconv, 
              injections = None):
    """Plots measured and modelled data from two Series objects."""
    #convert units of modelled data
    model_data = conc_to_units(model_data, units, cconv)
    
    #get a datetime index for better plotting
    model_data_df_index = [pd.to_datetime(str(timedelta(seconds=x))) for x in model_data.index]
    if not measured_data.empty:
        measured_data_df_index = [pd.to_datetime(str(timedelta(seconds=x))) for x in measured_data.index]
        ax.plot(measured_data_df_index, measured_data, ".", c="rebeccapurple", 
                label=f"Measured {name}")
        
    ax.plot(model_data_df_index, model_data, c="firebrick", 
            label=f"Modelled {name}")
    
    #plot injections if specified
    if injections:      
        for t in injections:
            ax.axvline(pd.to_datetime(str(timedelta(seconds=t))),ls="dotted",c="k")
        
    
    time_form = DateFormatter("%H:%M")
    ax.xaxis.set_major_formatter(time_form)
    
    ax.set_xlabel("Time (UTC)")
    ax.set_ylabel(f"{name} ({units})")
    
    ax.legend()
    
def plot_species(model_df, df_dict, trans_dict, specs_to_plot, cconv, title=None, 
                 injections = None):
    """Makes a plt figure with the requested species plotted"""
    nspecs = len(specs_to_plot)    
    fig = plt.Figure(figsize = (10, math.ceil(nspecs/2)*5))
        
    #get start and end of model to trim the measured data
    tstart = model_df.index[0]
    tend = model_df.index[-1]
    
    for i,spec in enumerate(specs_to_plot):
        try:
            #read measured data
            name = trans_dict[spec][0]
            meas_series = df_dict[trans_dict[spec][1]][name]
            
            closest_tstart = min(meas_series.index, key=lambda x:abs(x-tstart))
            closest_tend = min(meas_series.index, key=lambda x:abs(x-tend))
            
            meas_series = df_dict[trans_dict[spec][1]][name].loc[closest_tstart:closest_tend].dropna()
            
            #get units
            units = trans_dict[spec][2]
    
        except KeyError:
            meas_series = pd.Series()
            units = "ppt"
            
        #get model data
        if spec != "NOx":
            model_series = model_df[spec]
        else:
            model_series = model_df[["NO","NO2"]].sum(axis=1)

        #plot        
        ax = fig.add_subplot(math.ceil(nspecs/2),2,i+1)
        if injections and (spec in injections.keys()):
            make_axis(model_series,meas_series, ax, spec, units, 
                      df_dict["Monitors_Env"], cconv, injections[spec])
        else:
            make_axis(model_series,meas_series, ax, spec, units, 
                      df_dict["Monitors_Env"], cconv)
    
    if title:
        fig.suptitle(title)
    
    fig.tight_layout()
        
    return fig
        
        
        
