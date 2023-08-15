#imports
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
import pandas as pd
from datetime import timedelta
from .reading_concentrations import conc_to_units
import math

def plot_species(conc_df, species, nrows = 1, ncols = None, units = None, 
                 cconv = 2.45E19, title=None, ax_size = 5, convert_xaxis = True,
                 xaxis_units = "UTC", **kwargs):
    """Creates a multi-panel figure with a selection of model time series 
    plotted on separate axes."""
    #check that ncols*nrows is less than the number of species requested
    nspecs = len(species)    
    if ncols and nrows:
        if ncols*nrows < nspecs:
            raise Exception(f"""Provided number of rows and columns 
                            ({nrows} and {ncols}) does not accommodate the 
                            requested number of species ({nspecs}).
                            You can leave either ncols or nrows unassigned
                            to automatically produce the required number of axes.""")
    #calculate the number of rows or cols (depending on which hasn't been provided)
    if nrows and not ncols:
        ncols = math.ceil(nspecs/nrows)
    elif ncols and not nrows:
        nrows = math.ceil(nspecs/ncols)
    elif not nrows and not ncols:
        raise Exception("""nrows and ncols cannot both be left unassigned.""")
    
    #check that the units list is equal in length to the species list (if 
    #assigned at all)
    if units:
        if not (len(species) == len(units)):
            raise Exception("""Number of units provided does not match the number
                            of species ({len(units)} vs. {len(species)}). Either
                            assign a unit for each species or leave 'units' 
                            unassigned to plot in the original model units.""")
    
    fig = plt.Figure(figsize = (5*ncols, 5*nrows))

    for i,spec in enumerate(species):            
        #get model data
        if spec != "NOx":
            model_series = conc_df[spec]
        else:
            model_series = conc_df[["NO","NO2"]].sum(axis=1)
        
        #convert the units
        if units:
            sel_units = units[i]
            model_series = conc_to_units(model_series, sel_units,
                                         concconversionfactor=cconv)
        else:
            sel_units = "molecules/cm3"
            
        if convert_xaxis:
            times = pd.to_datetime(model_series.index, unit="s")
        else:
            times = model_series.index.to_list()
        
        #plot        
        ax = fig.add_subplot(nrows,ncols,i+1)
        ax.plot(times, model_series, **kwargs)
        
        if convert_xaxis:
            ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
                
        ax.set_xlabel(f"Time ({xaxis_units})")
        ax.set_ylabel(f"{spec} ({sel_units})")
    
    if title:
        fig.suptitle(title)

    fig.tight_layout()
        
    return fig
        
        
        
