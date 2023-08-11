#imports
import warnings

def convert_units(in_data, current_units, concconversionfactor=2.45E+19):
    """Converts a dataframe, series, or value from specified units to the units 
    required for AtChem2 model runs"""
    if (current_units.casefold() == "molecules/cm3".casefold() or 
                            current_units.casefold() == "K".casefold() or 
                            current_units.casefold() == "kelvin".casefold() or 
                            current_units.casefold() == "hPa".casefold() or 
                            current_units.casefold() == "s-1".casefold() or 
                            current_units.casefold() == "mbar".casefold() or 
                            current_units.casefold() == "%".casefold()):
        pass #these are the units required by atchem2, so no conversion required
    
    elif ((current_units.casefold() =="ppb".casefold()) or 
          (current_units.casefold() =="ppbv".casefold())): #convert ppb to molecules cm-3
        in_data = (in_data*1E-9)*float(concconversionfactor)
    
    elif ((current_units.casefold() =="ppt".casefold()) or 
          (current_units.casefold() =="pptv".casefold())): #convert ppt to molecules cm-3
        in_data = (in_data*1E-12)*float(concconversionfactor)
    
    elif ((current_units.casefold() == "celcius".casefold()) or 
          (current_units.casefold() == "C".casefold())): #convert celcius to K
        in_data += 273.15
    
    elif current_units.casefold() == "Pa".casefold():
        in_data /= 100

    elif current_units.casefold() == "h-1".casefold(): #convert h-1 to s-1
        in_data /= (60*60)

    else:
        warnings.warn(f"""Only recognised input units are ppb, ppt, celcius, 
                      mbar, Pa, and %. You provided units of {current_units}.
                      Your constraint file will be output, 
                      but the units may not be correct for AtChem2""")
    return in_data

def conc_to_units(in_data, target_units, concconversionfactor=2.45E19):
    """Converts a dataframe, series, or value from molecules/cm3 to the units 
    specified"""
    if (target_units.casefold() == "molecules/cm3".casefold()):
        pass #no conversion required
    
    elif ((target_units.casefold() =="ppb".casefold()) or 
          (target_units.casefold() =="ppbv".casefold())): #convert ppb to molecules cm-3
        in_data = ((in_data)/float(concconversionfactor))*1E9
    
    elif ((target_units.casefold() =="ppt".casefold()) or 
          (target_units.casefold() =="pptv".casefold())): #convert ppt to molecules cm-3
        in_data = ((in_data)/float(concconversionfactor))*1E12
    else:
        warnings.warn(f"""Only recognised target units for plotting conversion 
                      are ppb, ppt, or molecules/cm3. You provided units of 
                      {target_units}. The plot will be constructed, but the 
                      units may be wrong""")
    return in_data

def closest_conc(series, tstart, units, concconversionfactor=2.45E+19, 
                  name=None):
    """selects the concentration at a given time (or closest point in time) 
    from a series"""
    na_dropped = series.dropna()
    #if tstart is present in the index, then return that value (converted)
    if tstart in na_dropped.index:
        output_val = na_dropped[tstart]
    else: #otherwise find the closest start time and notify the user
        try:
            closest_time = min(na_dropped.index, key=lambda x:abs(x-tstart))
            output_val = na_dropped[closest_time]
            warnings.warn(f"""\nTime of {tstart} not present in measurement for {name}. Closest time used of {closest_time}.""")
        except ValueError:   
            output_val = 0
            warnings.warn(f"""\nMeasurement for {name} is empty. Using a value of 0.""")
            
    
    output_val = convert_units(output_val, units, concconversionfactor)
    
    return output_val

def peak_conc(series, time, units, concconversionfactor=2.45E+19, tol_range=600):
    """selects the maximum concentration across a range for a given time
    from a series"""
    na_dropped = series.dropna()
    
    range_start = time - (tol_range/2)
    range_end = time + (tol_range/2)
    
    peak_range = na_dropped.loc[range_start:range_end]
    
    peak_conc = peak_range.max()
    
    output_val = convert_units(peak_conc, units, concconversionfactor)
    
    return output_val

def initial_conc_dict(df_dict, trans_dict, initialised_specs, tstart,
                      concconversionfactor=2.45E+19):
    """Function to make dictionary of species initial concentrations from 
    a provided dataframes and compound names"""
    out_dict = {}
    for spec in initialised_specs:
        name = trans_dict[spec][0]
        series = df_dict[trans_dict[spec][1]][name]
        units = trans_dict[spec][2]
        
        out_dict[spec] = closest_conc(series,tstart,units,concconversionfactor,
                                      name)
    
    return out_dict

