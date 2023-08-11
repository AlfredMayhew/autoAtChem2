#imports
import sys

args=sys.argv

def calc_cfactor(env_df,tstart,tend):
    """Function to calculate conversion factor from mixing ratio to molecules cm-3
    using ideal gas law and average the Temp and Press in a EUPHORE 2021 experiment"""
      
    #calculate average pressure over exp period
    avp=env_df["PRESSURE"].loc[tstart:tend].mean()*100 #convert from hPa to Pa
    #calculate average temperature over exp period
    avt=env_df["TEMPERATURE"].loc[tstart:tend].mean()+273.15 #convert from c to K
    #use ideal gas law to calculate conversion factor
    conv=((avp*1E-6)/(8.314*avt))*6.022E23
    
    return(conv)
