from pysolar.solar import get_altitude
import math
from datetime import datetime,time, timedelta, timezone
import pandas as pd

def J_Calc(x,l,m,n):
    """Calculate a J value for a given sza"""
    cosx = math.cos(x)
    
    if cosx < 0:
        return 0
    
    j = l * (cosx**m)*math.exp(-n*(1/cosx))
    
    return j.real

def JNO2_Calc(lat, long, dt):
    """Calculate JNO2 value for a given lat, long, and datetime"""
    sza = math.radians(float(90)-get_altitude(lat, long, dt))    
    
    j_val = J_Calc(sza, 1.165e-02, 0.244, 0.267)
    
    return j_val

def calcJFAC_list(jNO2_series, date, lat, long, cutoff = 0.005):
    """Produces a series of JFAC values for a given series of measured JNO2 values"""  
    #get index as datetimes for calculation of JNO2
    dt_index = [datetime.combine(date,time()) + timedelta(seconds = x) for x in jNO2_series.index]
    dt_index = [d.replace(tzinfo = timezone.utc) for d in dt_index]
    
    calc_jno2s = [JNO2_Calc(lat,long, x) for x in dt_index]
    
    jfacs = jNO2_series.values/calc_jno2s
    
    test_cutoff = lambda x : x if x>cutoff else 0
    
    jfacs = [test_cutoff(x) for x in jfacs]
    
    return pd.Series(jfacs,index = jNO2_series.index)
    
    
    