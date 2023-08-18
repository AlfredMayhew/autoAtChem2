import re
from datetime import datetime

def is_number(s):
    """Tests whether the provided argument is a number"""
    try:
        float(s)
        return True
    except ValueError:
        return False

def convert_time_to_seconds(time):
    """Converts a time string provided in the format HH:MM:SS to seconds since 
    midnight"""
    if is_number(time):
        out_time = int(time)
    elif re.match("^\d\d:\d\d:\d\d$",time):
        start_dt = datetime.strptime(time,"%H:%M:%S")
        out_time = (start_dt.hour*60*60)+(start_dt.minute*60)+(start_dt.second)
    else:
        raise Exception(f"""Time must be numeric or in the form HH:MM:SS, provided
                        time was {time}.""")
    
    return out_time

def round_to_tstep(time, start_time, time_step):
    """Rounds a given time (in seconds) to the nearest whole time step for a 
    given start time and step length"""
    t_diff = time - start_time
    
    rounded_diff = time_step * round(t_diff/time_step)
    
    return start_time+rounded_diff

def string_to_bool(string):
    """Evaluates strings of 'True' and 'False' into bools"""
    if string.casefold() == "True".casefold():
        return True
    elif string.casefold() == "False".casefold():
        return False
    else:
        raise Exception("Provided bool does not match True or False")

def datetime_to_secs_since_midnight(datetime):
    """Converts a datetime object to the time in seconds since midnight"""
    return (datetime.hour*60*60) + (datetime.minute*60) + (datetime.second)