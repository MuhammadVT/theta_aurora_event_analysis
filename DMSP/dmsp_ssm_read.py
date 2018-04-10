import pandas as pd
import numpy as np
import datetime as dt
from dateutil.relativedelta import relativedelta

def dmsp_ssm_read_mfr(stime, sat_num, file_dir='./data/', drop_microsecond=True):
    '''reads DMSP SSM data of a single satellite for a given date

    NOTE: Spacecraft  coordinates  are  defined  where  x  is  along  the  local
          vertical  measured  positive  in  the  downward  direction,
          z  is  perpendicular  to  both  the  local  vertical  and  
          the  spacecraft  velocity  vector measured positive in 
          the anti -orbit normal direction (in approximately the anti -sunward direction)
          , and y completes a right hand coordinate system with positive y
          in the same general direction as the spacecraft velocity direction. 
          Note that the positive y direction is not the same as the spacecraft velocity direction 
          because in general the spacecraft velocity vector is not perpendicular to
          the local vertical; however, the local vertical direction, 
          the spacecraft velocity vector, and the y axis are always coplanar. 

    '''

    fname = "mfr_" + stime.strftime("%Y%m%d") + "_F" + str(sat_num) + ".dat"
    columns = ["Date-time", "Sec-of-day", "MissnID", "Avg/Smp", "Lat",
               "Lon", "Alt", "EphSrc", "1st/Dif", "TotX", "TotY", "TotZ", 
               "CalDate", "Meas-ModX", "Meas-ModY", "Meas-ModZ",
               "ModDate", "ModTyp"]

    df = pd.read_csv(file_dir+fname, skiprows=6, delim_whitespace=True,
                     header=None, names=columns, parse_dates=[0],
                     date_parser=lambda x: pd.datetime.strptime(x, "%Y%j%H%M%S.%f"))

    # Sample the data at 1 second
    if drop_microsecond:
        df.loc[:, "Date-time"] = df.loc[:, "Date-time"].apply(lambda x: x.replace(microsecond=0))
    return df

if __name__ == "__main__":
    import datetime as dt
    import matplotlib.pyplot as plt
    stime = dt.datetime(2002,3,18,16,12)
    sat_num = 13
    file_dir = '../data/ssm/' + stime.strftime("%Y%m%d") + "/"
    df=dmsp_ssm_read_mfr(stime, sat_num, file_dir=file_dir)

