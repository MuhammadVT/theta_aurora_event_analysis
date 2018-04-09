#import matplotlib
#matplotlib.use("Agg")

def overlay_magnetometer_stations(mobj, ax, stime, station_label=False,
                                  file_dir="../data/magnetometer/"):
    """Overlays magnetometer stations on a map"""

    import numpy as np
    import pandas as pd

    # read supermag data
    fname = stime.strftime("%Y%m%d") + "_supermag.csv"
    df = pd.read_csv(file_dir + fname, parse_dates=[0]) 

    # Select the data of interest
    dfn = df.loc[df.Date_UTC == stime, :] 
    stations = [x for x in dfn.IAGA]

    # Plot the stations on a map
    X, Y = mobj(dfn.MLT.as_matrix()*15.,
                dfn.MLAT.as_matrix(),
                coords="mlt")
    mobj.scatter(X, Y, s=5.0, facecolors='r',edgecolors='r', zorder = 15)

    # Label the stations
    if station_label:
        for i in range(len(stations)):
            mobj.ax.text(X[i], Y[i], stations[i], ha='left',va='center', 
                         fontsize=10, color="r")

    return

def overlay_magnetometer_hvecs(mobj, ax, stime, rot_90_clockwise=True,
                               vecscl=1000., zorder=5,
                               file_dir="../data/magnetometer/"):
    """Overlays magnetometer stations on a map"""

    import numpy as np
    import pandas as pd
    from matplotlib.collections import LineCollection

    # read supermag data
    fname = stime.strftime("%Y%m%d") + "_supermag.csv"
    df = pd.read_csv(file_dir + fname, parse_dates=[0]) 

    # Select the data of interest
    dfn = df.loc[df.Date_UTC == stime, :] 
    dfn.loc[:, "H"] = np.sqrt(np.square(dfn.N.as_matrix()) +\
                              np.square(dfn.E.as_matrix()))
    x1,y1 = mobj(dfn.MLT.as_matrix()*15., dfn.MLAT.as_matrix(),
                 coords='mlt')
    xp,yp = mobj(0, 90, coords='mlt')
    theta_p = np.arctan2(yp-y1,xp-x1)
    theta_H = np.arctan2(dfn.E.as_matrix(), dfn.N.as_matrix())
    theta = theta_p - theta_H
    if rot_90_clockwise:
        # Rotate by 90 deg clcwise to make vecs the same as the current direction
        theta = theta - np.deg2rad(90)
    x1 = np.array(x1)
    y1 = np.array(y1)
    x2 = x1+np.array(dfn.H.tolist())*vecscl*(+1.0)*np.cos(theta)
    y2 = y1+np.array(dfn.H.tolist())*vecscl*(+1.0)*np.sin(theta)
    lines = []
    lines.extend(zip(zip(x1,y1),zip(x2,y2)))

    # Plot the vectors
    coll = LineCollection(np.array(lines), linewidths=2.0,
                          color="g", zorder=zorder)

    mobj.ax.add_collection(coll)

    return


if __name__ == "__main__":

    from davitpy import utils
    import matplotlib.pyplot as plt
    import datetime as dt
    import numpy as np
    import pandas as pd
    
    # Create a map
    coords = "mlt"
    map_lat0 = 90
    map_lon0= 0
    map_width=40*111e3
    map_height=40*111e3
    stime = dt.datetime(2002, 3, 18, 16, 0)
   
    fig, ax = plt.subplots()
    map_obj = utils.mapObj(coords=coords, projection='stere',
                           width=map_width, height=map_height,
                           lat_0=map_lat0, lon_0=map_lon0,
                           resolution="l", gridLatRes=10.,
                           datetime=stime, showCoords=True)

    overlay_magnetometer_stations(map_obj, ax, stime, station_label=False,
                                  file_dir="../data/magnetometer/")

    overlay_magnetometer_hvecs(map_obj, ax, stime, rot_90_clockwise=True,
                               vecscl=2.e3, zorder=5,
                               file_dir="../data/magnetometer/")

    plt.show()
