import matplotlib
matplotlib.use("Agg")

def overlay_guvi_JE_E0_JNE(mobj, stime, etime, orbit,
                           color_norm=None, cmap="jet", 
                           data_level="L3", overlay_param="JEe",
                           hemi="north", coords="mlt"):
    """Plots energy flux, average energy, and number flux from
       TIMED-GUVI data on a map"""

    import numpy as np
    import pandas as pd
    from read_timed_guvi_aurora_data import read_guvi_aurora_L3_data

    # read GUVI data
    if data_level == "L3":
        df = read_guvi_aurora_L3_data(stime, etime, orbit, file_name=None)

    # Plot the data on a map
    X, Y = mobj(df.loc[:, "mlt_"+hemi].as_matrix()*15.,
                df.loc[:, "mlat_"+hemi].as_matrix(),
                coords=coords)
    
    var = df.loc[:, overlay_param+"_"+hemi].as_matrix()
    mappable = mobj.scatter(X, Y, c=var, s=5.0, marker="s", norm=color_norm, cmap=cmap)
    
    return mappable


def overlay_ssj_JE_E0_JNE(mobj, stime, etime, sat_num,
                          color_norm=None, cmap="jet", 
                          overlay_param="JEe", aacgm=True):
    """Plots energy flux, average energy, and number flux from
       DMSP SSJ data on a map"""
    
    import numpy as np
    import pandas as pd
    from dmsp_ssj_read import dmsp_ssj_read
    import math
    import matplotlib
    from matplotlib.colors import ListedColormap as lcm
    from matplotlib.collections import PolyCollection,LineCollection

    # read ssj data
    stm = stime - dt.timedelta(seconds=20*60)
    etm = etime + dt.timedelta(seconds=20*60)
    df = dmsp_ssj_read(stm, etm, sat_num, file_name=None)

    # Overlay the SSJ data
    df_pos  = df[['GLAT', 'GLON', 'MLAT', 'MLT']]
    if aacgm:
        df_lats = df_pos.loc[:, 'GLAT'].as_matrix()
        df_lons = df_pos.loc[:, 'GLON'].as_matrix()
        df_lat = df_pos.loc[(df.datetime >= stime) & (df.datetime <= etime), 'GLAT'].as_matrix()
        df_lon = df_pos.loc[(df.datetime >= stime) & (df.datetime <= etime), 'GLON'].as_matrix()
        tmp_coords = "geo"
    else:
        df_lats = df_pos.loc[:, 'MLAT'].as_matrix()
        df_lons = df_pos.loc[:, 'MLT'].as_matrix() * 15.
        df_lat = df_pos.loc[(df.datetime >= stime) & (df.datetime <= etime), 'MLAT'].as_matrix()
        df_lon = df_pos.loc[(df.datetime >= stime) & (df.datetime <= etime), 'MLT'].as_matrix() * 15.
        tmp_coords = mobj.coords

    df_var = df.loc[(df.datetime >= stime) & (df.datetime <= etime), overlay_param]

    # plot the DMSP path
    x1s, y1s = mobj(df_lons, df_lats, coords=tmp_coords)
    mobj.scatter(x1s, y1s,
                 s=1.0, zorder=5, marker='o', color='gray',
                 edgecolors='face', linewidths=.5)


    # plot values at the time interval of interest
    x1, y1 = mobj(df_lon, df_lat, coords=tmp_coords)
    mappable = mobj.scatter(x1, y1,
                            s=9.0,zorder=5,marker='o',
                            c=df_var.as_matrix(),
                            edgecolors='face', linewidths=.5,
                            norm=color_norm, cmap=cmap)
    return mappable


if __name__ == "__main__":

    import datetime as dt
    import matplotlib.pyplot as plt
    from davitpy.utils import plotUtils
    import numpy as np
    from matplotlib import colors

    #stime = dt.datetime(2002,3,18,17,24)
    #etime = dt.datetime(2002,3,18,17,40)
    stime = dt.datetime(2002,3,18,16,15)
    etime = dt.datetime(2002,3,18,16,58)

    dmsp_sat_num = [13, 15]
    guvi_orbit = 1495

    cmap = "jet"
    coords = "mlt"
    overlay_param="JNEe"
    dmsp_aacgm=False
    if overlay_param == "JEe":
        color_norm = colors.LogNorm(vmin=0.1, vmax=1.e2)
        guvi_cbar_label = "GUVI Flux " + r"(ergs/s/cm$^{2}$)"
    if overlay_param == "AvgEe":
        color_norm = colors.LogNorm(vmin=0.5, vmax=15)
        guvi_cbar_label = "GUVI E0 (KeV)"
    if overlay_param == "JNEe":
        color_norm = colors.LogNorm(vmin=1.e7, vmax=1.e10)
        guvi_cbar_label = "GUVI Number Flux " + "r(#/s/cm$^{2}$)"

    # Plot a map
    fig, ax = plt.subplots(figsize=(15,15))
    mobj = plotUtils.mapObj(ax=ax, datetime=stime, lon_0=0., boundinglat=50.,
                            gridLabels=False, gridLatRes=10., coords=coords,
                            fillContinents='None')
    # Overlay DMSP SSJ Data
    for sat in dmsp_sat_num:
        ssj_mappable = overlay_ssj_JE_E0_JNE(mobj, stime, etime, sat,
                                             color_norm=color_norm, cmap=cmap, 
                                             overlay_param=overlay_param, aacgm=dmsp_aacgm)

    # Overlay TIMED-GUVI data
    guvi_mappable = overlay_guvi_JE_E0_JNE(mobj, stime, etime, guvi_orbit,
                                          color_norm=color_norm, cmap=cmap, 
                                          data_level="L3", overlay_param=overlay_param, 
                                          hemi="north", coords=coords)
    # Add a colorbar
    cbar = fig.colorbar(mappable=guvi_mappable, ax=ax, shrink=0.6)
    cbar.set_label(guvi_cbar_label, size=15.)
    
    title = stime.strftime("%b %d, %Y")  +   "       Orbit: " + str(guvi_orbit)
    ax.set_title(title)

    #title = stime.strftime("%b%d, %Y") + "    DMSP F" + str(sat_num)

    # save fig
    fig_dir = "../plots/multiinstrument_data_overlay/JE_E0_JNE/"
    fig_name = "APEX_" + overlay_param + "_" + stime.strftime("%Y%m%d.%H%M%S") + "_" +\
               etime.strftime("%Y%m%d.%H%M%S") + ".png"

    fig.savefig(fig_dir + fig_name, dpi=200, bbox_inches="tight")


