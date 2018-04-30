class DMSP_sat():
    """Reads and overlays data from DMSP satellite onto a map"""

    def __init__(self, dtm, sat_num, data_type="ssies", file_dir="../data/"):

        import sys
        sys.path.append('../DMSP/')
        from dmsp_ssies_read_utdallas import dmsp_ssies_read_utdallas
        import datetime as dt

        self.datetime = dtm
        self.sat_num = sat_num

        file_dir = file_dir + data_type + "/" + dtm.strftime("%Y%m%d") + "/"
        self.file_dir = file_dir
        if data_type == "ssies":
            # Read DMSP ion drift data
            df=dmsp_ssies_read_utdallas(dt.datetime(dtm.year, dtm.month, dtm.day),
                                        sat_num, file_dir=file_dir)
        if data_type == "ssj":
            pass
        if data_type == "ssm":
            # Read SSM data from MFR data files obtained from NGDC NOAA
            from dmsp_ssm_read import dmsp_ssm_read_mfr
            df = dmsp_ssm_read_mfr(dtm, sat_num, file_dir=file_dir)

        self.data = df

    def overlay_ssies_data(self, mobj, ax,
                           interval=60, velscl=500.e3,
                           vec_cmap=None,
                           vel_scale=[0, 1000.],
                           plot_path=True,
                           plot_vecs_on_full_path=False,
                           quality_flag=False):

        """Overlays data onto a map"""
    
        import matplotlib
        from matplotlib.colors import ListedColormap as lcm
        import sys
        import datetime as dt
        import math
        import numpy as np
        from matplotlib.collections import LineCollection
    
        # create a colormap for the quality flags of Vx and Vy 
        cmj = matplotlib.cm.jet
        cmap_flag = lcm(['k', 'y', 'r', cmj(.27)])
        bounds_flag = np.round(np.linspace(0.5, 4.5, 5))
        norm_flag = matplotlib.colors.BoundaryNorm(bounds_flag, cmap_flag.N)
    
        # drop NaN values for 'Vy', 'GLAT', and 'GLONG'
        df = self.data[['Vy', 'GLAT', 'GLONG', 'I']].dropna()
        df_vys_tmp = df['Vy']
        df_pos  = df[['GLAT', 'GLONG']]
        df_Is = df['I']
        df_vy = df_vys_tmp.loc[self.datetime.strftime("%Y%m%d%H%M%S"):\
                           (self.datetime+dt.timedelta(seconds=interval-1)).\
                           strftime("%Y%m%d%H%M%S")]
        df_lat = df_pos.loc[self.datetime.strftime("%Y%m%d%H%M%S"):\
                            (self.datetime+dt.timedelta(seconds=interval-1)).\
                            strftime("%Y%m%d%H%M%S"), 'GLAT']
        df_lon = df_pos.loc[self.datetime.strftime("%Y%m%d%H%M%S"):\
                            (self.datetime+dt.timedelta(seconds=interval-1)).\
                            strftime("%Y%m%d%H%M%S"), 'GLONG']
        df_I = df_Is.loc[self.datetime.strftime("%Y%m%d%H%M%S"):\
                         (self.datetime+dt.timedelta(seconds=interval-1)).\
                         strftime("%Y%m%d%H%M%S")]
        try:
            xxs, yys = mobj(df_lon[0], df_lat[0], coords='geo')
            xxe, yye = mobj(df_lon[-1], df_lat[-1], coords='geo')
        except IndexError:
            return
        the_x = math.atan2(yye-yys, xxe-xxs)
        the_vel = the_x - np.deg2rad(90)

        # Get the locs along the full path over polar region
        df_lats = df_pos.loc[(self.datetime-dt.\
                             timedelta(seconds=20*interval)).strftime("%Y%m%d%H%M%S"):\
                             (self.datetime+dt.timedelta(seconds=20*interval-1)).\
                             strftime("%Y%m%d%H%M%S"), 'GLAT']
        df_lons = df_pos.loc[(self.datetime-dt.\
                             timedelta(seconds=20*interval)).strftime("%Y%m%d%H%M%S"):\
                             (self.datetime+dt.timedelta(seconds=20*interval-1)).\
                             strftime("%Y%m%d%H%M%S"), 'GLONG']
        df_vys = df_vys_tmp.loc[(self.datetime-dt.\
                                 timedelta(seconds=20*interval)).strftime("%Y%m%d%H%M%S"):\
                                 (self.datetime+dt.timedelta(seconds=20*interval-1)).\
                                 strftime("%Y%m%d%H%M%S")]

        # plot the DMSP path
        if plot_path:
            for k in range(len(df_lats)):
                x1s, y1s = mobj(df_lons[k], df_lats[k], coords='geo')
                mobj.scatter(np.array(x1s), np.array(y1s),
                             s=1.0, zorder=5, marker='o', color='gray',
                             edgecolors='face', linewidths=.5)

        if plot_vecs_on_full_path:
            verts = [[],[]]
            tails_dmsp = []

            # Plot vectors along the full path
            for k in range(len(df_lats)):
                x1, y1 = mobj(df_lons[k], df_lats[k], coords='geo')
                verts[0].append(x1)
                verts[1].append(y1)
                x2 = x1+df_vys[k]*velscl*(-1.0)*math.cos(the_vel)
                y2 = y1+df_vys[k]*velscl*(-1.0)*math.sin(the_vel)
                tails_dmsp.append(((x1,y1),(x2,y2)))
        
            lcoll = LineCollection(np.array(tails_dmsp),
                                   linewidths=.6,
                                   zorder=12, alpha=0.6, color="b")
            ax.add_collection(lcoll)

        verts = [[],[]]
        tails_dmsp = []
        # plot the measurement points and velocity vectors at a specified time
        for k in range(len(df_lat)):
            x1, y1 = mobj(df_lon[k], df_lat[k], coords='geo')
            verts[0].append(x1)
            verts[1].append(y1)
            x2 = x1+df_vy[k]*velscl*(-1.0)*math.cos(the_vel)
            y2 = y1+df_vy[k]*velscl*(-1.0)*math.sin(the_vel)
            tails_dmsp.append(((x1,y1),(x2,y2)))
    
        xx = mobj.scatter(np.array(verts[0]),np.array(verts[1]),
                          s=2.5,zorder=5,marker='o',
                          c=np.abs(df_vy.as_matrix()),
                          vmin=vel_scale[0], vmax=vel_scale[1],
                          edgecolors='face', linewidths=.5, 
                          cmap=vec_cmap)
        if quality_flag:
            lcoll = LineCollection(np.array(tails_dmsp), linewidths=.6,
                                   zorder=12, cmap=cmap_flag,
                                   norm=norm_flag, alpha=1)
            lcoll.set_array(df_I)
        else:
            lcoll = LineCollection(np.array(tails_dmsp),
                                   linewidths=.6,
                                   zorder=12, alpha=1,
                                   cmap=vec_cmap)
            lcoll.set_array(np.abs(df_vy.as_matrix()))
            lcoll.set_clim(vmin=vel_scale[0], vmax=vel_scale[1])
        ax.add_collection(lcoll)
    
        return

    def overlay_ssm_data(self, mobj, ax, rot_clockwise=90,
                         interval=60, velscl=500.e3,
                         plot_path=False,
                         plot_vecs_on_full_path=False,
                         correct_bias=False):

        """Overlays SSM data onto a map"""
    
        import matplotlib
        import datetime as dt
        import math
        import numpy as np
        import pandas as pd
        from matplotlib.collections import LineCollection
     
        df = self.data[['Date-time', 'Lat', 'Lon', 'Meas-ModY', 'Meas-ModZ']].dropna()
        df.set_index('Date-time', inplace=True)

        # Remove duplicated indices
        df = df[~df.index.duplicated(keep="first")]

        # Correct the bias in SSM data
        if correct_bias:
            fln = "bias_corrected_mfr_" + self.datetime.strftime("%Y%m%d") +\
                  "_F" + str(self.sat_num) + ".dat" 
            dfc = pd.read_csv(self.file_dir+fln, index_col=0)
            # Remove duplicated indices
            dfc = dfc[~dfc.index.duplicated(keep="first")]
            df.loc[pd.to_datetime(dfc.index), 'Meas-ModY'] = dfc.loc[:, 'Meas-ModY'].as_matrix()
            df.loc[pd.to_datetime(dfc.index), 'Meas-ModZ'] = dfc.loc[:, 'Meas-ModZ'].as_matrix()

        df.loc[:, "H"] = np.sqrt(np.square(df['Meas-ModY'].as_matrix()) +\
                            np.square(df['Meas-ModZ'].as_matrix()))
        df.loc[:,"theta_H_to_Zdir"] = np.arctan2(df['Meas-ModZ'].as_matrix(),
                                                 df['Meas-ModY'].as_matrix())

        # Select data for the time of interval
        df_H = df.loc[self.datetime.strftime("%Y%m%d%H%M%S"):\
                      (self.datetime+dt.timedelta(seconds=interval-1)).\
                      strftime("%Y%m%d%H%M%S"), "H"]
        df_Htheta = df.loc[self.datetime.strftime("%Y%m%d%H%M%S"):\
                           (self.datetime+dt.timedelta(seconds=interval-1)).\
                           strftime("%Y%m%d%H%M%S"), "theta_H_to_Zdir"]
        df_lat = df.loc[self.datetime.strftime("%Y%m%d%H%M%S"):\
                        (self.datetime+dt.timedelta(seconds=interval-1)).\
                        strftime("%Y%m%d%H%M%S"), 'Lat']
        df_lon = df.loc[self.datetime.strftime("%Y%m%d%H%M%S"):\
                        (self.datetime+dt.timedelta(seconds=interval-1)).\
                        strftime("%Y%m%d%H%M%S"), 'Lon']
        try:
            xxs, yys = mobj(df_lon[0], df_lat[0], coords='geo')
            xxe, yye = mobj(df_lon[-1], df_lat[-1], coords='geo')
        except IndexError:
            return
        theta_y = math.atan2(yye-yys, xxe-xxs)
        theta_H = theta_y - df_Htheta.as_matrix()

        # Rotate clockwise
        theta_H = theta_H - np.deg2rad(rot_clockwise)

        # Get the locs along the full path over polar region
        df_lats = df.loc[(self.datetime-dt.\
                         timedelta(seconds=20*interval)).strftime("%Y%m%d%H%M%S"):\
                         (self.datetime+dt.timedelta(seconds=20*interval-1)).\
                         strftime("%Y%m%d%H%M%S"), 'Lat']
        df_lons = df.loc[(self.datetime-dt.\
                         timedelta(seconds=20*interval)).strftime("%Y%m%d%H%M%S"):\
                         (self.datetime+dt.timedelta(seconds=20*interval-1)).\
                         strftime("%Y%m%d%H%M%S"), 'Lon']
        df_Hs = df.loc[(self.datetime-dt.\
                       timedelta(seconds=20*interval)).strftime("%Y%m%d%H%M%S"):\
                       (self.datetime+dt.timedelta(seconds=20*interval-1)).\
                       strftime("%Y%m%d%H%M%S"), "H"]
        df_Hthetas = df.loc[(self.datetime-dt.\
                            timedelta(seconds=20*interval)).strftime("%Y%m%d%H%M%S"):\
                            (self.datetime+dt.timedelta(seconds=20*interval-1)).\
                            strftime("%Y%m%d%H%M%S"), "theta_H_to_Zdir"]

        # plot the DMSP path
        if plot_path:
            for k in range(len(df_lats)):
                x1s, y1s = mobj(df_lons[k], df_lats[k], coords='geo')
                mobj.scatter(np.array(x1s), np.array(y1s),
                             s=1.0, zorder=5, marker='o', color='gray',
                             edgecolors='face', linewidths=.5)

        if plot_vecs_on_full_path:
            verts = [[],[]]
            tails_dmsp = []
            theta_Hs = theta_y - df_Hthetas.as_matrix()
            # Plot vectors along the full path
            for k in range(len(df_lats)):
                x1, y1 = mobj(df_lons[k], df_lats[k], coords='geo')
                verts[0].append(x1)
                verts[1].append(y1)
                x2 = x1+df_Hs[k]*velscl*(+1.0)*math.cos(theta_Hs[k])
                y2 = y1+df_Hs[k]*velscl*(+1.0)*math.sin(theta_Hs[k])
                tails_dmsp.append(((x1,y1),(x2,y2)))
        
            lcoll = LineCollection(np.array(tails_dmsp),
                                   linewidths=0.6, color="g",
                                   zorder=12, alpha=0.6)
            ax.add_collection(lcoll)

        verts = [[],[]]
        tails_dmsp = []
        # plot the measurement points and velocity vectors at a specified time
        for k in range(len(df_lat)):
            x1, y1 = mobj(df_lon[k], df_lat[k], coords='geo')
            verts[0].append(x1)
            verts[1].append(y1)
            x2 = x1+df_H[k]*velscl*(+1.0)*math.cos(theta_H[k])
            y2 = y1+df_H[k]*velscl*(+1.0)*math.sin(theta_H[k])
            tails_dmsp.append(((x1,y1),(x2,y2)))
    
        xx = mobj.scatter(np.array(verts[0]),np.array(verts[1]),
                          s=2.5,zorder=5,marker='o',
                          color="r", edgecolors='face', linewidths=.5)
        lcoll = LineCollection(np.array(tails_dmsp),
                               linewidths=0.6, color="g",
                               zorder=12, alpha=1)
        ax.add_collection(lcoll)
    
        return

if __name__ == "__main__":

    import datetime as dt
    import matplotlib.pyplot as plt
    from davitpy.utils import plotUtils
    
    dtm = dt.datetime(2002,3,18,16,50)
    DMSP_sat_num = 13
    data_type = "ssm"
    #data_type = "ssies"
    vmin=200.; vmax=1000.
    coords = "mlt"
    
    fig, ax = plt.subplots(figsize=(12,8))
    #ax.set_facecolor('black')

    # Plot a map
    mobj = plotUtils.mapObj(ax=ax, datetime=dtm, lon_0=0., boundinglat=60.,
                            gridLabels=False, gridLatRes=10., coords=coords,
                            fillContinents='None')
    DMSP_obj = DMSP_sat(dtm, DMSP_sat_num, data_type=data_type)

    if data_type == "ssies":
        DMSP_obj.overlay_ssies_data(mobj, ax,
                                    interval=2*60, velscl=3*111.,
                                    vec_cmap="jet",
                                    vel_scale=[0, 1000.],
                                    quality_flag=False)

    if data_type == "ssm":
        DMSP_obj.overlay_ssm_data(mobj, ax, interval=60, velscl=1.e3,
                                  rot_clockwise=90)

    plt.show()




