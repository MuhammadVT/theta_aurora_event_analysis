class DMSP_sat():
    """Reads and overlays data from DMSP satellite onto a map"""

    def __init__(self, dtm, sat_num, data_type="ssies", file_dir="../data/"):

        import sys
        sys.path.append('../DMSP/')
        from dmsp_ssies_read_utdallas import dmsp_ssies_read_utdallas
        import datetime as dt

        self.datetime = dtm
        self.sat_num = sat_num
        self.file_dir = file_dir

        file_dir = file_dir + data_type + "/" + dtm.strftime("%Y%m%d") + "/"
        if data_type == "ssies":
            # Read DMSP ion drift data
            df=dmsp_ssies_read_utdallas(dt.datetime(dtm.year, dtm.month, dtm.day),
                                        sat_num, file_dir=file_dir)
        if data_type == "ssj":
            pass
        self.data = df

    def overlay_ssies_data(self, mobj, ax,
                           interval=60, velscl=500.e3,
                           vec_cmap=None,
                           vel_scale=[0, 1000.],
                           quality_flag=False):

        """Overlays data onto a map"""
    
        import matplotlib
        from matplotlib.colors import ListedColormap as lcm
        import sys
        import datetime as dt
        import math
        import numpy as np
        from matplotlib.collections import PolyCollection,LineCollection
    
        # create a colormap for the quality flags of Vx and Vy 
        cmj = matplotlib.cm.jet
        cmap_flag = lcm(['k', 'y', 'r', cmj(.27)])
        bounds_flag = np.round(np.linspace(0.5, 4.5, 5))
        norm_flag = matplotlib.colors.BoundaryNorm(bounds_flag, cmap_flag.N)
    
        verts = [[],[]]
        tails_dmsp = []
        # drop NaN values for 'Vy', 'GLAT', and 'GLONG'
        df = self.data[['Vy', 'GLAT', 'GLONG', 'I']].dropna()
        df_vys = df['Vy']
        df_pos  = df[['GLAT', 'GLONG']]
        df_Is = df['I']
        df_vy = df_vys.loc[self.datetime.strftime("%Y%m%d%H%M%S"):\
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
        xxs, yys = mobj(df_lon[0], df_lat[0], coords='geo')
        xxe, yye = mobj(df_lon[-1], df_lat[-1], coords='geo')
        the_x = math.atan2(yye-yys, xxe-xxs)
        the_vel = the_x - np.deg2rad(90)

        # plot the DMSP path
        df_lats = df_pos.loc[(self.datetime-dt.\
                             timedelta(seconds=20*interval)).strftime("%Y%m%d%H%M%S"):\
                             (self.datetime+dt.timedelta(seconds=20*interval-1)).\
                             strftime("%Y%m%d%H%M%S"), 'GLAT']
        df_lons = df_pos.loc[(self.datetime-dt.\
                             timedelta(seconds=20*interval)).strftime("%Y%m%d%H%M%S"):\
                             (self.datetime+dt.timedelta(seconds=20*interval-1)).\
                             strftime("%Y%m%d%H%M%S"), 'GLONG']
        for k in range(len(df_lats)):
            x1s, y1s = mobj(df_lons[k], df_lats[k], coords='geo')
            mobj.scatter(np.array(x1s), np.array(y1s),
                         s=1.0, zorder=5, marker='o', color='gray',
                         edgecolors='face', linewidths=.5)

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

