class IMAGE_sat():
    """Reads and overlays WIC and SI data from IMAGE satellite onto a map"""

    def __init__(self, dtm, datatype="WIC"):

        import sys
        sys.path.append("../track_wic_si_spot/")
        from track_max_point import read_wic_si_data

        self.datatype = datatype
        self.datetime = dtm

        # Read WIC or SI data
        read_wic=False; read_si=False
        if datatype == "WIC":
            read_wic=True
        if datatype == "SI":
            read_si=True
        fname, data = read_wic_si_data(dtm, read_si=read_si, read_wic=read_wic)
        self.fname = fname
        self.data = data

    def overlay_data(self, mobj, param="image",
                     interpolate_data=True,
                     overlay_spot_only=False,
                     mark_max_intensity_point=False,
                     max_point_color='r',
                     max_point_size=5,
                     gauss_filter=False,
                     vmin=0, vmax=1000,
                     alpha=1.0,
                     zorder=None,
                     cmap='gist_gray'):
        """Overlays data onto a map"""

        import sys
        sys.path.append("../track_wic_si_spot/")
        from track_max_point import find_max_intensity_point
        import numpy as np
        from scipy.interpolate import griddata
        from scipy.ndimage.filters import gaussian_filter

        self.interpolate_data = interpolate_data 
        self.mark_max_intensity_point=mark_max_intensity_point
        self.param = param
        if gauss_filter:
            self.gauss_filter = gauss_filter

        X =self.data["mlt_aacgm"] * 15.
        Y =self.data["mlat_aacgm"]
        Z =self.data[param]
        X, Y = mobj(X, Y, coords=mobj.coords)

        if overlay_spot_only:
            if param == "image": 
                (maxval, maxloc) = find_max_intensity_point(self.data)
                maxloc_X = X[maxloc[1], maxloc[0]]
                maxloc_Y = Y[maxloc[1], maxloc[0]]
                mobj.scatter(maxloc_X, maxloc_Y, color=max_point_color,
                             s=max_point_size, zorder=zorder)
            else:
                pass

        else:
            if interpolate_data:
                # Interpolate the 'image' data to make it 256 by 256 
                Xnew = np.linspace(0, 60*111.e3, 256)   
                Ynew = np.linspace(0, 60*111.e3, 256)   
                Xnew, Ynew = np.meshgrid(Xnew, Ynew) 
                 
                # Replace the missing values
                Xtmp = X; Ytmp = Y
                Xtmp[np.where(Xtmp > 100*111.e3)] = 30*111.e3
                Ytmp[np.where(Ytmp > 100*111.e3)] = 0*111.e3
                Znew = griddata((Xtmp.flatten(), Ytmp.flatten()), Z.flatten(),
                                (Xnew.flatten(), Ynew.flatten()), method='linear',
                                fill_value=0.0)
                Znew = np.reshape(Znew, (256, 256))
                if gauss_filter:
                    Znew = gaussian_filter(Znew, 1)
                mappable = mobj.pcolor(Xnew, Ynew, Znew, vmin=vmin, vmax=vmax,\
                                       alpha=1.0, zorder=zorder, cmap=cmap)
            else:
                #mobj.pcolormesh(X, Y, Z, vmin=vmin, vmax=vmax,
                #                alpha=1.0, cmap=cmap)
                if gauss_filter:
                    Z = gaussian_filter(Z, 1)
                mappable = mobj.pcolor(X, Y, Z, vmin=vmin, vmax=vmax,
                                       alpha=1.0, zorder=zorder, cmap=cmap)
            
            self.mappable = mappable

            #mobj.scatter(333e4, 333e4, color="r")
            #mobj.scatter(1e3, 1e3, color="r")
            #mobj.scatter(X.flatten(), Y.flatten(), color="g", s=0.2)
            #mobj.scatter(Xnew.flatten(), Ynew.flatten(), color="r", s=0.2)
            
            # Mark the max intensity point
            if mark_max_intensity_point:
                if param == "image": 
                    (maxval, maxloc) = find_max_intensity_point(self.data)
                    maxloc_X = X[maxloc[1], maxloc[0]]
                    maxloc_Y = Y[maxloc[1], maxloc[0]]
                    mobj.scatter(maxloc_X, maxloc_Y, color=max_point_color,
                                 s=max_point_size, zorder=zorder)
                else:
                    pass

        return

    def add_cbar(self, fig, ax, shrink=0.6,
                 label="Rayleigh", label_fontsize=15):
        cbar = fig.colorbar(mappable=self.mappable,
                            ax=ax, shrink=shrink)
        cbar.set_label(label, size=label_fontsize)
        self.cbar = cbar
        return

