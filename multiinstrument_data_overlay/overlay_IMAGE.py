class IMAGE_sat():
    """Reads and overlays WIC and SI data from IMAGE satellite onto a map"""

    def __init__(self, dtm, datatype="WIC", wic_si_ratio=False, si13_si_ratio=False):

        import sys
        sys.path.append("../track_wic_si_spot/")
        from track_max_point import read_wic_si_data

        """
        Note: if wic_si_ratio=True, datatype will be ignored
        """

        self.datatype = datatype
        self.datetime = dtm

        # Read WIC or SI data
        read_wic=False; read_si=False; read_si13=False
        if datatype == "WIC":
            read_wic=True
        if datatype == "SI":
            read_si=True
        if datatype == "SI13":
            read_si13=True
        if wic_si_ratio:
            fname_wic, data_wic = read_wic_si_data(dtm, read_si=False, read_wic=True, read_si13=False)
            fname_si, data_si = read_wic_si_data(dtm, read_si=True, read_wic=False, read_si13=False)
            self.data_wic = data_wic
            self.data_si = data_si
            self.fname_wic = fname_wic
            self.fname_si = fname_si
        elif si13_si_ratio:
            fname_wic, data_wic = read_wic_si_data(dtm, read_si=False, read_si13=True, read_wic=False)
            fname_si, data_si = read_wic_si_data(dtm, read_si=True, read_si13=False, read_wic=False)
            self.data_wic = data_wic
            self.data_si = data_si
            self.fname_wic = fname_wic
            self.fname_si = fname_si

        else:
            fname, data = read_wic_si_data(dtm, read_si=read_si, read_si13=read_si13, read_wic=read_wic)
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
                                fill_value=np.nan)
                Znew = np.reshape(Znew, (256, 256))

                if gauss_filter:
                    Znew = gaussian_filter(Znew, 1)

                # Remove points whose intensities are below vmin
                Znew[np.where(Znew < vmin)] = np.nan

                #mappable = mobj.pcolor(Xnew, Ynew, Znew, vmin=vmin, vmax=vmax,\
                #                       alpha=alpha, zorder=zorder, cmap=cmap)
                mappable = mobj.scatter(Xnew, Ynew, c=Znew, marker="s", s=1,
                                        vmin=vmin, vmax=vmax,
                                        alpha=alpha, zorder=zorder, cmap=cmap)
            else:
                #mobj.pcolormesh(X, Y, Z, vmin=vmin, vmax=vmax,
                #                alpha=1.0, cmap=cmap)
                if gauss_filter:
                    Z = gaussian_filter(Z, 1)
                # Remove points whose intensities are below vmin
                Z[np.where(Z < vmin)] = np.nan

                mappable = mobj.pcolor(X, Y, Z, vmin=vmin, vmax=vmax,
                                       alpha=alpha, zorder=zorder, cmap=cmap)
            
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

    def overlay_wic_si_ratio(self, mobj, param="image",
                             interpolate_data=False,
                             gauss_filter=False,
                             mark_max_intensity_point=False,
                             wic_si_spot_max_intensity_ratio_only=True,
                             normalize=False,
                             vmin=0, vmax=1000.,
                             vmin_ratio=0, vmax_ratio=1.,
                             alpha=1.0,
                             zorder=None,
                             norm_factor_fname=None,
                             cmap='jet'):
        """Overlays WIC/SI or SI13/SI12 ratio onto a map"""

        import sys
        import numpy as np
        import pandas as pd
        from scipy.interpolate import griddata
        from scipy.ndimage.filters import gaussian_filter
        import sys
        sys.path.append("../track_wic_si_spot/")
        from track_max_point import find_max_intensity_point

        self.param = param
        if gauss_filter:
            self.gauss_filter = gauss_filter


        X_wic =self.data_wic["mlt_aacgm"] * 15.
        Y_wic =self.data_wic["mlat_aacgm"]
        Z_wic =self.data_wic[param]
        X_wic, Y_wic = mobj(X_wic, Y_wic, coords=mobj.coords)

        X_si =self.data_si["mlt_aacgm"] * 15.
        Y_si =self.data_si["mlat_aacgm"]
        Z_si =self.data_si[param]
        X_si, Y_si = mobj(X_si, Y_si, coords=mobj.coords)


        if not wic_si_spot_max_intensity_ratio_only:
            if interpolate_data:
                # Interpolate the 'image' data to make it 256 by 256 
                Xnew = np.linspace(0, 60*111.e3, 256)   
                Ynew = np.linspace(0, 60*111.e3, 256)   
                Xnew, Ynew = np.meshgrid(Xnew, Ynew) 
                 
                # Replace the missing values
                Xtmp_wic = X_wic; Ytmp_wic = Y_wic
                Xtmp_wic[np.where(Xtmp_wic > 100*111.e3)] = 30*111.e3
                Ytmp_wic[np.where(Ytmp_wic > 100*111.e3)] = 0*111.e3
                Znew_wic = griddata((Xtmp_wic.flatten(), Ytmp_wic.flatten()), Z_wic.flatten(),
                                (Xnew.flatten(), Ynew.flatten()), method='linear',
                                fill_value=np.nan)
                Znew_wic = np.reshape(Znew_wic, (256, 256))

                Xtmp_si = X_si; Ytmp_si = Y_si
                Xtmp_si[np.where(Xtmp_si > 100*111.e3)] = 30*111.e3
                Ytmp_si[np.where(Ytmp_si > 100*111.e3)] = 0*111.e3
                Znew_si = griddata((Xtmp_si.flatten(), Ytmp_si.flatten()), Z_si.flatten(),
                                (Xnew.flatten(), Ynew.flatten()), method='linear',
                                fill_value=np.nan)
                Znew_si = np.reshape(Znew_si, (256, 256))
            else:
                # No need to interpolate is SI13/SI12 ratio is needed
                # Because they both have the same shape of (128, 128)
                Xnew = X_si
                Ynew = Y_si
                Znew_wic = Z_wic
                Znew_si = Z_si

            if gauss_filter:
                Znew_wic = gaussian_filter(Znew_wic, 1)
                Znew_si = gaussian_filter(Znew_si, 1)

#            # Remove points whose intensities are very low
#            Znew_wic[np.where(Znew_wic < 0.1)] = np.nan
#            Znew_si[np.where(Znew_si < 0.1)] = np.nan

            # normalize
            if normalize:
                #Znew_wic = Znew_wic/Znew_wic[~np.isnan(Znew_wic)].max()
                #Znew_si = Znew_si/Znew_si[~np.isnan(Znew_si)].max()
                if norm_factor_fname:
                    # Read from a file
                    df_factor = pd.read_csv(norm_factor_fname, index_col=0)
                    df_factor.set_index("datetime", inplace=True)
                    factor = df_factor.wic.sum()/df_factor.si.sum()
                else:
                    factor = Znew_wic[~np.isnan(Znew_wic)].sum() / Znew_si[~np.isnan(Znew_si)].sum()
                Znew_si = Znew_si * factor

            # Remove points whose intensities are below vmin
            Znew_wic[np.where(Znew_wic < vmin)] = np.nan
            Znew_si[np.where(Znew_si < vmin)] = np.nan

            # Take the ratio of WIC/SI
            Znew = Znew_wic / Znew_si

            #mappable = mobj.pcolor(Xnew, Ynew, Znew, vmin=vmin_ratio, vmax=vmax_ratio,\
            #                       alpha=alpha, zorder=zorder, cmap=cmap)
            mappable = mobj.scatter(Xnew, Ynew, c=Znew, marker="s", s=1,
                                    vmin=vmin_ratio, vmax=vmax_ratio,
                                    alpha=alpha, zorder=zorder, cmap=cmap)
                    
            self.mappable = mappable

            #mobj.scatter(X.flatten(), Y.flatten(), color="g", s=0.2)
            #mobj.scatter(Xnew.flatten(), Ynew.flatten(), color="r", s=0.2)

            # Mark the max intensity point
            if mark_max_intensity_point:
                if param == "image": 
                    (maxval_wic, maxloc) = find_max_intensity_point(self.data_wic)
                    maxloc_X = X_wic[maxloc[1], maxloc[0]]
                    maxloc_Y = Y_wic[maxloc[1], maxloc[0]]
                    mobj.scatter(maxloc_X, maxloc_Y, color="r",
                                 s=30., zorder=zorder)

                    (maxval_si, maxloc) = find_max_intensity_point(self.data_si)
                    maxloc_X = X_si[maxloc[1], maxloc[0]]
                    maxloc_Y = Y_si[maxloc[1], maxloc[0]]
                    mobj.scatter(maxloc_X, maxloc_Y, color="b",
                                 s=30., zorder=zorder)

                else:
                    pass

        if wic_si_spot_max_intensity_ratio_only:
            if param == "image": 
                (maxval, maxloc_wic) = find_max_intensity_point(self.data_wic)
                (maxval, maxloc_si) = find_max_intensity_point(self.data_si)
                # Assign the center location to WIC or SI13 spot
                maxloc_X = X_wic[maxloc_wic[1], maxloc_wic[0]]
                maxloc_Y = Y_wic[maxloc_wic[1], maxloc_wic[0]]
                # Find WIC/SI12 or SI13/SI12 spot ratio
                maxval_wic = self.data_wic['image'][maxloc_wic[1], maxloc_wic[0]]
                maxval_si = self.data_si['image'][maxloc_si[1], maxloc_si[0]]
                wic_si_factor = maxval_wic/maxval_si
                # Plot on a map
                mappable = mobj.scatter(maxloc_X, maxloc_Y, c=wic_si_factor,
                                        marker="o", s=50,
                                        vmin=vmin_ratio, vmax=vmax_ratio,
                                        alpha=alpha, zorder=zorder, cmap=cmap)
                self.mappable = mappable
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

