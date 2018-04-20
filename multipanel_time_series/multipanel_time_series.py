from matplotlib.dates import DateFormatter 
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import datetime as dt
import numpy as np 
import pandas as pd

class multipanel_plot():
    """ Makes multipanel time series plots of an event """

    def __init__(self, stime, etime, figsize=None, nrows=6, ncols=1,
                 hspace=None, sharex=True):

        self.stime = stime
        self.etime = etime

        # Draw a figure and axes
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize,
                                 sharex=sharex)
        if hspace is not None:
            fig.subplots_adjust(hspace=hspace)

        self.figsize = figsize
        self.nrows=nrows
        self.ncols = ncols
        self.fig = fig
        self.axes = axes

    def plot_ace(self, ax_Vx=None, ax_Pdynx=None, ax_IMF=None,ax_theta=None,
            ylim_V=[0,800],ylim_Pdynx=[10,25], ylim_IMF=[-20, 20],ylim_theta=[-180, 180],
            marker='.', linestyle='--', markersize=2, ylabel_fontsize=9, zero_line=True,
            padding_method=None, drop_na=False):

        # read ace data that is minute-averaged
        import sys
        sys.path.append('/home/muhammad/Documents/Important/RISRpy/RISR_06222015/codes')
        from read_sw_imf import ace_read

        # read data
        date_lag_dict = {"20020318" : 48,
                         "20030214" : 45,
                         "20040325" : 75,
                         "20060306" : 70}
        imf_delay = date_lag_dict[self.stime.strftime("%Y%m%d")]
        ylim_imf = [-40, 40]
        df_imf, df_sw =  ace_read(self.stime, self.etime, res=0, how='mean', coord='GSM', delay=imf_delay)

        df_imf.loc[:, 'Bt'] = np.sqrt((df_imf.By ** 2 + df_imf.Bz ** 2))
        # clock angle
        df_imf.loc[:, 'theta_Bt'] = np.degrees(np.arctan2(df_imf.By, df_imf.Bz))
        # dynamic pressure
        mp = 1.6726219 * 1e-27  # kg
        df_sw.loc[:, 'Pdynx'] = (df_sw.Np * mp * 1e6) * (df_sw.Vx * 1e3)**2 * 1e+9
        if drop_na:
            df_imf.dropna(how='any', inplace=True)
            df_sw.dropna(how='any', inplace=True)
        if padding_method is not None:
            df_imf.fillna(method=padding_method, inplace=True)
            df_sw.fillna(method=padding_method, inplace=True)

        # plot solar wind Vx
        if ax_Vx is not None:
            df_sw.loc[df_sw.Vx < -1e4, 'Vx'] = np.nan 
            ax_Vx.plot_date(df_sw.index.to_pydatetime(), -df_sw.Vx, color='k',
                    marker=marker, linestyle=linestyle, markersize=markersize)
            ax_Vx.set_ylabel('Vx [km/s]', fontsize=ylabel_fontsize)
            ax_Vx.set_ylim([ylim_V[0], ylim_V[1]])
            ax_Vx.locator_params(axis='y', nbins=4)
            ax_Vx.annotate('ACE', xy=(0.90, 0.8), xycoords='axes fraction',
                           ha='center', fontsize=9)
#            # add text to indicate sheath region
#            dt1 = dt.datself.etime(2014,9,12,15,26)
#            dt2 = dt.datself.etime(2014,9,12,16,10)
#            dt3 = dt.datself.etime(2014,9,12,16,57)
#            ax_Vx.annotate('', xy=(dt1, 400), xycoords='data', xytext=(dt3, 400),
#                    arrowprops={'arrowstyle': '<->'})
#            ax_Vx.annotate('Sheath', xy=(dt2, 400), xycoords='data', xytext=(0, 5), ha='center',
#                    textcoords='offset points')

        # plot Pdynx [nPa]
        if ax_Pdynx is not None:
            df_sw.loc[df_sw.Pdynx < 0, 'Pdynx'] = np.nan
            ax_Pdynx.plot_date(df_sw.index.to_pydatetime(), df_sw.Pdynx, color='k',
                    marker=marker, linestyle=linestyle, markersize=markersize)
            ax_Pdynx.set_ylabel('Pdynx [nPa]\n', fontsize=ylabel_fontsize)
            ax_Pdynx.set_ylim([ylim_Pdynx[0], ylim_Pdynx[1]])
            ax_Pdynx.locator_params(axis='y', nbins=4)
            ax_Pdynx.annotate('ACE', xy=(0.90, 0.3), xycoords='axes fraction',
                               ha='center', fontsize=9)

        # plot solar wind IMF
        if ax_IMF is not None:
            ax_IMF.plot_date(df_imf.index.to_pydatetime(), df_imf.Bx, color='k',
                    marker=marker, linestyle=linestyle, markersize=markersize)
            ax_IMF.plot_date(df_imf.index.to_pydatetime(), df_imf.By, color='g',
                    marker=marker, linestyle=linestyle, markersize=markersize)
            ax_IMF.plot_date(df_imf.index.to_pydatetime(), df_imf.Bz, color='r',
                    marker=marker, linestyle=linestyle, markersize=markersize)
            lns = ax_IMF.get_lines()
            ax_IMF.legend(lns,['Bx', 'By', 'Bz'], frameon=False, bbox_to_anchor=(0.98, 0.5),
                    loc='center left', fontsize='medium')
            #ax_IMF.legend(lns,['Bx', 'By', 'Bz'], frameon=False, fontsize='small', mode='expand')
            ax_IMF.set_ylabel('IMF_gsm\n' + "(nT)", fontsize=ylabel_fontsize)
            ax_IMF.set_ylim([ylim_IMF[0], ylim_IMF[1]])
            ax_IMF.locator_params(axis='y', nbins=4)
            if zero_line:
                ax_IMF.axhline(y=0, color='r', linewidth=0.20, linestyle='--')
            ax_IMF.annotate('ACE', xy=(0.90, 0.89), xycoords='axes fraction', ha='center', fontsize=9)

        # plot solar wind clock angle
        date_fname_dict = {"20020318" : "C3_CP_FGM_SPIN_31086.txt",
                           "20030214" : "C1_CP_FGM_SPIN_21206.txt",
                           "20040325" : "C3_CP_FGM_SPIN_3152599.txt",
                           "20060306" : "C3_CP_FGM_SPIN_31759.txt"}
        def my_func(row):
            clmn = "theta_Bt"
            if clmn in row.keys():
                row[clmn] = row[clmn] % 360
                #if row[clmn] > 270:
                #    row[clmn] = row[clmn] - 360
            return row
        if self.stime.strftime("%Y%m%d") == "20060306":
            df_imf.apply(my_func, axis=1)
        if ax_theta is not None:
            ax_theta.plot_date(df_imf.index.to_pydatetime(), df_imf.theta_Bt, color='k',
                    marker=marker, linestyle=linestyle, markersize=markersize)
            #ax_theta.set_ylabel(r'$\theta$ [deg]', fontsize=ylabel_fontsize)
            ax_theta.set_ylabel('IMF Clock Angle\n'+r"$(^\circ)$", fontsize=ylabel_fontsize)
            ax_theta.set_ylim([ylim_theta[0], ylim_theta[1]])
            ax_theta.locator_params(axis='y', nbins=4)

            # add extra ytick labels
            #poss1 = np.append(ax_theta.yaxis.get_majorticklocs(), ax_theta.yaxis.get_minorticklocs())
            poss1 = ax_theta.yaxis.get_majorticklocs()
            poss1 = np.append(poss1, [-180, -90, 90, 180])
            poss1 = np.delete(poss1, [0, 1, 3, 4])
            ax_theta.yaxis.set_ticks(poss1)

            # draw horizontal lines
            if zero_line:
                ax_theta.axhline(y=0, color='r', linewidth=0.10, linestyle='--')
            ninety_line = False 
            if ninety_line:
                ax_theta.axhline(y=90, color='r', linewidth=0.10, linestyle='--')
                ax_theta.axhline(y=-90, color='r', linewidth=0.10, linestyle='--')

        return

    def plot_cluster(self, ax_c_B=None,ax_c_Btheta=None,
                     ax_c_Pdynx=None, time_lag=0,
                     ylim_c_B=[-60, 60],ylim_c_theta=[-180, 180],
                     ylim_c_Pdynx=[0, 20],
                     marker='.', linestyle='--', markersize=2,
                     ylabel_fontsize=9, zero_line=True,
                     interpolation_method=None, drop_na=False):

        # read data
        date_parser = lambda x: pd.datetime.strptime(x, '%d-%m-%Y %H:%M:%S.%f')
        #date_parser = lambda x, y: pd.datself.etime.strptime(" ".join([str(x), str(y)]), '%d-%m-%Y %H:%M:%S.%f')

        # reading the data into pandas' DataFrame object
        date_fname_dict = {"20020318" : "C3_CP_FGM_SPIN_31086.txt",
                           "20030214" : "C1_CP_FGM_SPIN_21206.txt",
                           "20040325" : "C3_CP_FGM_SPIN_3152599.txt",
                           "20060306" : "C3_CP_FGM_SPIN_31759.txt"}
        fname = "../data/cluster/" + date_fname_dict[self.stime.strftime("%Y%m%d")]
        date_srows_dict = {"20020318" : 114,
                           "20030214" : 127,
                           "20040325" : 114,
                           "20060306" : 114}

        srows = range(date_srows_dict[self.stime.strftime("%Y%m%d")])
        df_c1 = pd.read_csv(fname, index_col=0, skipinitialspace=True,
                delim_whitespace=True, skiprows=srows,
                header=None,
                parse_dates={'datself.etime': [0,1]},
                date_parser=date_parser,
                error_bad_lines=False)

        df_c1.columns = ["B", "BX_GSE", "BY_GSE", "BZ_GSE", "X_GSE", "Y_GSE", "Z_GSE"]

        # select the period of interest
        df_c1 = df_c1.loc[self.stime:self.etime,]

        # convert Bfield data from GSE to GSM
        import spacepy.coordinates as spc
        import spacepy.time as spt
        import numpy as np
        GSE_array = np.array([df_c1.BX_GSE.tolist(), df_c1.BY_GSE.tolist(),
                                  df_c1.BZ_GSE.tolist()])
        GSE_array = GSE_array.transpose()
        GSE_position = spc.Coords(GSE_array, "GSE", "car")
        GSE_position.ticks = spt.Ticktock(df_c1.index.to_pydatetime(), "UTC")
        GSM_position = GSE_position.convert("GSM", "car")
        df_c1.loc[:, "BX_GSM"] = GSM_position.x
        df_c1.loc[:, "BY_GSM"] = GSM_position.y
        df_c1.loc[:, "BZ_GSM"] = GSM_position.z

        df_c1.loc[:, 'Bt'] = np.sqrt((df_c1.BY_GSM ** 2 + df_c1.BZ_GSM ** 2))
        # clock angle
        df_c1.loc[:, 'theta_Bt'] = np.degrees(np.arctan2(df_c1.BY_GSM, df_c1.BZ_GSM))
        def my_func(row):
            clmn = "theta_Bt"
            if clmn in row.keys():
                row[clmn] = row[clmn] % 360
                #if row[clmn] > 270:
                #    row[clmn] = row[clmn] - 360
            return row
        if self.stime.strftime("%Y%m%d") == "20060306":
            df_c1.apply(my_func, axis=1)

        if drop_na:
            df_c1.dropna(how='any', inplace=True)
        if interpolation_method is not None:
            df_c1.fillna(method=interpolation_method, inplace=True)

        # plot cluster Bfield data 
        if ax_c_B is not None:
            ax_c_B.plot_date(df_c1.index.to_pydatetime(), df_c1.BX_GSM, color='k',
                    marker=marker, linestyle=linestyle, markersize=markersize)
            ax_c_B.plot_date(df_c1.index.to_pydatetime(), df_c1.BY_GSM, color='g',
                    marker=marker, linestyle=linestyle, markersize=markersize)
            ax_c_B.plot_date(df_c1.index.to_pydatetime(), df_c1.BZ_GSM, color='r',
                    marker=marker, linestyle=linestyle, markersize=markersize)
            lns = ax_c_B.get_lines()
            ax_c_B.legend(lns,['Bx', 'By', 'Bz'], frameon=False, bbox_to_anchor=(0.98, 0.5),
                    loc='center left', fontsize='medium')
            #ax_c_B.legend(lns,['Bx', 'By', 'Bz'], frameon=False, fontsize='small', mode='expand')
            ax_c_B.set_ylabel('B_gsm\n'+"(nT)", fontsize=ylabel_fontsize)
            ax_c_B.set_ylim([ylim_c_B[0], ylim_c_B[1]])
            ax_c_B.locator_params(axis='y', nbins=4)
            if zero_line:
                ax_c_B.axhline(y=0, color='r', linewidth=0.30, linestyle='--')
            ax_c_B.annotate("Cluster C3", xy=(0.90, 0.89),
                            xycoords='axes fraction', ha='center', fontsize=9)


        # plot cluster Bfield clock angle
        if ax_c_Btheta is not None:
            ax_c_Btheta.plot_date(df_c1.index.to_pydatetime(), df_c1.theta_Bt, color='r',
                    marker=marker, linestyle=linestyle, markersize=markersize)
            #ax_c_Btheta.set_ylabel(r'$\theta$ [deg]', fontsize=ylabel_fontsize)
            ax_c_Btheta.set_ylabel('IMF Clock Angle\n'+r"$(^\circ)$", fontsize=ylabel_fontsize)
            ax_c_Btheta.set_ylim([ylim_c_theta[0], ylim_c_theta[1]])
            ax_c_Btheta.locator_params(axis='y', nbins=4)
            # add extra ytick labels
            #poss1 = np.append(ax_c_Btheta.yaxis.get_majorticklocs(), ax_c_Btheta.yaxis.get_minorticklocs())
            poss1 = ax_c_Btheta.yaxis.get_majorticklocs()
            poss1 = np.append(poss1, [-180, -90, 90, 180])
            #poss1 = np.delete(poss1, [0, 1, 3, 4])
            ax_c_Btheta.yaxis.set_ticks(poss1)
            # draw horizontal lines
            if zero_line:
                ax_c_Btheta.axhline(y=0, color='r', linewidth=0.30, linestyle='--')
            ninety_line = True 
            if ninety_line:
                ax_c_Btheta.axhline(y=90, color='r', linewidth=0.30, linestyle='--')
                ax_c_Btheta.axhline(y=-90, color='r', linewidth=0.30, linestyle='--')

            # Annotate some text
            xy_coord_ace = (0.90, 0.25)
            xy_coord_cname = (0.90, 0.15)
            xy_coord_lagtxt = (0.90, 0.03)
            ax_c_Btheta.annotate("ACE", xy=xy_coord_ace,
                    xycoords='axes fraction', ha='center', fontsize=9)
            ax_c_Btheta.annotate("Cluster C3", xy=xy_coord_cname,
                    xycoords='axes fraction', ha='center', fontsize=9, color="r")
            ax_c_Btheta.annotate("Tp=48min", xy=xy_coord_lagtxt,
                    xycoords='axes fraction', ha='center', fontsize=9)

        # Plot cluster Pdyn
        if ax_c_Pdynx is not None:
            # read cluster CIS data
            import sys
            sys.path.append('../')
            from read_cluster_data import read_cis_data
            
            fnm = "../data//cluster/C3_PP_CIS_28071.txt" 
            # read the ace data again with proper lag time
            delay_in_response = time_lag
            df_sw = read_cis_data(self.stime, self.etime, fnm, start_row=77,
                                   res=0, delay=delay_in_response,
                                   coords="GSM")
            # dynamic pressure
            mp = 1.6726219 * 1e-27  # kg
            df_sw.loc[:, 'Pdynx'] = (df_sw.N_HIA * mp * 1e6) * (df_sw.VX_HIA_GSM * 1e3)**2 * 1e+9
            #df_sw.loc[:, 'Pdyn'] = (df_sw.N_HIA * mp * 1e6) * (df_sw.V_HIA * 1e3)**2 * 1e+9
            df_sw.loc[df_sw.Pdynx < 0, 'Pdynx'] = np.nan
            #df_sw.dropna(how='any', inplace=True)
    
            df_sw.loc[df_sw.N_HIA < 0, 'N_HIA'] = np.nan
    
    #        # clock angle
    #        df_sw.loc[:, 'theta_Bt'] = np.degrees(np.arctan2(df_sw.By, df_sw.Bz))
    
            # plot IMF data
            ax_c_Pdynx.plot_date(df_sw.index.to_pydatetime(), df_sw.Pdynx, color='k',
                                 marker=marker, linestyle=linestyle, markersize=markersize)

            ax_c_Pdynx.set_ylim([ylim_c_Pdynx[0], ylim_c_Pdynx[1]])
            ax_c_Pdynx.annotate("Cluster C3", xy=(0.90, 0.89),
                                xycoords='axes fraction', ha='center', fontsize=9, color="k")
            ax_c_Pdynx.set_ylabel('Pdynx [nPa]\n', fontsize=ylabel_fontsize)    


    def plot_symh(self, ax_symh, ylim_symh=[-100, 100], ylabel_fontsize=9,
            marker='.', linestyle='--', markersize=2, zero_line=True):

        import symasy

        # read SYMH data
        sym_list = symasy.readSymAsyWeb(sTime=self.stime, eTime=self.etime)
        #sym_list = gme.ind.symasy.readSymAsy(sTime=self.stime, eTime=self.etime)
        symh = []
        symh_time = []
        for k in range(len(sym_list)):
            symh.append(sym_list[k].symh)
            symh_time.append(sym_list[k].time)
        
        # plot symh
        indx = [symh_time.index(x) for x in symh_time if (x>= self.stime) and (x<=self.etime)]
        symh_time = [symh_time[i] for i in indx]
        symh = [symh[i] for i in indx]

        ax_symh.plot_date(symh_time, symh, 'k', marker=marker,
                          linestyle=linestyle, markersize=markersize)
        ax_symh.set_ylabel('SYM-H', fontsize=9)
        ax_symh.set_ylim([ylim_symh[0], ylim_symh[1]])
        ax_symh.locator_params(axis='y', nbins=4)
        if zero_line:
            ax_symh.axhline(y=0, color='r', linewidth=0.30, linestyle='--')


    def plot_ae(self, ax_ae, ylim_ae=[0, 500], ylabel_fontsize=9,
            marker='.', linestyle='--', markersize=2):

        import ae

        # read AE data
        #AE_list = gme.ind.readAeWeb(sTime=self.stime, eTime=self.etime,res=1)
        AE_list = ae.readAeWeb(sTime=self.stime, eTime=self.etime, res=1)
        AE = []
        AE_time = []
        for m in range(len(AE_list)):
            AE.append(AE_list[m].ae)
            AE_time.append(AE_list[m].time)

        # plot AE 
        indx = [AE_time.index(x) for x in AE_time if (x>= self.stime) and (x<=self.etime)]
        AE_time = [AE_time[i] for i in indx]
        AE = [AE[i] for i in indx]
        ax_ae.plot_date(AE_time, AE, 'k', marker=marker,
                        linestyle=linestyle, markersize=markersize)
        ax_ae.set_ylabel('AE', fontsize=9)
        ax_ae.set_ylim([ylim_ae[0], ylim_ae[1]])
        ax_ae.locator_params(axis='y', nbins=4)

    def plot_kp(self, ax_kp, ylim_kp=[0, 9], ylabel_fontsize=9,
            marker='.', linestyle='--', markersize=2):

        from davitpy import gme
        # Kp data
        Kp_list = gme.ind.readKpFtp(stime=self.stime, etime=self.etime)
        Kp = []
        Kp_time = []
        for n in range((self.etime-self.stime).days):
            kp_tmp = Kp_list[n].kp
            time_tmp = Kp_list[n].time
            for l in range(len(kp_tmp)):
                if len(kp_tmp[l])== 2:
                    if kp_tmp[l][1] == '+':
                        Kp.append(int(kp_tmp[l][0])+0.3)
                    elif kp_tmp[l][1] == '-':
                        Kp.append(int(kp_tmp[l][0])-0.3)
                else:
                    Kp.append(int(kp_tmp[l][0]))
                Kp_time.append(time_tmp + dt.timedelta(hours=3*l))

        # plot Kp 
        ax_kp.plot_date(Kp_time, Kp, 'k', marker=marker, linestyle=linestyle, markersize=markersize)
        ax_kp.set_ylabel('Kp', fontsize=9)
        ax_kp.set_ylim([ylim_kp[0], ylim_kp[1]])
        ax_kp.locator_params(axis='y', nbins=4)

    def plot_wic_si_max_point(self, ax_mlt=None, ax_mlat=None, 
                              ax_intensity=None,
                              delay_in_response=10,
                              ace_delay_at_cluster=48,
                              read_si=True, read_wic=False,
                              mlt_range=[9, 15], mlat_range=[70, 89],
                              plot_both_wic_si=True, overlay_ace=False,
                              overlay_cluster=False,
                              xcorr_pair={"mlt":"By", "mlat":"Bz",
                                          "intensity":"Pdynx"},
                              spot_patch_size=(0,0),
                              marker='.', linestyle='-', markersize=2):
        """ Plots the max intensity point location of the cusp spot with IMF parameters"""
        import sys
        sys.path.append("../track_wic_si_spot")
        from track_max_point import find_filenames, read_wic_si_data
        from track_max_point import find_max_intensity_point
        from matplotlib.dates import DateFormatter
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        from lag_time_finder import find_best_lag_time
    
#        df_tmp = []
#        color_list = ["r", "k"]
#        if plot_both_wic_si:
#            niter = 2
#            si_bool = [True, False]
#            wic_bool = [False, True]
#            legend_text = ["SI", "WIC"]
#        else:
#            niter = 1
#        for i in range(niter):
#            if niter == 1:
#                read_si = read_si; read_wic = read_wic
#                if read_si:
#                    label_txt = "SI"
#                if read_wic:
#                    label_txt = "WIC"
#            else:
#                read_si = si_bool[i]; read_wic = wic_bool[i]
#                label_txt = legend_text[i]
#            dtms, fnames = find_filenames(self.stime, self.etime, read_si=read_si, read_wic=read_wic)
#            mlats = []
#            mlts = []
#            intensities = []
#            for dtm in dtms:
#                fname, data = read_wic_si_data(dtm, read_si=read_si, read_wic=read_wic)
#                (maxval, maxloc) = find_max_intensity_point(data, mlt_range=mlt_range,
#                                                            mlat_range=mlat_range)
#                mlats.append(data['mlat'][maxloc[1], maxloc[0]])
#                mlts.append(data['mlt'][maxloc[1], maxloc[0]])
#                # image intensity values in a given patch
#                patch_intensity = data['image'][maxloc[1]-spot_patch_size[1]:maxloc[1]+spot_patch_size[1]+1,
#                                                maxloc[0]-spot_patch_size[0]:maxloc[0]+spot_patch_size[0]+1]
#                intensities.append(np.mean(patch_intensity))
#
#            # construct a df 
#            df_var = pd.DataFrame(data={"mlat":mlats, "mlt":mlts, 
#                                        "intensity":intensities}, index=dtms)
#            df_tmp.append(df_var)


        df_tmp = []
        color_list = ["r", "k"]
        if plot_both_wic_si:
            niter = 2
            si_bool = [True, False]
            wic_bool = [False, True]
            legend_text = ["SI", "WIC"]
        else:
            niter = 1
        for i in range(niter):
            if niter == 1:
                if read_si:
                    label_txt = "SI"
                    clr = color_list[0]
                if read_wic:
                    label_txt = "WIC"
                    clr = color_list[1]
            else:
                read_si = si_bool[i]; read_wic = wic_bool[i]
                label_txt = legend_text[i]
                clr = color_list[i]
            fltmp = "./" + label_txt + "_spot_max_point_" +\
                    str(2*spot_patch_size[0]+1) + "_by_" + str(2*spot_patch_size[1]+1) + ".csv" 
            df_var = pd.read_csv(fltmp, parse_dates=True, index_col=0)
            # Select data for the time interval of interest
            df_var = df_var.loc[self.stime:self.etime, :]
            mlats = df_var.mlat.as_matrix()
            mlts = df_var.mlt.as_matrix()
            dtms = pd.to_datetime(df_var.index)
            intensities = df_var.intensity.as_matrix()
            df_tmp.append(df_var)

            if ax_mlt:
                ax_mlt.plot_date(dtms, mlts, '.-', color=clr, label=label_txt)
                ax_mlt.xaxis.set_major_formatter(DateFormatter('%H:%M'))
                ax_mlt.set_ylabel("MLT")
                ax_mlt.set_ylim([9, 15])
                ax_mlt.legend(loc="best", frameon=False, fontsize='medium')

            if ax_mlat:
                ax_mlat.plot_date(dtms, mlats, '.-', color=clr, label=label_txt)
                ax_mlat.xaxis.set_major_formatter(DateFormatter('%H:%M'))
                ax_mlat.set_ylabel("MLAT")
                ax_mlat.set_ylim([75, 85])
                ax_mlat.legend(loc="best", frameon=False, fontsize='medium')

            if ax_intensity:
                ax_intensity.plot_date(dtms, intensities, '.-', color=clr, label=label_txt)
                ax_intensity.xaxis.set_major_formatter(DateFormatter('%H:%M'))
                ax_intensity.set_ylabel("Intensity")
                ax_intensity.set_ylim([0, 5000])
                ax_intensity.legend(loc="best", frameon=False, fontsize='medium')
            
        if plot_both_wic_si:
            df_wic_si = df_tmp[0].join(df_tmp[1], rsuffix="_wic")
            # Annotate the correlation
            if ax_mlat:
                # calculate the best xcor #df_tmp[0]["mlat_wic"] = df_tmp[1]["mlat"]
                xcor_max = df_wic_si[["mlat_wic", "mlat"]].corr().iloc[0,1]
                # annotate the correlation coefficient
                ax_mlat.annotate('CorrCoef = ' + '{0:.2f}'.format(xcor_max),\
                                xy=(0.5, 0.05), xycoords='axes fraction', ha="center",
                                fontsize='medium', color='r')
            if ax_mlt:
                # calculate the best xcor
                #df_tmp[0]["mlt_wic"] = df_tmp[1]["mlt"]
                xcor_max = df_wic_si[["mlt_wic", "mlt"]].corr().iloc[0,1]
                # annotate the correlation coefficient
                ax_mlt.annotate('CorrCoef = ' + '{0:.2f}'.format(xcor_max),\
                            xy=(0.5, 0.05), xycoords='axes fraction', ha="center",
                            fontsize='medium', color='r')
            if ax_intensity:
                # calculate the best xcor
                #df_tmp[0]["intensity_wic"] = df_tmp[1]["intensity"]
                xcor_max = df_wic_si[["intensity_wic", "intensity"]].corr().iloc[0,1]
                # annotate the correlation coefficient
                ax_intensity.annotate('CorrCoef = ' + '{0:.2f}'.format(xcor_max),\
                            xy=(0.5, 0.05), xycoords='axes fraction', ha="center",
                            fontsize='medium', color='r')
        
        if overlay_ace:
            # read ace data that is minute-averaged
            import sys
            sys.path.append('/home/muhammad/Documents/Important/RISRpy/RISR_06222015/codes')
            from read_sw_imf import ace_read
    
            # read the ace data again with proper lag time
            df_imf, df_sw =  ace_read(self.stime, self.etime, res=0, how='mean', coord='GSM',
                                      delay=ace_delay_at_cluster + delay_in_response)
            # clock angle
            df_imf.loc[:, 'theta_Bt'] = np.degrees(np.arctan2(df_imf.By, df_imf.Bz))
    
            # dynamic pressure
            mp = 1.6726219 * 1e-27  # kg
            df_sw.loc[:, 'Pdynx'] = (df_sw.Np * mp * 1e6) * (df_sw.Vx * 1e3)**2 * 1e+9
            df_sw.loc[df_sw.Pdynx < 0, 'Pdynx'] = np.nan
            df_sw.dropna(how='any', inplace=True)
    
            # calculate the best xcor
            def plot_xcor_pair(ax, xcorr_pair, spot_param, df_sw, df_imf, df_var,
                               marker='.', linestyle='-', markersize=2):
                IMF_param=xcorr_pair[spot_param]
                if IMF_param=="Pdynx":
                    dfn_tmp = df_var.join(df_sw.resample('1Min').mean(), how='left').dropna()
                    dfn = df_sw
                else:
                    dfn_tmp = df_var.join(df_imf.resample('1Min').mean(), how='left').dropna()
                    dfn = df_imf

                # calculate the best xcor
                xcor_max = dfn_tmp[[spot_param, IMF_param]].corr().iloc[0,1]

                # Plot the IMF & SW params
                if IMF_param=="Bx":
                    ax_IMF.plot_date(dfn.index.to_pydatetime(), dfn[IMF_param].as_matrix(), color='k',
                            marker=marker, linestyle=linestyle, markersize=markersize)
                if IMF_param=="By":
                    ax_IMF.plot_date(dfn.index.to_pydatetime(), dfn[IMF_param].as_matrix(), color='g',
                            marker=marker, linestyle=linestyle, markersize=markersize, label="By")
                if IMF_param=="Bz":
                    ax_IMF.plot_date(dfn.index.to_pydatetime(), dfn[IMF_param].as_matrix(), color='b',
                            marker=marker, linestyle=linestyle, markersize=markersize, label="Bz")
                if IMF_param=="theta_Bt":
                    ax_IMF.plot_date(dfn.index.to_pydatetime(), dfn[IMF_param].as_matrix(), color='k',
                            marker=marker, linestyle=linestyle, markersize=markersize)
                if IMF_param=="Pdynx":
                    ax_IMF.plot_date(dfn.index.to_pydatetime(), dfn[IMF_param].as_matrix(), color='b',
                            marker=marker, linestyle=linestyle, markersize=markersize)
                lns = ax_IMF.get_lines()
                #ax_IMF.legend(lns,['Bx', 'By', 'Bz'], frameon=False, fontsize='small', mode='expand')
                if IMF_param=="theta_Bt":
                    ax_IMF.set_ylim([-150, 150])
                    ax_IMF.set_ylabel('IMF Theta ' + "[Degree]")
                    ax_IMF.legend(lns,[IMF_param], frameon=False, 
                                  loc='best', fontsize='medium')
                elif IMF_param=="Pdynx":
                    ax_IMF.set_ylim([10, 25])
                    ax_IMF.set_ylabel('Pdyn ' + "[nPa]")
                    ax_IMF.legend(lns,[IMF_param], frameon=False,
                                  loc='best', fontsize='medium')
                else:
                    ax_IMF.set_ylim([-25, 25])
                    ax_IMF.set_ylabel('IMF ' + "[nT]")
                    ax_IMF.locator_params(axis='y', nbins=4)
                    #ax_IMF.legend(lns,[IMF_param], frameon=False, bbox_to_anchor=(0.98, 0.7),
                    #        loc='center left', fontsize='medium')
                    ax_IMF.legend(frameon=False, mode='expand',
                                  loc='best', fontsize='medium')
        
                ax_IMF.axhline(y=0, color='k', linewidth=0.40, linestyle='--')

                return xcor_max
    
            # plot IMF data

            if ax_mlt:
                ax_IMF = ax_mlt.twinx()
                xcor_max = plot_xcor_pair(ax_IMF, xcorr_pair, "mlt", df_sw, df_imf, df_var,
                                          marker=marker, linestyle=linestyle, markersize=markersize)
                # annotate the correlation coefficient
                ax_mlt.annotate('CorrCoef = ' + '{0:.2f}'.format(xcor_max),\
                            xy=(0.5, 0.17), xycoords='axes fraction', ha="center",
                            fontsize='medium', color='r')
                ax_mlt.annotate('Lag Time = ' + '{0:.0f}'.format(ace_delay_at_cluster + delay_in_response) + ' min',\
                            xy=(0.5, 0.03), xycoords='axes fraction', ha="center",
                            fontsize='medium', color='r')

            if ax_mlat:
                ax_IMF = ax_mlat.twinx()
                xcor_max = plot_xcor_pair(ax_IMF, xcorr_pair, "mlat", df_sw, df_imf, df_var,
                                          marker=marker, linestyle=linestyle, markersize=markersize)
                # annotate the correlation coefficient
                ax_mlat.annotate('CorrCoef = ' + '{0:.2f}'.format(xcor_max),\
                            xy=(0.5, 0.17), xycoords='axes fraction', ha="center",
                            fontsize='medium', color='r')
                ax_mlat.annotate('Lag Time = ' + '{0:.0f}'.format(ace_delay_at_cluster + delay_in_response) + ' min',\
                            xy=(0.5, 0.03), xycoords='axes fraction', ha="center",
                            fontsize='medium', color='r')

            if ax_intensity:
                ax_IMF = ax_intensity.twinx()
                xcor_max = plot_xcor_pair(ax_IMF, xcorr_pair, "intensity", df_sw, df_imf, df_var,
                                          marker=marker, linestyle=linestyle, markersize=markersize)
                # annotate the correlation coefficient
                ax_intensity.annotate('CorrCoef = ' + '{0:.2f}'.format(xcor_max),\
                            xy=(0.5, 0.17), xycoords='axes fraction', ha="center",
                            fontsize='medium', color='r')
                ax_intensity.annotate('Lag Time = ' + '{0:.0f}'.format(ace_delay_at_cluster + delay_in_response) + ' min',\
                            xy=(0.5, 0.03), xycoords='axes fraction', ha="center",
                            fontsize='medium', color='r')

    
        if overlay_cluster:
            # read cluster data
            import sys
            sys.path.append('../')
            from read_cluster_data import read_cis_data

            fnm = "../data//cluster/C3_PP_CIS_28071.txt"
            # read the ace data again with proper lag time
            df_sw = read_cis_data(self.stime, self.etime, fnm, start_row=77,
                                   res=0, delay=delay_in_response,
                                   coords="GSM")
            # dynamic pressure
            mp = 1.6726219 * 1e-27  # kg
            df_sw.loc[:, 'Pdynx'] = (df_sw.N_HIA * mp * 1e6) * (df_sw.VX_HIA_GSM * 1e3)**2 * 1e+9
            #df_sw.loc[:, 'Pdyn'] = (df_sw.N_HIA * mp * 1e6) * (df_sw.V_HIA * 1e3)**2 * 1e+9
            df_sw.loc[df_sw.Pdynx < 0, 'Pdynx'] = np.nan
            #df_sw.dropna(how='any', inplace=True)
   
            df_sw.loc[df_sw.N_HIA < 0, 'N_HIA'] = np.nan

    #        # clock angle
    #        df_sw.loc[:, 'theta_Bt'] = np.degrees(np.arctan2(df_sw.By, df_sw.Bz))
    
            # Create a DF to calculate correlation
            dfn = df_var.join(df_sw.resample('1Min').mean(), how='left').dropna()
    
            # plot IMF data
            if ax_mlt:
                ax_IMF = ax_mlt.twinx()
                IMF_param=xcorr_pair["mlt"]
                ax_IMF.plot_date(df_sw.index.to_pydatetime(), df_sw[IMF_param].as_matrix(), color='k',
                        marker=marker, linestyle=linestyle, markersize=markersize, label=IMF_param)
                # calculate the best xcor
                xcor_max = dfn[["mlt", IMF_param]].corr().iloc[0,1]
                # annotate the correlation coefficient
                ax_mlt.annotate('CorrCoef = ' + '{0:.2f}'.format(xcor_max),\
                            xy=(0.5, 0.17), xycoords='axes fraction', ha="center",
                            fontsize='medium', color='r')
                ax_mlt.annotate('Lag Time = ' + '{0:.0f}'.format(delay_in_response) + ' min',\
                            xy=(0.5, 0.14), xycoords='axes fraction', ha="center",
                            fontsize='medium', color='r')

            if ax_mlat:
                ax_IMF = ax_mlat.twinx()
                IMF_param=xcorr_pair["mlat"]
                ax_IMF.plot_date(df_sw.index.to_pydatetime(), df_sw[IMF_param].as_matrix(), color='k',
                        marker=marker, linestyle=linestyle, markersize=markersize, label=IMF_param)
                # calculate the best xcor
                xcor_max = dfn[["mlat", IMF_param]].corr().iloc[0,1]
                # annotate the correlation coefficient
                ax_mlat.annotate('CorrCoef = ' + '{0:.2f}'.format(xcor_max),\
                            xy=(0.5, 0.17), xycoords='axes fraction', ha="center",
                            fontsize='medium', color='r')
                ax_mlat.annotate('Lag Time = ' + '{0:.0f}'.format(delay_in_response) + ' min',\
                            xy=(0.5, 0.14), xycoords='axes fraction', ha="center",
                            fontsize='medium', color='r')

            if ax_intensity:
                ax_IMF = ax_intensity.twinx()
                IMF_param=xcorr_pair["intensity"]
                ax_IMF.plot_date(df_sw.index.to_pydatetime(), df_sw[IMF_param].as_matrix(), color='b',
                        marker=marker, linestyle=linestyle, markersize=markersize, label=IMF_param)

                if IMF_param=="Pdynx":
                    ax_IMF.set_ylim([0, 15])
                    ax_IMF.set_ylabel(IMF_param + " [nPa]")
                elif IMF_param=="N_HIA":
                    ax_IMF.set_ylim([0, 250])
                    ax_IMF.set_ylabel(IMF_param + " [/cm3]")

                # calculate the best xcor
                xcor_max = dfn[["intensity", IMF_param]].corr().iloc[0,1]
                # annotate the correlation coefficient
                ax_intensity.annotate('CorrCoef = ' + '{0:.2f}'.format(xcor_max),\
                            xy=(0.5, 0.17), xycoords='axes fraction', ha="center",
                            fontsize='medium', color='r')
                ax_intensity.annotate('Lag Time = ' + '{0:.0f}'.format(delay_in_response) + ' min',\
                            xy=(0.5, 0.03), xycoords='axes fraction', ha="center",
                            fontsize='medium', color='r')


#            if IMF_param=="theta_Bt":
#                ax_IMF.plot_date(df_sw.index.to_pydatetime(), df_sw.theta_Bt, color='k',
#                        marker=marker, linestyle=linestyle, markersize=markersize)
#            if IMF_param=="Pdyn":
#                ax_IMF.plot_date(df_sw.index.to_pydatetime(), df_sw.Pdyn, color='k',
#                        marker=marker, linestyle=linestyle, markersize=markersize)
#            if IMF_param=="N_HIA":
#                ax_IMF.plot_date(df_sw.index.to_pydatetime(), df_sw.N_HIA, color='k',
#                        marker=marker, linestyle=linestyle, markersize=markersize)
#    
#            lns = ax_IMF.get_lines()
#            #ax_IMF.legend(lns,['Bx', 'By', 'Bz'], frameon=False, fontsize='small', mode='expand')
#            if IMF_param=="theta_Bt":
#                ax_IMF.set_ylim([-150, 150])
#                ax_IMF.set_ylabel('IMF Theta ' + "[Degree]")
#                ax_IMF.legend(lns,[IMF_param], frameon=False, bbox_to_anchor=(0.80, 0.2),
#                        loc='center left', fontsize='medium')
#            elif IMF_param=="Pdyn":
#                ax_IMF.set_ylim([0, 15])
#                ax_IMF.set_ylabel('Pdyn ' + "[nPa]")
#                ax_IMF.legend(lns,[IMF_param], frameon=False, bbox_to_anchor=(0.80, 0.2),
#                        loc='center left', fontsize='medium')
#            elif IMF_param=="N_HIA":
#                ax_IMF.set_ylim([0, 250])
#                ax_IMF.set_ylabel('N_HIA ' + "[/cm3]")
#                ax_IMF.legend(lns,[IMF_param], frameon=False, bbox_to_anchor=(0.80, 0.2),
#                        loc='center left', fontsize='medium')
#    
#            else:
#                ax_IMF.set_ylim([-25, 25])
#                ax_IMF.set_ylabel('IMF ' + "[nT]")
#                ax_IMF.locator_params(axis='y', nbins=4)
#                ax_IMF.legend(lns,[IMF_param], frameon=False, bbox_to_anchor=(0.98, 0.7),
#                        loc='center left', fontsize='medium')
    
            ax_IMF.axhline(y=0, color='k', linewidth=0.40, linestyle='--')
    
        return


# testing
if __name__ == "__main__":
    import datetime as dt

    stime = dt.datetime(2002,3,18,15,0) 
    etime = dt.datetime(2002,3,18,17,0) 
    
    IMF_geomag_indices = False
    WIC_SI_spot = True

    ###########################################################################
    if IMF_geomag_indices:
        # Plot IMF and geomagnetic indices
        #figsize = None
        figsize=(8,12)
        nrows = 8
        ncols = 1
        hspace=None
        obj = multipanel_plot(stime, etime, figsize=figsize, nrows=nrows, ncols=ncols,
                         hspace=hspace, sharex=True)
        ax_Vx = obj.axes[0]
        ax_Pdynx = obj.axes[1]
        ax_c_Pdynx = obj.axes[2]
        ax_IMF = obj.axes[3]
        ax_c_B = obj.axes[4]
        ax_theta = obj.axes[5]
        ax_c_Btheta = obj.axes[5]
        ax_ae = obj.axes[6]
        ax_symh = obj.axes[7]

        # Plot ACE IMF
        obj.plot_ace(ax_Vx=ax_Vx, ax_Pdynx=ax_Pdynx, ax_IMF=ax_IMF, ax_theta=ax_theta,
                     ylim_V=[300,700],ylim_Pdynx=[0,25], ylim_IMF=[-25, 25],ylim_theta=[-180, 180],
                     marker='.', linestyle='--', markersize=2, ylabel_fontsize=9, zero_line=True,
                     padding_method=None, drop_na=False)

        # Plot Cluster IMF
        obj.plot_cluster(ax_c_B=ax_c_B, ax_c_Btheta=ax_c_Btheta,
                         ax_c_Pdynx=ax_c_Pdynx, time_lag=0,
                         ylim_c_Pdynx=[0, 15],
                         ylim_c_B=[-100, 100], ylim_c_theta=[-180, 180],
                         marker='.', linestyle='--', markersize=2,
                         ylabel_fontsize=9, zero_line=True,
                         interpolation_method=None, drop_na=False)

        # Plot AE
        obj.plot_ae(ax_ae, ylim_ae=[0, 300], ylabel_fontsize=9,
                    marker='.', linestyle='--', markersize=2)

        # Plot SymH
        obj.plot_symh(ax_symh, ylim_symh=[-100, 100], ylabel_fontsize=9,
                      marker='.', linestyle='--', markersize=2, zero_line=True)

        fig_name =  stime.strftime("%Y%m%d.%H%M.")+etime.strftime("%Y%m%d.%H%M.")+\
                    'imf_geoindices_info.png' 

#        # Save the figure
#        dpi = 300
#        obj.fig.savefig('../plots/multipanel_time_series/' + fig_name ,dpi=dpi)
#        plt.close(obj.fig)
#        #plt.show()


    ###########################################################################
    # Plot WIC/SI related parameters
    if WIC_SI_spot:
        #figsize = None
        figsize=(6,12)
        nrows = 7
        ncols = 1
        hspace=None
        obj = multipanel_plot(stime, etime, figsize=figsize, nrows=nrows, ncols=ncols,
                         hspace=hspace, sharex=True)
        # Plot WIC/SI spot
        ax_mlat= obj.axes[0]
        ax_mlt= obj.axes[1]
        ax_intensity = obj.axes[2]
        xcorr_pair=None
        obj.plot_wic_si_max_point(ax_mlt=ax_mlt, ax_mlat=ax_mlat, 
                                  ax_intensity=ax_intensity,
                                  delay_in_response=10,
                                  ace_delay_at_cluster=48,
                                  read_si=True, read_wic=False,
                                  mlt_range=[9, 15], mlat_range=[70, 89],
                                  plot_both_wic_si=True, overlay_ace=False,
                                  overlay_cluster=False,
                                  xcorr_pair=xcorr_pair,
                                  spot_patch_size=(0,0))

        # Plot WIC/SI spot together with ACE params
        ax_mlt= obj.axes[3]
        ax_mlat= obj.axes[4]
        ax_intensity = None 
        xcorr_pair={"mlt":"By", "mlat":"Bz", "intensity":"Pdynx"}
        obj.plot_wic_si_max_point(ax_mlt=ax_mlt, ax_mlat=ax_mlat, 
                                  ax_intensity=ax_intensity,
                                  delay_in_response=10,
                                  ace_delay_at_cluster=48,
                                  read_si=False, read_wic=True,
                                  mlt_range=[9, 15], mlat_range=[70, 89],
                                  plot_both_wic_si=False, overlay_ace=True,
                                  overlay_cluster=False,
                                  xcorr_pair=xcorr_pair,
                                  spot_patch_size=(0,0))

        # Plot WIC/SI spot together with Cluster params
        ax_mlt=None 
        ax_mlat=None
        ax_intensity = obj.axes[5]
        xcorr_pair={"mlt":"By", "mlat":"Bz", "intensity":"Pdynx"}
        obj.plot_wic_si_max_point(ax_mlt=ax_mlt, ax_mlat=ax_mlat, 
                                  ax_intensity=ax_intensity,
                                  delay_in_response=10,
                                  ace_delay_at_cluster=48,
                                  read_si=False, read_wic=True,
                                  mlt_range=[9, 15], mlat_range=[70, 89],
                                  plot_both_wic_si=False, overlay_ace=False,
                                  overlay_cluster=True,
                                  xcorr_pair=xcorr_pair,
                                  spot_patch_size=(0,0))
        ax_mlt=None 
        ax_mlat=None
        ax_intensity = obj.axes[6]
        xcorr_pair={"mlt":"By", "mlat":"Bz", "intensity":"N_HIA"}
        obj.plot_wic_si_max_point(ax_mlt=ax_mlt, ax_mlat=ax_mlat, 
                                  ax_intensity=ax_intensity,
                                  delay_in_response=10,
                                  ace_delay_at_cluster=48,
                                  read_si=False, read_wic=True,
                                  mlt_range=[9, 15], mlat_range=[70, 89],
                                  plot_both_wic_si=False, overlay_ace=False,
                                  overlay_cluster=True,
                                  xcorr_pair=xcorr_pair,
                                  spot_patch_size=(0,0))


        fig_name =  stime.strftime("%Y%m%d.%H%M.")+etime.strftime("%Y%m%d.%H%M.")+\
                    'wic_si_spot.png' 

    ###########################################################################
    # Plot all the parameters



    # Adjust the tick params
    # format the datetime xlabels
    if (etime-stime).days >= 2:
        obj.axes[-1].xaxis.set_major_formatter(DateFormatter('%m/%d'))
        locs = obj.axes[-1].xaxis.get_majorticklocs()
        locs = locs[::2]
        locs = np.append(locs, locs[-1]+1)
        obj.axes[-1].xaxis.set_ticks(locs)
    if (etime-stime).days == 0:
        obj.axes[-1].xaxis.set_major_formatter(DateFormatter('%H:%M'))

    # adding extra xtick labels
#    if vline:
#        from matplotlib import dates
#        poss = np.append(axes[-1].xaxis.get_majorticklocs(), axes[-1].xaxis.get_minorticklocs())
#        poss = np.append(poss, dates.date2num([dt1, dt2]))
#        poss = np.delete(poss, 2)
#        axes[-1].xaxis.set_ticks(poss)

    # rotate xtick labels
    plt.setp(obj.axes[-1].get_xticklabels(), rotation=30)

    # set axis label and title
    obj.axes[-1].set_xlabel('Time UT')
    obj.axes[-1].xaxis.set_tick_params(labelsize=11)
    if (etime-stime).days > 0:
        obj.axes[0].set_title('  ' + 'Date: ' +\
                stime.strftime("%m/%d/%Y") + ' - ' + (etime-dt.timedelta(days=1)).strftime("%m/%d/%Y")\
                + '    ACE SW data')
    else:
        obj.axes[0].set_title(stime.strftime("%m/%d/%Y"))


    # Save the figure
    dpi = 300
    obj.fig.savefig('../plots/multipanel_time_series/' + fig_name ,dpi=dpi)
    plt.close(obj.fig)
    #plt.show()


    ###########################################################################

