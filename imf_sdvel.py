from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.dates import DateFormatter 
from matplotlib.ticker import MultipleLocator
from spacepy import pycdf
import matplotlib.pyplot as plt
import matplotlib as mpl
import datetime as dt
import bisect
import numpy as np
from davitpy import gme
from glob import glob
import pandas as pd

def params_panel_plot(stime, etime, marker='.', linestyle='--', markersize=2, vline=False,
             panel_num=3, fig_size=None, hspace=None):

    fig, axes = draw_axes(fig_size=fig_size, panel_num=panel_num, hspace=hspace)
    ax_IMF, ax_theta = axes[0], axes[2]
    ax_Pdyn = axes[3]
    #ax_V = axes[4]
    ax_V = None
    
    # plot ace solar wind and IMF data
    date_ylim_IMF_dict = {"20020318" : [-20, 25],
                          "20030214" : [-15, 15],
                          "20040325" : [-9, 9],
                          "20060306" : [-12, 9]}
    plot_ace(stime, etime, ax_V=ax_V, ax_Pdyn=ax_Pdyn, ax_IMF=ax_IMF, ax_theta=ax_theta,
             ylim_IMF=date_ylim_IMF_dict[stime.strftime("%Y%m%d")], marker=marker,
             linestyle=linestyle, zero_line=True,
             markersize=markersize, drop_na=False)


    date_ylim_c_B_dict = {"20020318" : [-80, 100],
                          "20030214" : [-60, 60],
                          "20040325" : [-40, 55],
                          "20060306" : [-70, 50]}

    date_ylim_theta_dict = {"20020318" : [-180, 100],
                            "20030214" : [-180, 180],
                            "20040325" : [-90, 100],
                            "20060306" : [60, 300]}

    ax_c_B, ax_c_Btheta = axes[1], axes[2]
    plot_cluster(stime, etime, ax_c_B=ax_c_B, ax_c_Btheta=ax_c_Btheta,
        ylim_c_B=date_ylim_c_B_dict[stime.strftime("%Y%m%d")],
        ylim_c_theta=date_ylim_theta_dict[stime.strftime("%Y%m%d")],
        marker=marker, linestyle=linestyle, markersize=markersize,
        ylabel_fontsize=9, zero_line=False,
        interpolation_method=None, drop_na=False)


#    ax_imf_theta = axes[3]
#    plot_ace(stime, etime, ax_V=None, ax_Pdyn=None, ax_IMF=None, ax_theta=ax_imf_theta,
#             ylim_IMF=date_ylim_theta_dict[stime.strftime("%Y%m%d")] ,zero_line=True,
#             marker=marker, linestyle=linestyle,
#             markersize=markersize, drop_na=False)
#
#    # plot sd fitted vels from actural measurement points 
#    #lon_range = [11, 13]   # mlt hours
#    lon_range = [11.6, 13.2]   # mlt hours
#    lon_del = 0.5  # degree
#    #lat_range = np.arange(79.5, 82.5, 0.5)  # mag lats
#    lat_range = np.arange(76.5, 84.5, 0.5)  # mag lats
#    #lat_range = np.arange(78.5, 83.5, 1)  # mag lats
#    color_list = ["r", "b", "g", "c", "y", "k"]
#    ax_vel_angle = axes[3]
#    stime = dt.datetime(2003,2,14,20,00) 
#    for j, lat_tmp in enumerate(lat_range):
#        plot_vel_actual(stime, etime, ftype="mapex", ax_vel_angle=ax_vel_angle, ylim_vel_angle=[-60, 200],
#                 lon_range=lon_range, lat_point=lat_tmp, lon_del = lon_del,
#                 #color=color_list[j], marker=marker, linestyle=linestyle,
#                 color='r', marker=marker, linestyle=linestyle,
#                 markersize=markersize)


    # plot sd fitted vels from arbiturary points where there may not have actual radar measurements 

    date_ylim_vel_angle_dict = {"20020318" : [-80, 270],
                                "20030214" : [-60, 200],
                                "20040325" : [-60, 270],
                                "20060306" : [60, 300]}
#    ax_vel_angle = axes[3]
#    plot_vel(stime, etime, ax_vel_mag=None, ax_vel_angle=ax_vel_angle,
#               ylim_vel_mag=[-1000, 1000],
#               ylim_vel_angle=date_ylim_vel_angle_dict[stime.strftime("%Y%m%d")],
#               ylabel_fontsize=9, marker='.', linestyle='-', markersize=2,
#               padding_method=None, zero_line=True)


#    # plot symh
#    ax_symh = axes[4]
#    plot_symh(stime, etime, ax_symh, ylim_symh=[-50, 50],
#            marker=marker, linestyle=linestyle, markersize=markersize, zero_line=True)
#
#    # plot AE
#    ax_ae = axes[5]
#    plot_ae(stime, etime, ax_ae, ylim_ae=[0, 1000],
#        marker=marker, linestyle=linestyle, markersize=markersize)

#    # plot Kp
#    ax_kp = axes[5]
#    plot_kp(stime, etime, ax_kp, ylim_kp=[0, 9],
#        marker=marker, linestyle=linestyle, markersize=markersize)


    # adjusting tick parameters
    for ll in range(panel_num): 
        axes[ll].yaxis.set_tick_params(labelsize=11)
        #axes[ll].xaxis.grid(True, which='major')
        axes[ll].yaxis.grid(True, which='major')

#        # lining up the ylabels and setting ylabel fontsize
#        labelx = -0.10  # axes coords
#        axes[ll].yaxis.set_label_coords(labelx, 0.5)
#        axes[ll].yaxis.label.set_size(11)
#
#        # add vlines at the points of interest
#        if vline:
#            dt1 = dt.datetime(2014,9,12,15,26)
#            dt2 = dt.datetime(2014,9,12,15,55)
#            dt3 = dt.datetime(2014,9,12,16,57)
#            #dt4 = dt.datetime(2014,9,12,17,51)
#            axes[ll].axvline(dt1, color='r', linewidth=0.50, linestyle='--')
#            axes[ll].axvline(dt2, color='r', linewidth=0.50, linestyle='--')
#            axes[ll].axvline(dt3, color='r', linewidth=0.50, linestyle='--')
#            #axes[ll].axvline(dt4, color='r', linewidth=0.50, linestyle='--')
#
    # format the datetime xlabels
    if (etime-stime).days >= 2:   
        axes[-1].xaxis.set_major_formatter(DateFormatter('%m/%d'))
        locs = axes[-1].xaxis.get_majorticklocs()
        locs = locs[::2]
        locs = np.append(locs, locs[-1]+1)
        axes[-1].xaxis.set_ticks(locs)
    if (etime-stime).days == 0:   
        axes[-1].xaxis.set_major_formatter(DateFormatter('%H:%M'))
#
#    # adding extra xtick labels
#    if vline:
#        from matplotlib import dates
#        poss = np.append(axes[-1].xaxis.get_majorticklocs(), axes[-1].xaxis.get_minorticklocs())
#        poss = np.append(poss, dates.date2num([dt1, dt2]))
#        poss = np.delete(poss, 2)
#        axes[-1].xaxis.set_ticks(poss)
#
#    # remove xtick labels of the subplots except for the last one
#    for i in range(panel_num):
#        if i < panel_num-1:
#            axes[i].tick_params(axis='x',which='both',labelbottom='off') 
#
    # rotate xtick labels
    plt.setp(axes[-1].get_xticklabels(), rotation=30)

    # set axis label and title
    axes[-1].set_xlabel('Time UT')
    axes[-1].xaxis.set_tick_params(labelsize=11)
    if (etime-stime).days > 0:
        axes[0].set_title('  ' + 'Date: ' +\
                stime.strftime("%Y/%m/%d") + ' - ' + (etime-dt.timedelta(days=1)).strftime("%Y/%m/%d")\
                + '    IMF and SuperDARN Noon HL Flow')
    else:
        #axes[0].set_title('  ' + 'Date: ' + stime.strftime("%Y/%m/%d") + '    IMF and SuperDARN Noon HL Flow')
        axes[0].set_title('  ' + 'Date: ' + stime.strftime("%Y/%m/%d") + '    IMF and Cluster Data')

    # add annotations 
    # dicts to be used for annotation
    date_cname_dict = {"20020318" : "Cluster C3",
                       "20030214" : "Cluster C1",
                       "20040325" : "Cluster C3",
                       "20060306" : "Cluster C3"}

    date_lagtxt_dict = {"20020318" : "Tp=48min",
                        "20030214" : "Tp=45min",
                        "20040325" : "Tp=75min",
                        "20060306" : "Tp=70min"}

    #ax_IMF.annotate('ACE', xy=(0.85, 0.9), xycoords='axes points', xytext=(0, 5), ha='center',
    #        textcoords='offset points')

    ax_IMF.annotate('ACE', xy=(0.90, 0.89), xycoords='axes fraction', ha='center', fontsize=9)
    ax_c_B.annotate(date_cname_dict[stime.strftime("%Y%m%d")], xy=(0.90, 0.89),
            xycoords='axes fraction', ha='center', fontsize=9)

    if stime.strftime("%Y%m%d") == "20060306":
        xy_coord_ace = (0.90, 0.85)
        xy_coord_cname = (0.90, 0.77)
        xy_coord_lagtxt = (0.90, 0.69)
    else:
        xy_coord_ace = (0.90, 0.25)
        xy_coord_cname = (0.90, 0.17)
        xy_coord_lagtxt = (0.90, 0.09)

    ax_c_Btheta.annotate("ACE", xy=xy_coord_ace,
            xycoords='axes fraction', ha='center', fontsize=9)
    ax_c_Btheta.annotate(date_cname_dict[stime.strftime("%Y%m%d")], xy=xy_coord_cname,
            xycoords='axes fraction', ha='center', fontsize=9, color="r")
    ax_c_Btheta.annotate(date_lagtxt_dict[stime.strftime("%Y%m%d")], xy=xy_coord_lagtxt,
            xycoords='axes fraction', ha='center', fontsize=9)

    return fig

def draw_axes(fig_size=None, panel_num=6, hspace=None):
    
    fig, axes = plt.subplots(nrows=panel_num, figsize=fig_size, sharex=True)
    if hspace is not None:
        fig.subplots_adjust(hspace=hspace)

    return fig, axes

def plot_ace(stime, etime, ax_V=None, ax_Pdyn=None, ax_IMF=None,ax_theta=None,
        ylim_V=[0,800],ylim_Pdyn=[10,25], ylim_IMF=[-30, 30],ylim_theta=[-180, 180],
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
    imf_delay = date_lag_dict[stime.strftime("%Y%m%d")]
    ylim_imf = [-40, 40]
    df_imf, df_sw =  ace_read(stime, etime, res=0, how='mean', coord='GSM', delay=imf_delay)

    df_imf.loc[:, 'Bt'] = np.sqrt((df_imf.By ** 2 + df_imf.Bz ** 2))
    # clock angle
    df_imf.loc[:, 'theta_Bt'] = np.degrees(np.arctan2(df_imf.By, df_imf.Bz))
    # dynamic pressure
    mp = 1.6726219 * 1e-27  # kg
    df_sw.loc[:, 'Pdyn'] = (df_sw.Np * mp * 1e6) * (df_sw.Vx * 1e3)**2 * 1e+9
    if drop_na:
        df_imf.dropna(how='any', inplace=True)
        df_sw.dropna(how='any', inplace=True)
    if padding_method is not None:
        df_imf.fillna(method=padding_method, inplace=True)
        df_sw.fillna(method=padding_method, inplace=True)

    # plot solar wind Vx
    if ax_V is not None:
        df_sw.loc[df_sw.Vx < -1e4, 'Vx'] = np.nan 
        ax_V.plot_date(df_sw.index.to_pydatetime(), -df_sw.Vx, color='k',
                marker=marker, linestyle=linestyle, markersize=markersize)
        ax_V.set_ylabel('V [km/s]', fontsize=ylabel_fontsize)
        ax_V.set_ylim([ylim_V[0], ylim_V[1]])
        ax_V.locator_params(axis='y', nbins=4)
        # add text to indicate sheath region
        dt1 = dt.datetime(2014,9,12,15,26)
        dt2 = dt.datetime(2014,9,12,16,10)
        dt3 = dt.datetime(2014,9,12,16,57)
        ax_V.annotate('', xy=(dt1, 400), xycoords='data', xytext=(dt3, 400),
                arrowprops={'arrowstyle': '<->'})
        ax_V.annotate('Sheath', xy=(dt2, 400), xycoords='data', xytext=(0, 5), ha='center',
                textcoords='offset points')

    # plot Pdyn [nPa]
    if ax_Pdyn is not None:
        df_sw.loc[df_sw.Pdyn < 0, 'Pdyn'] = np.nan
        ax_Pdyn.plot_date(df_sw.index.to_pydatetime(), df_sw.Pdyn, color='k',
                marker=marker, linestyle=linestyle, markersize=markersize)
        ax_Pdyn.set_ylabel('Pdyn [nPa]\n', fontsize=ylabel_fontsize)
        ax_Pdyn.set_ylim([ylim_Pdyn[0], ylim_Pdyn[1]])
        ax_Pdyn.locator_params(axis='y', nbins=4)

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
    if stime.strftime("%Y%m%d") == "20060306":
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

#        # add Bt component
#        ax_theta_twinx = ax_theta.twinx()
#        ax_theta_twinx.plot_date(df_imf.index.to_pydatetime(), df_imf.Bt, color='b',
#                marker=marker, linestyle=linestyle, markersize=markersize)
#        ax_theta_twinx.locator_params(axis='y', nbins=4)
#        ax_theta_twinx.set_ylim([ylim_IMF[0], ylim_IMF[1]])
#        ax_theta_twinx.set_ylabel('IMF |Bt| [nT]', color='b')
#        ax_theta_twinx.yaxis.label.set_size(11)
#        ax_theta_twinx.yaxis.set_tick_params(labelsize=11, labelcolor='b')



def plot_symh(stime, etime, ax_symh, ylim_symh=[-100, 100], ylabel_fontsize=9,
        marker='.', linestyle='--', markersize=2, zero_line=True):
    # read SYMH data
    #sym_list = gme.ind.symasy.readSymAsyWeb(stime=stime,etime=etime)
    sym_list = gme.ind.symasy.readSymAsy(sTime=stime,eTime=etime)
    symh = []
    symh_time = []
    for k in range(len(sym_list)):
        symh.append(sym_list[k].symh)
        symh_time.append(sym_list[k].time)
    
    # plot symh
    ax_symh.plot_date(symh_time, symh, 'k', marker=marker, linestyle=linestyle, markersize=markersize)
    ax_symh.set_ylabel('SYM-H', fontsize=9)
    ax_symh.set_ylim([ylim_symh[0], ylim_symh[1]])
    ax_symh.locator_params(axis='y', nbins=4)
    if zero_line:
        ax_symh.axhline(y=0, color='r', linewidth=0.30, linestyle='--')

def plot_ae(stime, etime, ax_ae, ylim_ae=[0, 500], ylabel_fontsize=9,
        marker='.', linestyle='--', markersize=2):
    # read AE data
    #AE_list = gme.ind.readAeWeb(stime=stime,etime=etime,res=1)
    AE_list = gme.ind.readAe(sTime=stime,eTime=etime,res=1)
    AE = []
    AE_time = []
    for m in range(len(AE_list)):
        AE.append(AE_list[m].ae)
        AE_time.append(AE_list[m].time)

    # plot AE 
    ax_ae.plot_date(AE_time, AE, 'k', marker=marker, linestyle=linestyle, markersize=markersize)
    ax_ae.set_ylabel('AE', fontsize=9)
    ax_ae.set_ylim([ylim_ae[0], ylim_ae[1]])
    ax_ae.locator_params(axis='y', nbins=4)

def plot_kp(stime, etime, ax_kp, ylim_kp=[0, 9], ylabel_fontsize=9,
        marker='.', linestyle='--', markersize=2):
    # Kp data
    Kp_list = gme.ind.readKp(stime=stime,etime=etime)
    Kp = []
    Kp_time = []
    for n in range((etime-stime).days):
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

def plot_cluster(stime, etime, ax_c_B=None,ax_c_Btheta=None,
        ylim_c_B=[-60, 60],ylim_c_theta=[-180, 180],
        marker='.', linestyle='--', markersize=2, ylabel_fontsize=9, zero_line=True,
        interpolation_method=None, drop_na=False):

    # read data
    date_parser = lambda x: pd.datetime.strptime(x, '%d-%m-%Y %H:%M:%S.%f')
    #date_parser = lambda x, y: pd.datetime.strptime(" ".join([str(x), str(y)]), '%d-%m-%Y %H:%M:%S.%f')

    # reading the data into pandas' DataFrame object
    date_fname_dict = {"20020318" : "C3_CP_FGM_SPIN_31086.txt",
                       "20030214" : "C1_CP_FGM_SPIN_21206.txt",
                       "20040325" : "C3_CP_FGM_SPIN_3152599.txt",
                       "20060306" : "C3_CP_FGM_SPIN_31759.txt"}
    fname = "./data/" + date_fname_dict[stime.strftime("%Y%m%d")]
    date_srows_dict = {"20020318" : 114,
                       "20030214" : 127,
                       "20040325" : 114,
                       "20060306" : 114}

    srows = range(date_srows_dict[stime.strftime("%Y%m%d")])
    df_c1 = pd.read_csv(fname, index_col=0, skipinitialspace=True,
            delim_whitespace=True, skiprows=srows,
            header=None,
            parse_dates={'datetime': [0,1]},
            date_parser=date_parser,
            error_bad_lines=False)

    df_c1.columns = ["B", "BX_GSE", "BY_GSE", "BZ_GSE", "X_GSE", "Y_GSE", "Z_GSE"]

    # select the period of interest
    df_c1 = df_c1.loc[stime:etime,]

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
    if stime.strftime("%Y%m%d") == "20060306":
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
        poss1 = np.delete(poss1, [0, 1, 3, 4])
        ax_c_Btheta.yaxis.set_ticks(poss1)
        # draw horizontal lines
        if zero_line:
            ax_c_Btheta.axhline(y=0, color='r', linewidth=0.30, linestyle='--')
        ninety_line = True 
        if ninety_line:
            ax_c_Btheta.axhline(y=90, color='r', linewidth=0.30, linestyle='--')
            ax_c_Btheta.axhline(y=-90, color='r', linewidth=0.30, linestyle='--')


def plot_vel_actual(stime, etime, ftype="grdex", ax_vel_mag=None, ax_vel_angle=None,
               lon_range=[11, 13], lat_point=80.5, lon_del = 0.5,
               ylim_vel_mag=[-1000, 1000], ylim_vel_angle=[-100, 100],
               ylabel_fontsize=9, marker='.', linestyle='--', markersize=2,
               padding_method=None, zero_line=True, color="r"):

    df_vel_mag, df_vel_angle = read_point_sdvel(stime, etime, ftype=ftype, coord="mlt",
                              lon_range=lon_range, lat_point=lat_point, lon_del=lon_del)
    if padding_method is not None:
        df_vel_mag.fillna(method=padding_method, inplace=True)
        df_vel_angle.fillna(method=padding_method, inplace=True)

    # plot velocity magnitude data
    color_list = ["r", "b", "g", "c", "y", "k"]
    if ax_vel_mag is not None:
        for clmn in df_vel_mag.columns:
        #for j, clmn in enumerate(df_vel_angle.columns):
            ax_vel_mag.plot_date(df_vel_mag.index.to_pydatetime(), df_vel_mag.loc[:,clmn],
                    color=color, marker=marker, linestyle=linestyle, markersize=markersize)
                    #color=color_list[j], marker=marker, linestyle=linestyle, markersize=markersize)
        ax_vel_mag.set_ylabel('Vel. [m/s]', fontsize=9)
        ax_vel_mag.set_ylim([ylim_vel_mag[0], ylim_vel_mag[1]])
        ax_vel_mag.locator_params(axis='y', nbins=4)
#        lns = ax_vel_mag.get_lines()
#        ax_vel_mag.legend(lns,['NS', 'EW'], frameon=False, bbox_to_anchor=(0.98, 0.5),
#                loc='center left', fontsize='medium')
        if zero_line:
            ax_vel_mag.axhline(y=0, color='r', linewidth=0.10, linestyle='--')

    # plot velocity angle data
    if ax_vel_angle is not None:
        for clmn in df_vel_angle.columns:
        #for j, clmn in enumerate(df_vel_angle.columns):
            ax_vel_angle.plot_date(df_vel_angle.index.to_pydatetime(), df_vel_angle.loc[:,clmn],
                                   color=color, marker=marker, linestyle=linestyle,
                                   #color=color_list[j], marker=marker, linestyle=linestyle,
                                   markersize=markersize)
        ax_vel_angle.set_ylabel('Flow Direction\n'+r"$(^\circ)$", fontsize=9)
        ax_vel_angle.set_ylim([ylim_vel_angle[0], ylim_vel_angle[1]])
        ax_vel_angle.locator_params(axis='y', nbins=4)
    #    lns = ax_vel_angle.get_lines()
    #    ax_vel_angle.legend(lns,['NS', 'EW'], frameon=False, bbox_to_anchor=(0.98, 0.5),
    #            loc='center left', fontsize='medium')
        if zero_line:
            ax_vel_angle.axhline(y=0, color='r', linewidth=0.30, linestyle='--')


def plot_vel(stime, etime, ax_vel_mag=None, ax_vel_angle=None,
               ylim_vel_mag=[-1000, 1000], ylim_vel_angle=[-100, 100],
               ylabel_fontsize=9, marker='.', linestyle='--', markersize=2,
               padding_method=None, zero_line=True):

    df_vel_mag, df_vel_angle = rad_map_read_fitvel(stime, etime)

    write_output = True
    #write_output = False
    ########## write the output into file ####################
    if write_output:
        vel_mag_output = df_vel_mag.copy()
        vel_mag_output.columns = ["vel_magnitude"]
        vel_mag_output.to_csv("./data/data_to_stefan/moving_points_20161223/" +
                #"20060306_FitVelMagnitude_Moving_point_low_lat.csv")
                "20060306_FitVelMagnitude_Moving_point_high_lat.csv")
        vel_angle_output = df_vel_angle.copy()
        vel_angle_output.columns = ["vel_angle"]
        vel_angle_output.to_csv("./data/data_to_stefan/moving_points_20161223/" + 
                #"20060306_FitVelAngle_Moving_point_low_lat.csv")
                "20060306_FitVelAngle_Moving_point_high_lat.csv")

    ##########################################################
    if padding_method is not None:
        df_vel_mag.fillna(method=padding_method, inplace=True)
        df_vel_angle.fillna(method=padding_method, inplace=True)

    color_list = ["r", "b", "g", "c", "y", "k"]
    # plot velocity magnitude data
    if ax_vel_mag is not None:
        #for clmn in df_vel_mag.columns:
        for j, clmn in enumerate(df_vel_angle.columns):
            ax_vel_mag.plot_date(df_vel_mag.index.to_pydatetime(), df_vel_mag.loc[:,clmn],
                    #color=color, marker=marker, linestyle=linestyle, markersize=markersize)
                    color=color_list[j], marker=marker, linestyle=linestyle, markersize=markersize)
        ax_vel_mag.set_ylabel('Vel. [m/s]', fontsize=9)
        ax_vel_mag.set_ylim([ylim_vel_mag[0], ylim_vel_mag[1]])
        ax_vel_mag.locator_params(axis='y', nbins=4)
        if zero_line:
            ax_vel_mag.axhline(y=0, color='r', linewidth=0.10, linestyle='--')

    # plot velocity angle data
    if ax_vel_angle is not None:
        #for clmn in df_vel_angle.columns:
        for j, clmn in enumerate(df_vel_angle.columns):
            ax_vel_angle.plot_date(df_vel_angle.index.to_pydatetime(), df_vel_angle.loc[:,clmn],
                                   #color=color, marker=marker, linestyle=linestyle,
                                   color=color_list[j], marker=marker, linestyle=linestyle,
                                   markersize=markersize)
        ax_vel_angle.set_ylabel('Flow Direction\n'+r"$(^\circ)$", fontsize=9)
        ax_vel_angle.set_ylim([ylim_vel_angle[0], ylim_vel_angle[1]])
        ax_vel_angle.locator_params(axis='y', nbins=4)
        #switch = True
        switch = False
        if switch:
            lns = ax_vel_angle.get_lines()[2:]
            ax_vel_angle.legend(lns,['78', '78.5', '79', '79.5', '80'], frameon=False,
            #ax_vel_angle.legend(lns,['84', '84.5', '85', '85.5'], frameon=False,
                    bbox_to_anchor=(0.98, 0.5), loc='center left', fontsize='small')
            ax_vel_angle.annotate('Flow at 12 MLT', xy=(0.82, 0.43), xycoords='axes fraction', ha='center', fontsize=9)
        if zero_line:
            ax_vel_angle.axhline(y=0, color='r', linewidth=0.30, linestyle='--')


def rad_map_read_fitvel(stime, etime):

    # read fitted velocity data generated by rad_map_fitvel.pro

    date_parser = lambda x: pd.datetime.strptime(x, '%Y%m%d-%H%M')
    # reading the data into pandas' DataFrame object
    df_vel_angle = pd.read_csv("./data/"+stime.strftime("%Y%m%d")+"_fitvel_angle.txt", index_col=0,
                               delim_whitespace=True, parse_dates=[0], date_parser=date_parser)
    df_vel_mag = pd.read_csv("./data/"+stime.strftime("%Y%m%d")+"_fitvel_mag.txt", index_col=0,
                             delim_whitespace=True, parse_dates=[0], date_parser=date_parser)

    df_vel_angle = df_vel_angle.loc[stime:etime, :]
    df_vel_mag = df_vel_mag.loc[stime:etime, :]

    def my_func(row):
        for clmn in row.keys():
            if row[clmn] < 0:
                row[clmn] = row[clmn] + 180
            else:
                row[clmn] = row[clmn] - 180
            row[clmn] = row[clmn] % 360
            if row[clmn] > 350:
            #if row[clmn] > 270:
                row[clmn] = row[clmn] - 360
        return row
    # set the velocity direction parametr, kvect, such that 0 deg is sunward, 90 is dawnward,
    # 180 is antisunward and -90 is duskward
    df_vel_angle.apply(my_func, axis=1)

#    # convert mag to mlt
#    xlon, xlat = coord_conv(sdrec.grid.vector.mlon, sdrec.grid.vector.mlat, "mag", coord, 
#                            altitude=100., date_time=sdrec.eTime)

    return df_vel_mag, df_vel_angle


def read_point_sdvel(stime, etime, hemi="north", ftype="grdex",
                     coord="mlt", lon_range=[11, 13], lat_point=80.5,
                     lon_del = 0.5):

    from davitpy.pydarn.sdio.sdDataRead import sdDataOpen, sdDataReadAll
    from davitpy.utils.coordUtils import coord_conv

    my_ptr = sdDataOpen(stime, hemi=hemi, eTime=etime, fileType=ftype)
    my_list = sdDataReadAll(my_ptr)

    # convert mlt_lon to mlon
    lon_range_tmp =np.arange(15*lon_range[0], 180, lon_del)
    np.append(lon_range_tmp, np.arange(-180, (15*lon_range[1] - 360), lon_del))
    lon_range = lon_range_tmp
    tms = []
    df_vel_mag = pd.DataFrame(index=range(len(my_list)), columns=lon_range) 
    df_vel_angle = pd.DataFrame(index=range(len(my_list)), columns=lon_range) 
    for k, sdrec in enumerate(my_list):
        # convert mag to mlt
        if ftype in ["grd", "grdex"]:
            xlon, xlat = coord_conv(sdrec.vector.mlon, sdrec.vector.mlat, "mag", coord, 
                                    altitude=100., date_time=sdrec.eTime)
        elif ftype in ["map", "mapex"]:
            xlon, xlat = coord_conv(sdrec.grid.vector.mlon, sdrec.grid.vector.mlat, "mag", coord, 
                                    altitude=100., date_time=sdrec.eTime)
        #lons_tmp = []
        vel_mag, vel_angle = [], []
        for i in range(len(lon_range)):
            diffs = np.array(xlon) - lon_range[i]
            if np.min(abs(diffs)) <= lon_del:
                min_indx = np.argmin(abs(diffs))
                if (xlat[min_indx] - lat_point) == 0:
                    #lons_tmp.append(xlon[min_indx])
                    if ftype in ["grd", "grdex"]:
                        vel_mag.append(sdrec.vector.velmedian[min_indx])
                        kvect = sdrec.vector.kvect[min_indx]
                    elif ftype in ["map", "mapex"]:
                        vel_mag.append(sdrec.grid.vector.velmedian[min_indx])
                        kvect = sdrec.grid.vector.kvect[min_indx]
                    # set the velocity direction parametr, kvect, such that 0 deg is sunward, 90 is dawnward,
                    # 180 is antisunward and -90 is duskward
                    if kvect < 0:
                        kvect = kvect + 180
                    else:
                        kvect = kvect - 180
                    vel_angle.append(kvect)
                else:
                    vel_mag.append(np.nan)
                    vel_angle.append(np.nan)
            else:
                vel_mag.append(np.nan)
                vel_angle.append(np.nan)


        # populate the empty dataframe
        vel_mag_dict = dict(zip(lon_range, vel_mag))
        vel_angle_dict = dict(zip(lon_range, vel_angle))
        df_vel_mag.loc[k] = vel_mag_dict 
        df_vel_angle.loc[k] = vel_angle_dict 

        tms.append(sdrec.sTime)

    # change the dataframe index into datetime index
    df_vel_mag.index = tms
    df_vel_angle.index = tms

    return df_vel_mag, df_vel_angle

def plot_fitvel_loc(stime, etime):

    fig_pos, ax_pos = plt.subplots(nrows=2, sharex=True)
    ax_mlat, ax_mltlon = ax_pos

    date_parser = lambda x: pd.datetime.strptime(x, '%Y%m%d-%H%M')
    # reading the data into pandas' DataFrame object
    df_mlat = pd.read_csv("./data/"+stime.strftime("%Y%m%d")+"_fitvel_mlatloc.txt", index_col=0,
                               delim_whitespace=True, parse_dates=[0], date_parser=date_parser)
    df_mltlon = pd.read_csv("./data/"+stime.strftime("%Y%m%d")+"_fitvel_mltlonloc.txt", index_col=0,
                             delim_whitespace=True, parse_dates=[0], date_parser=date_parser)

    df_mlat = df_mlat.loc[stime:etime, :]
    df_mltlon = df_mltlon.loc[stime:etime, :]

    write_output = True
    #write_output = False
    ########## write the output into file ####################
    if write_output:
        vel_loc_output = df_mlat.copy()
        vel_loc_output.columns = ["mlat"]
        vel_loc_output["mlt"] = df_mltlon.as_matrix()/15.
        vel_loc_output.to_csv("./data/data_to_stefan/moving_points_20161223/" +
                #"20060306_fitvel_low_latitude_mlt_mlat.csv")
                "20060306_fitvel_high_latitude_mlt_mlat.csv")


    #################################################################
    # plot  data
    marker=''; linestyle='-'; markersize=1
    color_list = ["r", "b", "g", "c", "y", "k"]
    for j, clmn in enumerate(df_mlat.columns):
        ax_mlat.plot_date(df_mlat.index.to_pydatetime(), df_mlat.loc[:,clmn],
                #color=color, marker=marker, linestyle=linestyle, markersize=markersize)
                color=color_list[j], marker=marker, linestyle=linestyle, markersize=markersize)
    for j, clmn in enumerate(df_mltlon.columns):
        ax_mltlon.plot_date(df_mltlon.index.to_pydatetime(), df_mltlon.loc[:,clmn]/15.,
                #color=color, marker=marker, linestyle=linestyle, markersize=markersize)
                color=color_list[j], marker=marker, linestyle=linestyle, markersize=markersize)

    ax_mlat.set_ylabel('MLAT [degree]', fontsize=9)
    ax_mltlon.set_ylabel('MLT [hr]', fontsize=9)
    ax_mlat.set_ylim([75, 90])
    ax_mltlon.set_ylim([9, 15])
    ax_mltlon.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    ax_mlat.set_title("Location (MLAT, MLT) of Extracted Velocity Vectors    Date:" + \
            stime.strftime("%Y/%m/%d"), fontsize=8)
    return fig_pos


# testing
import datetime as dt

stime = dt.datetime(2002,3,18,15,00) 
etime = dt.datetime(2002,3,18,17,20) 

#stime = dt.datetime(2003,2,14,19,30) 
#etime = dt.datetime(2003,2,14,22,30) 

#stime = dt.datetime(2004,3,25,8,00) 
#etime = dt.datetime(2004,3,25,11,00) 


#stime = dt.datetime(2006,3,06,16,59) 
#etime = dt.datetime(2006,3,06,19,30) 
#stime = dt.datetime(2006,3,06,18,19) 
#etime = dt.datetime(2006,3,06,18,37) 

panel_num = 4; fig_size=(8,10); hspace=None
dpi = 300
fig = params_panel_plot(stime, etime, marker='', linestyle='-', markersize=1, vline=True,
             panel_num=panel_num, fig_size=fig_size, hspace=0.15)
fig.savefig('./plots/' +\
            stime.strftime("%Y%m%d.%H%M.")+etime.strftime("%Y%m%d.%H%M.")+\
            'Pdyn_imf_info.png' ,dpi=dpi)

#fig.savefig('./' +\
        #stime.strftime("%Y%m%d.%H%M.")+etime.strftime("%Y%m%d.%H%M.")+'imf_sdfitvel_moving_points_high_lat.pdf', format='pdf')
#        stime.strftime("%Y%m%d.%H%M.")+etime.strftime("%Y%m%d.%H%M.")+'imf_sdfitvel_moving_points_high_lat.png', dpi=dpi)
#plt.close(fig)
#fig.show()

#fig_pos = plot_fitvel_loc(stime, etime)
#fig_pos.savefig('./' +\
        #stime.strftime("%Y%m%d.%H%M.")+etime.strftime("%Y%m%d.%H%M.")+'fitvel_loc_fixedpoints.png' ,dpi=dpi)
#        stime.strftime("%Y%m%d.%H%M.")+etime.strftime("%Y%m%d.%H%M.")+'fitvel_high_latitude_mlt_mlat_loc.png' ,dpi=dpi)
#plt.show()

