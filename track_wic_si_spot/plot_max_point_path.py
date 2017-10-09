def plot_max_point_path(stm, etm, read_si=True, read_wic=False, param="mlt",
                        mlt_range=[8, 16], mlat_range=[70, 85],
                        plot_both_wic_si=True, overlay_IMF=False):
    """ Plots the path of the max intensity point motion"""
    from track_max_point import find_filenames, read_wic_si_data
    from track_max_point import find_max_intensity_point
    from matplotlib.dates import DateFormatter
    import matplotlib.pyplot as plt
    import numpy as np        

    fig, ax = plt.subplots()

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
            read_si = read_si; read_wic = read_wic
            if read_si:
                label_txt = "SI"
            if read_wic:
                label_txt = "WIC"
        else:
            read_si = si_bool[i]; read_wic = wic_bool[i] 
            label_txt = legend_text[i]
        dtms, fnames = find_filenames(stm, etm, read_si=read_si, read_wic=read_wic)
        mlats = []
        mlts = []
        vals = [] 
        for dtm in dtms:            
            fname, data = read_wic_si_data(dtm, read_si=read_si, read_wic=read_wic)
            (maxval, maxloc) = find_max_intensity_point(data, mlt_range=mlt_range,
                                                        mlat_range=mlat_range)
            mlats.append(data['mlat'][maxloc[1], maxloc[0]]) 
            mlts.append(data['mlt'][maxloc[1], maxloc[0]])
            vals.append(maxval)
            
        if param == "mlat":
            var = mlats
        if param == "mlt":
            var = mlts
        if param == "intensity":
            var = vals
        ax.plot_date(dtms, var, '.-', color=color_list[i], label=label_txt)
        ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
        ax.set_xlabel("Time [UT]")
        ax.set_ylabel(param.upper())

    if overlay_IMF:
        # read ace data that is minute-averaged
        import sys
        sys.path.append('/home/muhammad/Documents/Important/RISRpy/RISR_06222015/codes')
        from read_sw_imf import ace_read
        df_imf, df_sw =  ace_read(stm, etm, res=0, how='mean', coord='GSM', delay=48)

        # clock angle
        df_imf.loc[:, 'theta_Bt'] = np.degrees(np.arctan2(df_imf.By, df_imf.Bz))

        # dynamic pressure
        mp = 1.6726219 * 1e-27  # kg
        df_sw.loc[:, 'Pdyn'] = (df_sw.Np * mp * 1e6) * (df_sw.Vx * 1e3)**2 * 1e+9
        df_sw.loc[df_sw.Pdyn < 0, 'Pdyn'] = np.nan
        df_sw.dropna(how='any', inplace=True)

        # plot IMF data
        ax_IMF = ax.twinx()
        #plot_param="Theta"    # change this accordingly
        plot_param="Pdyn"
        #plot_param="Bx"
        marker='.'; linestyle='-'; markersize=2;

        if plot_param=="Bx":
            ax_IMF.plot_date(df_imf.index.to_pydatetime(), df_imf.Bx, color='k',
                    marker=marker, linestyle=linestyle, markersize=markersize)
        if plot_param=="By":
            ax_IMF.plot_date(df_imf.index.to_pydatetime(), df_imf.By, color='g',
                    marker=marker, linestyle=linestyle, markersize=markersize, label="By")
        if plot_param=="Bz":
            ax_IMF.plot_date(df_imf.index.to_pydatetime(), df_imf.Bz, color='b',
                    marker=marker, linestyle=linestyle, markersize=markersize, label="Bz")
        if plot_param=="Theta":
            ax_IMF.plot_date(df_imf.index.to_pydatetime(), df_imf.theta_Bt, color='k',
                    marker=marker, linestyle=linestyle, markersize=markersize)
        if plot_param=="Pdyn":
            ax_IMF.plot_date(df_sw.index.to_pydatetime(), df_sw.Pdyn, color='k',
                    marker=marker, linestyle=linestyle, markersize=markersize)
        lns = ax_IMF.get_lines()
        #ax_IMF.legend(lns,['Bx', 'By', 'Bz'], frameon=False, fontsize='small', mode='expand')
        if plot_param=="Theta":
            ax_IMF.set_ylim([-150, 150])
            ax_IMF.set_ylabel('IMF Theta ' + "[Degree]")
            ax_IMF.legend(lns,[plot_param], frameon=False, bbox_to_anchor=(0.80, 0.2),
                    loc='center left', fontsize='medium')
        elif plot_param=="Pdyn":
            ax_IMF.set_ylim([10, 25])
            ax_IMF.set_ylabel('Pdyn ' + "[nPa]")
            ax_IMF.legend(lns,[plot_param], frameon=False, bbox_to_anchor=(0.80, 0.2),
                    loc='center left', fontsize='medium')
        else:
            ax_IMF.set_ylim([-25, 25])
            ax_IMF.set_ylabel('IMF ' + "[nT]")
            ax_IMF.locator_params(axis='y', nbins=4)
            ax_IMF.legend(lns,[plot_param], frameon=False, bbox_to_anchor=(0.98, 0.7),
                    loc='center left', fontsize='medium')

        ax_IMF.axhline(y=0, color='k', linewidth=0.40, linestyle='--')
     
    if plot_both_wic_si:
        title = dtm.strftime("%Y/%m/%d")
    else:
        if read_si:
            title = dtm.strftime("%Y/%m/%d") + "   SI Data"
        if read_wic:
            title = dtm.strftime("%Y/%m/%d") + "   WIC Data"
    ax.set_title(title)
    ax.legend(loc="best")

    return fig, ax

# run the code
def main():
    import datetime as dt
    import matplotlib.pyplot as plt

    read_si=True; read_wic=False
    #read_si=False; read_wic=True
    plot_both_wic_si=False
    overlay_IMF = True        
    mlt_range = [9, 15]
    mlat_range = [70, 89]
    param = "mlt"     # parameter to plot
    #param = "mlat"
    #param = "intensity"
    fig_path = "../plots/track_spot/"
    #fig_name = "max_point_" + param + ".png"
    #fig_name = "max_point_" + param + "_with_theta.png"
    #fig_name = "max_point_" + param + "_with_Bx.png"
    fig_name = "max_point_" + param + "_with_Pdyn.png"
    stm = dt.datetime(2002, 3, 18, 15, 00)
    etm = dt.datetime(2002, 3, 18, 17, 20)
    fig, ax = plot_max_point_path(stm, etm, read_si=read_si, read_wic=read_wic,
                                  plot_both_wic_si=plot_both_wic_si,
                                  overlay_IMF=overlay_IMF,
                                  mlt_range=mlt_range, mlat_range=mlat_range,
                                  param=param)
    fig.savefig(fig_path+fig_name, dpi=200)
    #plt.show()

if __name__ == "__main__":
    main()

