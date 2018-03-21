import matplotlib.pyplot as plt
plt.style.use("ggplot")

def plot_ssj_JE_E0_JNE(df, title):
    """Plots energy flux, average energy, and number flux from DMSP SSJ data"""
    
    import numpy as np
    from matplotlib.dates import DateFormatter
        
    fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True) 

    # Plot Energy Flux
    axes[0].plot_date(df.datetime, np.log10(df.JEe), ".k", ms=1.5, label="JEe")
    axes[0].plot_date(df.datetime, np.log10(df.JEi), ".r", ms=1.5, label="JEi")
    axes[0].set_ylabel("Log(JE)\n" + r"[eV/cm$^{2}$ s sr]")
    axes[0].set_ylim([8, 13])

    # Plot Average Energy
    axes[1].plot_date(df.datetime, np.log10(df.AvgEe), ".k", ms=1.5, label="AvgEe")
    axes[1].plot_date(df.datetime, np.log10(df.AvgEi), ".r", ms=1.5, label="AvgEi")
    axes[1].set_ylabel("Log(AvgE)\n" + "[eV]")
    axes[1].set_ylim([1, 5])

    # Plot Number Energy
    axes[2].plot_date(df.datetime, np.log10(df.JNEe), ".k", ms=1.5, label="JNEe")
    axes[2].plot_date(df.datetime, np.log10(df.JNEi), ".r", ms=1.5, label="JNEi")
    axes[2].set_ylabel("Log(JNE)\n" + r"[#/cm$^{2}$ s sr]")
    axes[2].set_ylim([4, 12])

    axes[2].xaxis.set_major_formatter(DateFormatter("%H:%M"))
    axes[0].set_title(title)

    axes[0].legend(frameon=True, fontsize=8, loc="upper right")
    axes[1].legend(frameon=True, fontsize=8, loc="upper right")
    axes[2].legend(frameon=True, fontsize=8, loc="upper right")

    return fig

if __name__ == "__main__":

    from dmsp_ssj_read import dmsp_ssj_read
    import datetime as dt
    #stime = dt.datetime(2002,3,18,17,24)
    #etime = dt.datetime(2002,3,18,17,40)
    stime = dt.datetime(2002,3,18,16,56)
    etime = dt.datetime(2002,3,18,17,12)

    sat_num = 15

    # read ssj data
    df = dmsp_ssj_read(stime, etime, sat_num, file_name=None)

    # Make panel plots of JE, E0 and JNE
    title = stime.strftime("%b%d, %Y") + "    DMSP F" + str(sat_num)
    fig = plot_ssj_JE_E0_JNE(df, title)

    # save fig
    fig_dir = "../plots/DMSP/JE_E0_JNE/"
    fig_name = "dmsp_ssj_F" + str(sat_num) + "_" + stime.strftime("%Y%m%d.%H%M%S") + "_" +\
               etime.strftime("%Y%m%d.%H%M%S") + ".png"

    fig.savefig(fig_dir + fig_name, dpi=200)

    #plt.show()

