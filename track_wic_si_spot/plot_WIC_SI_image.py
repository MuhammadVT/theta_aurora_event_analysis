from track_max_point import find_filenames, read_wic_si_data
from track_max_point import find_max_intensity_point

def plot_WIC_SI_image(dtm, read_si=True, read_wic=False,
                      param="image", mark_max_intensity_point=True,
                      mlt_range=[8, 16], mlat_range=[70, 85]):

    import scipy
    import cv2
    import matplotlib.pyplot as plt


    fname, data = read_wic_si_data(dtm, read_si=read_si, read_wic=read_wic)
    ss_txt = ":" + fname.split("/")[-1].split("_")[-1][4:6]

    fig, ax = plt.subplots()
    if param == "image":
        vmin = 0; vmax = 1000; cmap = 'gist_gray'
    if param == "mlt_img":
        vmin = 0; vmax = 1000; cmap = 'gist_gray'
    if param == "mlat":
        vmin = 60; cmap=None
    if param == "mlon":
        vmin=-180; vmax=180; cmap=None
    if param == "mlt":
        vmin=0; vmax=24; cmap=None

    #img = scipy.misc.bytescale(data[param], cmin=0)
    img = data[param]
    cax = ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax)
    if mark_max_intensity_point:
        if param == "image":
            (maxval, maxloc) = find_max_intensity_point(data, mlt_range=mlt_range,
                                                        mlat_range=mlat_range)
            #cv2.circle(img, maxloc, 9, (255, 0, 0), 2)
            ax.plot(maxloc[0], maxloc[1], 'or', markersize=4)
        else:
            pass
    ax.set_axis_off()
    if read_si:
        title = dtm.strftime("%Y/%m/%d %H:%M") + ss_txt + "   SI " + param
    if read_wic:
        title = dtm.strftime("%Y/%m/%d %H:%M") + ss_txt + "   WIC " + param
    ax.set_title(title)

    cbar = fig.colorbar(cax, shrink=0.7)
    return fig

def iteratively_plot_WIC_SI_image(stm, etm, read_si=True, read_wic=False,
                                  param="image", mark_max_intensity_point=True,
                                  mlt_range=[8, 16], mlat_range=[70, 85]):

    dtms, fnames = find_filenames(stm, etm, read_si=read_si, read_wic=read_wic)
    for dtm in dtms:
        fig = plot_WIC_SI_image(dtm, read_si=read_si, read_wic=read_wic,
                                mark_max_intensity_point=mark_max_intensity_point,
                                mlt_range=mlt_range, mlat_range=mlat_range,
                                param=param)
        if read_si:
            fdir ="../plots/track_spot/SI/"
            fname = dtm.strftime("%H%M") + "_SI_" + param
        if read_wic:
            fdir ="../plots/track_spot/WIC/"
            fname = dtm.strftime("%H%M") + "_WIC_" + param
        fig.savefig(fdir + fname, dpi=300)
    return



# run the code
def main():
    import datetime as dt
    import matplotlib.pyplot as plt

    read_si=True; read_wic=False
    #read_si=False; read_wic=True
    param = "mlt_img"        
    mark_max_intensity_point=True
    mlt_range = [9, 15]
    mlat_range = [70, 89]
    dtm = dt.datetime(2002, 3, 18, 16, 59)
    stm = dt.datetime(2002, 3, 18, 15, 00)
    etm = dt.datetime(2002, 3, 18, 17, 22)
    fig = plot_WIC_SI_image(dtm, read_si=read_si, read_wic=read_wic,
                            mark_max_intensity_point=mark_max_intensity_point,
                            mlt_range=mlt_range, mlat_range=mlat_range,
                            param=param)
    #iteratively_plot_WIC_SI_image(stm, etm, read_si=read_si, read_wic=read_wic,
    #                              param="image", mark_max_intensity_point=True,
    #                              mlt_range=mlt_range, mlat_range=mlat_range)

    plt.show()

if __name__ == "__main__":
    main()

