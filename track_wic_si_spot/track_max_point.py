import datetime as dt
import matplotlib.pyplot as plt
from skimage import feature
import numpy as np


def read_wic_si_data(dtm, read_si=True, read_wic=False):

    """ Reads IMAGE SI and WIC image data
    """
    from scipy.io import readsav
    import glob
    import numpy as np

    if read_si:
        fname_image = "../from_harald_frey/IMF12LSI_2002_0318_" + dtm.strftime("%H%M") + "*.sav"
    if read_wic:
        fname_image = "../from_harald_frey/IMFHWIC_2002_0318_" + dtm.strftime("%H%M") + "*.sav"
    fname_image = glob.glob(fname_image)[0]
    dat = readsav(fname_image, python_dict=True, verbose=False)
    data_dict = {}

    # Note: the data is originally saved in IDL
    # and IDL considers the first dimension to be the column.
    # Therefore, we do the transpose to fix it. Also we need to
    # rotate the image so that the dayside is on the top

#    data_dict["mlt_img"] = dat['imageinfo'].mlt_img[0]
#    data_dict["image"] = dat['imageinfo'].image[0]
#    data_dict["mlat"] = dat['imageinfo'].mlat[0]
#    data_dict["mlon"] = dat['imageinfo'].mlon[0]
#    data_dict["mlt"] = dat['imageinfo'].mlt[0]


    data_dict["mlt_img"] = np.flip(dat['imageinfo'].mlt_img[0], axis=0)
    data_dict["image"] = np.flip(dat['imageinfo'].image[0], axis=0)
    data_dict["mlat"] = np.flip(dat['imageinfo'].mlat[0], axis=0)
    data_dict["mlon"] = np.flip(dat['imageinfo'].mlon[0], axis=0)
    data_dict["mlt"] = np.flip(dat['imageinfo'].mlt[0], axis=0)

    return fname_image, data_dict

def find_filenames(stm, etm, read_si=True, read_wic=False):
    import glob
    import datetime as dt

    dtm = stm
    fnames = []
    dtms = []
    while(dtm<=etm):
        if read_si:
            fname_image = "../from_harald_frey/IMF12LSI_2002_0318_" + dtm.strftime("%H%M") + "*.sav"
        if read_wic:
            fname_image = "../from_harald_frey/IMFHWIC_2002_0318_" + dtm.strftime("%H%M") + "*.sav"
        try:
            fname = glob.glob(fname_image)[0]
            fnames.append(fname)
            dtms.append(dtm)
        except:
            pass
        dtm = dtm + dt.timedelta(minutes=1)
        
    return dtms, fnames

def find_max_intensity_point(data, mlt_range=[8, 16],
                             mlat_range=[75, 85], plot_image=False):
    import cv2
    import scipy

    # Crop the image for 
    func_mlt = lambda x: 1 if (x>mlt_range[0] and x<mlt_range[1]) else 0
    func_mlat = lambda x: 1 if (x>mlat_range[0] and x<mlat_range[1]) else 0
    f_mlt = np.vectorize(func_mlt)
    f_mlat = np.vectorize(func_mlat)
    arr_mlt = f_mlt(data['mlt'])
    arr_mlat = f_mlat(data['mlat'])
    arr_mask = np.multiply(arr_mlt, arr_mlat)
    croped_image = np.multiply(data['image'], arr_mask)

    # find the max point
    #image = cv2.imread("../from_harald_frey/IMF12LSI_2002_0318_160223m.jpg")
    orig = scipy.misc.bytescale(data['image'], cmin=0)
    cropped = scipy.misc.bytescale(croped_image, cmin=0)

    # Blur the cropped image
    cropped = cv2.GaussianBlur(cropped, (5,5), 0)

    # find the max point
    (minval, maxval, minloc, maxloc) = cv2.minMaxLoc(cropped)

    if plot_image:
        cv2.circle(orig, maxloc, 9, (255, 0, 0), 2)
        #cv2.imshow("Naive", image)
        fig, ax = plt.subplots()
        #ax.imshow(image)
        cax = ax.imshow(orig)
        ax.set_axis_off()
        ax.set_title(dtm.strftime("%Y/%m/%d %H:%M"))
        cbar = fig.colorbar(cax)

    return (maxval, maxloc)

def iteratively_plot_WIC_SI_image(stm, etm, read_si=True, read_wic=False,
                                  param="image", mark_max_intensity_point=True,
                                  mlt_range=[8, 16], mlat_range=[70, 85]):

    dtms, fnames = find_filenames(stm, etm, read_si=read_si, read_wic=read_wic)
    for dtm in dtms:
        fig = plot_WIC_SI_image(dtm, read_si=read_si, read_wic=read_wic, 
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
        
def plot_WIC_SI_image(dtm, read_si=True, read_wic=False,
                      param="image", mark_max_intensity_point=True,
                      mlt_range=[8, 16], mlat_range=[70, 85]):

    import scipy
    import cv2

    fname, data = read_wic_si_data(dtm, read_si=read_si, read_wic=read_wic)
    ss_txt = ":" + fname.split("/")[-1].split("_")[-1][4:6]

    fig, ax = plt.subplots()
    if param == "image":
        vmin = 0; vmax = 1000; cmap = 'gist_gray'
    if param == "mlat":
        vmin = 60; cmap=None
    if param == "mlon":
        vmin=-180; vmax=180; cmap=None
    if param == "mlt":
        vmin=0; vmax=24; cmap=None

    #img = scipy.misc.bytescale(data[param], cmin=0)
    img = data[param]
    if mark_max_intensity_point:
        (maxval, maxloc) = find_max_intensity_point(data, mlt_range=mlt_range,
                                                    mlat_range=mlat_range)
        #cv2.circle(img, maxloc, 9, (255, 0, 0), 2)
    cax = ax.imshow(img, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.plot(maxloc[0], maxloc[1], 'or', markersize=4)
    ax.set_axis_off()
    if read_si:
        title = dtm.strftime("%Y/%m/%d %H:%M") + ss_txt + "   SI " + param
    if read_wic:
        title = dtm.strftime("%Y/%m/%d %H:%M") + ss_txt + "   WIC " + param
    ax.set_title(title)

    cbar = fig.colorbar(cax, shrink=0.7)
    return fig

# run the code
def main():
    #read_si=True; read_wic=False
    read_si=False; read_wic=True
    mlt_range = [9, 15]
    mlat_range = [70, 85]
    dtm = dt.datetime(2002, 3, 18, 16, 02) 
    stm = dt.datetime(2002, 3, 18, 15, 00) 
    etm = dt.datetime(2002, 3, 18, 17, 22) 
    fname, data = read_wic_si_data(dtm, read_si=read_si, read_wic=read_wic)
    #fig = plot_WIC_SI_image(dtm, read_si=read_si, read_wic=read_wic, param="image")
    #iteratively_plot_WIC_SI_image(stm, etm, read_si=read_si, read_wic=read_wic,
    #                              param="image", mark_max_intensity_point=True,
    #                              mlt_range=mlt_range, mlat_range=mlat_range)
    #find_max_intensity_point(data, mlt_range=mlt_range, mlat_range=mlat_range, plot_image=True)
    #fnames = find_filenames(stm, etm, read_si=True, read_wic=False)

    #plt.show()

if __name__ == "__main__":
    main()

