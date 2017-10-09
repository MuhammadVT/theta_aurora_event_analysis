from track_max_point import find_filenames, read_wic_si_data
from track_max_point import find_max_intensity_point
from matplotlib.dates import DateFormatter
import matplotlib.pyplot as plt
import numpy as np        
import datetime as dt
import skimage
import scipy

def overlay_edges_on_WIC_SI_image(dtm, read_si=True, read_wic=False,
                                  param="mlt_img"):

    # read the image data
    fname, data = read_wic_si_data(dtm, read_si=read_si, read_wic=read_wic)
    orig_img = data[param]
    img = orig_img 

    # gaussina filter the image
    img = scipy.ndimage.gaussian_filter(img, 3)

    # find the edges
    edges = skimage.feature.canny(img, sigma=7)

    # blend the original image and the edge image
    new_img = np.multiply(orig_img, edges+0)
    new_img[np.where(new_img > 0)] = np.max(new_img)
    new_img = new_img + orig_img

    # plot the images
    fig, ax = plt.subplots()
    #cax = ax.imshow(new_img, cmap="gray", vmin=0, vmax=np.max(orig_img))
    cax = ax.imshow(new_img, cmap="gray", vmin=0, vmax=1000)
    #cax = ax.imshow(new_img, vmin=0)
    #cax = ax.imshow(edges + 0, cmap="gray")
    fig.colorbar(cax)

    if read_si:
        title = dtm.strftime("%Y/%m/%d %H:%M") + "   SI Data"
    if read_wic:
        title = dtm.strftime("%Y/%m/%d %H:%M") + "   WIC Data"
    ax.set_title(title)
    
    return fig
       
def iteratively_overlay_edges_on_WIC_SI_image(stm, etm, read_si=True, read_wic=False,
                                              param="mlt_img"):

    dtms, fnames = find_filenames(stm, etm, read_si=read_si, read_wic=read_wic)
    for dtm in dtms:
        fig = overlay_edges_on_WIC_SI_image(dtm, read_si=read_si, read_wic=read_wic,
                                            param=param)
        if read_si:
            fdir ="../plots/track_spot/SI_edges/"
            fname = dtm.strftime("%H%M") + "_SI_" + param
        if read_wic:
            fdir ="../plots/track_spot/WIC_edges/"
            fname = dtm.strftime("%H%M") + "_WIC_" + param
        fig.savefig(fdir + fname, dpi=200)

    return

# run the code
def main():
    read_si=True; read_wic=False
    #read_si=False; read_wic=True
    param="mlt_img"
    dtm = dt.datetime(2002, 3, 18, 16, 02)
    stm = dt.datetime(2002, 3, 18, 15, 00)
    etm = dt.datetime(2002, 3, 18, 17, 20)
#    fig = overlay_edges_on_WIC_SI_image(dtm, read_si=read_si, read_wic=read_wic,
#                                        param=param)
#    plt.show()
    iteratively_overlay_edges_on_WIC_SI_image(stm, etm, read_si=read_si,
                                              read_wic=read_wic, param=param)

if __name__ == "__main__":
    main()
