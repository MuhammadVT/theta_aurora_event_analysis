def read_wic_si_data(dtm, read_si=True, read_wic=False):
    """
        This procedure reads IMAGE SI and WIC image data
    """
    from scipy.io import readsav
    import glob
    import numpy as np

    if read_si:
        fname_image = "./from_harald_frey/IMF12LSI_2002_0318_" + dtm.strftime("%H%M") + "*.sav"
    if read_wic:
        fname_image = "./from_harald_frey/IMFHWIC_2002_0318_" + dtm.strftime("%H%M") + "*.sav"
    fname_image = glob.glob(fname_image)[0]
    dat = readsav(fname_image, python_dict=True, verbose=False)
    data_dict = {}

    # Note: the data is originally saved in IDL
    # and IDL considers the first dimension to be the column.
    # Therefore, we do the transpose to fix it. Also we need to
    # rotate the image so that the dayside is on the top

    data_dict["img"] = np.rot90(np.transpose(dat['imageinfo'].mlt_img[0]))
    data_dict["mlat"] = np.rot90(np.transpose(dat['imageinfo'].mlat[0]))
    data_dict["mlon"] = np.rot90(np.transpose(dat['imageinfo'].mlon[0]))
    data_dict["mlt"] = np.rot90(np.transpose(dat['imageinfo'].mlt[0]))

    return fname_image, data_dict

def find_max_intensity_point():
    pass

    return

# run the code
import datetime as dt
import matplotlib.pyplot as plt
from skimage import feature
import numpy as np

dtm = dt.datetime(2002, 3, 18, 16, 02) 
fname, data = read_wic_si_data(dtm, read_si=True, read_wic=False)
ss_txt = ":" + fname.split("/")[-1].split("_")[-1][4:6]

fig, ax = plt.subplots()
#ax.imshow(np.rot90(data['img']))
ax.imshow(data['img'])
ax.set_axis_off()
ax.set_title(dtm.strftime("%Y/%m/%d %H:%M") + ss_txt)

#fig, ax = plt.subplots()
#edge = feature.canny(data['img'], sigma=5)
#ax.imshow(edge)
#ax.set_axis_off()

plt.show()


