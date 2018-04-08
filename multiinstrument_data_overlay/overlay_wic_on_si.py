# Overlay IMAGE data

from overlay_IMAGE import IMAGE_sat
import datetime as dt
import matplotlib.pyplot as plt
from davitpy.utils import plotUtils
import sys
sys.path.append("../track_wic_si_spot/")
from track_max_point import find_filenames

stm = dt.datetime(2002,3,18,15,50)
etm = dt.datetime(2002,3,18,17,0)

vmin=200.; vmax=1000.
coords = "mlt"
interpolate_IMAGE_data=True
read_wic=True; read_si=False
wic_si_ratio = True
dtms, fnames = find_filenames(stm, etm, read_si=read_si, read_wic=read_wic)
#dtms = [dt.datetime(2002,3,18,16,41)]
for dtm in dtms:

    fig, ax = plt.subplots(figsize=(12,8))
    ax.set_facecolor('black')

    # Plot a map
    mobj = plotUtils.mapObj(ax=ax, datetime=dtm, lon_0=0., boundinglat=60.,
                            gridLabels=False, gridLatRes=10., coords=coords,
                            fillContinents='None')

    if wic_si_ratio:
        IMAGE_wic_si_ratio = IMAGE_sat(dtm, wic_si_ratio=True)
        IMAGE_wic_si_ratio.overlay_wic_si_ratio(mobj, param="image",
                                                gauss_filter=False,
                                                vmin=vmin, vmax=vmax,   # vmax is ignored here
                                                vmin_ratio=0., vmax_ratio=2.,
                                                alpha=1.0,
                                                zorder=3,
                                                cmap='bwr')
        IMAGE_wic_si_ratio.add_cbar(fig, ax, shrink=0.6,
                                    label="WIC/SI",
                                    label_fontsize=15)
        title = dtm.strftime("%m/%d/%Y    %H:%M") + "     IMAGE WIC/SI"
        fname = "wic_si_ratio_" + dtm.strftime("%H%M")
    else:
        IMAGE_wic = IMAGE_sat(dtm, datatype="WIC")
        IMAGE_wic.overlay_data(mobj, param="image",
                               interpolate_data=interpolate_IMAGE_data,
                               overlay_spot_only=False,
                               mark_max_intensity_point=True,
                               max_point_color='b',
                               max_point_size=30,
                               gauss_filter=False,
                               vmin=vmin, vmax=vmax,
                               alpha=0.5,
                               zorder=6,
                               cmap='gist_gray')
        IMAGE_wic.add_cbar(fig, ax, shrink=0.6,
                           label="Rayleigh (for WIC data)",
                           label_fontsize=15)

        IMAGE_si = IMAGE_sat(dtm, datatype="SI")
        IMAGE_si.overlay_data(mobj, param="image",
                               interpolate_data=interpolate_IMAGE_data,
                               overlay_spot_only=False,
                               mark_max_intensity_point=True,
                               max_point_color='r',
                               max_point_size=30,
                               gauss_filter=False,
                               vmin=vmin, vmax=vmax,
                               alpha=1.0,
                               zorder=3,
                               cmap='jet')
        IMAGE_si.add_cbar(fig, ax, shrink=0.6,
                          label="Rayleigh (for SI data)",
                          label_fontsize=15)

        title = dtm.strftime("%m/%d/%Y    %H:%M") + "     IMAGE WIC + SI"
        fname = "wic_si_overlayed_" + dtm.strftime("%H%M")

    ax.set_title(title, fontsize=15)

    fdir ="../plots/multiinstrument_data_overlay/WIC_SI/"
    fig.savefig(fdir + fname, dpi=200, bbox_inches="tight")

    #plt.show()
