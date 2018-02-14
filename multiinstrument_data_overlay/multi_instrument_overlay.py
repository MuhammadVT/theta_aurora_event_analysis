import matplotlib
matplotlib.use("Agg")

def main():

    import datetime as dt
    import matplotlib.pyplot as plt
    from davitpy.utils import plotUtils
    import plotMapGrd
    from imagers.timed import timed_utils
    from overlay_IMAGE import IMAGE_sat
    from overlay_DMSP import DMSP_sat
    from overlay_TIMEDGUVI import TIMEDGUVI_sat
    import matplotlib.cm as cm
    import sys
    sys.path.append("../track_wic_si_spot/")
    from track_max_point import find_filenames

    overlay_TIMEDGUVI_data=False
    overlay_IMAGE_data=True
    IMAGE_datatype="WIC"
    interpolate_IMAGE_data=True
    overlay_SuperDARN_data=True
    overlay_MapFitVel=True
    overlay_CnvCntrs=False
    overlay_HMB=False
    overlay_DMSP_data=False
    DMSP_sat_nums=[13]
    #DMSP_sat_nums=[14, 15]
    coords = "mlt"
    vec_cmap = cm.jet

    # Plot iteratively
    read_wic=False; read_si=False
    if IMAGE_datatype == "WIC":
        read_wic=True
    if IMAGE_datatype == "SI":
        read_si=True
    
#################################################################
# Times when DMSP data is available for 03182002 event
#    # Time interval where DMSP F13 is available
#    stm = dt.datetime(2002,3,18,15,2)
#    etm = dt.datetime(2002,3,18,15,20)

#    # Time interval where DMSP F14 and F15 is available
#    stm = dt.datetime(2002,3,18,16,2)
#    etm = dt.datetime(2002,3,18,16,6)

#    # Time interval where DMSP F15 is available
#    stm = dt.datetime(2002,3,18,16,8)
#    etm = dt.datetime(2002,3,18,16,20)

#    # Time interval where DMSP F13 is available
#    stm = dt.datetime(2002,3,18,16,44)
#    etm = dt.datetime(2002,3,18,17,0)

#################################################################

    stm = dt.datetime(2002,3,18,16,20)
    etm = dt.datetime(2002,3,18,16,42)

    dtms, fnames = find_filenames(stm, etm, read_si=read_si, read_wic=read_wic)
    #dtms = [dt.datetime(2002,3,18,16,41)]
    for dtm in dtms:

        fig, ax = plt.subplots(figsize=(12,8))
        ax.set_facecolor('black')

        # Plot a map
        mobj = plotUtils.mapObj(ax=ax, datetime=dtm, lon_0=0., boundinglat=60.,
                                gridLabels=False, gridLatRes=10., coords=coords,
                                fillContinents='None')

        # Overlay IMAGE data
        if overlay_IMAGE_data:

            IMAGE_obj = IMAGE_sat(dtm, datatype=IMAGE_datatype)
            IMAGE_obj.overlay_data(mobj, param="image",
                                   interpolate_data=interpolate_IMAGE_data,
                                   overlay_spot_only=False,
                                   mark_max_intensity_point=True,
                                   max_point_color='b',
                                   max_point_size=30,
                                   gauss_filter=False,
                                   vmin=0, vmax=1000,
                                   alpha=1.0, 
                                   zorder=3,
                                   cmap='gist_gray')
            IMAGE_obj.add_cbar(fig, ax, shrink=0.6,
                               label="Rayleigh (for " + IMAGE_datatype + " data)",
                               label_fontsize=15)

            IMAGE_obj = IMAGE_sat(dtm, datatype="SI")
            IMAGE_obj.overlay_data(mobj, param="image",
                                   interpolate_data=interpolate_IMAGE_data,
                                   overlay_spot_only=True,
                                   mark_max_intensity_point=True,
                                   max_point_color='r',
                                   max_point_size=30,
                                   gauss_filter=False,
                                   vmin=0, vmax=1000,
                                   alpha=1.0, 
                                   zorder=6,
                                   cmap='brg')
#            IMAGE_obj.add_cbar(fig, ax, shrink=0.6,
#                               label="Rayleigh (for WIC data)",
#                               label_fontsize=15)


        # Overlay SuperDARN convection flows
        if overlay_SuperDARN_data:
            mapDatObj = plotMapGrd.MapConv(dtm, mobj, ax, vec_len_factor=3)
            if overlay_MapFitVel:
                mapDatObj.overlayMapFitVel(pltColBar=True, overlayRadNames=True,
                             annotateTime=True, colorBarLabelSize=15.0,
                             marker_size=10.0, alpha=0.7, zorder=5.0,
                             edgecolor='none', cbar_shrink=0.6,
                             colMap=vec_cmap, label_style="web")
            if overlay_CnvCntrs:
                mapDatObj.overlayCnvCntrs()
            if overlay_HMB:
                mapDatObj.overlayHMB()

        # Overlay DMSP
        if overlay_DMSP_data:
            for DMSP_sat_num in DMSP_sat_nums:
                DMSP_obj = DMSP_sat(dtm, DMSP_sat_num)
                DMSP_obj.overlay_data(mobj, ax,
                                      interval=2*60, velscl=3*111.,
                                      vec_cmap=vec_cmap,
                                      vel_scale=[0, 1000.],
                                      quality_flag=False)
        
        # Overlay TIMED-GUVI
        if overlay_TIMEDGUVI_data:
            currDate = dt.datetime(dtm.year, dtm.month, dtm.day)
            TIMEDGUVI_obj = TIMEDGUVI_sat(currDate) 
            TIMEDGUVI_obj.overlay_data(mobj, ax, inpTime=dtm,
                         timeDelta=40, vmin=0., vmax=3000.,
                         alpha=1.0, zorder=1, timeZorder=7.,
                         timeColor="white", timeTextColor="r",
                         timeMarkerSize=2., timeFontSize=8., 
                         plotCBar=True, autoScale=False,
                         plotType='d135', overlayTime=True,
                         overlayTimeInterval=5, timeMarker='o',
                         plotTitle=False, cbar_shrink=0.6, titleString=None,
                         coords=coords, timedguviCmap='gist_gray')

        fdir ="../plots/multiinstrument_data_overlay/" + IMAGE_datatype + "/"
        fname = "map_overlay_" + dtm.strftime("%H%M") +\
                "_with_" + IMAGE_datatype
                #"_with_TIMEDGUVI" 
        fig.savefig(fdir + fname, dpi=200)

if __name__ == "__main__": 
    main()
