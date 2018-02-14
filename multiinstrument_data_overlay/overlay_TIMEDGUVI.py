class TIMEDGUVI_sat():

    def __init__(self, currDate, 
                 inpDir="../data/timed_guvi/processed/"): 

        from imagers.timed import timed_utils
        
        self.currDate = currDate
        self.inpDir = inpDir
        self.tgObj = timed_utils.UtilsTimedGuvi(inpDir, currDate)

        return 


    def download_data(self, endDate,
                      tempFileDir = "../data/timed_guvi/"):
        """ Download TIMED GUVI data """
        #import os
        #import datetime
        #from imagers.timed import dwnld_timed_guvi
        #
        #tgDwnldObj = dwnld_timed_guvi.TimedGuviDownload(\
        #                        outBaseDir = tempFileDir)
        #tDelta = datetime.timedelta(days=1)
        #currDate = self.currDate
        #while currDate <= endDate:
        #        print "currently downloading files for --> ",\
        #                        currDate.strftime("%Y-%m-%d")
        #        tgDwnldObj.download_files(currDate)
        #        currDate = currDate + tDelta
        pass
        
        return 

    def process_data(self):
        """ Processed TIMED GUVI data """
        #import os
        #import datetime
        #from imagers.ssusi import read_ssusi
        #from imagers.timed import read_timed_guvi
        #
        #rawFileDir = "../data/timed_guvi/" # Make sure you have this dir or create it
        #prcsdFileDir = "../data/timed_guvi/processed/" # Make sure you have this dir or create it
        #currDate = datetime.datetime( 2002, 3, 18 )
        #endDate = datetime.datetime( 2002, 3, 18 )
        #tDelta = datetime.timedelta(days=1)
        #while currDate <= endDate:
        #    for root, dirs, files in os.walk(rawFileDir):
        #        for nd, dd in enumerate(dirs):
        #            if currDate.strftime("%Y%m%d") not in dd:
        #                continue
        #            print "processing data --> ",\
        #                     currDate.strftime("%Y-%m-%d")
        #            tgRdObj = read_timed_guvi.ProcessTGData( [root + dd + "/"],\
        #                         prcsdFileDir, currDate )
        #            tgRdObj.processed_data_to_file()
        #    currDate += tDelta
        pass

        return 

    def overlay_data(self, mobj, ax, inpTime=None,
                     timeDelta=40, vmin=0., vmax=3000.,
                     alpha=1.0, zorder=1, timeZorder=7.,
                     timeColor="white", timeTextColor="r",
                     timeMarkerSize=2., timeFontSize=8., 
                     plotCBar=True, autoScale=True,
                     plotType='d135', overlayTime=True,
                     overlayTimeInterval=5, timeMarker='o',
                     plotTitle=True, cbar_shrink=0.7, titleString=None,
                     coords="mag", timedguviCmap='gist_gray'):

        """ overlay TIMED GUVI data """        
        self.fDict = self.tgObj.filter_data_by_time(inpTime,
                                                    timeDelta=timeDelta)
        self.tgObj.overlay_sat_data(self.fDict, mobj, ax,
                                    inpTime=inpTime,
                                    vmin=vmin, vmax=vmax,
                                    alpha=alpha, zorder=zorder,
                                    timeZorder=timeZorder,
                                    timeColor=timeColor, timeTextColor=timeTextColor,
                                    timeMarkerSize=timeMarkerSize,
                                    timeFontSize=timeFontSize, 
                                    plotCBar=plotCBar, autoScale=autoScale,
                                    plotType=plotType, overlayTime=overlayTime,
                                    overlayTimeInterval=overlayTimeInterval,
                                    timeMarker=timeMarker,
                                    plotTitle=plotTitle, cbar_shrink=cbar_shrink,
                                    titleString=titleString,
                                    coords=coords, timedguviCmap=timedguviCmap)
        return



