import os
import pandas
import datetime
import numpy
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.interpolate import interp1d
import aacgmv2
import seaborn as sns
import netCDF4

class PlotUtils(object):
    """
    A class to plot data from POES satellites.
    There are three plotting options available:
    1) given a time range and satellite number,
    plot the satellite pass
    1) given a time plot all nearby poes satellites 
    3) plot auroral boundary estimated using POES
    NOTE : We expect processed data to be present
    in the given folder!
    """
    def __init__(self, inpDate, pltCoords="mag"):
        # set up a few constants
        self.inpDate = inpDate
        if pltCoords not in ["mag", "mlt"]:
            print "Currently only mag, mlt coords are supported!"
            return None
        self.pltCoords = pltCoords


    def overlay_sat_pass( self, timeRange,  mapHandle, ax,\
                     rawSatDir, overlayElecFlux=True,\
                     satList=["m01", "m02", "n15", "n19", "n18"], markerSize=15,\
                     overlayTime=True, inpCmap=ListedColormap(sns.color_palette("BuPu")),\
                     timeMarkerSize=2., zorder=5., alpha=0.6, vmin=0, vmax=5,\
                      timeZorder=7., timeFontSize=6., plotCBar=True,overlayTimeInterval=5,\
                      autoScale=True, plotTitle=True,titleString=None ):
        # Given a timeRange and satellite list plot 
        # corresponding satellite passes
        # If overlayElecFlux is true! Electron flux is overlayed
        # else proton flux is overlayed!
        fileList = []
        fileCnt = 0
        for root, dirs, files in os.walk(rawSatDir):
            for fNum, fName in enumerate(files):
                currFile = root + fName    
                if ( (currFile.endswith(".nc")) &\
                        ("poes" in currFile.lower()) &\
                         (self.inpDate.strftime("%Y%m%d") in currFile) ):
                    fileList.append( currFile )
                    fileCnt += 1
                    if fileCnt >= 7:
                        print "got 7 satellites! skipping"
                        break
        # We now have the files
        # get only the required sats
        # This is more efficient instead of
        # checking earlier
        rawPoesFileList = []
        for csat in satList:
            for ff in fileList:
                if csat in ff:
                    rawPoesFileList.append( ff )
        ( poesAllEleDataDF, poesAllProDataDF ) = self.read_poes_data(rawPoesFileList)
        if overlayElecFlux:
            poesRawPlotDF = poesAllEleDataDF[ (poesAllEleDataDF["date"] >= timeRange[0]) &\
                                            (poesAllEleDataDF["date"] <= timeRange[1]) ]
            fluxVals = poesRawPlotDF["log_ele_flux"]
        else:
            poesRawPlotDF = poesAllProDataDF[ (poesAllProDataDF["date"] >= timeRange[0]) &\
                                            (poesAllProDataDF["date"] <= timeRange[1]) ]
            fluxVals = poesRawPlotDF["log_pro_flux"].values
        # get color scales
        if autoScale:
            vmin = 0.
            vmax = numpy.round( numpy.max( fluxVals )/5 )*5.
        # get required data for plotting
        poesLats = poesRawPlotDF["aacgm_lat_foot"].values
        if self.pltCoords == "mag":
            poesLons = poesRawPlotDF["aacgm_lon_foot"].values
            xVecs, yVecs = mapHandle(poesLons, poesLats, coords=self.pltCoords)
        else:
            poesLons = poesRawPlotDF["MLT"].values
            # NOTE in davitpy MLT ranges between 0 and 360.
            xVecs, yVecs = mapHandle(poesLons*15., poesLats, coords=self.pltCoords)
        poesPlot = mapHandle.scatter(xVecs, yVecs,\
                            c=fluxVals, s=markerSize, zorder=zorder,
                            vmin=vmin, vmax=vmax, ax=ax, alpha=alpha, cmap=inpCmap)
        # plot colorbar
        if plotCBar:
            cbar = plt.colorbar(poesPlot, orientation='vertical')
            cbar.set_label(r"Log Elec. Flux $ [ergs\ cm^{-2}\ s^{-2}]$", size=10)
        # overlay time
        if overlayTime:
            for ss in satList:
                uniqueTimeList = poesRawPlotDF[ poesRawPlotDF["sat"] == ss[1:] ]["date"].unique()
                timeDiff = ( uniqueTimeList.max() -\
                                 uniqueTimeList.min()\
                                  ).astype('timedelta64[m]')
                delRange = overlayTimeInterval*\
                            uniqueTimeList.shape[0]/timeDiff.astype('int')
                minDate = datetime.datetime.utcfromtimestamp(\
                                uniqueTimeList.min().tolist()/1e9)
                maxDate = datetime.datetime.utcfromtimestamp(\
                                uniqueTimeList.max().tolist()/1e9)
                # get lists of datetimes to plot 
                allDayDatesList = []
                allDayTSList = []
                currDt = minDate
                while currDt <= maxDate:
                    ts = (numpy.datetime64(currDt) - \
                                numpy.datetime64('1970-01-01T00:00:00Z')\
                                ) / numpy.timedelta64(1, 's')
                    allDayTSList.append( ts )
                    allDayDatesList.append( currDt )
                    currDt += datetime.timedelta(minutes=overlayTimeInterval)
                    if self.pltCoords == "mlt":
                        poesLats = poesRawPlotDF[ poesRawPlotDF["sat"] == ss[1:] ]["aacgm_lat_foot"].values
                        # Remember in davitpy mlt ranges between 0,360
                        poesLons = poesRawPlotDF[ poesRawPlotDF["sat"] == ss[1:] ]["MLT"].values * 15.
                    else:
                        poesLats = poesRawPlotDF[ poesRawPlotDF["sat"] == ss[1:] ]["aacgm_lat_foot"].values
                        poesLons = poesRawPlotDF[ poesRawPlotDF["sat"] == ss[1:] ]["aacgm_lon_foot"].values
                        
                    poesTimes = poesRawPlotDF[ poesRawPlotDF["sat"] == ss[1:] ]["date"].values
                    satTSArr = (poesTimes - \
                            numpy.datetime64('1970-01-01T00:00:00Z')\
                            ) / numpy.timedelta64(1, 's')
                    # Interpolate the values to get times
                    for dd in range( len(allDayDatesList) ):
                        (x,y) = self.pol2cart( poesLats, poesLons )
                        fXIn = interp1d(satTSArr, x)
                        fYIn = interp1d(satTSArr, y)
                        xArr = fXIn(allDayTSList[dd])
                        yArr = fYIn(allDayTSList[dd])
                        (timePlotLatArr, timePlotLonArr) = self.cart2pol( xArr, yArr )
                        xTVecs, yTVecs = mapHandle(timePlotLonArr,\
                                         timePlotLatArr, coords=self.pltCoords)
                        timeStr = allDayDatesList[dd].strftime("%H:%M") + "(" + str(ss) + ")"
                        timeXVecs, timeYVecs = mapHandle(timePlotLonArr,\
                                 timePlotLatArr, coords=self.pltCoords)
                        ax.text(timeXVecs, timeYVecs, timeStr,\
                                fontsize=timeFontSize,
                                ha='left',va='center',color='grey',\
                                 clip_on=True, zorder=timeZorder)

        # Title
        if plotTitle:
            if titleString is not None:
                plt.title(titleString)
            else:
                inpTimeStr = timeRange[0].strftime("%Y-%m-%d  %H:%M") +" -- " +\
                                timeRange[1].strftime("%Y-%m-%d  %H:%M")
                plt.title( inpTimeStr )


    def read_poes_data(self, fileList):
        # read data from given POES files.
        poesAllEleDataDF = pandas.DataFrame( columns =  ["timestamp", "date", "aacgm_lat_foot",\
                         "aacgm_lon_foot", "MLT", "log_ele_flux", "sat"] )
        poesAllProDataDF = pandas.DataFrame( columns =  ["timestamp", "date", "aacgm_lat_foot",\
                                 "aacgm_lon_foot", "MLT", "log_pro_flux", "sat"] )
        try:
            for f in fileList:
                # print "reading file-->", f
                # read variable from the netCDF files
                # Check size of file and if its not in Mbs skip
                if os.path.getsize(f) < 1096.:
                    continue
                poesRawData = netCDF4.Dataset(f)
                poesDF = pandas.DataFrame( poesRawData.variables['time'][:], columns=[ "timestamp" ] )
                poesDF['date'] = pandas.to_datetime(poesDF['timestamp'], unit='ms')
                poesDF["alt"] = poesRawData.variables['alt'][:]
                poesDF["aacgm_lat_foot"] = poesRawData.variables['aacgm_lat_foot'][:]

                poesDF["aacgm_lon_foot"] = poesRawData.variables['aacgm_lon_foot'][:]
                poesDF["MLT"] = poesRawData.variables['MLT'][:]
                # round of to 2 decimals
                poesDF['alt'] = [ round( x, 2 ) for x in poesDF['alt']]
                poesDF['aacgm_lat_foot'] = [ round( x, 2 ) for x in poesDF['aacgm_lat_foot']]
                poesDF['aacgm_lon_foot'] = [ round( x, 2 ) for x in poesDF['aacgm_lon_foot']]
                poesDF['MLT'] = [ round( x, 2 ) for x in poesDF['MLT']]
                # Add up the fluxes
                poesDF["ted_ele_total_flux"] = poesRawData.variables['ted_ele_tel0_flux_4'][:] +\
                        poesRawData.variables['ted_ele_tel0_flux_8'][:] + \
                        poesRawData.variables['ted_ele_tel0_flux_11'][:] + \
                        poesRawData.variables['ted_ele_tel0_flux_14'][:] + \
                        poesRawData.variables['ted_ele_tel30_flux_4'][:] +\
                        poesRawData.variables['ted_ele_tel30_flux_8'][:] + \
                        poesRawData.variables['ted_ele_tel30_flux_11'][:] + \
                        poesRawData.variables['ted_ele_tel30_flux_14'][:]
                poesDF["ted_pro_total_flux"] = poesRawData.variables['ted_pro_tel0_flux_4'][:] +\
                        poesRawData.variables['ted_pro_tel0_flux_8'][:] + \
                        poesRawData.variables['ted_pro_tel0_flux_11'][:] + \
                        poesRawData.variables['ted_pro_tel0_flux_14'][:] + \
                        poesRawData.variables['ted_pro_tel30_flux_4'][:] +\
                        poesRawData.variables['ted_pro_tel30_flux_8'][:] + \
                        poesRawData.variables['ted_pro_tel30_flux_11'][:] + \
                        poesRawData.variables['ted_pro_tel30_flux_14'][:]
                poesDF['log_ele_flux'] = [0. if x <= 0. else round( numpy.log10(x), 2 )\
                             for x in poesDF['ted_ele_total_flux']]
                poesDF['log_pro_flux'] = [0. if x <= 0. else round( numpy.log10(x), 2 )\
                             for x in poesDF['ted_pro_total_flux']]
                # the current satellite number
                poesDF["sat"] = f[-19:-17]
                # seperate out electron and proton flux and discard all zeros
                currPoesEleFluxDF = poesDF[poesDF["log_ele_flux"] > 0.][ ["timestamp",\
                                 "date", "aacgm_lat_foot", "aacgm_lon_foot", "MLT",\
                                 "log_ele_flux", "sat"] ].reset_index(drop=True)
                currPoesProFluxDF = poesDF[poesDF["log_pro_flux"] > 0.][ ["timestamp",\
                                 "date", "aacgm_lat_foot", "aacgm_lon_foot", "MLT",\
                                 "log_pro_flux", "sat"] ].reset_index(drop=True)
                poesAllEleDataDF = poesAllEleDataDF.append( currPoesEleFluxDF )
                poesAllProDataDF = poesAllProDataDF.append( currPoesProFluxDF )
                # now delete all the rows for prev DFs
                # we don't want to duplicate data
                poesDF = poesDF.drop( poesDF.index )
                currPoesEleFluxDF = currPoesEleFluxDF.drop( currPoesEleFluxDF.index )
                currPoesProFluxDF = currPoesProFluxDF.drop( currPoesProFluxDF.index )
            # create a date and time columns
            poesAllEleDataDF["dateStr"] = poesAllEleDataDF["date"].map(lambda x: x.strftime('%Y%m%d'))
            poesAllEleDataDF["time"] = poesAllEleDataDF["date"].map(lambda x: x.strftime('%H%M'))
            poesAllProDataDF["dateStr"] = poesAllProDataDF["date"].map(lambda x: x.strftime('%Y%m%d'))
            poesAllProDataDF["time"] = poesAllProDataDF["date"].map(lambda x: x.strftime('%H%M'))
            return ( poesAllEleDataDF, poesAllProDataDF )
        except:
            print "data read failed-->" + fileList
            return ( None, None )


    def cart2pol(self, x, y):
        """
        convert from cartesian to polar coords
        """
        colat = numpy.sqrt(x**2 + y**2)
        lat = 90. - colat
        lon = numpy.rad2deg( numpy.arctan2(y, x) )
        return (lat, lon)

    def pol2cart(self, lat, lon):
        """
        convert from polar to cartesian coords
        """
        colat = 90. - lat
        x = colat * numpy.cos(numpy.deg2rad(lon))
        y = colat * numpy.sin(numpy.deg2rad(lon))
        return (x, y)