import pandas
import os
import numpy
from scipy.interpolate import interp1d
import datetime
import matplotlib.pyplot as plt
import aacgmv2
from davitpy import utils
import matplotlib.pyplot as plt

class UtilsSsusi(object):
    """
    A class to Download SSUSI data
    given a date and datatype!
    """
    def __init__(self, inpDir, fileDate, satList=[ "F18", "F17", "F16" ]):
        """
        Given a input dir read data from the satellites.
        """
        # Read data from the files and store them in a
        # dict of dataframes
        self.frames = {}
        dirList = []
        self.fileDate = fileDate
        currDtStr = self.fileDate.strftime("%Y%m%d")
        # Loop through the satList and read data 
        # from each satellite.
        # NOTE we expect the sub directory in the
        # parent directory to be the name of the
        # satellite listed in satList
        for sat in satList:
            currFname = inpDir + sat + "/" + currDtStr + ".txt"
            # check if file exists
            if os.path.exists( currFname ):
                print "reading data from--->", currFname
                self.frames[ "ssusi" + sat ] = pandas.read_csv(\
                                 currFname, delim_whitespace=True,\
                                infer_datetime_format=True,\
                                parse_dates=["date"] )
            else:
                print "file not found-->", currFname

    def get_pole_times(self, hemi="north"):
        """
        From the files, get times during each orbit
        when the satellite was closest to the pole.
        """
        # Return a dict of poleTimes
        poleTimesDict = {}
        for key in self.frames.keys():
            ssusiDF = self.frames[key]
            if "mlat.1" in ssusiDF.columns:
                latCols = [col for col in ssusiDF if col.startswith('mlat')]
                     
            else:
                latCols = [col for col in ssusiDF if col.startswith('glat')]
            selCols = [ "orbitNum", "date" ] + latCols
            # filter out values where no lats greater than 85 are found
            cutOffLat = 85.
            if hemi == "north":
                poleLat = 90.
                evalStr = "(ssusiDF['{0}'] >" + str( int(cutOffLat) ) + ".)" #
            else:
                poleLat = -90.
                evalStr = "(ssusiDF['{0}'] <" + str( int(-1*cutOffLat) ) + ".)" #
            ssusiDF = ssusiDF[selCols][eval(" | ".join([\
                            evalStr.format(col) 
                            for col in latCols]))].reset_index(drop=True)
            # We'll use a simple and naive way to estimate swath closest to
            # the pole, from each glat/mlat col we'll subtract 90(or -90)
            # and sum up the differences and get row with min value.
            for l in latCols:
                ssusiDF[l] = abs(ssusiDF[l] - poleLat)
            ssusiDF['coLatSum'] = ssusiDF[latCols].sum(axis=1)
            # groupby orbit to get min colatsum
            poleTimesDF = ssusiDF[ ["orbitNum", "coLatSum"]\
                         ].groupby( "orbitNum" ).min().reset_index()
            poleTimesDF = pandas.merge( poleTimesDF, ssusiDF,\
                            on=["orbitNum", "coLatSum"] )
            # sometimes we get multiple rows for same orbit,
            # so we'll just take the first row
            poleTimesDF = poleTimesDF.groupby('orbitNum').first()
            poleTimesDict[key] = poleTimesDF["date"].values
        return poleTimesDict

    def filter_data_by_time(self, inpTime, hemi="north",\
             timeDelta=20, filterLat=40.):
        """
        Filter the processed data for 
        the desired time (+/- timeDel) of
        a given time and hemisphere.
        NOTE - timeDel is in minutes
        """
        # We'll output the results in a dict
        filteredDict = {}
        for key in self.frames.keys():
            ssusiDF = self.frames[key]
            ssusiDF = ssusiDF.fillna(0.)
            # get the time ranges that confine
            # inptime and timedelta
            timeMin = inpTime - datetime.timedelta(minutes=timeDelta)
            timeMax = inpTime + datetime.timedelta(minutes=timeDelta)
            # Choose DF rows which lie between timeMin and timeMax
            ssusiDF = ssusiDF[ (ssusiDF["date"] >= timeMin) &\
                                (ssusiDF["date"] <= timeMax) ]
            # select data based on hemi
            if hemi == "north":
                evalStr = "(ssusiDF['{0}'] >" + str( int(filterLat) ) + ".)" #
                # select all rows where lats are positive
                # we'll use the eval func for this purpose
                # First check if we have MLAT or GLAT coords
                # we have in the file!
                filterCol = [col for col in ssusiDF if col.startswith('mlat')]
                if len(filterCol) > 0:
                    ssusiDF = ssusiDF[eval(" & ".join([\
                            evalStr.format(col) 
                            for col in filterCol]))].reset_index(drop=True)
                else:
                    filterCol = [col for col in ssusiDF if col.startswith('glat')]
                    ssusiDF = ssusiDF[eval(" & ".join([\
                            evalStr.format(col) 
                            for col in filterCol]))].reset_index(drop=True)
            else:
                evalStr = "(ssusiDF['{0}'] <" + str( int(-1*filterLat) ) + ".)" #
                filterCol = [col for col in ssusiDF if col.startswith('mlat')]
                if len(filterCol) > 0:
                    ssusiDF = ssusiDF[eval(" & ".join([\
                            evalStr.format(col) 
                            for col in filterCol]))].reset_index(drop=True)
                else:
                    filterCol = [col for col in ssusiDF if col.startswith('glat')]
                    ssusiDF = ssusiDF[eval(" & ".join([\
                            evalStr.format(col) 
                            for col in filterCol]))].reset_index(drop=True)
            if ssusiDF.shape[0] == 0:
                print "********NO DATA FOUND, CHECK FOR A " +\
                         "DIFFERENT TIME OR INCREASE TIMEDEL********"
            # Sort the DF by time, since orbits mess it up
            ssusiDF = ssusiDF.sort_values('date')
            filteredDict[key] = ssusiDF
        return filteredDict


    def filter_data_by_orbit(self, inpTime, hemi="north", filterLat=0.):
        """
        Filter the processed data for 
        the desired closest orbit in time
        and hemisphere
        """
        # We'll output the results in a dict
        filteredDict = {}
        for key in self.frames.keys():
            ssusiDF = self.frames[key]
            ssusiDF = ssusiDF.fillna(0.)
            # get min and max times in each orbit
            orbitMin = ssusiDF[ ["date", "sat", "orbitNum"] \
                        ].groupby(["orbitNum", "sat"]).min().reset_index()
            orbitMin.columns = [ "orbitNum", "sat", "date_min" ]
            orbitMax = ssusiDF[ ["date", "sat", "orbitNum"] \
                        ].groupby(["orbitNum", "sat"]).max().reset_index()
            orbitMax.columns = [ "orbitNum", "sat", "date_max" ]
            orbitDF = pandas.merge( orbitMin, orbitMax,\
                         on=["orbitNum", "sat"] )
            selOrbit = orbitDF[ (orbitDF["date_min"] <= inpTime) &\
                 (orbitDF["date_max"] >= inpTime)\
                  ].reset_index(drop=True)
            # Only select the required orbit
            ssusiDF = ssusiDF.merge( selOrbit, on=[ "orbitNum", "sat" ] )
            # select data based on hemi
            if hemi == "north":
                evalStr = "(ssusiDF['{0}'] >" + str( int(filterLat) ) + ".)" #
                # select all rows where lats are positive
                # we'll use the eval func for this purpose
                filterCol = [col for col in ssusiDF if col.startswith('mlat')]
                if len(filterCol) > 0:
                    ssusiDF = ssusiDF[eval(" & ".join([\
                            evalStr.format(col) 
                            for col in filterCol]))].reset_index(drop=True)
                else:
                    filterCol = [col for col in ssusiDF if col.startswith('glat')]
                    ssusiDF = ssusiDF[eval(" & ".join([\
                            evalStr.format(col) 
                            for col in filterCol]))].reset_index(drop=True)
            else:
                evalStr = "(ssusiDF['{0}'] <" + str( int(-1*filterLat) ) + ".)" #
                filterCol = [col for col in df if col.startswith('mlat')]
                if len(filterCol) > 0:
                    ssusiDF = ssusiDF[eval(" & ".join([\
                            evalStr.format(col) 
                            for col in filterCol]))].reset_index(drop=True)
                else:
                    filterCol = [col for col in ssusiDF if col.startswith('glat')]
                    ssusiDF = ssusiDF[eval(" & ".join([\
                            evalStr.format(col) 
                            for col in filterCol]))].reset_index(drop=True)
            filteredDict[key] = ssusiDF
        return filteredDict

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
        
    def overlay_sat_data(self, filteredDict, mapHandle, ax,\
                        satList=["F18", "F17", "F16"], plotType='d135',\
                        overlayTime=True, overlayTimeInterval=5, timeLinestyle=':',\
                        timeColor="black", timeFontSize=8.,\
                         plotCBar=True, autoScale=True, vmin=0., vmax=1000.,\
                         plotTitle=True, titleString=None, inpTime=None,alpha=0.5,\
                         markSatName=True, coords="mag", ssusiCmap="Greens"):
        """
        Plot SSUSI data on a map
        # overlayTimeInterval is in minutes
        """
        # Loop through and read data
        for key in filteredDict.keys():
            ssusiDF = filteredDict[key]
            satNameKey= key[-3:]
            if satNameKey not in satList:
                continue
            # plot according to coords
            # also check if we have MLAT or GLAT coords
            # we have in the file!
            if coords != "geo":
                if "mlat.1" not in ssusiDF.columns:
                    print "converting from geo to aacgm coordinates"
                    # ssusiDF = self.convert_aacgm_geo(ssusiDF, a2g=False)    
                    ssusiDF = ssusiDF.apply(self.convert_aacgm_geo, args=(False,), axis=1)
                ssusiLats = ssusiDF\
                                [ssusiDF.columns[pandas.Series(\
                                ssusiDF.columns).str.startswith('mlat')\
                                ]].values
                if coords == "mag":
                    ssusiLons = ssusiDF\
                                    [ssusiDF.columns[pandas.Series(\
                                    ssusiDF.columns).str.startswith('mlon')\
                                    ]].values
                else:
                    ssusiLons = ssusiDF\
                                    [ssusiDF.columns[pandas.Series(\
                                    ssusiDF.columns).str.startswith('mlt')\
                                    ]].values * 15.
                ssusiDisk = ssusiDF\
                                [ssusiDF.columns[pandas.Series(\
                                ssusiDF.columns).str.startswith(plotType)\
                                ]].values
            else:
                if "glat.1" not in ssusiDF.columns:
                    print "converting from geo to aacgm coordinates"
                    # ssusiDisk = self.convert_aacgm_geo(ssusiDisk, a2g=True)   
                    ssusiDF = ssusiDF.apply(self.convert_aacgm_geo, args=(True,), axis=1) 
                ssusiLats = ssusiDF\
                                [ssusiDF.columns[pandas.Series(\
                                ssusiDF.columns).str.startswith('glat')\
                                ]].values
                ssusiLons = ssusiDF\
                                [ssusiDF.columns[pandas.Series(\
                                ssusiDF.columns).str.startswith('glon')\
                                ]].values
                ssusiDisk = ssusiDF\
                                [ssusiDF.columns[pandas.Series(\
                                ssusiDF.columns).str.startswith(plotType)\
                                ]].values
            # Need to get max and min for the plot
            # we'll round off to the nearest 500
            # and keep a cap of 1000.
            if autoScale:
                vmin = 0.
                vmax = numpy.round( numpy.max( ssusiDisk )/500. )*500.
            xVecs, yVecs = mapHandle(ssusiLons, ssusiLats, coords=coords)
            # ssusiPlot = mapHandle.scatter(xVecs, yVecs, c=ssusiDisk, s=10.,\
            #            cmap=ssusiCmap, alpha=0.7, zorder=5., \
            #                      edgecolor='none', marker="s",\
            #                       vmin=vmin, vmax=vmax)
            ssusiPlot = mapHandle.pcolormesh(xVecs, yVecs,\
                            ssusiDisk, zorder=5.,
                            vmin=0, vmax=vmax,
                            ax=ax, alpha=alpha, cmap=ssusiCmap)
            ssusiPlot.set_rasterized(True)
            # overlay time
            if overlayTime:
                uniqueTimeList = ssusiDF["date"].unique()
                timeDiff = ( uniqueTimeList.max() -\
                             uniqueTimeList.min()\
                              ).astype('timedelta64[m]')
                delRange = overlayTimeInterval*\
                            uniqueTimeList.shape[0]/timeDiff.astype('int')
                # get a list of times every timeinterval
                # for the given day
                nextDayTime = self.fileDate + datetime.timedelta(days=1)
                # for some early orbits, there might be values from
                # the previous day! so our date starts a day earlier
                # than the date of the file.
                currDt = self.fileDate - datetime.timedelta(days=1)
                allDayDatesList = []
                allDayTSList = []
                minDate = datetime.datetime.utcfromtimestamp(\
                                uniqueTimeList.min().tolist()/1e9)
                maxDate = datetime.datetime.utcfromtimestamp(\
                                uniqueTimeList.max().tolist()/1e9)
                while currDt <= nextDayTime:
                    if ( currDt >= minDate ) & ( currDt <= maxDate ):
                        ts = (numpy.datetime64(currDt) - \
                                numpy.datetime64('1970-01-01T00:00:00Z')\
                                ) / numpy.timedelta64(1, 's')
                        allDayTSList.append( ts )
                        allDayDatesList.append( currDt )
                    currDt += datetime.timedelta(minutes=overlayTimeInterval)
                if coords == "mag":
                    timessusiLats = ssusiDF\
                            [ssusiDF.columns[pandas.Series(\
                            ssusiDF.columns).str.startswith('mlat')\
                            ]].values
                    timessusiLons = ssusiDF\
                            [ssusiDF.columns[pandas.Series(\
                            ssusiDF.columns).str.startswith('mlon')\
                            ]].values
                elif coords == "mlt":
                    timessusiLats = ssusiDF\
                            [ssusiDF.columns[pandas.Series(\
                            ssusiDF.columns).str.startswith('mlat')\
                            ]].values
                    timessusiLons = ssusiDF\
                            [ssusiDF.columns[pandas.Series(\
                            ssusiDF.columns).str.startswith('mlt')\
                            ]].values * 15.
                else:
                    timessusiLats = ssusiDF\
                            [ssusiDF.columns[pandas.Series(\
                            ssusiDF.columns).str.startswith('glat')\
                            ]].values
                    timessusiLons = ssusiDF\
                            [ssusiDF.columns[pandas.Series(\
                            ssusiDF.columns).str.startswith('glon')\
                            ]].values
                timeSSusiTimes = ssusiDF["date"].values
                satTSArr = (timeSSusiTimes - \
                            numpy.datetime64('1970-01-01T00:00:00Z')\
                            ) / numpy.timedelta64(1, 's')
                # Interpolate the values to get times
                for dd in range( len(allDayDatesList) ):
                    for pixel in range(timessusiLons.shape[1]):
                        currPixelMlons = timessusiLons[:,pixel]
                        currPixelMlats = timessusiLats[:,pixel]
                        (x,y) = self.pol2cart( currPixelMlats, currPixelMlons )
                        # xArr = numpy.interp(allDayTSList[dd], satTSArr, x)
                        # yArr = numpy.interp(allDayTSList[dd], satTSArr, y)
                        fXIn = interp1d(satTSArr, x)
                        fYIn = interp1d(satTSArr, y)
                        xArr = fXIn(allDayTSList[dd])
                        yArr = fYIn(allDayTSList[dd])
                        (timePlotLatArr, timePlotLonArr) = self.cart2pol( xArr, yArr )
                        xTVecs, yTVecs = mapHandle(timePlotLonArr,\
                                         timePlotLatArr, coords=coords)
                        timeStr = " " + allDayDatesList[dd].strftime("%H:%M")
                        # Write Sat names used in plotting
                        if pixel == 0:
                            xTVecFirst = xTVecs
                            yTVecFirst = yTVecs
                        if pixel == timessusiLons.shape[1]-1:
                            xTVecLast = xTVecs
                            yTVecLast = yTVecs
                            if markSatName:
                                timeStr = timeStr + " (" + satNameKey + ")"
                            timeXVecs, timeYVecs = mapHandle(timePlotLonArr,\
                                 timePlotLatArr, coords=coords)
                            ax.text(timeXVecs, timeYVecs, timeStr,\
                                fontsize=timeFontSize,fontweight='bold',
                                ha='left',va='center',color='k',\
                                 clip_on=True, zorder=7.)
                    mapHandle.plot([xTVecFirst, xTVecLast], [yTVecFirst, yTVecLast],\
                                 timeColor, zorder=7., linestyle=timeLinestyle)
            # plot colorbar
            if plotCBar:
                cbar = plt.colorbar(ssusiPlot, orientation='vertical')
                cbar.set_label('Rayleighs', size=14)
            # Title
            if plotTitle:
                if titleString is not None:
                    plt.title(titleString)
                else:
                    if inpTime is not None:
                        inpTimeStr = inpTime.strftime("%Y-%m-%d  %H:%M")
                        plt.title( inpTimeStr + " UT" )
                    else:
                        print "***********NEED INPTIME FOR TITLE***********"
            

    def convert_aacgm_geo(self, row, a2g=False):
        """
        For the SSUSI DF convert all the 42
        Given glat, glon and date return
        mlat, mlon and mlt
        """
        for i in range( row["shapeArr"] ):
            indStr = str(i+1)
            if a2g:
                glat, glon = aacgmv2.convert(row["mlat." + indStr], row["mlon." + indStr],\
                               300, row["date"], a2g=a2g)
                row["glat." + indStr] = numpy.round( glat, 2)
                row["glon." + indStr] = numpy.round( glon, 2)
            else:
                mlat, mlon = aacgmv2.convert(row["glat." + indStr], row["glon." + indStr],\
                                   300, row["date"])
                mlt = aacgmv2.convert_mlt(mlon, row["date"], m2a=False)
                row["mlat." + indStr] = numpy.round( mlat, 2)
                row["mlon." + indStr] = numpy.round( mlon, 2)
                row["mlt." + indStr] = numpy.round( mlt, 2)
        return row