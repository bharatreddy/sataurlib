import pandas
import os
import numpy
from scipy.interpolate import interp1d
import datetime
import matplotlib.pyplot as plt
import aacgmv2
from davitpy import utils
import matplotlib.pyplot as plt

class UtilsTimedGuvi(object):
    """
    A class to Download TIMED GUVI data
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
        # Very similar to what we do for TIMED GUVI
        # except we just have one sat.
        currFname = inpDir + currDtStr + ".txt"
        # check if file exists
        if os.path.exists( currFname ):
            print "reading data from--->", currFname
            self.frames[ "timed" ] = pandas.read_csv(\
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
            timedDF = self.frames[key]
            if "mlat.1" in timedDF.columns:
                latCols = [col for col in timedDF if col.startswith('mlat')]
                     
            else:
                latCols = [col for col in timedDF if col.startswith('glat')]
            selCols = [ "orbitNum", "date" ] + latCols
            # filter out values where no lats greater than 85 are found
            cutOffLat = 60.
            if hemi == "north":
                poleLat = 90.
                evalStr = "(timedDF['{0}'] >" + str( int(cutOffLat) ) + ".)" #
            else:
                poleLat = -90.
                evalStr = "(timedDF['{0}'] <" + str( int(-1*cutOffLat) ) + ".)" #

            timedDF = timedDF[selCols][eval(" | ".join([\
                            evalStr.format(col) 
                            for col in latCols]))].reset_index(drop=True)
            # We'll use a simple and naive way to estimate swath closest to
            # the pole, from each glat/mlat col we'll subtract 90(or -90)
            # and sum up the differences and get row with min value.
            for l in latCols:
                timedDF[l] = abs(timedDF[l] - poleLat)
            timedDF['coLatSum'] = timedDF[latCols].sum(axis=1)
            # groupby orbit to get min colatsum
            poleTimesDF = timedDF[ ["orbitNum", "coLatSum"]\
                         ].groupby( "orbitNum" ).min().reset_index()
            poleTimesDF = pandas.merge( poleTimesDF, timedDF,\
                            on=["orbitNum", "coLatSum"], how="inner" )
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
            timedDF = self.frames[key]
            timedDF = timedDF.fillna(0.)
            # get the time ranges that confine
            # inptime and timedelta
            timeMin = inpTime - datetime.timedelta(minutes=timeDelta)
            timeMax = inpTime + datetime.timedelta(minutes=timeDelta)
            # Choose DF rows which lie between timeMin and timeMax
            timedDF = timedDF[ (timedDF["date"] >= timeMin) &\
                                (timedDF["date"] <= timeMax) ]
            # select data based on hemi
            if hemi == "north":
                evalStr = "(timedDF['{0}'] >" + str( int(filterLat) ) + ".)" #
                # select all rows where lats are positive
                # we'll use the eval func for this purpose
                # First check if we have MLAT or GLAT coords
                # we have in the file!
                filterCol = [col for col in timedDF if col.startswith('mlat')]
                if len(filterCol) > 0:
                    timedDF = timedDF[eval(" & ".join([\
                            evalStr.format(col) 
                            for col in filterCol]))].reset_index(drop=True)
                else:
                    filterCol = [col for col in timedDF if col.startswith('glat')]
                    timedDF = timedDF[eval(" & ".join([\
                            evalStr.format(col) 
                            for col in filterCol]))].reset_index(drop=True)
            else:
                evalStr = "(timedDF['{0}'] <" + str( int(-1*filterLat) ) + ".)" #
                filterCol = [col for col in timedDF if col.startswith('mlat')]
                if len(filterCol) > 0:
                    timedDF = timedDF[eval(" & ".join([\
                            evalStr.format(col) 
                            for col in filterCol]))].reset_index(drop=True)
                else:
                    filterCol = [col for col in timedDF if col.startswith('glat')]
                    timedDF = timedDF[eval(" & ".join([\
                            evalStr.format(col) 
                            for col in filterCol]))].reset_index(drop=True)
            if timedDF.shape[0] == 0:
                print "********NO DATA FOUND, CHECK FOR A " +\
                         "DIFFERENT TIME OR INCREASE TIMEDEL********"
            # Sort the DF by time, since orbits mess up plotting
            timedDF = timedDF.sort_values('date')
            filteredDict[key] = timedDF
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
            timedDF = self.frames[key]
            timedDF = timedDF.fillna(0.)
            # get min and max times in each orbit
            orbitMin = timedDF[ ["date", "sat", "orbitNum"] \
                        ].groupby(["orbitNum", "sat"]).min().reset_index()
            orbitMin.columns = [ "orbitNum", "sat", "date_min" ]
            orbitMax = timedDF[ ["date", "sat", "orbitNum"] \
                        ].groupby(["orbitNum", "sat"]).max().reset_index()
            orbitMax.columns = [ "orbitNum", "sat", "date_max" ]
            orbitDF = pandas.merge( orbitMin, orbitMax,\
                         on=["orbitNum", "sat"] )
            selOrbit = orbitDF[ (orbitDF["date_min"] <= inpTime) &\
                 (orbitDF["date_max"] >= inpTime)\
                  ].reset_index(drop=True)
            # Only select the required orbit
            timedDF = timedDF.merge( selOrbit, on=[ "orbitNum", "sat" ] )
            # select data based on hemi
            if hemi == "north":
                evalStr = "(timedDF['{0}'] >" + str( int(filterLat) ) + ".)" #
                # select all rows where lats are positive
                # we'll use the eval func for this purpose
                filterCol = [col for col in timedDF if col.startswith('mlat')]
                if len(filterCol) > 0:
                    timedDF = timedDF[eval(" & ".join([\
                            evalStr.format(col) 
                            for col in filterCol]))].reset_index(drop=True)
                else:
                    filterCol = [col for col in timedDF if col.startswith('glat')]
                    timedDF = timedDF[eval(" & ".join([\
                            evalStr.format(col) 
                            for col in filterCol]))].reset_index(drop=True)
            else:
                evalStr = "(timedDF['{0}'] <" + str( int(-1*filterLat) ) + ".)" #
                filterCol = [col for col in df if col.startswith('mlat')]
                if len(filterCol) > 0:
                    timedDF = timedDF[eval(" & ".join([\
                            evalStr.format(col) 
                            for col in filterCol]))].reset_index(drop=True)
                else:
                    filterCol = [col for col in timedDF if col.startswith('glat')]
                    timedDF = timedDF[eval(" & ".join([\
                            evalStr.format(col) 
                            for col in filterCol]))].reset_index(drop=True)
            filteredDict[key] = timedDF
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
                         plotType='d135', overlayTime=True,\
                         overlayTimeInterval=5, timeMarker='o',\
                         timeMarkerSize=2., timeColor="grey", timeTextColor="k",
                         zorder=5., timeZorder=7.,\
                         timeFontSize=8., plotCBar=True, autoScale=True,\
                         vmin=0., vmax=1000., plotTitle=True,\
                         titleString=None, inpTime=None,alpha=0.6,\
                         coords="mag", timedguviCmap="Greens"):
        """
        Plot TIMED GUVI data on a map
        # overlayTimeInterval is in minutes
        """
        # Loop through and read data
        for key in filteredDict.keys():
            timedDF = filteredDict[key]
            # plot according to coords
            # also check if we have MLAT or GLAT coords
            # we have in the file!
            if coords != "geo":
                if "mlat.1" not in timedDF.columns:
                    print "converting from geo to aacgm coordinates"
                    # timedDF = self.convert_aacgm_geo(timedDF, a2g=False)    
                    timedDF = timedDF.apply(self.convert_aacgm_geo, args=(False,), axis=1)
                timedLats = timedDF\
                                [timedDF.columns[pandas.Series(\
                                timedDF.columns).str.startswith('mlat')\
                                ]].values
                if coords == "mag":
                    timedLons = timedDF\
                                    [timedDF.columns[pandas.Series(\
                                    timedDF.columns).str.startswith('mlon')\
                                    ]].values
                else:
                    # Some changes for MLT plotting
                    # NOTE in davitpy MLT is in degrees 
                    # (between 0 and 360.). So we multiply
                    # MLT values by 15 to convert to degrees.
                    timedLons = timedDF\
                                    [timedDF.columns[pandas.Series(\
                                    timedDF.columns).str.startswith('mlt')\
                                    ]].values * 15.
                timedDisk = timedDF\
                                [timedDF.columns[pandas.Series(\
                                timedDF.columns).str.startswith(plotType)\
                                ]].values
            else:
                if "glat.1" not in timedDF.columns:
                    print "converting from geo to aacgm coordinates"
                    # timedDisk = self.convert_aacgm_geo(timedDisk, a2g=True)   
                    timedDF = timedDF.apply(self.convert_aacgm_geo, args=(True,), axis=1) 
                timedLats = timedDF\
                                [timedDF.columns[pandas.Series(\
                                timedDF.columns).str.startswith('glat')\
                                ]].values
                timedLons = timedDF\
                                [timedDF.columns[pandas.Series(\
                                timedDF.columns).str.startswith('glon')\
                                ]].values
                timedDisk = timedDF\
                                [timedDF.columns[pandas.Series(\
                                timedDF.columns).str.startswith(plotType)\
                                ]].values
            # Need to get max and min for the plot
            # we'll round off to the nearest 500
            # and keep a cap of 1000.
            if autoScale:
                vmin = 0.
                vmax = numpy.round( numpy.max( timedDisk )/500. )*500.
            xVecs, yVecs = mapHandle(timedLons, timedLats, coords=coords)
            # timedPlot = mapHandle.scatter(xVecs, yVecs, c=timedDisk, s=75.,\
            #            cmap=timedguviCmap, alpha=0.7, zorder=zorder, \
            #                      edgecolor='none', marker="s",\
            #                       vmin=vmin, vmax=vmax)
            timedPlot = mapHandle.pcolormesh(xVecs, yVecs,\
                            timedDisk, zorder=zorder,
                            vmin=0, vmax=vmax,
                            ax=ax, alpha=alpha, cmap=timedguviCmap)
            timedPlot.set_rasterized(True)
            # overlay time
            if overlayTime:
                uniqueTimeList = timedDF["date"].unique()
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
                    timetimedLats = timedDF\
                            [timedDF.columns[pandas.Series(\
                            timedDF.columns).str.startswith('mlat')\
                            ]].values
                    timetimedLons = timedDF\
                            [timedDF.columns[pandas.Series(\
                            timedDF.columns).str.startswith('mlon')\
                            ]].values
                elif coords == "mlt":
                    timetimedLats = timedDF\
                            [timedDF.columns[pandas.Series(\
                            timedDF.columns).str.startswith('mlat')\
                            ]].values
                    timetimedLons = timedDF\
                            [timedDF.columns[pandas.Series(\
                            timedDF.columns).str.startswith('mlt')\
                            ]].values * 15.
                else:
                    timetimedLats = timedDF\
                            [timedDF.columns[pandas.Series(\
                            timedDF.columns).str.startswith('glat')\
                            ]].values
                    timetimedLons = timedDF\
                            [timedDF.columns[pandas.Series(\
                            timedDF.columns).str.startswith('glon')\
                            ]].values
                timeTimedGuviTimes = timedDF["date"].values
                satTSArr = (timeTimedGuviTimes - \
                            numpy.datetime64('1970-01-01T00:00:00Z')\
                            ) / numpy.timedelta64(1, 's')
                # Interpolate the values to get times
                for dd in range( len(allDayDatesList) ):
                    for pixel in range(timetimedLons.shape[1]):
                        currPixelMlons = timetimedLons[:,pixel]
                        currPixelMlats = timetimedLats[:,pixel]
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
                        mapHandle.plot(xTVecs, yTVecs,\
                             marker=timeMarker,color=timeColor,\
                              markersize=timeMarkerSize, zorder=timeZorder)
                        timeStr = allDayDatesList[dd].strftime("%H:%M")
                        # Write Sat names used in plotting
                        if pixel == 0:
                            timeXVecs, timeYVecs = mapHandle(timePlotLonArr,\
                                 timePlotLatArr, coords=coords)
                            ax.text(timeXVecs, timeYVecs, timeStr,\
                                fontsize=timeFontSize,fontweight='bold',
                                ha='left',va='center',color=timeTextColor,\
                                 clip_on=True, zorder=timeZorder)
            # plot colorbar
            if plotCBar:
                cbar = plt.colorbar(timedPlot, orientation='vertical')
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
        For the TIMED GUVI DF convert all the 42
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
