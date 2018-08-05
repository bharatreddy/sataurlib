import urllib
import bs4
import ssl
import shutil
import os
import netCDF4
import pandas
import datetime
import numpy
import math
from scipy import signal, ndimage, optimize

class PoesAur(object):
    """
    A class to read data from all POES satellites
    and estimate the equatorward auroral boundary 
    location by fitting a circle.
    """
    def __init__(self):
        # set up a few constants
        self.minCutoffFitLat = 45.
        self.delTimeCutOffNrstPass = 50 # min
        self.mlonDiffOtrEndCutoff = 50.
        self.delLatCutoff = 2.
        self.delCtimeCutoff = 60. #min
        # Set some parameters for gaussian fitting!
        self.gauss_smooth_sigma = 1#2 
        self.diffElctrCutoffBnd = 0.1#0.15
        # More than an order of magnitude, remember its a log scale
        self.filtEleFluxCutoffMagn = 1.25 

    def read_poes_data_files(self, fileList=None, poesRawDate=None, poesRawDir=None):
        # read data from given POES files, process it for 
        # details such as the auroral boundary loc..
        # Note the input here can be a list of files or
        # date and poes raw data directory.
        # If both are given fileList is chosen
        if ( ( poesRawDate is None) & (poesRawDir is None) & (fileList is None) ):
            print "none of the input options set, use either\
                     fileList or poesRawDate & poesRawDir"
            return None
        if ( ( poesRawDate is None) & (poesRawDir is None) ):
            print "poesRawDate & poesRawDir not set! working with file list"
            if not isinstance(fileList, list):
                print "input fileList should be a list type"
                return None
            if len(fileList) == 0:
                print "fileList input is empty"
                return None
        if fileList is None:
            print "fileList not set! Working with poesRawDate & poesRawDir"
            # We'll loop through all the files in the directory and get the required files
            fileList = []
            fileCnt = 0 # we can get no more than 7 files
            if not isinstance(poesRawDate, datetime.datetime):
                print "poesRawDate should be a datetime.datetime type"
                return None 
            for root, dirs, files in os.walk(poesRawDir):
                for fNum, fName in enumerate(files):
                    currFile = root + fName    
                    if ( (currFile.endswith(".nc")) &\
                            ("poes" in currFile.lower()) &\
                             (poesRawDate.strftime("%Y%m%d") in currFile) ):
                        fileList.append( currFile )
                        fileCnt += 1
                        if fileCnt >= 7:
                            print "got 7 satellites! skipping"
                            break
        # We'll store all the data into two dataframes
        # one for electron flux and the other for protons
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
                poesDF = pandas.DataFrame( poesRawData.variables['time'][:].data, columns=[ "timestamp" ] )
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
            print "data read failed-->" + str(fileList)
            return ( None, None )

    def get_closest_sat_passes( self, poesAllEleDataDF, poesAllProDataDF, timeRange,\
         timeInterval=datetime.timedelta(minutes=30) ):
        # given a timeRange, timestep
        # get the closest 45 MLAT - 45 MLAT passes
        # for each of the satellites.
        try:
            outDFList = []
            ctime = timeRange[0]
            while ctime <= timeRange[1]:
                # We only need those times when POES was above self.minCutoffFitLat(45) MLAT
                poesAllEleDataDF = poesAllEleDataDF[ \
                                ( abs( poesAllEleDataDF["aacgm_lat_foot"] ) >= self.minCutoffFitLat )\
                                ].reset_index(drop=True)
                # We only need a few columns, discard the rest
                poesAllEleDataDF = poesAllEleDataDF[ [ 'sat', 'date',\
                                        'aacgm_lat_foot', 'aacgm_lon_foot',\
                                            'MLT', 'log_ele_flux' ] ]
                poesAllEleDataDF["delCtime"] = abs(poesAllEleDataDF["date"] - ctime)
                poesAllEleDataDF["delLatFit"] = abs( abs( poesAllEleDataDF["aacgm_lat_foot"] ) -\
                                                    abs( self.minCutoffFitLat ) )
                # We are sorting by sats, dates and lats to pick the nearest time
                # when the satellite is between two 45 MLATs
                poesAllEleDataDFNth = poesAllEleDataDF[ poesAllEleDataDF["aacgm_lat_foot"]\
                                        >= 0. ].sort_values( ['sat', 'date', 'aacgm_lat_foot'],\
                                                ascending=True ).reset_index(drop=True).drop_duplicates()
                poesAllEleDataDFSth = poesAllEleDataDF[ poesAllEleDataDF["aacgm_lat_foot"]\
                                        < 0. ].sort_values( ['sat', 'date', 'aacgm_lat_foot'],\
                                                ascending=True ).reset_index(drop=True).drop_duplicates()
                # Now we need to pick the satellite path
                # which is closest to the selected time.!
                # Northern Hemisphere
                poesAllEleDataDFNthST = poesAllEleDataDFNth[ \
                                            (poesAllEleDataDFNth["date"] \
                                                >= ctime-datetime.timedelta(\
                                                minutes=self.delTimeCutOffNrstPass)) &\
                                                (poesAllEleDataDFNth["date"] <=\
                                                 ctime+datetime.timedelta(\
                                                minutes=self.delTimeCutOffNrstPass))].reset_index(drop=True)
                poesAllEleDataDFNthST = poesAllEleDataDFNthST.sort_values(\
                                            ["sat","date"], ascending=[True, True]\
                                            ).reset_index(drop=True)
                # We'll get the the satellite pass which is moving polewards
                # Basically percent change in latitudes should be positive
                # for a satellite moving polewards (percent change would help
                # with the southern hemisphere lcoations.)
                poesAllEleDataDFNthST["latRowDiffs"] = poesAllEleDataDFNthST.groupby("sat")[[\
                                "aacgm_lat_foot" ] ].pct_change()
                poesAllEleDataDFNthST = poesAllEleDataDFNthST[\
                                        poesAllEleDataDFNthST["latRowDiffs"] > 0.\
                                        ].reset_index(drop=True)
                poesAllEleDataDFNthST = poesAllEleDataDFNthST.sort_values(\
                                            ["sat", "aacgm_lat_foot","delCtime"]\
                                            ).reset_index(drop=True)
                # get the start time
                selTimeRangeNthDF = poesAllEleDataDFNthST.groupby("sat").first().reset_index()
                # Now if the time difference is too large, discard the satellite data
                selTimeRangeNthDF = selTimeRangeNthDF[ selTimeRangeNthDF["delCtime"] <= \
                                    datetime.timedelta(minutes=self.delTimeCutOffNrstPass)\
                                    ].reset_index()
                selTimeRangeNthDF = selTimeRangeNthDF[ ["sat", "date"] ]
                selTimeRangeNthDF.columns = [ "sat", "start_time" ] 

                # Now get the end times, simply get all times that are
                # greater than start time, sort them by date and get
                # lowest deLatFit
                poesAllEleDataDFNthET = pandas.merge( poesAllEleDataDFNth,\
                                            selTimeRangeNthDF, on="sat" )
                poesAllEleDataDFNthET = poesAllEleDataDFNthET[ (\
                                                poesAllEleDataDFNthET["date"] >=\
                                                poesAllEleDataDFNthET["start_time"] ) &\
                                                (poesAllEleDataDFNthET["date"] <=\
                                                poesAllEleDataDFNthET["start_time"]+datetime.timedelta(\
                                                minutes=self.delTimeCutOffNrstPass)) ].reset_index(drop=True)
                poesAllEleDataDFNthET = poesAllEleDataDFNthET.sort_values(\
                                            ["sat","date"], ascending=[True, True]\
                                            ).reset_index(drop=True)
                # We'll get the the satellite pass which is moving equatorwards
                # Basically percent change in latitudes should be negative
                # for a satellite moving polewards (percent change would help
                # with the southern hemisphere lcoations.)
                poesAllEleDataDFNthET["latRowDiffs"] = poesAllEleDataDFNthET.groupby("sat")[[\
                                "aacgm_lat_foot" ] ].pct_change()
                poesAllEleDataDFNthET = poesAllEleDataDFNthET[\
                                        poesAllEleDataDFNthET["latRowDiffs"] < 0.\
                                        ].reset_index(drop=True)
                poesAllEleDataDFNthET = poesAllEleDataDFNthET.sort_values(\
                                            ["sat", "aacgm_lat_foot","delCtime"]\
                                            ).reset_index(drop=True)
                # get the start time
                eTimeNthDF = poesAllEleDataDFNthET.groupby("sat").first().reset_index()
                eTimeNthDF = eTimeNthDF[ ["sat", "date"] ]
                eTimeNthDF.columns = [ "sat", "end_time" ] 
                selTimeRangeNthDF = pandas.merge( selTimeRangeNthDF, eTimeNthDF, on="sat" )
                selTimeRangeNthDF["selTime"] = ctime
                # Now we need to pick the satellite path
                # which is closest to the selected time.!
                # Northern Hemisphere
                poesAllEleDataDFSthST = poesAllEleDataDFSth[ \
                                            (poesAllEleDataDFSth["date"] \
                                                >= ctime-datetime.timedelta(\
                                                minutes=self.delTimeCutOffNrstPass)) &\
                                                (poesAllEleDataDFSth["date"] <=\
                                                 ctime+datetime.timedelta(\
                                                minutes=self.delTimeCutOffNrstPass))].reset_index(drop=True)
                poesAllEleDataDFSthST = poesAllEleDataDFSthST.sort_values(\
                                            ["sat","date"], ascending=[True, True]\
                                            ).reset_index(drop=True)
                # We'll get the the satellite pass which is moving polewards
                # Basically percent change in latitudes should be positive
                # for a satellite moving polewards (percent change would help
                # with the southern hemisphere lcoations.)
                poesAllEleDataDFSthST["latRowDiffs"] = poesAllEleDataDFSthST.groupby("sat")[[\
                                "aacgm_lat_foot" ] ].pct_change()
                poesAllEleDataDFSthST = poesAllEleDataDFSthST[\
                                        poesAllEleDataDFSthST["latRowDiffs"] > 0.\
                                        ].reset_index(drop=True)
                poesAllEleDataDFSthST = poesAllEleDataDFSthST.sort_values(\
                                            ["sat", "aacgm_lat_foot","delCtime"],\
                                            ascending=[True, False, True]\
                                            ).reset_index(drop=True)
                # # get the start time
                selTimeRangeSthDF = poesAllEleDataDFSthST.groupby("sat").first().reset_index()
                # Now if the time difference is too large, discard the satellite data
                selTimeRangeSthDF = selTimeRangeSthDF[ selTimeRangeSthDF["delCtime"] <= \
                                    datetime.timedelta(minutes=self.delTimeCutOffNrstPass)\
                                    ].reset_index()
                selTimeRangeSthDF = selTimeRangeSthDF[ ["sat", "date"] ]
                selTimeRangeSthDF.columns = [ "sat", "start_time" ] 


                # # Now get the end times, simply get all times that are
                # # greater than start time, sort them by date and get
                # # lowest deLatFit
                poesAllEleDataDFSthET = pandas.merge( poesAllEleDataDFSth,\
                                            selTimeRangeSthDF, on="sat" )
                poesAllEleDataDFSthET = poesAllEleDataDFSthET[ (\
                                                poesAllEleDataDFSthET["date"] >=\
                                                poesAllEleDataDFSthET["start_time"] ) &\
                                                (poesAllEleDataDFSthET["date"] <=\
                                                poesAllEleDataDFSthET["start_time"]+datetime.timedelta(\
                                                minutes=self.delTimeCutOffNrstPass)) ].reset_index(drop=True)
                poesAllEleDataDFSthET = poesAllEleDataDFSthET.sort_values(\
                                            ["sat","date"], ascending=[True, True]\
                                            ).reset_index(drop=True)
                # We'll get the the satellite pass which is moving equatorwards
                # Basically percent change in latitudes should be negative
                # for a satellite moving polewards (percent change would help
                # with the southern hemisphere lcoations.)
                poesAllEleDataDFSthET["latRowDiffs"] = poesAllEleDataDFSthET.groupby("sat")[[\
                                "aacgm_lat_foot" ] ].pct_change()
                poesAllEleDataDFSthET = poesAllEleDataDFSthET[\
                                        poesAllEleDataDFSthET["latRowDiffs"] < 0.\
                                        ].reset_index(drop=True)
                poesAllEleDataDFSthET = poesAllEleDataDFSthET.sort_values(\
                                            ["sat", "aacgm_lat_foot","delCtime"],\
                                            ascending=[True, False, True]\
                                            ).reset_index(drop=True)
                # get the start time
                eTimeSthDF = poesAllEleDataDFSthET.groupby("sat").first().reset_index()
                eTimeSthDF = eTimeSthDF[ ["sat", "date"] ]
                eTimeSthDF.columns = [ "sat", "end_time" ] 
                selTimeRangeSthDF = pandas.merge( selTimeRangeSthDF, eTimeSthDF, on="sat" )
                selTimeRangeSthDF["selTime"] = ctime
                # Merge the two time range DFs to one
                currselTimeRangeDF = pandas.merge( selTimeRangeNthDF, selTimeRangeSthDF,\
                                 on=["sat", "selTime"], how="outer", suffixes=( '_nth', '_sth' ) )
                outDFList.append( currselTimeRangeDF )
                ctime += timeInterval
            # Concat all the DFs for differnt time ranges
            selTimeRangeDF = pandas.concat( outDFList )
            return selTimeRangeDF
        except:
            print "closest pass failed!!"
            return None

    def get_nth_ele_eq_bnd_locs( self, poesDataDF, poesAllEleDataDF,\
             remove_outliers=False, cutoff_iqr_prop=1.5 ):
        # given a dataframe, loop through times and
        # get the locations of auroral boundaries
        # for each of the satellites.
        try:
            aurEqBndList = []
            for currTime in poesDataDF["selTime"].unique():
                # For each unique time, get the pass times
                passTimeRange = poesDataDF[ \
                    poesDataDF["selTime"] == currTime ][ [\
                    "sat","start_time_nth","end_time_nth"] ].dropna()
                currPOESDF = pandas.merge( poesAllEleDataDF, passTimeRange, on="sat")
                # get data from poesSatellites
                currPOESDF = currPOESDF[ \
                             (currPOESDF["date"] >= currPOESDF["start_time_nth"]) &\
                             (currPOESDF["date"] <= currPOESDF["end_time_nth"])\
                             ].reset_index(drop=True)
                # NOTE we are only calculating values for northern hemi
                # This will have a lot of impact on future calcs
                currPOESDF = currPOESDF[ currPOESDF["aacgm_lat_foot"] > 0. ]\
                                            .reset_index(drop=True)
                # Divide satellite data to two passes
                # we'll get boundary data from each pass
                # In the first pass, sat is moving from 
                # low to high latitudes and in the second one
                # We'll get the opposite case
                currPOESDF = currPOESDF.sort_values( ["sat","date"],\
                                ascending=[True, True] \
                                ).reset_index(drop=True)
                # We'll get the the satellite pass which is moving polewards
                # Basically percent change in latitudes should be positive
                # for a satellite moving polewards.
                currPOESDF["latRowDiffs"] = currPOESDF.groupby("sat")[[\
                                "aacgm_lat_foot" ] ].pct_change()
                currPOESDFPolewards = currPOESDF[\
                                        currPOESDF["latRowDiffs"] > 0.\
                                        ].reset_index(drop=True)
                currPOESDFEquatorwards = currPOESDF[\
                                        currPOESDF["latRowDiffs"] < 0.\
                                        ].reset_index(drop=True)
                # Remove outliers if the option if chosen
                # we'll use the box-plot technique here.
                # Remove all those values which are more than
                # certain proportion of inter quartile range 
                # this proportion is set by "cutoff_iqr_prop"
                if remove_outliers:
                    # poleward pass
                    poleDFDesc = currPOESDFPolewards["log_ele_flux"].describe()
                    poleIQR = poleDFDesc["75%"] - poleDFDesc["25%"]
                    minOutlierCutoffPolePass = poleDFDesc["25%"] - cutoff_iqr_prop*poleIQR
                    maxOutlierCutoffPolePass = poleDFDesc["75%"] + cutoff_iqr_prop*poleIQR
                    currPOESDFPolewards = currPOESDFPolewards[\
                                 ( currPOESDFPolewards["log_ele_flux"] >=\
                                  minOutlierCutoffPolePass ) &\
                                   ( currPOESDFPolewards["log_ele_flux"] <=\
                                    maxOutlierCutoffPolePass ) ]
                    # equatorward pass
                    equDFDesc = currPOESDFEquatorwards["log_ele_flux"].describe()
                    equIQR = equDFDesc["75%"] - equDFDesc["25%"]
                    minOutlierCutoffEquPass = equDFDesc["25%"] - cutoff_iqr_prop*equIQR
                    maxOutlierCutoffEquPass = equDFDesc["75%"] + cutoff_iqr_prop*equIQR
                    currPOESDFEquatorwards = currPOESDFEquatorwards[\
                                ( currPOESDFEquatorwards["log_ele_flux"] >=\
                                minOutlierCutoffEquPass ) &\
                                ( currPOESDFEquatorwards["log_ele_flux"] <=\
                                 maxOutlierCutoffEquPass ) ]

                # Apply Gaussian filter to smooth the data
                currPOESDFPolewards["filtEleFluxPoleArr"] = ndimage.filters.gaussian_filter1d(\
                                        currPOESDFPolewards["log_ele_flux"],self.gauss_smooth_sigma) 
                currPOESDFEquatorwards["filtEleFluxEquatorArr"] = \
                                    ndimage.filters.gaussian_filter1d(\
                                    currPOESDFEquatorwards["log_ele_flux"],self.gauss_smooth_sigma)
                # Now we'll use simple method to calculate the boudaries
                # We'll calculate the max, min and std values of electron 
                # flux and estimate the cutoff as the closest (in time) point
                # where e flux is with in (min value + 1 std. dev) to the max value.
                # There are 5 sats we'll use groupby
                polePassMinFlux = currPOESDFPolewards[["filtEleFluxPoleArr",\
                                     "sat"]].groupby("sat").min().reset_index()
                polePassMinFlux.columns = [ "sat", "min_flux" ]
                polePassMaxFlux = currPOESDFPolewards[["filtEleFluxPoleArr",\
                                     "sat"]].groupby("sat").max().reset_index()
                polePassMaxFlux.columns = [ "sat", "max_flux" ]
                polePassStdFlux = currPOESDFPolewards[["filtEleFluxPoleArr",\
                                     "sat"]].groupby("sat").std().reset_index()
                polePassStdFlux.columns = [ "sat", "std_flux" ]
                # equ
                equPassMinFlux = currPOESDFEquatorwards[["filtEleFluxEquatorArr",\
                                     "sat"]].groupby("sat").min().reset_index()
                equPassMinFlux.columns = [ "sat", "min_flux" ]
                equPassMaxFlux = currPOESDFEquatorwards[["filtEleFluxEquatorArr",\
                                         "sat"]].groupby("sat").max().reset_index()
                equPassMaxFlux.columns = [ "sat", "max_flux" ]
                equPassStdFlux = currPOESDFEquatorwards[["filtEleFluxEquatorArr",\
                                         "sat"]].groupby("sat").std().reset_index()
                equPassStdFlux.columns = [ "sat", "std_flux" ]
                polePassStatSummDF = pandas.merge( polePassMinFlux,\
                                                     polePassMaxFlux, on="sat" )
                polePassStatSummDF = pandas.merge( polePassStatSummDF,\
                                                     polePassStdFlux, on="sat" )
                equPassStatSummDF = pandas.merge( equPassMinFlux,\
                                                     equPassMaxFlux, on="sat" )
                equPassStatSummDF = pandas.merge( equPassStatSummDF,\
                                                     equPassStdFlux, on="sat" )
                # get the time instance where the flux is max
                # for each satellite
                maxDTPolInds = currPOESDFPolewards.groupby(['sat'])\
                                    ['filtEleFluxPoleArr'].transform(max) ==\
                                     currPOESDFPolewards['filtEleFluxPoleArr']
                maxDTPolDF = currPOESDFPolewards[ maxDTPolInds ][["date", "sat"]]
                maxDTPolDF.columns = [ "max_date", "sat" ]
                # equatorward pass
                maxDTEquInds = currPOESDFEquatorwards.groupby(['sat'])\
                                    ['filtEleFluxEquatorArr'].transform(max) ==\
                                     currPOESDFEquatorwards['filtEleFluxEquatorArr']
                maxDTEquDF = currPOESDFEquatorwards[ maxDTEquInds ][["date", "sat"]]
                maxDTEquDF.columns = [ "max_date", "sat" ]
                polePassStatSummDF["cut_off"] = polePassStatSummDF["min_flux"] +\
                                                polePassStatSummDF["std_flux"]
                equPassStatSummDF["cut_off"] = equPassStatSummDF["min_flux"] +\
                                                equPassStatSummDF["std_flux"]
                # Now get the time instances of boundary
                cutOffPoleDF = pandas.merge( currPOESDFPolewards,\
                                                 maxDTPolDF, on="sat" )
                cutOffPoleDF = pandas.merge( cutOffPoleDF,\
                                         polePassStatSummDF, on="sat" )
                cutOffPoleDF = cutOffPoleDF[ ( cutOffPoleDF["date"] <=\
                                     cutOffPoleDF["max_date"] ) &\
                                   ( cutOffPoleDF["filtEleFluxPoleArr"] <=\
                                         cutOffPoleDF["cut_off"] )
                                   ]
                cutOffPoleDF = cutOffPoleDF[["sat", "date"]].groupby(\
                                "sat").max().reset_index()
                # equatorward pass
                cutOffEquDF = pandas.merge( currPOESDFEquatorwards,\
                                                 maxDTEquDF, on="sat" )
                cutOffEquDF = pandas.merge( cutOffEquDF,\
                                         equPassStatSummDF, on="sat" )

                cutOffEquDF = cutOffEquDF[ ( cutOffEquDF["date"] >=\
                                     cutOffEquDF["max_date"] ) &\
                                   ( cutOffEquDF["filtEleFluxEquatorArr"] <=\
                                         cutOffEquDF["cut_off"] )
                                   ]
                cutOffEquDF = cutOffEquDF[["sat", "date"]].groupby(\
                                "sat").min().reset_index()
                # Now get sat, lat, mlon, mlt of the bounadry
                polePassEqBndDF = pandas.merge( currPOESDFPolewards,\
                             cutOffPoleDF, on=["sat", "date"] )
                polePassEqBndDF = polePassEqBndDF[ ["sat", "aacgm_lat_foot",\
                                    "aacgm_lon_foot", "MLT"] ]
                polePassEqBndDF.columns = [ "sat", "pole_mlat",\
                                             "pole_mlon", "pole_mlt"  ]
                equatorPassEqBndDF = pandas.merge( currPOESDFEquatorwards,\
                             cutOffEquDF, on=["sat", "date"] )
                equatorPassEqBndDF = equatorPassEqBndDF[ ["sat", "aacgm_lat_foot",\
                                    "aacgm_lon_foot", "MLT"] ]
                equatorPassEqBndDF.columns = [ "sat", "equator_mlat",\
                                            "equator_mlon", "equator_mlt"  ]
                currAurEqBndDF = pandas.merge( polePassEqBndDF, equatorPassEqBndDF,\
                                 on="sat", how="outer" )
                currAurEqBndDF["time"] = currTime
                aurEqBndList.append( currAurEqBndDF )
            aurEqBndDF = pandas.concat( aurEqBndList )
            return aurEqBndDF
        except:
            print "closest pass failed!!"
            return None

    def fit_circle_aurbnd( self, bndLocDF, save_to_file=True,\
             fileFormat="txt",outDir="./", save_fit_coeff=False):
        # Given the boundary locations obtained
        # from different satellites, estimate the
        # auroral oval boundary by fitting a circle!
        # make a list of DFs to return at the end
        fitDFList = []
        firstWrite = True
        # check if fileformat to save are in the right format!
        if save_to_file:
            if fileFormat not in ["txt", "csv"]:
                print "only txt and csv file\
                         formats allowed, TRY AGAIN!"
                return None
        for currTime in bndLocDF["time"].unique():
            # try:
            currBndDF = bndLocDF[ bndLocDF["time"] == currTime ]
            if currBndDF.shape[0] <= 3:
                continue
            # Convert to numpy arrays 
            poleMlatArr = currBndDF["pole_mlat"].values
            poleMlonArr = currBndDF["pole_mlon"].values
            poleMltArr = currBndDF["pole_mlt"].values
            equMlatArr = currBndDF["equator_mlat"].values
            equMlonArr = currBndDF["equator_mlon"].values
            equMltArr = currBndDF["equator_mlt"].values
            # discard nan values
            poleMlatArr = poleMlatArr[~numpy.isnan(poleMlatArr)]
            poleMlonArr = poleMlonArr[~numpy.isnan(poleMlonArr)]
            poleMltArr = poleMltArr[~numpy.isnan(poleMltArr)]
            equMlatArr = equMlatArr[~numpy.isnan(equMlatArr)]
            equMlonArr = equMlonArr[~numpy.isnan(equMlonArr)]
            equMltArr = equMltArr[~numpy.isnan(equMltArr)]
            # Concat the arrays together
            latPoesAll = numpy.append( poleMlatArr, equMlatArr )
            lonPoesAll = numpy.append( poleMlonArr, equMlonArr )
            # Drop na's again
            lonPoesAll = lonPoesAll[~numpy.isnan(lonPoesAll)]
            latPoesAll = latPoesAll[~numpy.isnan(lonPoesAll)]
            # Now we do the fitting part...
            # Target function
            fitfunc = lambda p, x: p[0] + \
                        p[1]*numpy.cos(\
                        2*math.pi*(x/360.)+p[2]) 
            # Distance to the target function
            errfunc = lambda p, x,\
                         y: fitfunc(p, x) - y 
            # get the fitting results
            # Initial guess
            p0Equ = [ 1., 1., 1.]
            p1Equ, successEqu = optimize.leastsq(errfunc,\
                         p0Equ[:], args=(lonPoesAll, latPoesAll))
            eqPlotLons = numpy.linspace(0., 360., 25.)
            eqPlotLons[-1] = 0.
            eqBndLocs = []
            for xx in eqPlotLons :
                currLatEst = p1Equ[0] +\
                        p1Equ[1]*numpy.cos(2*math.pi*(xx/360.)+p1Equ[2] )
                eqBndLocs.append( ( round(currLatEst,1), xx ) )
            # Convert to DF
            if save_fit_coeff:
                # Save the fitted coefficients of the Aur Bnd Shape
                aurFitDF = pandas.DataFrame( [tuple(np.rount(p1Equ, 3))],
                            columns=["p_0", "p_1", "p_2"] )
                aurFitDF["sat"] = [currBndDF["sat"].values]
            else:
                # Save the fitted Aur Bnd points
                aurFitDF = pandas.DataFrame( eqBndLocs, \
                            columns=["MLAT", "MLON"] )
            cnvrtTime = pandas.to_datetime(str(currTime)) 
            aurFitDF["date"] = cnvrtTime.strftime( "%Y%m%d" )
            aurFitDF["time"] = cnvrtTime.strftime( "%H%M" )
            if save_to_file:
                outFitResFil = outDir + "poes-fit-" +\
                        cnvrtTime.strftime( "%Y%m%d" ) + "." + fileFormat
                if firstWrite:
                    with open(outFitResFil, 'w') as fra:
                        aurFitDF.to_csv(fra, header=True,\
                                          index=False, sep=' ' )
                        print "saving to file--->", outFitResFil
                    firstWrite = False
                else:
                    with open(outFitResFil, 'a') as fra:
                        aurFitDF.to_csv(fra, header=False,\
                                          index=False, sep=' ' )                
                fitDFList.append( aurFitDF )
            # except:
            #     print "couldnt get a fit! skipping!"
            #     continue
        if len(fitDFList) > 0:
            return pandas.concat(fitDFList)
