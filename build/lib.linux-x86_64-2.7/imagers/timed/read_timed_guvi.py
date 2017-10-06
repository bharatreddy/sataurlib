import netCDF4
from cdf.internal import EPOCHbreakdown
import pandas
import os
import numpy
import datetime
import matplotlib.pyplot as plt
import aacgmv2
import shutil

class ProcessTGData(object):
    """
    A class to Download SSUSI data
    given a date and datatype!
    """
    def __init__(self, inpDirs, outDir, inpDate):
        """
        Given a list of dirs (SSUSI has multiple files per date
        for the same satellite). Read all files from it.
        """
        # loop through the dir and get a list of the files
        self.fileList = []
        for currDir in inpDirs:
            for root, dirs, files in os.walk(currDir):
                for fNum, fName in enumerate(files):
                    self.fileList.append( root + fName )
        self.outDir = outDir
        self.inpDate = inpDate
        print self.fileList

    def processed_data_to_file(self, coords="geo", keepRawFiles=False):
        """
        read the required data into a dataframe
        select only required columns, if aacgm
        coordinates are selected convert geo to
        AACGM coords and save data to file!
        """
        # We'll delete raw File dirs at the end, keep alist of them
        delDirList = []
        for fileInd, currFile in enumerate(self.fileList):
            # if selFname not in currFile:
            #     continue
            # Get Sat name
            print "currently working with file-->", currFile
            print "processing--->", fileInd+1, "/", len(self.fileList), "files"
            if ( (".nc" not in currFile) & (".NC" not in currFile) ):
                print "Not a valid netcdf file!!"
                continue
            currDataSet = netCDF4.Dataset(currFile)
            # get datetime from epoch list
            dtList = numpy.array( [ datetime.datetime( *EPOCHbreakdown( e ) ) \
                        for e in currDataSet.variables["TIME_EPOCH_DAY"][:] ] )
            currDate = self.inpDate.strftime("%Y%m%d")
            # get peircepoints
            prpntLats = currDataSet.variables['PIERCEPOINT_DAY_LATITUDE'][:]
            prpntLons = currDataSet.variables['PIERCEPOINT_DAY_LONGITUDE'][:]
            prpntAlts = currDataSet.variables['PIERCEPOINT_DAY_ALTITUDE'][:]
            # GET DISK intensity data - waveband/color radiance data
            # 5 colors are - 121.6, 130.4, 135.6 nm and LBH short and LBH long
            dskInt121 = currDataSet.variables['DISK_INTENSITY_DAY'][:, :, 0]
            dskInt130 = currDataSet.variables['DISK_INTENSITY_DAY'][:, :, 1]
            dskInt135 = currDataSet.variables['DISK_INTENSITY_DAY'][:, :, 2]
            dskIntLBHS = currDataSet.variables['DISK_INTENSITY_DAY'][:, :, 3]
            dskIntLBHL = currDataSet.variables['DISK_INTENSITY_DAY'][:, :, 4]
            # We'll store the data in a DF. Now we need to be a little cautious
            # when storing the data in a DF. SSUSI measures flux as swaths, so
            # at each time instance we have multiple lats and lons and disk data.
            # I'm taking a simple approach where I take each lat (lon and other
            #  data) at a time instance as a column and time as rows. So if 
            # the array ishaving a dimention of 42x1632, each of the 42 
            # elements becomes a column and the 1632 time instances become rows.
            latColList = [ "glat." + str(cNum+1) for cNum in range(prpntLats.shape[0]) ]
            lonColList = [ "glon." + str(cNum+1) for cNum in range(prpntLats.shape[0]) ]
            d121ColList = [ "d121." + str(cNum+1) for cNum in range(prpntLats.shape[0]) ]
            d130ColList = [ "d130." + str(cNum+1) for cNum in range(prpntLats.shape[0]) ]
            d135ColList = [ "d135." + str(cNum+1) for cNum in range(prpntLats.shape[0]) ]
            dLBHSColList = [ "dlbhs." + str(cNum+1) for cNum in range(prpntLats.shape[0]) ]
            dLBHLColList = [ "dlbhl." + str(cNum+1) for cNum in range(prpntLats.shape[0]) ]
            # create dataframes with
            dfLat = pandas.DataFrame(prpntLats.T,columns=latColList, index=dtList)
            dfLon = pandas.DataFrame(prpntLons.T,columns=lonColList, index=dtList)
            dfD121 = pandas.DataFrame(dskInt121.T,columns=d121ColList, index=dtList)
            dfD130 = pandas.DataFrame(dskInt130.T,columns=d130ColList, index=dtList)
            dfD135 = pandas.DataFrame(dskInt135.T,columns=d135ColList, index=dtList)
            dfDLBHS = pandas.DataFrame(dskIntLBHS.T,columns=dLBHSColList, index=dtList)
            dfDLBHL = pandas.DataFrame(dskIntLBHL.T,columns=dLBHLColList, index=dtList)
            # Merge the dataframes
            timedDF = reduce(lambda left,right: pandas.merge(left,right,\
                         left_index=True, right_index=True), [ dfLat, \
                        dfLon, dfD121, dfD130, dfD135, dfDLBHL, dfDLBHS ])
            timedDF["orbitNum"] = currDataSet.variables['ORBIT_DAY'][:]
            # Lets also keep track of the sat name and shape of arrays
            timedDF["sat"] = "TIMED"
            timedDF["shapeArr"] = prpntLats.shape[0]
            # # reset index, we need datetime as a col
            timedDF = timedDF.reset_index()
            timedDF = timedDF.rename(columns = {'index':'date'})
            if coords != "geo":
                # Now we need to convert the GLAT, GLON into MLAT, MLON and MLT
                timedDF = timedDF.apply(self.convert_to_aacgm, axis=1)
                timedDF = timedDF.round(2)
                # We'll only need aacgm coords, discard all geo coords
                mlatColList = [ "mlat." + str(cNum+1) for cNum in range(prpntLats.shape[0]) ]
                mlonColList = [ "mlon." + str(cNum+1) for cNum in range(prpntLats.shape[0]) ]
                mltColList = [ "mlt." + str(cNum+1) for cNum in range(prpntLats.shape[0]) ]
                outCols = ["date", "sat", "orbitNum"] + mlatColList + mlonColList + mltColList + d121ColList + \
                            d130ColList + d135ColList + dLBHSColList + dLBHLColList
            else:
                outCols = ["date", "sat", "orbitNum", "shapeArr"] + latColList + lonColList + d121ColList + \
                            d130ColList + d135ColList + dLBHSColList + dLBHLColList
            timedDF = timedDF[ outCols ]
            # We now need to write the processed data to a file
            if not os.path.exists(self.outDir):
                os.makedirs(self.outDir)
            # if the file for the date exists append data
            # else create the file and write data!!!
            outFileName = self.outDir + "/" + currDate + ".txt"
            if not os.path.exists( outFileName ):
                # NOTE we only need header when writing data for the first time!
                with open(outFileName, 'w') as ftB:
                    timedDF.to_csv(ftB, header=True,\
                                      index=False, sep=' ' )
            else:
                # sometimes the file is already present!
                # from a previous run.
                # we'll simply append to existing data
                # in that case!. So delete the existing file
                # if that is the case.
                if fileInd == 0:
                    print "FILE Exists already! deleting and overwriting"
                    os.remove(outFileName)
                    with open(outFileName, 'w') as ftB:
                        timedDF.to_csv(ftB, header=True,\
                                          index=False, sep=' ' )
                else:
                    with open(outFileName, 'a') as ftB:
                        timedDF.to_csv(ftB, header=False,\
                                          index=False, sep=' ' )
            if not keepRawFiles:
                os.remove(currFile)
                currDelDir = "/".join(currFile.split("/")[:-1]) + "/"
                if currDelDir not in delDirList:
                    delDirList.append( currDelDir )
        if not keepRawFiles:
            for dd in delDirList:
                shutil.rmtree(dd)


    def convert_to_aacgm(self, row):
        """
        For the SSUSI DF convert all the 42
        Given glat, glon and date return
        mlat, mlon and mlt
        """
        for i in range( row["shapeArr"] ):
            indStr = str(i+1)
            mlat, mlon = aacgmv2.convert(row["glat." + indStr], row["glon." + indStr],\
                               300, row["date"])
            # mlon, mlat = utils.coord_conv( row["glon." + indStr], row["glat." + indStr], \
            #                      "geo", "mag", altitude=300., \
            #                      date_time=row["date"] )
            mlt = aacgmv2.convert_mlt(mlon, row["date"], m2a=False)
            row["mlat." + indStr] = numpy.round( mlat, 2)
            row["mlon." + indStr] = numpy.round( mlon, 2)
            row["mlt." + indStr] = numpy.round( mlt, 2)
        return row