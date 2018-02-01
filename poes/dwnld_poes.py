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

class PoesDwnld(object):
    """
    A class to download poes data from noaa website.
    """
    def __init__(self, inpDate):
        # set up urls and dates
        self.homepage = "http://satdat.ngdc.noaa.gov/" +\
             "sem/poes/data/processed/ngdc/uncorrected/full/"
        self.inpDate = inpDate
        self.minCutoffFitLat = 45.
        self.delTimeCutOffNrstPass = 50 # min
        self.mlonDiffOtrEndCutoff = 50.
        self.delLatCutoff = 2.
        self.delCtimeCutoff = 60. #min
        # Roughly corresponds to 1 deg in MLAT
        self.gauss_smooth_sigma = 5. 
        self.diffElctrCutoffBnd = 0.1
        # More than an order of magnitude, remember its a log scale
        self.filtEleFluxCutoffMagn = 1.25 

    def get_all_sat_urls(self, dataFolder="./"):
        # ctx = ssl.create_default_context()
        # ctx.check_hostname = False
        # ctx.verify_mode = ssl.CERT_NONE
        # get a list of satellites avaiable for the date
        yearUrl = self.homepage + str( self.inpDate.year )
        try:
            conn = urllib.urlopen(yearUrl)
            htmlSource = conn.read()
            soup = bs4.BeautifulSoup(htmlSource, 'html.parser')
            # Get all the urls
            urlDict = {}
            for a in soup.find_all('a', href=True):
                if ( "metop" in a.contents[0] or "noaa" in a.contents[0] ):
                    urlDict[str(a['href'])] = yearUrl + "/" + a['href']
        except:
            print "data download from url failed-->" + yearUrl
            return None
        return urlDict

    def get_all_sat_data(self,outDir="./"):
        # generate urls to download POES satellite data from 
        # all files for a given date
        urlDict = self.get_all_sat_urls()
        if urlDict is None:
            print "url retreival failed!"
            return None
        if len(urlDict.keys()) == 0.:
            print "no urls/sats found!"
            return None
        try:
            # Get data from all the urls
            fileList = []
            for currSat in urlDict.keys():
                currFileName = "poes_" + currSat[0] + currSat[-3:-1] + "_"+\
                    self.inpDate.strftime("%Y%m%d") + "_proc.nc"
                print "downloading file from url-->" + \
                    urlDict[currSat] + currFileName
                self.get_file_from_url(urlDict[currSat], currFileName)
                # List of files to return
                fileList.append( outDir + "/" + currFileName )
                # Move the files to destination folder
                if outDir != "./":
                    # check if file exists and then transfer!
                    if os.path.isfile(outDir + "/" + currFileName):
                        print "file exists! check again..."
                    else:
                        print "moving file to destination folder", currFileName
                        print "outDir-->", outDir
                        shutil.move("./" + currFileName, outDir)
            return fileList
        except:
            print "download failed!!"
            return None

    def get_file_from_url(self, url, fileName):
        # Download a given poes file
        urllib.urlretrieve(url + fileName, fileName)