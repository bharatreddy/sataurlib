import pandas
import datetime
import requests
import bs4
import os
import zlib

class SSUSIDownload(object):
    """
    A class to Download SSUSI data
    given a date and datatype!
    """
    def __init__(self, outBaseDir = "../data/"):
        """
        Set up some constants.
        Used for fitting.
        """
        self.baseUrl = "https://ssusi.jhuapl.edu/"
        self.outBaseDir = outBaseDir

    def download_files(self, inpDate, dataTypeList,\
             satList = [ "f18", "f17", "f16" ]):
        """
        Get a list of the urls from input date
        and datatype and download the files 
        and also move them to the corresponding
        folders.!!!
        """
        # construct day of year from date
        inpDoY = inpDate.timetuple().tm_yday
        strDoY = str(inpDoY)
        if inpDoY < 10:
            strDoY = "00" + str(inpDoY)
        if ( inpDoY > 10) & (inpDoY < 100):
            strDoY = "0" + str(inpDoY)
        for sat in satList:
            # construct url to get the list of files for the day
            for dataType in dataTypeList:
                payload = { "spc":sat, "type":dataType,\
                           "Doy":strDoY,\
                           "year":str(inpDate.year) }
                # get a list of the files from dmsp ssusi website
                # based on the data type and date
                r = requests.get(self.baseUrl + "data_retriver/",\
                                 params=payload, verify=False)
                soup = bs4.BeautifulSoup(r.text, 'html.parser')
                divFList = soup.find("div", {"id": "filelist"})
                hrefList = divFList.find_all(href=True)
                urlList = [ self.baseUrl + hL["href"] for hL in hrefList ]
                for fUrl in urlList:
                    # we only need data files which have .NC
                    if ".NC" not in fUrl:
                        continue
                    # If working with sdr data use only
                    # sdr-disk files
                    if dataType == "sdr":
                        if "SDR-DISK" not in fUrl:
                            continue

                    print "currently downloading-->", fUrl
                    rf = requests.get( fUrl, verify=False )
                    currFName = rf.url.split("/")[-1]
                    if not os.path.exists(self.outBaseDir + "/" + sat + "/"):
                        os.makedirs(self.outBaseDir + "/" + sat + "/")
                    outDir = self.outBaseDir + "/" + sat +\
                            "/" + inpDate.strftime( "%Y%m%d" ) + "/"
                    if not os.path.exists(outDir):
                        os.makedirs(outDir)
                    with open( outDir + currFName, "wb" ) as ssusiData:
                        ssusiData.write( rf.content )
