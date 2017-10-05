import pandas
import datetime
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import bs4
import os
import requests

class TimedGuviDownload(object):
    """
    A class to Download SSUSI data
    given a date and datatype!
    """
    def __init__(self, outBaseDir = "../data/"):
        """
        Set up some constants.
        Used for fitting.
        """
        self.baseUrl = "http://guvitimed.jhuapl.edu/"
        self.outBaseDir = outBaseDir

    def download_files(self, inpDate):
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

        dwnldUrl = self.baseUrl +\
             "data_fetch_l1c_imaging_v013?y="+\
              str(inpDate.year) + "&d="+strDoY
        driver = webdriver.Chrome()
        driver.get(dwnldUrl)

        try:
            element = WebDriverWait(driver, 10).until(EC.element_to_be_clickable((By.ID, 'output')))
            filesDiv = driver.find_element_by_id("output")
            fileLinks = filesDiv.find_elements_by_css_selector('a')
            for uEl in fileLinks:
                fUrl = uEl.get_attribute('href')
                if "L1C-2-disk" not in fUrl:
                    continue
                print "currently downloading-->", fUrl
                rf = requests.get( fUrl, verify=False )
                currFName = rf.url.split("/")[-1]
                outDir = self.outBaseDir + inpDate.strftime( "%Y%m%d" ) + "/"
                if not os.path.exists(outDir):
                    os.makedirs(outDir)
                with open( outDir + currFName, "wb" ) as ssusiData:
                    ssusiData.write( rf.content )
        finally:
            driver.quit()

        