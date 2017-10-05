import os
from setuptools import setup, find_packages

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "sataurlib",
    version = "1.0.0b1",
    author = "Bharat Kunduri",
    author_email = "bharatr@vt.edu",
    description = ("Tools to download, process and plot data from auroral imagers. "
                                   "currently the library can work with DMSP SSUSI and TIMED GUVI satellites."),
    license = "MIT",
    keywords = "imagers, dmsp, timed guvi, aurora",
    url = "",
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
    ],
)