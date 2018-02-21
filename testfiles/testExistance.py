#!/usr/local/bin python2

# @author Jake Runyan
# @notes
#   Intended for use with Python 2.7

#**** IMPORTS ****#
import rpy2
from rpy2 import robjects #how r objects are defined in python, functions callable this way.
from rpy2.rinterface import R_VERSION_BUILD #simply for printing version info
import sys #For python system things, namely exit or raising exceptions

R_FILE_NAME = "HELLOWORLD.R"

print('1) R Global Environment holds:\n')
print(robjects.r.ls(robjects.globalenv))

rfilepath = "/Users/runyanjake/Desktop/Academia/Research/Brookslab Folder/R Paper Things/Development/HELLOWORLD.R"
try:
    robjects.r['source'](rfilepath) #runs but also adds to r env? how to add to r env alone?
except (rpy2.rinterface.RRuntimeError, rpy2.rinterface.RRuntimeWarning) as err:
    exitmsg = '\n\nA RRuntimeWarning/Error was caught.\nPlease verify that the file at ' + rfilepath + ' exists.\nExiting...'
    sys.exit(exitmsg)

print('2) R Global Environment holds:\n')
print(robjects.r.ls(robjects.globalenv))