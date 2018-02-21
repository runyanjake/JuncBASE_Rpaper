#!/usr/local/bin python2

# @author Jake Runyan
# @notes
#   Intended for use with Python 2.7

#**** IMPORTS ****#
import rpy2
from rpy2 import robjects #how r objects are defined in python, functions callable this way.
from rpy2.rinterface import R_VERSION_BUILD #simply for printing version info

#Print Version Information
print('\n**************VERSION INFORMATION**************')
print('rpy2 version: ' + rpy2.__version__)
print('This version of rpy2 was built on R version ')
print(R_VERSION_BUILD)
print('**************END VERSION INFORM***************\n')

#execution of .R file from rpy2 using R's !!! "source(<filename>)" !!! function
#NOTE: this assumes that the .R file can be runnable on its own. 
#      Other methods may be utilized to define it in R namespace (see test.py) and then call.abs
#      ALSO: .R file must be specified by direct path or 
rfile = "HELLOWORLD.R"
robjects.r['source'](rfile)

#For dealing with the R global environment, use http://rpy.sourceforge.net/rpy2/doc-dev/html/introduction.html 
#   for info on defining R functions from Python to deal with the global environment here.
#   Potentially can do a different method of calling rfile using the global environment to check if it exists, etc. Safer?
if (robjects.r['exists'](rfile)): 
    print(rfile + " Exists.\n")
else:
    print(rfile + " Does Not Exist.\n")




#TODO: determine input type to DBGLM1.
#TODO: parse relevant JuncBASE table info into R object form (someething within robjects, possibly StrVector or equiv?)
