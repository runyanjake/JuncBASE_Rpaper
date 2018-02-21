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

R_FILE_NAME = "HELLOWORLD.R"

# Check existence of desired file.
#   For dealing with the R global environment, use http://rpy.sourceforge.net/rpy2/doc-dev/html/introduction.html 
#   for info on defining R functions from Python to deal with the global environment here.
#   Potentially can do a different method of calling rfile using the global environment to check if it exists, etc. Safer?

exists = robjects.r['exists'](R_FILE_NAME)
print("1st value of exists: ") 
print(exists)

if (robjects.r['exists'](R_FILE_NAME)): 
    print(R_FILE_NAME + " Exists.\n")
else:
    print(R_FILE_NAME + " Does Not Exist.\n")

#execution of .R file from rpy2 using R's !!! "source(<filename>)" !!! function
# NOTE: Other methods may be utilized to define it in R namespace (see test.py) and then call.
# METHOD 1: Calling R file as script, therefore executes only what's in the file.
rfile = "HELLOWORLD.R"
robjects.r['source'](rfile)

exists = robjects.r['exists'](R_FILE_NAME)
print("2nd value of exists: ") 
print(exists)

# METHOD 2: Calling function itself (this reads the file specified by string + .R ext)
helloworld = robjects.r['HELLOWORLD']
helloworld()
#dbglm1 = robjects.r['DBGLM1'] #need full path or local file!

# Print R global environment
print('R Global Environment holds:\n')
print(robjects.r.ls(robjects.globalenv))

#TODO: check file existence etc
#TODO: determine input type to DBGLM1.
#TODO: parse relevant JuncBASE table info into R object form (someething within robjects, possibly StrVector or equiv?)
