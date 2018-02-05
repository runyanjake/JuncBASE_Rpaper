#!/usr/local/bin python3 

# @author Jake Runyan
# @notes
#   Using Python3
#   R version 3.3.3 installed for MacOS 10.9 and higher (have 10.10.5)
#   rpy2 installation failed via pip3, succeeded after R installed.

import subprocess #for calling a R function from Python

import rpy2
from rpy2 import robjects #how r objects are defined in python (and how we'll recognize functions)
from rpy2.rinterface import R_VERSION_BUILD #simply for printing version info

#import the R package into Python, and "expose" them as Python objects.
#generally this means changing .'s in R to _'s (b/c '.' shows ownership in Python)
from rpy2.robjects.packages import importr
base = importr('base')
utils = importr('utils')

#Print Version Information
print('\n**************VERSION INFORMATION**************')
print('rpy2 version: ' + rpy2.__version__)
print('This version of rpy2 was built on R version ')
print(R_VERSION_BUILD)
print('**************END VERSION INFORM***************\n')

#creating an R vector
print('Creating an R string vector:')
tmp = robjects.StrVector(['hello', 'world', 'this', 'is', 'Jake']) #create an R vector
print(tmp.r_repr() + '\n') #print said R vector

robjects.r("""sdfg <- function (){
    print("hello world")
    }""")

#printing out R globalenv
print(robjects.r.ls(robjects.globalenv))

#a note on robjects.r vs robjects.globalenv: the r version is meant for use when R is running 
#as an embedded process (i.e. r running in python) and the globalenv version is working 
#directly with R's namespace (this is what's searched first when you make a function call from a R console.

#use rpy2's importr function for calling R functions?
#use rpy2's robjects.r function for calling R functions? 
#use python's subprocess module to run the R file?


#use Python's subprocess module to run the R file. (https://docs.python.org/2/library/subprocess.html)
#           or https://www.mango-solutions.com/blog/integrating-python-and-r-part-ii-executing-r-from-python-and-vice-versa
# command = 'Rscript'
# script = 'HELLOWORLD.R'
# args = []
# cmd = [command, script] + args #fully assembled command
# print('\nRunning Command: ')
# for x in range(len(cmd)):
#     print(cmd[x])
# output = subprocess.check_output(cmd, universal_newlines = True)
# print('Done.')

# command = 'HELLOWORLD.R'
# print('\nRunning Command: ' + command)
# subprocess.call(command)
# print('Done.')
