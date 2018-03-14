#!/usr/local/bin python2

# @author Jake Runyan
# @notes
#   Intended for use with Python 2.7

#**** IMPORTS ****#
import rpy2
from rpy2 import robjects #how r objects are defined in python.
from rpy2.rinterface import R_VERSION_BUILD #print version info
import optparse #OptionParser
import sys #sys.exit etc
import os #os.path.exists
import csv #for input of jb_table



#######################################################################
####################### CONSTANT DEFINITIONS ##########################
#######################################################################
#OptionParser Defaults
DEF_THRESH = 10
DEF_DPSI_THRESH = 5.0

#######################################################################
######################### CLASS DEFINITIONS ###########################
#######################################################################
#Borrowed from Angela Brooks' JuncBASE/comparesamplesets.py
class OptionParser(optparse.OptionParser):
    """
    Adding a method for required arguments.
    Taken from:
    http://www.python.org/doc/2.3/lib/optparse-extending-examples.html
    """
    def check_required(self, opt):
        option = self.get_option(opt)
        # Assumes the option's 'default' is set to None!
        if getattr(self.values, option.dest) is None:
            print ("%s option not supplied. See usage below:" % option)
            self.print_help()
            sys.exit(1)

#######################################################################
###################### Main Loop Definition ###########################
#######################################################################

def main():
        #greet with message and version info
        print('doubleExpSeq.py: Starting main loop...')
        printVersionInfo()

        #initialize an OptionParser
        ##**** Options ****##
        # --jb_table       | full path + name of juncBASE table
        # --all_psi_output | full path + name of output file
        # --mt_correction  | is this necessary for WEB/DEB seq?
        # --thresh         | see notes for usage
        # --delta_thresh   | see notes for usage
        # --sample_set1    | prefix for one set of samples
        # --sample_set2    | prefix for the other set of samples
        # NOTE: prefixes come from JBase table entries (last 2n cols)
        optionParser = OptionParser()
        optionParser.add_option("--jb_table",
                          dest="jb_table",
                          type="string",
                          help="""The full path and filename of the 
                                JuncBASE table that will be used
                                for all calculations.""",
                          default=None)
        optionParser.add_option("--all_psi_output",
                          dest="all_psi_output",
                          type="string",
                          help="""Output file that will contain the PSI
                                values for all events and samples. The 
                                last two columns will correspond to the 
                                raw-pvalue and corrected p-value. If a 
                                generic file is used, this will be the 
                                output file""",
                          default=None)
        optionParser.add_option("--mt_correction",
                          dest="mt_method",
                          type="string",
                          help="""Multiple testing correction Method: 
                                "BH" - Benjamini & Hochberg, 
                                "bonferroni".  Must select these strings 
                                as the option""",
                          default=None)
        optionParser.add_option("--thresh",
                          dest="threshold",
                          type="float",
                          help="""Threshold for minimum abundance
                                  in an event. Default=%d""" 
                                  % DEF_THRESH,
                          default=DEF_THRESH)
        optionParser.add_option("--delta_thresh",
                          dest="delta_thresh",
                          type="float",
                          help="""Minimum PSI(or generic value) 
                                difference between the maximum and 
                                minimum values for a given event to be
                                considered a change. This
                                  should probably be less than the delta
                                  threshold used to filter significantly
                                  associated events. Default=%s""" 
                                  % DEF_DPSI_THRESH,
                          default=DEF_DPSI_THRESH)
        optionParser.add_option("--sample_set1",
                          dest="sample_set1",
                          type="string",
                          help="""Comma delimited list of samples in set 
                                1 or a file with a list of names, one 
                                per line. Names must be in header 
                                columns of input files.""",
                          default=None)
        optionParser.add_option("--sample_set2",
                          dest="sample_set2",
                          type="string",
                          help="""Comma delimited list of samples in set 
                                2 or a file with a list of names, one 
                                per line. Names must be in header 
                                columns of input files.""",
                          default=None)

        #grab arguments & put into list
        (options, args) = optionParser.parse_args()

        #validate supplied arguments
        #NOTE: still need to check if files passed in exist etc.
        optionParser.check_required("--jb_table")
        # optionParser.check_required("--all_psi_output")
        # optionParser.check_required("--mt_correction")
        optionParser.check_required("--thresh")
        # optionParser.check_required("--delta_thresh")
        # optionParser.check_required("--sample_set1")
        # optionParser.check_required("--sample_set2")

        #ensure jb table and r file exist.
        R_FUNC_PATH = './DoubleExpSeq/R/DBGLM1.R'
        checkImportantFiles(options.jb_table, R_FUNC_PATH)

        #ensure that file size is adequate.
        fileLength = numLines(options.jb_table)
        checkSampleSizeThresh(fileLength, options.threshold)

        #Read JB tables into R-objects compatible with DBGLM1
        #NOTE: I'm using python lists for matrices becausethey're mutable
        #and can be passed by reference/changed in method
        y = [] #"numeric matrix of inclusion counts"
        x = [] #"numeric matrix of total counts (incl+excl)"
        groups = False #"vector/factor w expr grp/cond for each sample"
        shrinkMethod = "WEB" #WEB- or DEB-Seq (default WEB)
        contrast = False #"size 2 vector specifying which group levels to compare"
        fdrLevel = 0.05 #"thresh for significant event"
        useAllGroups = True #"use contrast's groups to est dispersion?"

        #alternatively, do the call in a helper function
        #NOTE: remember this can only be done if y and x are mutable
        #(lists are mutable so use these?)
        parseJBTable(options.jb_table, y, x)
        


#######################################################################
################# Auxiliary Function Definitions ######################
#######################################################################

# Version Information printing function.
def printVersionInfo():
    print('\n**************VERSION INFORMATION**************')
    print('rpy2 version: ' + rpy2.__version__)
    print('This version of rpy2 was built on R version ')
    print(R_VERSION_BUILD)
    print('*************END VERSION INFORMATION*************\n')

def numLines(filepath):
    ctr = 0
    with open(filepath, 'rt') as ctrfile:
        reader = csv.reader(ctrfile, delimiter='\t')
        for line in reader: 
            ctr += 1
    return ctr-1 #since jb_table's first line is a key.

#Attempts to read in a JuncBASE output table.
# NOTE: requires knowledge of the size of sample_set1 and sample_set2 (maybe)
# REF: https://docs.python.org/2/library/csv.html
# @param filePath Takes in an output juncBASE table (.txt ext)
# @return: RETURNS SOME ITEM CONTAINING list of inclusion/exclusion counts
def parseJBTable(filepath, y, x):
    linenr = 1
    print('Reading from ' + filepath + "...\n")
    with open(filepath, 'rt') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for line in reader: 
            #line is indexable w/ array indices 0 - n-1
            #Its an (11+n)-tuple, with the last n being total samples.
            #order is specified in Step 6A. organized by --sample_setX
            for itor in range(11, len(line)):
                print(itor)
                print(': ')
                print(line[itor])
            print('\n')
            # print(len(line))
            # print(linenr)
            # print(': ')
            # print(line)
            # print('\n')
            # linenr += 1

#checks to make sure the juncbase table and R file exist.
def checkImportantFiles(jb_table, r_file):
    #Error if the Juncbase table is an existing file.
    if(not os.path.exists(jb_table)):
        print('doubleExpSeq.py: ERROR: the --jb_table option (' 
            + options.jb_table + ") does not exist.\n")
        sys.exit(1)
    #Error if the R function's file cannot be found.
    # (for now is just relative to pwd)
    if(not os.path.exists(r_file)):
        print('doubleExpSeq.py: ERROR: expected DBGLM1.R file (' 
            + r_file + ") does not exist.\n")
        sys.exit(1)

#checks to make sure file contains more than thresh events 
def checkSampleSizeThresh(num_lines, thresh):
    if(thresh > num_lines):
        print('doubleExpSeq.py: ERROR: file does not contain more than '
            + str(thresh) + 'values. Update --thresh param or try with '+ 
            ' different file.')
        sys.exit(1)
    

#TODO: check file existence etc & read in inclusion/exclusion counts ^
#TODO: parse relevant JuncBASE table info into R object form 
#   (someething within robjects, possibly StrVector or equiv?)
#TODO: determine input type to DBGLM1.

#######################################################################
###################### Main Loop Execution ############################
#######################################################################
if __name__ == "__main__": main()

