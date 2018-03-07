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
        # --in_prefix      | full path + name prefix of input file
        # --all_psi_output | full path + name of output file
        # --mt_correction  | is this necessary for WEB/DEB seq?
        # --thresh         | see notes for usage
        # --delta_thresh   | see notes for usage
        # --sample_set1    | prefix for one set of samples
        # --sample_set2    | prefix for the other set of samples
        # NOTE: prefixes come from JBase table entries (last 2n cols)
        optionParser = OptionParser()
        optionParser.add_option("--in_prefix",
                          dest="in_prefix",
                          type="string",
                          help="""Prefix of output files created from
                                  createAS_CountTables. In 
                                  createAS_CountTables this was the -o 
                                  option""",
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
        optionParser.check_required("--in_prefix")
        optionParser.check_required("--all_psi_output")
        optionParser.check_required("--mt_correction")
        optionParser.check_required("--thresh")
        optionParser.check_required("--delta_thresh")
        optionParser.check_required("--sample_set1")
        optionParser.check_required("--sample_set2")

        #Read JB tables into an R-object compatible with DBGLM1
        y = False #"numeric matrix of inclusion counts"
        x = False #"numeric matrix of total counts (incl+excl)"
        groups = False #"vector/factor w expr grp/cond for each sample"
        shrinkMethod = "WEB" #WEB- or DEB-Seq (default WEB)
        contrast = False #"size 2 vector specifying which group levels to compare"
        fdrLevel = 0.05 #"thresh for significant event"
        useAllGroups = True #"use contrast's groups to est dispersion?"

        #alternatively, do the call in a helper function
        #someVal = parseJBTable("some filename passed in. or hardcoded")



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

# Check existence of desired file.
#   For dealing with the R global environment, use 
#       http://rpy.sourceforge.net/rpy2/doc-dev/html/introduction.html 
#   For info on defining R functions from Python to deal with the global 
#       environment here.
#   Potentially can do a different method of calling rfile using the 
#       global environment to check if it exists, etc. Safer?
# @param r_filename: A string containing the filepath to R file.
# @return True: Returns true iff r_filename points to a existing file.
# @return False: Returns false if no such file exists.
def RFileExists(r_filename):
    exists = robjects.r['exists'](r_filename)
    print(exists)
    if (robjects.r['exists'](r_filename)): 
        return True
    else:
        return False

#Attempts to read in a JuncBASE output table.
# NOTE: requires knowledge of the size of sample_set1 and sample_set2
# @param filePath Takes in an output juncBASE table (.txt ext)
# @return: RETURNS SOME ITEM CONTAINING list of inclusion/exclusion counts
def parseJBTable(filepath):
    print('TODO')
    

#TODO: check file existence etc & read in inclusion/exclusion counts ^
#TODO: parse relevant JuncBASE table info into R object form 
#   (someething within robjects, possibly StrVector or equiv?)
#TODO: determine input type to DBGLM1.

#######################################################################
###################### Main Loop Execution ############################
#######################################################################
if __name__ == "__main__": main()
