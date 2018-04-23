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

#messing with R matrices
from rpy2.robjects.packages import importr #messing with R matrices
from rpy2.robjects import FloatVector #messing with R matrices

#######################################################################
####################### CONSTANT DEFINITIONS ##########################
#######################################################################
#OptionParser Defaults
DEF_THRESH = 10         #min number of lines    
DEF_DPSI_THRESH = 5.0   #min numerical diff needed to consider change
#debug statement toggle
DEBUG_STMTS = False

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
        optionParser.add_option("--debug",
                          dest="debug",
                          type="string",
                          help="""A value of 0 disables debug statements.
                                Any other value enables them. This arg
                                is optional.""",
                          default="0")
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

        #setup debug statments
        if(options.debug != "0"):
            DEBUG_STMTS = True
        else:
            DEBUG_STMTS = False

        #ensure jb table and r file exist.
        R_FUNC_PATH = './DoubleExpSeq/R/DBGLM1.R'
        checkImportantFiles(options.jb_table, R_FUNC_PATH)

        #ensure that file size is adequate.
        fileLength = getNumLinesNoKey(options.jb_table)
        checkSampleSizeThresh(fileLength, options.threshold)

        #Read JB tables into R-objects compatible with DBGLM1
        #and can be passed by reference/changed in method
        yvalues = [] #python list to be y
        mvalues = [] #python list to be m
        groups = None #"vector/factor w expr grp/cond for each sample"
        shrinkMethod = "WEB" #WEB- or DEB-Seq (default WEB)
        contrast = None #"size 2 vector specifying which group levels to compare"
        fdrLevel = 0.05 #"thresh for significant event (not same as delta_thresh?)"
        useAllGroups = None #"use contrast's groups to est dispersion?"

        parseJBTable(options.jb_table, yvalues, mvalues, options.delta_thresh, getArity(options.jb_table)-11.0, fileLength)

        #convert python lists to matrices for R 
        rmatrix = robjects.r['matrix']
        yFloatVec = FloatVector(yvalues) #"numeric matrix of inclusion counts"
        mFloatVec = FloatVector(mvalues) #"numeric matrix of total counts (incl+excl)"
        #matrix rows is arity, number of cols is number of file rows-1
        y = rmatrix(yFloatVec, nrow=fileLength, ncol=getArity(options.jb_table)-11.0,byrow=True)
        m = rmatrix(mFloatVec, nrow=fileLength, ncol=getArity(options.jb_table)-11.0,byrow=True)

        print('Y inputs as rpy2 FloatVector: ')
        print(yFloatVec)
        print('M inputs as rpy2 FloatVector: ')
        print(mFloatVec)

        print('==================================================================')
        print('==================================================================')

        print('Y inputs as R matrix: ')
        print(y)
        print('M inputs as R matrix: ')
        print(m)

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

#returns number of lines in file, not including the key.
#equivalent to the number of recorded events in the data file.
def getNumLinesNoKey(filepath):
    ctr = 0
    with open(filepath, 'rt') as ctrfile:
        reader = csv.reader(ctrfile, delimiter='\t')
        for line in reader: 
            ctr += 1
    return ctr-1 #jb_table first line is a key.

#returns the arity of the jb table. (num cols per entry)
#can be used with static offset to obtain number of samples
def getArity(filepath):
    with open(filepath, 'rt') as readfile:
        reader = csv.reader(readfile, delimiter='\t')
        for line in reader: 
            return len(line) #break on first

#checks to make sure the juncbase table and R file exist.
def checkImportantFiles(jb_table, r_file):
    if(not os.path.exists(jb_table)):
        print('doubleExpSeq.py: ERROR: the --jb_table option (' 
            + options.jb_table + ") does not exist.\n")
        sys.exit(1)
    # (for now is just relative to pwd)
    if(not os.path.exists(r_file)):
        print('doubleExpSeq.py: ERROR: expected DBGLM1.R file (' 
            + r_file + ") does not exist.\n")
        sys.exit(1)

#checks to make sure file contains more than --thresh events 
def checkSampleSizeThresh(num_lines, thresh):
    if(thresh > num_lines):
        print('doubleExpSeq.py: ERROR: file does not contain more than '
            + str(thresh) + ' values. Update --thresh param or try with '+ 
            'different file.')
        sys.exit(1)

#verify a line to ensure it satisfies the delta-thresh condition
#return TRUE if line is ok to use.
#return FALSE if line does not satisfy delta-thresh
#FOR RIGHT NOW THIS IS TAKEN TO BE PERCENT DIFFERENCE FROM MEAN VALUE ###############
def checkDeltaThresh(line, linenr, numSamples, dthresh):
    total = 0.0
    avg = 0.0
    if(linenr > 1):
        for itor in range(11, len(line)):
            inclexcl = line[itor].split(';')
            total += float(inclexcl[-1])
            #print(total)
        avg = total / numSamples
        #print('Average: ' + str(avg))
        #test all vals
        for itor in range(11, len(line)):
            inclexcl = line[itor].split(';')
            inclAndExclCount = float(inclexcl[-1])
            confidenceRange = avg * (dthresh / 100.0)
            if(inclAndExclCount > (avg + confidenceRange) or inclAndExclCount < (avg - confidenceRange)):
                #print('The line satisfies the delta_thresh condition. At least one value (' + str(inclAndExclCount) + ') appeared outside ' + str(avg) + ' +/- ' + str(confidenceRange) + '.')
                #print('Returning True value for line ' + str(linenr))
                return True
        #print('The line did not satisfy the delta_thresh condition. No value appeared outside ' + str(avg) + ' +/- ' + str(avg + confidenceRange) + '.')
        #print('Returning False value for line ' + str(linenr))
        return False

# Reads a JuncBASE table's values into yvalues and mvalues.
# @param filepath The filepath to the table.
# @param yvalues A python list to hold inclusion counts.
# @param mvalues A python list to hold inclusion+exclusion counts.
# @param dthresh The --delta-thresh value.
# @param numSamples The number of samples that appear in the table.
# @param numLines The number of recorded events in the table (number of lines -1)
def parseJBTable(filepath, yvalues, mvalues, dthresh, numSamples, numLines):
    linenr = 1
    print('Reading from ' + filepath + "...\n")
    with open(filepath, 'rt') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')

        for line in reader: 
            if(linenr > 1):
                keep = checkDeltaThresh(line, linenr, numSamples, dthresh)
                if(keep):
                    #DO THESE IF LINE PASSES DELTA_THRESH
                    #  line is indexable w/ array indices 0 - n-1
                    #  Its an (11+n)-tuple, with the last n being total samples.
                    #  order is specified in Step 6A. organized by --sample_setX
                    # print('Line '+str(linenr)+':')
                    for itor in range(11, len(line)):
                        inclexcl = line[itor].split(';')
                        yvalues.append(inclexcl[0])
                        sum = float(inclexcl[0]) + float(inclexcl[1])
                        mvalues.append(sum)
                        # print(str(itor)+': '+line[itor] + '  (Inc: ' + inclexcl[0] + ',Exc: ' + inclexcl[-1] + ')')
                    # print('\n')
            linenr += 1
    
#TODO: parse relevant JuncBASE table info into R object form 
#   (someething within robjects, possibly StrVector or equiv?)
#TODO: determine input type to DBGLM1.

#######################################################################
###################### Main Loop Execution ############################
#######################################################################
if __name__ == "__main__": main()

