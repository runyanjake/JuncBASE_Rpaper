# @author Jake Runyan
# @notes
#   Intended for use with Python 2.7

#**** IMPORTS ****#
import rpy2
from rpy2.robjects.packages import importr #importing R packages.
from rpy2 import robjects #how r objects are defined in python.
from rpy2.rinterface import R_VERSION_BUILD #print version info
import optparse #OptionParser
import sys #sys.exit etc
import os #os.path.exists
import csv #for input of jb_table
import datetime #for naming generated files
import logging #for debug statements

#messing with R matrices
from rpy2.robjects.packages import importr #messing with R matrices
from rpy2.robjects import FloatVector #messing with R matrices
from rpy2.robjects import StrVector #messing with R matrices

#######################################################################
####################### CONSTANT DEFINITIONS ##########################
#######################################################################
#OptionParser Defaults
DEF_THRESH = 10         #min number of lines    
DEF_DPSI_THRESH = 5.0   #min numerical diff needed to consider change
#debug statement toggle
DEBUG_STMTS = []

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
    #initialize an OptionParser
    ##**** Options ****##
    # --initialize     | facilitates R import of DoubleExpSeq package
    # --debug          | debug statements
    # --jb_table       | full path + name of juncBASE table
    # --all_psi_output | full path + name of output file
    # --mt_correction  | is this necessary for WEB/DEB seq?
    # --thresh         | see notes for usage
    # --delta_thresh   | see notes for usage
    # --sample_set1    | prefix for one set of samples
    # --sample_set2    | prefix for the other set of samples
    # NOTE: prefixes come from JBase table entries (last 2n cols)
    optionParser = OptionParser()
    optionParser.add_option("--initialize",
                        action="store_true", 
                        dest="is_first_run", 
                        default=False,
                        help="""Before running, we must import the R 
                            package. First run this script with 
                            '--initialize' to import the package.
                            You should see about 8 runtime warnings
                            and a message noting where files have been
                            installed. """)
    optionParser.add_option("--debug",
                        dest="debug",
                        type="string",
                        help="""Optional debugging statements.
                            '--debug 1' to enable. """,
                        default="0")
    optionParser.add_option("--jb_table",
                        dest="jb_table",
                        type="string",
                        help="""The full path and filename of the 
                            JuncBASE table that will be used
                            for all calculations.""",
                        default=None)
    optionParser.add_option("--thresh",
                        dest="thresh",
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

    #setup debug statments
    if(options.debug != "0"):
        DEBUG_STMTS.append(True)
    else:
        DEBUG_STMTS.append(False)

    #greet with version info
    printVersionInfo()

    #### FIRST TIME INSTALLATION ####
    if(options.is_first_run):
        #Import R package
        #trying to mirror from CRAN (Stack overflow and rpy2 documentation point to this way)
        utils = importr('utils')
        log('Setting CRAN as default to install DoubleExpSeq from...')
        utils.chooseCRANmirror(ind=1) #default the CRAN repo
        log('Done.')
        #First time installation of R Package
        log('Performing first time installation of the DoubleExpSeq package...')
        utils.install_packages('DoubleExpSeq') #is this really required?
        log('Done.')
        sys.exit()

    #### NORMAL PROGRAM EXECUTION ####
    else:

        #### CHECK REQUIRED ARGUMENTS ARE SUPPLIED ####
        optionParser.check_required("--jb_table")
        # optionParser.check_required("--thresh")
        # optionParser.check_required("--delta_thresh")
        # optionParser.check_required("--sample_set1")
        # optionParser.check_required("--sample_set2")

        #### ENSURE SUPPLIED TABLE EXISTS #### 
        checkImportantFiles(options.jb_table)

        #Read JB tables into R-objects compatible with DBGLM1
        #and can be passed by reference/changed in method
        y = None #R numeric matrix holding incl only.
        m = None #R numeric matrix holding incl + excl
        groups = None #"vector/factor w expr grp/cond for each sample"
        shrinkMethod = None #WEB- or DEB-Seq (default WEB)
        contrast = None #"size 2 vector specifying which group levels to compare"
        fdrLevel = None #"thresh for significant event (not same as delta_thresh?)"
        useAllGroups = None #"use contrast's groups to est dispersion?"
        numRetainedLines = [0] #required for making array. This is clunky python that has to be done this way apparently.

        yvalues = [] #python list to be y
        mvalues = [] #python list to be m
        rnames = [] #python list holding the exon #'s that are kept (holds ints)
        fileLength = getNumLinesNoKey(options.jb_table)
        parseJBTable(options.jb_table, yvalues, mvalues, options.thresh, (options.delta_thresh / 100.0), getArity(options.jb_table)-11.0, fileLength, numRetainedLines, rnames)

        log("THIS MANY WERE KEPT: "  + str(numRetainedLines[0]))

        # ###################################################################################################################
        # sys.exit()
        # ###################################################################################################################



        #convert python lists to matrices for R 
        rmatrix = robjects.r['matrix'] #matrix creation
        rc = robjects.r['c'] #R vector creation
        yFloatVec = FloatVector(yvalues) #"numeric matrix of inclusion counts"
        mFloatVec = FloatVector(mvalues) #"numeric matrix of total counts (incl+excl)"
        #matrix cols is arity, number of rows is number of rows kept
        y = rmatrix(yFloatVec, nrow=numRetainedLines[0], ncol=getArity(options.jb_table)-11.0,byrow=True)
        m = rmatrix(mFloatVec, nrow=numRetainedLines[0], ncol=getArity(options.jb_table)-11.0,byrow=True)

        #Add row and colnames so the R function can read things
        #Exons are now labelled by 
        cnames = rc("G1_1", "G1_2", "G1_3", "G2_1", "G2_2", "G2_3", "G3_1", "G3_2", "G3_3", "G4_1", "G4_2", "G4_3")
        m.colnames = StrVector(cnames) 
        m.rownames = StrVector(rnames)
        y.colnames = StrVector(cnames)
        y.rownames = StrVector(rnames)

        log("HERE ARE THE LABELLED MATRICEs")
        log(m)
        log(y)
        log("THIS MANY WERE KEPT: "  + str(numRetainedLines[0]))

        f = open("Matrices.txt", 'w')
        f.write(str(y))
        f.write(str(m))
        f.close()

        # ###################################################################################################################
        # sys.exit()
        # ###################################################################################################################

        #set other params
        groups = rc("CTRL", "CTRL", "CTRL",     #1 #NOTE: Groups should probably be determined based on the cnames variable that must be passed in.
                    "E7107", "E7107", "E7107",  #2
                    "MELPH", "MELPH", "MELPH",  #3  
                    "CFZ", "CFZ", "CFZ")        #4
        shrinkMethod = rc("WEB")
        contrast = rc(2,3) #the INDICES we compare ^
        fdrLevel = 0.05
        useAllGroups = True

        log('Loading the DoubleExpSeq package into this R environment...')
        DoubleExpSeq = importr('DoubleExpSeq') 
        log('Done.')

        log('Running DBGLM1...')
        # resultsG1G2WEB = DoubleExpSeq.DBGLM1(y,m,groups,shrinkMethod,contrast,fdrLevel,useAllGroups)
        resultsG1G2WEB = DoubleExpSeq.DBGLM1(y,m,groups)
        WEBsig = resultsG1G2WEB.rx2("Sig") #grab just the $Sig matrix
        rrownames = robjects.r['rownames']
        rownames = rrownames(WEBsig) #grab the rownames
        log('Done.')

        #Dump R script output to a text file if it is required.
        #assumes there's a .txt ending. otherwise filesize must be > 4chars.
        log('Dumping DBGLM1 output to file...')
        now = datetime.datetime.now()
        datafile = options.jb_table[(options.jb_table.rfind("/")+1):(len(options.jb_table) - 4)] + '_doubleExpSeqOutdata_' + str(now.month) + '-' + str(now.day) + '.' + str(now.hour) + ':' + str(now.minute) + '_' + '.txt'
        f = open(datafile, 'w')
        f.write(str(resultsG1G2WEB))
        f.close()
        log('Done.')

        #Generate an M-A Plot
        log('Creating an M-A plot...')
        MAtitle = "M-A Plot Sample " + str(contrast[0]) + " vs " + str(contrast[1]) + " on " + str(now.month) + '/' + str(now.day) + ' ' + str(now.hour) + ':' + str(now.minute)
        tmp = DoubleExpSeq.DB_MAPlot(y,m,groups,contrast=contrast, de_tags=rownames,main=MAtitle,xlab="XLABEL",ylab="YLABEL")
        log(tmp)
        #(this returns 2 lists of points as hidden data (give var for return val) along with the pdf)
        log('Done.')

#######################################################################
################# Auxiliary Function Definitions ######################
#######################################################################

# Version Information printing function.
def printVersionInfo():
    log('\n**************VERSION INFORMATION**************')
    log('rpy2 version: ' + rpy2.__version__)
    log('This version of rpy2 was built on R version ')
    log(R_VERSION_BUILD)
    log('************END VERSION INFORMATION************\n')

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
def checkImportantFiles(jb_table):
    if(not os.path.exists(jb_table)):
        log('doubleExpSeq.py: ERROR: the --jb_table option (' 
            + jb_table + ") does not exist.\n")
        sys.exit(1)

#checks to make sure each ASEvent from each sample was looked at more than thresh times
#return TRUE if line is ok to use.
#return FALSE if line does not satisfy thresh
def checkThresh(line, linenr, thresh):
    for itor in range(11, len(line)):
        inclexcl = line[itor].split(';')
        if(float(inclexcl[0]) + float(inclexcl[1]) < thresh):
            # print("Line " + str(linenr-2) + " failed thresh test. total read count: " + str(float(inclexcl[0]) + float(inclexcl[1])) + " < " + str(thresh))
            return False
        # print("Line " + str(linenr-2) + " passed thresh test. total read count: " + str(float(inclexcl[0]) + float(inclexcl[1])) + " >= " + str(thresh))
    return True

#verify a line to ensure it satisfies the delta-thresh condition
#return TRUE if line is ok to use.
#return FALSE if line does not satisfy delta-thresh
def checkDeltaThresh(line, linenr, numSamples, dthresh):
    first = line[11].split(';') #Default to first in line.
    max_psi = 0.0
    min_psi = 0.0
    if(not(float(first[1]) == 0 and float(first[1]) == 0)):
        max_psi = float(first[1]) / (float(first[0]) + float(first[1]))
        min_psi = float(first[1]) / (float(first[0]) + float(first[1]))

    if(linenr > 1):
        for itor in range(12, len(line)): #since we default to first in line, one less per line to check
            inclexcl = line[itor].split(';')
            this_psi = 0.0
            if(not(float(inclexcl[1]) == 0 and float(inclexcl[1]) == 0)):
                this_psi = float(inclexcl[1]) / (float(inclexcl[0]) + float(inclexcl[1]))
            if(this_psi > max_psi):
                max_psi = this_psi
            if(this_psi < min_psi):
                min_psi = this_psi
        if(max_psi - min_psi < dthresh):
            # print("Line " + str(linenr-2) + " failed delta_thresh test. delta_psi: " + str(max_psi - min_psi) + " < " + str(dthresh))
            return False
        # print("Line " + str(linenr-2) + " passed delta_thresh test. delta_psi: " + str(max_psi - min_psi) + " >= " + str(dthresh))
        return True

# Reads a JuncBASE table's values into yvalues and mvalues.
# @param filepath The filepath to the table.
# @param yvalues A python list to hold inclusion counts.
# @param mvalues A python list to hold inclusion+exclusion counts.
# @param dthresh The --delta-thresh value.
# @param numSamples The number of samples that appear in the table.
# @param numLines The number of recorded events in the table (number of lines -1)
def parseJBTable(filepath, yvalues, mvalues, thresh, dthresh, numSamples, numLines, numRetained, rnames):
    linenr = 1

    numnotkept = 0
    numfailthresh = 0
    numfaildthresh = 0
    numfailboth = 0

    log('Reading from ' + filepath + "...\n")
    with open(filepath, 'rt') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')

        for line in reader: 
            if(linenr > 1):
                passthresh = checkThresh(line, linenr, thresh)
                passdthresh = checkDeltaThresh(line, linenr, numSamples, dthresh)
                keep = passthresh and passdthresh
                if(keep):
                    rnames.append("exon_" + str(linenr-1))
                    numRetained[0] +=1 #####
                    #### INCLUDE LINE IF PASSES THRESH TESTS ####
                    for itor in range(11, len(line)):
                        inclexcl = line[itor].split(';')
                        yvalues.append(inclexcl[0])
                        sum = float(inclexcl[0]) + float(inclexcl[1])
                        mvalues.append(sum)
                else:
                    numnotkept += 1 #####
                    if(not(passthresh)):
                        numfailthresh += 1
                    if(not(passdthresh)):
                        numfaildthresh += 1
                    if(not(passthresh) and not(passdthresh)):
                        numfailboth += 1
            linenr += 1
        log("Out of " + str(linenr-2) + " lines, " + str(numRetained[0]) + " were kept and " + str(numnotkept) + " were not kept.")
        log("Number of lines failing thresh test: " + str(numfailthresh))
        log("Number of lines failing delta_thresh test: " + str(numfaildthresh))
        log("Number of lines failing both tests: " + str(numfailboth))

def log(s):
    if DEBUG_STMTS is None or DEBUG_STMTS[0]:
        print(s)

#######################################################################
###################### Main Loop Execution ############################
#######################################################################
if __name__ == "__main__": main()
