# @author Jake Runyan
# Written for BrooksLab 2018

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
import numpy as numpy

#R imports
from rpy2.robjects.packages import importr #importing R functions
from rpy2.robjects import FloatVector #messing with R matrices
from rpy2.robjects import StrVector #messing with R matrices
from rpy2.robjects import BoolVector #reading R matrix (setting columns to read from)

#######################################################################
####################### CONSTANT DEFINITIONS ##########################
#######################################################################
#OptionParser Defaults
DEF_THRESH = 10
DEF_DPSI_THRESH = 5.0 
DEF_FDR = 0.05
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

def tokenizeargs(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split(','))

#######################################################################
###################### Main Loop Definition ###########################
#######################################################################
def main():

    #held for creating files
    now = datetime.datetime.now()

    #initialize an OptionParser
    optionParser = OptionParser()
    optionParser.add_option("--initialize",
                        action="store_true", 
                        dest="is_first_run", 
                        default=False,
                        help="""(REQUIRED ON FIRST RUN ONLY)
                            Before running, we must import the R 
                            package. First run this script with 
                            '--initialize' to import the package.
                            You should see about 8 runtime warnings
                            and a message noting where files have been
                            installed. """)
    optionParser.add_option("--jb_table",
                        dest="jb_table",
                        type="string",
                        help="""(REQUIRED) The full path and filename of the 
                            JuncBASE table that will be used
                            for all calculations.""",
                        default=None)
    optionParser.add_option("--col_labels",
                        dest="col_labels",
                        type="string",
                        action='callback',
                        callback=tokenizeargs,
                        help="""(REQUIRED) Comma delimited list of column
                         identifiers for replicates. Section before '_' 
                         is used to group samples. MUST be 
                         of form 'S1_1' where LHS is sample 
                         identifier and RHS is replicate identifier.
                         Example: --col_labels "G1_1, G1_2, 
                         G1_3, G2_1, G2_2, G2_3, G3_1, 
                         G3_2, G3_3, G4_1, G4_2, G4_3" """,
                        default=None)
    optionParser.add_option("--contrast",
                        dest="contrast",
                        type="string",
                        action='callback',
                        callback=tokenizeargs,
                        help="""(REQUIRED) Length 2 list of the indexes
                         in --col_labels to compare.
                         Example: --contrast "2, 3" will compare the 
                         G2 and G3 samples as described in the above 
                         example by looking at all replicates with a 
                         G2_ or G3_ prefix in the supplied table.
                         NOTE: behavior undefined if parameters do not
                         match to a group to contrast. """,
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
    optionParser.add_option("--fdrlevel",
                        dest="fdrlevel",
                        type="float",
                        help="""False discovery rate value, directly 
                        plugged into DBGLM1. Default=%d""" 
                                % DEF_FDR,
                        default=DEF_FDR)
    optionParser.add_option("--debug",
                        action="store_true", 
                        dest="debug", 
                        default=False,
                        help="""Optional debugging statements.
                            '--debug' to enable. """)
    optionParser.add_option("--store_dbglm1_output",
                        action="store_true", 
                        dest="store_dbglm1_output", 
                        default=False,
                        help="""Stores the return value of the 
                        DBGLM1 call into the file specified by 
                        --dbglm1_output_filename.""")
    optionParser.add_option("--dbglm1_output_filename",
                        dest="dbglm1_output_filename",
                        type="string",
                        help="""The name of the file to store 
                        DBGLM1 output if --store_dbglm1_output 
                        option is supplied. Default is dbglm1out.txt""",
                        default="dbglm1out.txt")
    optionParser.add_option("--shrinkmethod",
                        dest="shrinkmethod",
                        type="string",
                        help="""The shrinkage method as described 
                        by the DoubleExpSeq documentation on the CRAN.
                        Can be either 'WEB' (default) or 'DEB'.""",
                        default="WEB")
    #NOTE THIS ARG IS NOT FULLY LINKED UP BUT CAN BE SET
    optionParser.add_option("--useallgroups",
                        action="store_true", 
                        dest="useallgroups", 
                        default=False,
                        help="""By default, we use only ASEvent 
                        inclusion/exclusion counts that pass thresh 
                        tests for ONLY the replicates in groups we're 
                        contrasting. First run this script with 
                        '--useallgroups' to instead use all groups' 
                        counts (will provide more data points at risk 
                        of modifying distribution). Default=FALSE NOTE: To be implemented fully.""")
    optionParser.add_option("--store_MAplot",
                        action="store_true", 
                        dest="store_MAplot", 
                        default=False,
                        help="""Creates an M-A plot comparing the original 
                        and computed distributions in the file specified 
                        by --store_MAplot_filename.""")
    optionParser.add_option("--store_MAplot_filename",
                        dest="store_MAplot_filename",
                        type="string",
                        help="""The name of the file to store 
                        MA plot if --store_MAplot 
                        option is supplied. Default is maplot.txt""",
                        default="maplot.txt")
    optionParser.add_option("--store_MAraw",
                        action="store_true", 
                        dest="store_MAraw", 
                        default=False,
                        help="""Stores the MA plot's raw values into 
                        the file specified by --store_MAraw_filename.""")
    optionParser.add_option("--store_MAraw_filename",
                        dest="store_MAraw_filename",
                        type="string",
                        help="""The name of the file to store 
                        MA plot's raw values if --store_MAraw 
                        option is supplied. Default is marawout.txt""",
                        default="marawout.txt")
    #grab arguments & put into list
    (options, args) = optionParser.parse_args()

    #setup debug statments
    if(options.debug):
        DEBUG_STMTS.append(True)
    else:
        DEBUG_STMTS.append(False)

    #greet with version info
    printVersionInfo()

    ##################################################################
    ######################## FIRST TIME RUN ##########################
    if(options.is_first_run):
        #Import R package
        #trying to mirror from CRAN (Stack overflow and rpy2 documentation
        #                            recommend doing it this way)
        utils = importr('utils')
        log('Setting CRAN as default to install DoubleExpSeq from...')
        utils.chooseCRANmirror(ind=1) #default the CRAN repo
        log('Done.')
        #First time installation of R Package
        log('Performing first time installation of the DoubleExpSeq package...')
        utils.install_packages('DoubleExpSeq') #is this really required?
        log('Done. Now you can run without the --initialize flag.')
        sys.exit()

    ##################################################################
    ####################### NORMAL  EXECUTION ########################
    else:
        #### CHECK REQUIRED ARGUMENTS ARE SUPPLIED ####
        optionParser.check_required("--jb_table")
        optionParser.check_required("--col_labels")
        optionParser.check_required("--contrast")

        #### ENSURE SUPPLIED TABLE EXISTS #### 
        checkImportantFiles(options.jb_table)

        #### CHECK VALIDITY OF PASSED ARGUMENTS BEFORE MUCH COMPUTATION HAPPENS ####
        if options.shrinkmethod != "WEB" and options.shrinkmethod != "DEB":
            print("\n\n")
            tb = sys.exc_info()[2]
            raise Exception("\nThe shrinkage method provided (\"" + options.shrinkmethod + "\") wasn't recognized. The only parameters accepted are WEB or DEB (caps only, ex: --shrinkmethod WEB).").with_traceback(tb)
            exit(1)
        #----
        tmp = 0
        for item in options.contrast:
            tmp = tmp + 1
        if not tmp == 2:
            print("\n\n")
            tb = sys.exc_info()[2]
            raise Exception("\nThe contrast option provided (\"" + str(options.contrast) + "\") wasn't the right size. Only 2 groups can be contrasted at a time.").with_traceback(tb)
            exit(1)
        #----
        if options.fdrlevel < 0 or options.fdrlevel > 1:
            print("\n\n")
            tb = sys.exc_info()[2]
            raise Exception("\nThe fdrlevel option provided (\"" + str(options.fdrlevel) + "\") cannot be less than 0 (0%) or greater than 1 (100%).").with_traceback(tb)
            exit(1)
        #----
        groups_pylist = [] #create the group labels.
        for label in options.col_labels:
            firstinstance = label.find('_')
            identifier = label[:firstinstance]
            groups_pylist.append(identifier)

        #### CHECK ARITY OF TABLE MATCHES SIZE OF COLNAMES AND THAT EACH LABEL IS FORMATTED CORRECTLY ####
        numcolsgiven = len(options.col_labels)
        filecolumns = getArity(options.jb_table)-11
        if numcolsgiven != filecolumns: #arity must match for labels and matrix
            print("\n\n")
            tb = sys.exc_info()[2]
            raise Exception("\nNumber of column labels (" + str(numcolsgiven) + ") did not equal the number of data columns inferred from the given JuncBASE file (" + str(filecolumns) + "). Review the sizes and try again.").with_traceback(tb)
            exit(1)
        for label in options.col_labels: #section before each sample's underscore should be a non-null substring
            firstinstance = label.find('_')
            if firstinstance == -1: #has no underscore NOTE: this might be ok without but for now will leave in.
                print("\n\n")
                tb = sys.exc_info()[2]
                raise Exception("Column labels must be of the form 'S1_1' where LHS is sample identifier and RHS is replicate identifier. Example: --col_labels 'G1_1, G1_2, G1_3, G2_1, G2_2, G2_3, G3_1, G3_2, G3_3, G4_1, G4_2, G4_3'").with_traceback(tb)
                exit(1)
            identifier = label[:firstinstance]
            if identifier=="": #identifiers of the form _1 are not ok
                print("\n\n")
                tb = sys.exc_info()[2]
                raise Exception("Samples cannot be identified with an empty string. (Put something to the left of the underscore that is a unique identifier for a sample.)").with_traceback(tb)
                exit(1)

        #the correct way to use splat operator to map array as function args
        # rc = robjects.r['c']
        # col_labels = rc(*options.col_labels)
        # log("Column Labels: " + str(col_labels))

        #Read JB tables into R-objects compatible with DBGLM1
        y = None #R numeric matrix holding incl only.
        m = None #R numeric matrix holding incl + excl
        groups = None #"vector/factor w expr grp/cond for each sample"
        shrinkMethod = None #WEB- or DEB-Seq (default WEB)
        contrast = None #"size 2 vector specifying which group levels to compare"
        fdrLevel = None #"thresh for significant event (not same as delta_thresh?)"
        useAllGroups = None #"use contrast's groups to est dispersion?"
        numRetainedLines = [0] #required for passing by reference. This is clunky python that has to be done this way apparently.

        yvalues = [] #python list to be y
        mvalues = [] #python list to be m
        rnames = [] #python list holding the exon #'s that are kept (holds ints)
        fileLength = getNumLinesNoKey(options.jb_table)
        parseJBTable(options.jb_table, yvalues, mvalues, options.thresh, (options.delta_thresh / 100.0), getArity(options.jb_table)-11.0, fileLength, numRetainedLines, rnames, options.useallgroups, options.contrast, groups_pylist)
        log("THIS MANY TABLE LINES WERE KEPT: "  + str(numRetainedLines[0]))

        #convert python lists to matrices for R 
        rmatrix = robjects.r['matrix'] #matrix creation
        rc = robjects.r['c'] #R vector creation
        yFloatVec = FloatVector(yvalues) #"numeric matrix of inclusion counts"
        mFloatVec = FloatVector(mvalues) #"numeric matrix of total counts (incl+excl)"
        #matrix cols is arity, number of rows is number of ASEvents kept
        y = rmatrix(yFloatVec, nrow=numRetainedLines[0], ncol=getArity(options.jb_table)-11.0,byrow=True)
        m = rmatrix(mFloatVec, nrow=numRetainedLines[0], ncol=getArity(options.jb_table)-11.0,byrow=True)

        #Add row and colnames so the R function can read things
        #Exons are now labelled by 
        cnames = rc(*options.col_labels) #rc("G1_1", "G1_2", "G1_3", "G2_1", "G2_2", "G2_3", "G3_1", "G3_2", "G3_3", "G4_1", "G4_2", "G4_3")
        m.colnames = StrVector(cnames) 
        m.rownames = StrVector(rnames)
        y.colnames = StrVector(cnames)
        y.rownames = StrVector(rnames)

        #set other params
        groups = rc(*groups_pylist)
        shrinkMethod = rc(options.shrinkmethod)
        contrast = rc(*options.contrast) #the INDICES we compare
        fdrLevel = options.fdrlevel
        useAllGroups = options.useallgroups

        log('Loading the DoubleExpSeq package into this R environment...')
        DoubleExpSeq = importr('DoubleExpSeq') 
        log('Done.')

        log('Running DBGLM1...')
        # resultsG1G2 = DoubleExpSeq.DBGLM1(y,m,groups,shrinkMethod,contrast,fdrLevel,useAllGroups)
        if useAllGroups: 
            #TODO: MODIFY WITH USEALLGROUPS=TRUE
            resultsG1G2 = DoubleExpSeq.DBGLM1(y,m,groups,use_all_groups=True)
        if not useAllGroups: 
            resultsG1G2 = DoubleExpSeq.DBGLM1(y,m,groups)
        WEBsig = resultsG1G2.rx2("Sig") #grab just the $Sig matrix
        rrownames = robjects.r['rownames']
        rownames = rrownames(WEBsig) #grab the rownames
        log('Done.')

        #Dump R script output to a text file if it is required.
        if options.store_dbglm1_output:
            log('Dumping DBGLM1 output to file...')
            f = open(options.dbglm1_output_filename, 'w')
            f.write(str(resultsG1G2))
            f.close()
            log('Done.')

        #Generate an M-A Plot and its raw data if either required
        if options.store_MAplot:
            log('Creating an M-A plot...')
            MAtitle = "M-A Plot Sample " + str(contrast[0]) + " vs " + str(contrast[1]) + " on " + str(now.month) + '/' + str(now.day) + ' ' + str(now.hour) + ':' + str(now.minute)
            MAraw = DoubleExpSeq.DB_MAPlot(y,m,groups,contrast=contrast, de_tags=rownames,main=MAtitle,xlab="XLABEL",ylab="YLABEL")
            #(this returns 2 lists of points as hidden data (give var for return val) along with the pdf)
            if options.store_MAraw:
                log('Dumping MA raw output to file...')
                f = open(options.store_MAraw_filename, 'w')
                f.write(str(MAraw))
                f.close()
                log('Done.')
            log('Done.')

        #Create the file that is the main output for the program.
        log('Generating main output file...')
        filename = "doubleExpSeqOut_"
        if options.useallgroups:
            filename = filename + "UAG_"
        if not options.useallgroups:
            filename = filename + "BG_"
        filename = filename + str(now.month) + '-' + str(now.day) + '_' + str(now.hour) + '.' + str(now.minute) + ".txt"
        f = open(filename, 'w')

        #Gets the indexes of the members of each group (used below)
        label_groups = [] #create the group labels. this will have a problem if 
        for label in groups_pylist:
            identifier = label[1:] #this assumes we have split by underscore leaving a G and a #
            label_groups.append(identifier)
        group1consider = [] #list of what parts of the line to take from jbase line (11 offset)
        group2consider = [] #list of what parts of the line to take from jbase line (11 offset)
        itor = 0
        for label in label_groups:
            if label == contrast[0]:
                group1consider.append(itor)
            elif label == contrast[1]:
                group2consider.append(itor)
            itor = itor + 1

        #itorate through dbglm1 output as the main control:
        routput = resultsG1G2.rx2("All") #All holds data, #sig held the column names for earlier use.
        routput_rows = numpy.asarray(routput) #breaks matrix into list of 9tuples
        #open jb table
        jbtablefile = open(options.jb_table, 'rt')
        jbtable = csv.reader(jbtablefile, delimiter='\t') #$#$!@!@#!

        #write the header
        # isPval>0.05  ref#inInputTable  TAGS_OF_COLS_FROM_INPUT_FILE  POSSIBLY_SPLICE_IN/OUT_COUNTS_FROM_INPUT_FILE  median_psi_group1  median_psi_group2  delta_psi\traw_pval  corrected_pval
        f.write("# isPval>0.05\tref#inInputTable\tTAGS_OF_COLS_FROM_INPUT_FILE\tPOSSIBLY_SPLICE_IN/OUT_COUNTS\tmedian_psi_group1\tmedian_psi_group2\tdelta_psi\traw_pval\tcorrected_pval\n")
        
        ritor = 0
        for row in routput_rows: #ritor: 1-n jbitor: 0-(n-1)
            group1psitotal = []
            group2psitotal = []
            jbitor = 0
            print("R file row " + str(ritor) + " refers to ASEvent number " + str(int(rnames[ritor][5:])) + ", looking for the matching jb file line.")
            while not jbitor == int(rnames[ritor][5:]):
                jbitor = jbitor + 1
            print("R output row " + str(ritor) + " matches with the jb file line number " + str(jbitor) + " with text " + str(jbline))
            f.write("N\t" + str(rnames[ritor][5:]) + "\tTAGS_OF_COLS_FROM_INPUT_FILE\tPOSSIBLY_SPLICE_IN/OUT_COUNTS\t" + "\n")
            ritor = ritor + 1

        jbtablefile.close()
        f.close()
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

#IF useallgroups=True
#checks to make sure each ASEvent from each sample was looked at more than thresh times
#return TRUE if line is ok to use.
#return FALSE if line does not satisfy thresh
def checkThresh_UAG(line, linenr, thresh):
    for itor in range(11, len(line)):
        inclexcl = line[itor].split(';')
        if(float(inclexcl[0]) + float(inclexcl[1]) < thresh):
            # print("Line " + str(linenr-2) + " failed thresh test. total read count: " + str(float(inclexcl[0]) + float(inclexcl[1])) + " < " + str(thresh))
            return False
        # print("Line " + str(linenr-2) + " passed thresh test. total read count: " + str(float(inclexcl[0]) + float(inclexcl[1])) + " >= " + str(thresh))
    return True

#IF useallgroups=True
#verify a line to ensure it satisfies the delta-thresh condition
#return TRUE if line is ok to use.
#return FALSE if line does not satisfy delta-thresh
def checkDeltaThresh_UAG(line, linenr, numSamples, dthresh):
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

#IF useallgroups=False ("By Group")
#Difference between this and checkThresh_BG is the addition of the check: "(itor-11) in toconsider and"
# @param toconsider A list containing the offsets from 11 (aka # of moves from first pair of incl/excl to a pair to consider in test)
def checkThresh_BG(line, linenr, thresh, toconsider):
    for itor in range(11, len(line)):
        inclexcl = line[itor].split(';')
        if((itor-11) in toconsider and float(inclexcl[0]) + float(inclexcl[1]) < thresh):
            # print("Line " + str(linenr-2) + " failed thresh test. total read count at " + str(itor-11) + ": " + str(float(inclexcl[0]) + float(inclexcl[1])) + " < " + str(thresh))
            return False
    # print("Line " + str(linenr-2) + " passed thresh test. total read count: " + str(float(inclexcl[0]) + float(inclexcl[1])) + " >= " + str(thresh))
    return True

#IF useallgroups=False ("By Group")
#Difference between this and checkDeltaThresh_UAG is the addition of line: "if (itor-11) in toconsider:"
def checkDeltaThresh_BG(line, linenr, numSamples, dthresh, toconsider):
    first = line[11].split(';') #Default to first in line.
    max_psi = 0.0
    min_psi = 0.0
    if(not(float(first[1]) == 0 and float(first[1]) == 0)):
        max_psi = float(first[1]) / (float(first[0]) + float(first[1]))
        min_psi = float(first[1]) / (float(first[0]) + float(first[1]))

    if(linenr > 1):
        for itor in range(12, len(line)): #since we default to first in line, one less per line to check
            if (itor-11) in toconsider:
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
def parseJBTable(filepath, yvalues, mvalues, thresh, dthresh, numSamples, numLines, numRetained, rnames, usingallgroups, contrast, labels):
    linenr = 1

    numnotkept = 0
    numfailthresh = 0
    numfaildthresh = 0
    numfailboth = 0

    log('Reading from ' + filepath + "...\n")
    with open(filepath, 'rt') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        #determine what parts of a row to consider for each line read (this is relevant only if not usingallgroups)
        label_groups = [] #create the group labels. this will have a problem if 
        for label in labels:
            identifier = label[1:] #this assumes we have split by underscore leaving a G and a #
            label_groups.append(identifier)
        log("Groups: " + str(contrast))
        log("Labels: " + str(label_groups))
        for line in reader: 
            if(linenr > 1):
                passthresh = None
                passdthresh = None
                if usingallgroups:
                    passthresh = checkThresh_UAG(line, linenr, thresh)
                    passdthresh = checkDeltaThresh_UAG(line, linenr, numSamples, dthresh)
                if not usingallgroups:
                    toconsider = [] #list of what parts of the line to take from jbase line (11 offset)
                    itor = 0
                    for label in label_groups:
                        if label in contrast:
                            toconsider.append(itor)
                        itor = itor + 1
                    passthresh = checkThresh_BG(line, linenr, thresh, toconsider)
                    passdthresh = checkDeltaThresh_BG(line, linenr, numSamples, dthresh, toconsider)

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
