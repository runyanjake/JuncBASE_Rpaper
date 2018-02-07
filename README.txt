Notes to Jake

Python version I had to use was 3.6.1 (will this work with the rest of JuncBASE?)
R version and rpy2 versions were just the latest versions.

Generally working function call: python3 test.py

See here for using rpy2 vectors with R: http://rpy2.readthedocs.io/en/version_2.8.x/introduction.html#calling-r-functions
Here is a good documentation for rpy2: http://rpy.sourceforge.net/rpy2/doc-2.1/html/robjects.html

Questions for Thursday 2/8 meeting with Angela: 
-I may have done something wrong in my juncBASE run, all tmp files in /scratch/jmrunyan/MM1s_wiita_lab/getASEventReadCounts_first 
        seem to have only headers. These seem to be combined to get the count tables output from JuncBASE
    -going into a chromosome it seems better, e.g. 
            /scratch/jmrunyan/MM1s_wiita_lab/getASEventReadCounts_first/cfz_1_S3/cfz_1_S3_chr1/cfz_1_S3_chr1_all_AS_event_info.txt 
            has items in it, though the numbers seem much lower than what's in the one Angela gave as an example