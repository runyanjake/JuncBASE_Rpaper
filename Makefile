all:
	clear
	make runpy3

help:
	python3 doubleExpSeq.py --help

run10lines:
	@echo Failing the --thresh parameter with a file having less that 10 recorded ASEvents with Python3.x.
	python3 doubleExpSeq.py --thresh 10 --jb_table ./ExampleCountTables/MM1s_juncBASE_171117_AS_exclusion_inclusion_counts_9Lines.txt --col_labels "G1_1,G1_2,G1_3,G2_1,G2_2,G2_3,G3_1,G3_2,G3_3,G4_1,G4_2,G4_3" --contrast 2,3

faildeltathresh:
	@echo Failing the --delta_thresh parameter with a file having not enough difference between samples  with Python3.x.
	python3 doubleExpSeq.py --thresh 10 --jb_table ./ExampleCountTables/MM1s_juncBASE_171117_AS_exclusion_inclusion_counts_BadLine.txt

runfullfile:
	@echo Running doubleExpSeq with Python3.x using the provided file in its entirety.
	python3 doubleExpSeq.py --thresh 10 --delta_thresh 5.0 --jb_table ./ExampleCountTables/MM1s_juncBASE_171117_AS_exclusion_inclusion_counts.txt --col_labels "G1_1,G1_2,G1_3,G2_1,G2_2,G2_3,G3_1,G3_2,G3_3,G4_1,G4_2,G4_3"

run10kfile:
	@echo Running doubleExpSeq with Python3.x using 10k of the lines from the provided file.
	python3 doubleExpSeq.py --thresh 10 --jb_table ./ExampleCountTables/MM1s_juncBASE_171117_AS_exclusion_inclusion_counts_10kLines.txt --col_labels "G1_1,G1_2,G1_3,G2_1,G2_2,G2_3,G3_1,G3_2,G3_3,G4_1,G4_2,G4_3"

run1kfile:
	@echo Running doubleExpSeq with Python3.x using 10k of the lines from the provided file.
	python3 doubleExpSeq.py --thresh 10 --jb_table ./ExampleCountTables/MM1s_juncBASE_171117_AS_exclusion_inclusion_counts_1kLines.txt --col_labels "G1_1,G1_2,G1_3,G2_1,G2_2,G2_3,G3_1,G3_2,G3_3,G4_1,G4_2,G4_3"

initialize:
	python3 doubleExpSeq.py --initialize --thresh 10 --jb_table ./ExampleCountTables/MM1s_juncBASE_171117_AS_exclusion_inclusion_counts_10kLines.txt

showhelp:
	@echo Printing help menu.
	python3 doubleExpSeq.py --help

cmdlinetester:
	@echo Testing command line arguments.
	python3 doubleExpSeq.py --thresh 10 --delta_thresh 5.0 --jb_table ./ExampleCountTables/MM1s_juncBASE_171117_AS_exclusion_inclusion_counts.txt --col_labels "G1_1,G1_2,G1_3,G2_1,G2_2,G2_3,G3_1,G3_2,G3_3,G4_1,G4_2,G4_3" --shrinkmethod "WEB" --contrast 2,4 --fdrlevel 0.05 --store_dbglm1_output

runpy2:
	@echo Running doubleExpSeq with Python2.x.
	python2 doubleExpSeq.py
runpy3:
	@echo Running doubleExpSeq with Python3.x.
	python3 doubleExpSeq.py --debug --thresh 10 --jb_table ./ExampleCountTables/MM1s_juncBASE_171117_AS_exclusion_inclusion_counts_SAMPLE.txt
