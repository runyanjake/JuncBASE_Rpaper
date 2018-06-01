all:
	clear
	make runpy3

run10lines:
	@echo Failing the --thresh parameter with a file having less that 10 recorded ASEvents with Python3.x.
	python3 doubleExpSeq.py --thresh 10 --jb_table ./ExampleCountTables/MM1s_juncBASE_171117_AS_exclusion_inclusion_counts_9Lines.txt

faildeltathresh:
	@echo Failing the --delta_thresh parameter with a file having not enough difference between samples  with Python3.x.
	python3 doubleExpSeq.py --thresh 10 --jb_table ./ExampleCountTables/MM1s_juncBASE_171117_AS_exclusion_inclusion_counts_BadLine.txt

runfullfile:
	@echo Running doubleExpSeq with Python3.x using the provided file in its entirety.
	python3 doubleExpSeq.py --thresh 4 --delta_thresh 5.0 --jb_table ./ExampleCountTables/MM1s_juncBASE_171117_AS_exclusion_inclusion_counts.txt

run10kfile:
	@echo Running doubleExpSeq with Python3.x using 10k of the lines from the provided file.
	python3 doubleExpSeq.py --thresh 10 --jb_table ./ExampleCountTables/MM1s_juncBASE_171117_AS_exclusion_inclusion_counts_10kLines.txt

initialize:
	python3 doubleExpSeq.py --initialize --thresh 10 --jb_table ./ExampleCountTables/MM1s_juncBASE_171117_AS_exclusion_inclusion_counts_10kLines.txt

runpy2:
	@echo Running doubleExpSeq with Python2.x.
	python2 doubleExpSeq.py
runpy3:
	@echo Running doubleExpSeq with Python3.x.
	python3 doubleExpSeq.py --debug 1 --thresh 10 --jb_table ./ExampleCountTables/MM1s_juncBASE_171117_AS_exclusion_inclusion_counts_SAMPLE.txt
