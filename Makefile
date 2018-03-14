all:
	clear
	make runpy3
runpy2:
	@echo Running doubleExpSeq with Python2.x.
	python2 doubleExpSeq.py
runpy3:
	@echo Running doubleExpSeq with Python3.x.
	python3 doubleExpSeq.py --thresh 10 --jb_table ./ExampleCountTables/MM1s_juncBASE_171117_AS_exclusion_inclusion_counts_SAMPLE.txt
