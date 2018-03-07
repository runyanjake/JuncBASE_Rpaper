all:
	clear
	make runpy3
runpy2:
	@echo Running doubleExpSeq with Python2.x.
	python2 doubleExpSeq.py
runpy3:
	@echo Running doubleExpSeq with Python3.x.
	python3 doubleExpSeq.py --jb_table ./notes.txt
