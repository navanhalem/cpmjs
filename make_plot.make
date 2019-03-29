all: 45_0

45_0: results/kt_45_be_0.pdf

results/kt_45_be_0.pdf: 45 data/sensitive/sensitive_summary.txt data/arresting/arresting_summary.txt data/motile/motile_summary.txt results output
	Rscript Rscripts/summarydata_plots.R $(time) 45 0 $@

data/sensitive/sensitive_summary.txt: data/sensitive/45/sensitive-lambdachem0-killingtime45-biasedentry0-sim1.txt
	bash data/sensitive/combine_sensitive.sh > $@

data/arresting/arresting_summary.txt: data/arresting/45/arresting-lambdachem0-killingtime45-biasedentry0-sim1.txt
	bash data/arresting/combine_arresting.sh > $@

data/motile/motile_summary.txt: data/motile/45/motile-lambdachem0-killingtime45-biasedentry0-sim1.txt
	bash data/motile/combine_motile.sh > $@

data/sensitive/45/sensitive-lambdachem0-killingtime45-biasedentry0-sim1.txt: makefiles/makefile_0_45.make
	make -f makefiles/makefile_0_45.make

data/arresting/45/arresting-lambdachem0-killingtime45-biasedentry0-sim1.txt: makefiles/makefile_0_45.make
	make -f makefiles/makefile_0_45.make

data/motile/45/motile-lambdachem0-killingtime45-biasedentry0-sim1.txt: makefiles/makefile_0_45.make
	make -f makefiles/makefile_0_45.make

makefiles/makefile_0_45.make:
	bash makefiles/make_makefiles.sh $(time) 0 45 $(no_sim) > makefiles/makefile_0_45.make

results: 
	mkdir results

output: 
	mkdir output

datafolders: data data/sensitive data/arresting data/motile

data:
	mkdir $@

data/sensitive:
	mkdir $@

data/arresting:
	mkdir $@

data/motile:
	mkdir $@

45: datafolders data/sensitive/45 data/arresting/45 data/motile/45

data/sensitive/45: 
	mkdir $@

data/arresting/45: 
	mkdir $@

data/motile/45:
	mkdir $@








