
  runner
  ├── bin
  │   ├── crop.pl        # you use it to merge histograms onse all jobs have finished
  │   └── run.pl         # used by Makefile
  ├── JetOrange2018.h
  ├── Makefile           # you use it to run jobs
  ├── MakeHist_true_v2.C # ROOT macro for hadron level distributions
  ├── MakeHist_v8.C      # ROOT macro for detector level distributions



####################################################################################
# Submit jobs:

Makefile MUST be evaluated in this directory ("runner"):
$ cd runner

Data:
$ make data out=<output-directory> macro=MakeHist_v8.C

Monte Carlo (detector level):
$ make MC_incl out=<output-directory> macro=MakeHist_v8.C

Monte Crlo (hadron level):
$ make MC_incl out=<output-directory> macro=MakeHist_true_v2.C

In order to produce control plots for different multiplicity and pt-range run:
$ make MC_incl out=<output-directory> macro=MakeHist_pt.C
You will need to edit a line with `Kt_njet_b`-cut to change multiplicity.


This will create <output-directory> (note: must not exist yet) with subdirectories
referring to the data-lists (i.e. "Sample_ari_incl_nc_DIS_lowQ2_040506e").

Supplied ROOT-macro is copied to the <output-directory> so feel free to modify
it after that `make ...` finished suubitting jobs.


You can do it with `hadd` or use a script `crop.pl`. (see example below)

Example including previous step with `make`:

$ make data out=data-blabla macro=MakeHist_v8.C
$ ./bin/crop.pl -o data-blabla-out data-blabla/*    # note a star at the end!!
And press 'y' to confirm.
And to merge all periouds:
$ cd data-blabla-out
$ hadd all.root *.root

