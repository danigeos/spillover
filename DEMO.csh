#!/bin/csh

#Demo for using spillover and produce a graphic output using GMT4

set prj = $0
set file = $prj.output.txt

set echo

#Bonneville flood:
spillover -V4 -w5/.5/t -M0/1 -k4e-3/1.5 -S5.2e10/1460/l -O1550/-120/1e4 -B1e13/0/l -z1553/0 -t0/5000/-10000  > $file 


#Mediterranean (Zanclean) flood:
spillover -S3.6e14/-5000/b -O-10/-1000/100000 -B2.51e12/-1509/l -w1.1/0.5 -M0/0 -k2.0e-4/1.5 -t0/80e3/-1 -z0/-1500 > $file


unset echo

$softdir/src/spillover/spillover_output2ps_linear_axes.csh $file 

