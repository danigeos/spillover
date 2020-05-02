#!/bin/csh

#Demo for using spillover and produce a graphic output using GMT4

set prj = $0
set file = $prj.output.txt

set echo

#Bonneville flood:
spillover -V4 -w5/.5/t -M0/1 -k4e-3/1.5 -S5.2e10/1460/l -O1550/-120/1e4 -B1e13/0/l -z1553/0 -t0/5000/-10000  > $file 

unset echo

$softdir/src/spillover/spillover_output2ps_linear_axes.csh $file 

