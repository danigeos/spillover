#!/bin/csh

#Demo script for using 'spillover' and producing a graphic output using GMT4

set prj = $0
set file = $prj.output.txt

set echo

#Bonneville flood:
#spillover -V4 -w5/.5/t -M0/2 -k5e-3/1.5 -S5.2e10/1460/l -O1550/-120/1e4 -B1e13/0/l -z1553/0 -t0/200/-10000  > $file 


#Pre-MSC Paratethys flood into the Mediterranean, 6.7 Ma:
#spillover -S3.0e12/-2500/l -O39/-40/10000 -B3.6e14/-0.1/b -w5/.5/t -M1/1 -k2e-4/1.5 -t0/580e3/-1 -z40/0./0 > $file

#Mediterranean (Zanclean) flood, 1 basin:
#spillover -S3.6e14/-5000/b -O-10/-1000/100000 -B2.51e12/-1509/l -w1.1/0.5 -M0/0 -k2.0e-4/1.5 -t0/80e3/-1 -z0/-1500 -V2 > $file
#Mediterranean (Zanclean) flood, 2 basins:
#spillover -w1.1/.5/t -M0/0 -k.1e-4/1.5 -S3.6e14/-5000/b -O-10/-1200/100000 -B1.8e12/-1509/l -o-430/-2000/100000 -b.5e12/-1900/b -t0/1580e3/-1 -z0/-1000./-1500 > $file

#Morella crater / Elaver Vallis, Mars (Coleman, 2013, JGR) (2.0 to 3.5e7 m3/s). Check mailing with NCOLEMAN@pitt.edu:
#Manning is better here to get the effect of the present flattening of the present outlet.
#spillover -g3.7 -V2 -w1.6/.5 -M0/2 -k2e-3/1.5 -S4.6e9/1250/b -O1776/-529/34e3 -B1e13/0/l -z1786/0 -t0/300/-10000  > $file 

#Eyre:
#spillover -V2 -w10/0/t -M0/2 -k9e1/1.5 -S6.5e9/10/b -O22/-12/25e3 -B1e13/0/l -z22.1/0 -t0/300/-10000  > $file 

#Sand experiment by Walder et al., Exp#4 (see Garcia-Castellanos & O'Connor, 2018): 
#spillover -V2 -w3/0/t -M1/0 -k5.1e1/1.5 -S2.e1/0/l -O.947/-1/2 -B1e5/0/l -z.95/0 -t0/.09/-.01 > $file 


#Reference runs for spillover Greg's app (simplified examples for the app, all consist of two lakes with bottom at z=0):
#Mediterranean (Zanclean) flood, 1 basin, triangular outlet's section:
#spillover -z1000/0 -S3.6e14/0/l -O990/-1000/100e3 -B2.51e12/0/l -M1/0 -k2.7e-4/1.5 -w5/0/t -t0/80e3/-1 -V2 > $file
#Bonneville flood, triangular outlet's section:
#spillover -z100/0 -S5.0e10/0/l -O98/-10/1e3 -B3.6e14/0/l -M1/0 -k2.5e-3/1.5 -w5/0/t -t0/24000/-1 -V2 > $file
#Tangjiashan Sichuan Quake Lake, Qp=6.50E+03 m3 s-1:
spillover -z26/0 -S6.9e6/0/l -O25/-10/1e3 -B3.6e14/0/l -M1/0 -k1.5e-0/1.5 -w5/0/t -t0/120/-1 -V2 > $file
#Sand experiment by Walder et al., Exp#4 (see Garcia-Castellanos & O'Connor, 2018): 
#spillover -z.95/0 -S2.e1/0/l -O.947/-1/2 -B1e5/0/l -M1/0 -k5.1e1/1.5 -w5/0/t -t0/.09/-.01 -V2  > $file 


unset echo

$softdir/src/spillover/spillover_output2ps_linear_axes.csh $file 

