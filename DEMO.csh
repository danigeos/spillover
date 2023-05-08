#!/bin/csh

#Demo for using spillover and produce a graphic output using GMT4

set prj = $0
set file = $prj.output.txt

set echo

#Bonneville flood:
#spillover -V4 -w5/.5/t -M0/2 -k5e-3/1.5 -S5.2e10/1460/l -O1550/-120/1e4 -B1e13/0/l -z1553/0 -t0/200/-10000  > $file 


#Pre-MSC Paratethys flood into the Mediterranean, 6.7 Ma:
spillover -S3.0e12/-2500/l -O39/-40/10000 -B3.6e14/-0.1/b -w5/.5/t -M1/1 -k2e-4/1.5 -t0/580e3/-1 -z40/0./0 > $file
#spillover -S3.0e12/-2500/l -O39/-40/10000 -B3.6e14/-0.1/b -w5/.5/t -M1/1 -k4e-4/1.5 -t0/580e3/-1 -z40/0./0 > $file
#spillover -S3.4e12/-2500/l -O79/-40/10000 -B3.6e14/-0.1/b -w5/.5/t -M1/1 -k2e-4/1.5 -t0/580e3/-1 -z80/0./0 > $file
#spillover -S3.0e12/-2500/l -O39/-40/10000 -B3.6e14/-0.1/b -w5/.5/t -M1/1 -k.05e-4/1.5 -t0/29580e3/-1 -z40/0./0 > $file

#Mediterranean (Zanclean) flood:
spillover -S3.6e14/-5000/b -O-10/-1000/100000 -B2.51e12/-1509/l -w1.1/0.5 -M0/0 -k2.0e-4/1.5 -t0/80e3/-1 -z0/-1500 -V2 > $file
#spillover -V2 -r0/4.5e3/12e3 -p0/.5 -e0/1.4/1.61 -w5/.5/t -M1/1 -k2.0e-4/1.5 -S3e12/-2500/l -O-10/-1200/100000 -B1.6e14/-0.1/b -o-430/-2000/100000 -b2.5e14/-.1/b -t0/580e3/-1 -z40/0./0 > $file

#Morella crater / Elaver Vallis, Mars (Coleman, 2013, JGR) (2.0 to 3.5e7 m3/s). Check mailing with NCOLEMAN@pitt.edu:
#Manning is better here to get the effect of the present flattening of the present outlet.
#spillover -g3.7 -V2 -w1.6/.5 -M0/2 -k2e-3/1.5 -S4.6e9/1250/b -O1776/-529/34e3 -B1e13/0/l -z1786/0 -t0/300/-10000  > $file 

#Eyre:
#.
#spillover -V2 -w10/0/t -M0/2 -k9e1/1.5 -S6.5e9/10/b -O22/-12/25e3 -B1e13/0/l -z22.1/0 -t0/300/-10000  > $file 


#Sand experiment by Walder et al., Exp#4 (see Garcia-Castellanos & O'Connor, 2018): 
#spillover -V0 -w14/0/t -M1/0 -k3.1e1/1.5 -S2.e1/0/l -O.9485/-1/2 -B1e5/0/l -z.95/0 -t0/.29/-.01 > $file 


unset echo

$softdir/src/spillover/spillover_output2ps_linear_axes.csh $file 

