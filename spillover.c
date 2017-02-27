/*
	Calculates the incision and water flow across an overtopping lake.
	Author: Daniel Garcia-Castellanos, 2009-2016, danielgc@ictja.csic.es
	License: Creative Commons 3.0 BY-SA. This is the original reference 
	to be cited when using this code: 
	Garcia-Castellanos et al., 2009, Catastrophic flood of the 
	Mediterranean after the Messinian salinity crisis. Nature, 462, 
	778-781, doi:10.1038/nature08555 
*/
/*
	Present-day Mediterranean Sea:
		3.77 10^6 km3 with an average salinity of 3.9%, 2.5e12 m2 
	Present-day world's oceans:
		1340 10^6 km3 (Gleick, 1993), 3.6e14 km2
*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
//#include <malloc.h>
#include <time.h>

#define	MAXLENLINE	1024			/*Max. length for character strings, input lines, ...*/
#define TAKE_LINE_2(x, y)	{ char auxstr[MAXLENLINE], *lin; int nfields=0; while (nfields<2) {lin=fgets(auxstr, MAXLENLINE-1, file); if (lin==NULL) break; nfields=sscanf(lin, "%f %f",       &x, &y);};      	if (lin==NULL) break;}
#define AUTHORSHIP		{ fprintf(stderr, "\n\t\t\t\t2009, Daniel Garcia-Castellanos, danielgc@ictja.csic.es\n");}
#define MIN_2(x,y)	(((x)<(y))? (x) : (y))	/*Gives minimum of two values*/
#define MAX_2(x,y)	(((x)>(y))? (x) : (y))	/*Gives maximum of two values*/

#define CAPTION "#time[h] zsill[m] Rh[m]\tSlope\tV[m/s]\tQ[m3/s] \tW[m]\te[m/yr] \tvol0[km3] \tvol1[km3] \tvol2[km3] \tvoltr[km3] \tz0[m]\tz1[m]\tz2[m]"
#define secsperyr 	(365.24  *24*3600)	/*Converts years to seconds*/
#define secsperday 	(24*3600)		/*Converts days to seconds*/

int level_and_area_from_volume (float *hypso_z, float *hypso_a, int np, double vol, double *z, double *area);
int volume_and_area_from_level (float *hypso_z, float *hypso_a, int np, double z, double *vol, double *area);
int read_basin_geometry (double z_sill, char *prm, float *hypso_z, float *hypso_a, int *np);
int syntax (int argc, char **argv);

int verbose_level=1;

int main(int argc, char **argv)
{
	int	i, iarg, niters=0, 
		model_eros=1, model_vel=0, 
		np0=0, np1=0, np2=0, npmax=2000, 
		nbasins=1, full_basins=0, water_conservative=0, decreasable_width=1, curvature_width=0, triangular_width=0, 
		switch_ps=1;
	double	area, 
		time, /*internally and externally in hours!*/
		timeini = 0, 
		timeend = -1e9, 
		dt, dt_default = 36, dtinput, /*[s]*/
		hl1 = -1000,	hl2 = -1000,
		z_sill1 = -1,	z_sill2 = -430,
		z_sill1_ini,
		z_bott0=0, z_bott1 = -2500, z_bott2 = -9999, 
		dist1 = 100e3, 	dist2 = 100e3, 
		Ke = .6e-3, 	expe = 1, 
		Kw = 1.5, 	expw = .5, 
		e0=0, e1=0, e2=0, 
		p0=0, p1=0, p2=0, 
		r0=0, r1=0, r2=0, 
		roughness = .035, Cz=40, /*Note this two should be related through Cz=Rh^(1/6) /n */
		erostotal=0,
		z0=0, z1=-2500, z2=-9999, 
		g = 9.81, denswater = 1020, 
		volini0, volini1, volini2, volini1full;
	float	*hypso0_z, *hypso0_a, 
		*hypso1_z, *hypso1_a, 
		*hypso2_z, *hypso2_a; 
	double	volumetr=0, vol0=0, vol1=0, vol2=0; /*1637338e9, vol2=1927865e9*/

	if (argc<=1) {syntax(argc, argv); exit(0);}
		
	hypso0_z = (float *) calloc(npmax, sizeof(float));
	hypso0_a = (float *) calloc(npmax, sizeof(float));
	hypso1_z = (float *) calloc(npmax, sizeof(float));
	hypso1_a = (float *) calloc(npmax, sizeof(float));
	hypso2_z = (float *) calloc(npmax, sizeof(float));
	hypso2_a = (float *) calloc(npmax, sizeof(float));
	
	/*Interpreting command line*/
	fprintf(stdout, "#");
	for (iarg=0; iarg<argc; iarg++) fprintf(stdout, "%s ", argv[iarg]);
	fprintf(stdout, "\n");
	for (iarg=1; iarg<argc; iarg++) {
		
		if (argv[iarg][0] == '-') {
			int 	ilet;
			double 	value, value2, area;
			char 	prm[MAXLENLINE], *ptr;

			for (ilet=2; ilet < strlen(argv[iarg])+2; ilet++) prm[ilet-2] = argv[iarg][ilet];
			value=atof(prm);
			if (verbose_level>=3) fprintf(stderr, "\n%s %e", argv[iarg], value);

			switch (argv[iarg][1]) {
				case 'B':
					read_basin_geometry(z_sill1, prm, hypso1_z, hypso1_a, &np1);
					break;
				case 'b':
					nbasins = 2;
					read_basin_geometry(z_sill2, prm, hypso2_z, hypso2_a, &np2);
					break;
				case 'e':
					e0 = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) e1 = atof(ptr);
					else e1=e0;
					if (ptr != NULL) e2 = atof(ptr);
					else e2=e1;
					break;
				case 'f':
					full_basins=1;
					break;
				case 'h':
					syntax(argc, argv);
					timeend = timeini+dt_default*5;
					break;
				case 'k':
					Ke = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) expe = atof(ptr);
					break;
				case 'm':
				case 'M':
					model_eros = atoi(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) model_vel = atoi(ptr);
					break;
				case 'O':
					z_sill1 = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) hl1 = atof(ptr);	else break;
					ptr=strtok(NULL, "/");
					if (ptr != NULL) dist1 = atof(ptr);	else break;
					break;
				case 'o':
					z_sill2 = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) hl2 = atof(ptr);	else break;
					ptr=strtok(NULL, "/");
					if (ptr != NULL) dist2 = atof(ptr);	else break;
					break;
				case 'P':
					switch_ps = 1;
					break;
				case 'p':
					p0 = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) p1 = atof(ptr);
					else p1=p2;
					if (ptr != NULL) p2 = atof(ptr);
					else p2=p1;
					break;
				case 'R':
					roughness = value;
					break;
				case 'r':
					r0 = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) r1 = atof(ptr);
					else r1=r0;
					ptr=strtok(NULL, "/");
					if (ptr != NULL) r2 = atof(ptr);
					else r2=r1;
					break;
				case 'S':
					read_basin_geometry(z_sill1, prm, hypso0_z, hypso0_a, &np0);
					break;
				case 't':
					timeini = atof(strtok(prm, "/"))*3600; /*[h] to [s]*/
					ptr=strtok(NULL, "/");
					if (ptr != NULL) timeend = atof(ptr)*3600; /*[h] to [s]*/
					ptr=strtok(NULL, "/");
					if (ptr != NULL) dtinput = atof(ptr)*3600; /*[h] to [s]*/
					else 		 dtinput=dt_default; /*[s]*/
					if (dtinput<0) dt = -dtinput;
					else 	       dt = dtinput;
					break;
				case 'V':
					verbose_level = value;
					break;
				case 'W':
					water_conservative=1;
					break;
				case 'w':
					Kw = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) expw = atof(ptr);
					else expw = 0;
					ptr=strtok(NULL, "/");
					if (ptr != NULL && *ptr == 'd') decreasable_width=1;
					if (ptr != NULL && *ptr == 'c') curvature_width=1;
					if (ptr != NULL && *ptr == 't') triangular_width=1;
					break;
				case 'z':
					z0 = atof(strtok(prm, "/"));
					ptr=strtok(NULL, "/");
					if (ptr != NULL) z1 = atof(ptr);
					else z1=hypso1_z[0];
					if (ptr != NULL) z2 = atof(ptr);
					else z2=hypso2_z[0];
					break;

				default:
					fprintf(stderr, "\nWarning in %s: incomprehensible parameter '%s' (argc=%d).", argv[0], argv[iarg], argc);
					exit(0);
					break;
			}
		}
	}

	z_sill1_ini = z_sill1;

	if (verbose_level>=4) for (i=0;i<np0;i++) fprintf(stderr, "\n#basin0  [%d]: z,a = %f m,  %e m2", i, hypso0_z[i], hypso0_a[i]);
	if (verbose_level>=4) for (i=0;i<np1;i++) fprintf(stderr, "\n#basin1  [%d]: z,a = %f m,  %e m2", i, hypso1_z[i], hypso1_a[i]);
	if (verbose_level>=4) for (i=0;i<np2;i++) fprintf(stderr, "\n#basin2  [%d]: z,a = %f m,  %e m2", i, hypso2_z[i], hypso2_a[i]);
	if (verbose_level>=1 && hypso0_a[0]!=0)      
		fprintf(stderr, "\nERROR: hypsometry of basin0 should have area 0 in first point: %.2f != 0 m", hypso0_a[0]);
	if (verbose_level>=1 && hypso1_a[0]!=0)      
		fprintf(stderr, "\nERROR: hypsometry of basin1 should have area 0 in first point: %.2f != 0 m", hypso1_a[0]);
	if (verbose_level>=1 && hypso2_a[0]!=0)      
		fprintf(stderr, "\nERROR: hypsometry of basin2 should have area 0 in first point: %.2f != 0 m", hypso2_a[0]);
	if (verbose_level>=2 && hypso0_z[np0-1]<z0)      
		fprintf(stderr, "\nERROR: basin0. Last hypsometry point should be > z0 altitude, but %.2f < %.2f m", hypso0_z[np1-1], z0);
	if (verbose_level>=2 && hypso1_z[np1-1]<z0)      
		fprintf(stderr, "\nERROR: basin1. Last hypsometry point should be at least z0 in altitude: %.2f < %.2f m", hypso1_z[np1-1], z0);
	if (verbose_level>=2 && hypso2_z[np2-1]!=z_sill2 && nbasins>=2) 
		fprintf(stderr, "\nERROR: basin2. Last point in hypsometry file should be at z_sill2 altitude: %.2f m", hypso2_z[np2-1]);
	
	/*Initial volume and area of basins*/
	volume_and_area_from_level(hypso0_z, hypso0_a, np0, z0, &volini0, &area);
	volume_and_area_from_level(hypso1_z, hypso1_a, np1, z1, &volini1, &area);
	volume_and_area_from_level(hypso2_z, hypso2_a, np2, z2, &volini2, &area);
	volume_and_area_from_level(hypso1_z, hypso1_a, np1, z0, &volini1full, &area);
	vol0=volini0; if (full_basins) {vol1=volini1*.995; vol2=volini2*.995;}

	if (verbose_level>=2) 
		fprintf(stderr, "\n======== %s's initial parameters ========"
			"\nHypsometry points: basin0= %d;  basin1= %d;  basin2= %d"
			"\nmodel eros,vel = %d (0 for shear stress; 1 for V in incision law) %d (0 for Manning; 1 for critical)"
			"\ntimeini,end,dt = %.1f , %.1f , %.2f h"
			"\nz0, z1, z2     = %.2f , %.2f , %.2f m"
			"\nz_sill1,2      = %.2f , %.2f m"
			"\nheadloss1,2    = %.2f , %.2f m"
			"\ndist1,dist2    = %.0f , %.0f m"
			"\nKe,expe        = %.3e m/y/Pa^-a (for m=0), %.2f"
			"\nnbasins        = %1d"
			"\nvolini0,1,2    = %.3e , %.3e , %.3e m3"
			"\ne0,e1,e2       = %.3f , %.3f , %.3f m/y"
			"\np0,p1,p2       = %.3f , %.3f , %.3f m/y"
			"\nr0,r1,r2       = %.2f , %.2f , %.2f m3/s"
			"\nKw,expw        = %.2f , %.2f"
			"\ndecr,curv,tr?  = %d   , %d, %d"
			"\nroughness      = %.2f"
			"\n", 
			argv[0], 
			np0, np1, np2, 
			model_eros, model_vel, timeini/3600, timeend/3600, dt/3600, z0, z1, z2, 
			z_sill1, z_sill2, hl1, hl2, dist1, dist2, Ke, expe, 
			nbasins, volini0, volini1, volini2, 
			e0, e1, e2, p0, p1, p2, r0, r1, r2, Kw, expw, 
			decreasable_width, curvature_width, triangular_width, roughness);
	
	time=timeini; 

	fprintf (stdout, CAPTION); 


	do {
		double dvolume, area0, area1, area2, widthlaw, 
			Rh1, slope1, areasill1, disch1, disch1ant, erosrate1, Dsill1, width1, shear1, vel1;
		double dvol0, dvol1, dvol2, dz;
		/*Cross-sectional area of the sill gate. Note z0-z_sill1 is the MEAN water depth of sill1*/
		Dsill1 = MAX_2(z0-z_sill1,0);
		if (time==timeini && Dsill1) width1 = 5*Dsill1; 
		areasill1 = width1*Dsill1;
		/*Hydraulic radius*/
		if (Dsill1>0) Rh1 = areasill1/(2*sqrt(4*Dsill1*Dsill1+width1*width1/4)); else Rh1=0;
		slope1 = (z1<z0)? (MAX_2(z1-z0,hl1))/dist1 : 0;
		/*Critical flow or Manning formula*/
		if (model_vel) 	vel1 = sqrt(g*Dsill1);
		else 		vel1 = 1/roughness*pow(Rh1,(double)2/3)*pow(-slope1,.5); 
		if (volini0 && vol0<=0) {vel1=0;}


		/*Water flow across sill1*/
		disch1 = vel1 * areasill1;
		if (dtinput<0) {
			/*variable dt*/
			float dt_disch, dt_eros;
			dt_disch = (MIN_2(volini0,volini1full)+1e0)/(disch1+1e-1)/1e4;
			dt_eros = fabs((z_sill1_ini-hypso0_z[0]+1e-2)/((erosrate1+1e-4)/secsperyr)/2e3);
//fprintf (stdout, "\n>>>>>>>%f %e %e %e %e %e", timeend, z_sill1, hypso0_z[0], erosrate1, dt_disch, dt_eros);
			if (z_sill1<=hypso0_z[0]) dt = dt_disch; else dt = MIN_2(dt_disch, dt_eros);
			dt *= 1/(1+1e2/(niters+1));
			if (time==timeini) dt=dt_default/*s*/;
//			dt = MAX_2(dt, 1e-5);
/*!!*/			if (time>timeini) dt = MIN_2(dt, (time-timeini)/(niters+1)*100);
		}
		dvolume = disch1*dt;

		/*TRANSFER WATER BETWEEN BASINS TO DEDUCE LEVELS z0,z1,z2, AND AREAS*/
		if (volini0) {
			dvolume = MAX_2(0,dvolume); 
			dvolume = MIN_2(vol0,dvolume); 
			disch1 = dvolume/dt;
			vol0 -= dvolume; 
		}
		volumetr += dvolume; 
		if (nbasins==1 || z1<=z_sill2 || z2>=z_sill2) {
			/*Keep water in basin1*/
			vol1 += dvolume;
		}
		else {
			/*Transfer water to basin2*/
			vol2 += dvolume;
		}
		/*Find level and area by filling vol*/
		if (volini0) level_and_area_from_volume(hypso0_z, hypso0_a, np0, vol0, &z0, &area0);
		level_and_area_from_volume(hypso1_z, hypso1_a, np1, vol1, &z1, &area1);
		level_and_area_from_volume(hypso2_z, hypso2_a, np2, vol2, &z2, &area2);

		if (z2>=z_sill2) z2 = z1;


		/*ADD RUNOFF+PRECIPITATION-EVAPORATION TO WATER FLOW*/
		dvol0 = r0*dt + p0*area0*dt/secsperyr - e0*area0*dt/secsperyr; 
		if (water_conservative) {dvol0 += e1*area1*dt;}
		dvol1 = r1*dt + p1*area1*dt/secsperyr - e1*area1*dt/secsperyr; 
		if (z1<hypso2_z[np2-1]) {
			dvol2 = r2*dt + p2*area2*dt/secsperyr - e2*area2*dt/secsperyr; 
			if (water_conservative) {dvol0 += e2*area2*dt/secsperyr;}
		} 
		else {dvol1+=r2*dt; dvol2=0;}

		vol0 += dvol0; vol0=MAX_2(vol0,0); 
		vol1 += dvol1; vol1=MAX_2(vol1,0); 
		vol2 += dvol2; vol2=MAX_2(vol2,0);

		if (volini0 && vol0<=0) {z0=z_sill1;}

		/*Erosion*/
		if (model_eros==2) 
			/*Darcy-Weisbach equation; Chezy's constant*/
			if (Rh1>0) shear1 = denswater*g / (Cz*Cz) /*g * pow(roughness,(double)2.)/pow(Rh1,(double)2./(double)6.)*/ * vel1*vel1; else shear1=0;
		else
			shear1 = denswater*g*Dsill1*(-slope1); 
		erosrate1 = Ke*pow((shear1)*((model_eros==1)?vel1:1), expe); 
		if (z_sill1<=hypso0_z[0]) erosrate1=0;
		//fprintf(stderr, "\nvel1=%.1e Rh1= %.1e  shear = %.1e erosrate=%.1e\n", vel1, Rh1, shear1, erosrate1);
//if (erostotal>75) erosrate1 /= 10;
//if (shear1<50) erosrate1 = 0;
		fprintf (stdout, 
			"\n%.5f\t%6.3f\t%6.1f\t%6.4f\t%6.2f\t%6.2e\t%6.3f\t%8.5f" 
			"\t%7.2e\t%7.2e\t%7.2e\t%7.2e\t%6.3f\t%6.3f\t%6.3f", 
			time/3600, z_sill1, Rh1, slope1, vel1, disch1, width1, erosrate1, 
			vol0/1e9, vol1/1e9, vol2/1e9, volumetr/1e9, z0, z1, z2);

		widthlaw = Kw*pow(fabs(disch1), expw); 
		if (curvature_width) widthlaw = pow(8*Kw*(z0-z_sill1), .5); 
		if (triangular_width) widthlaw = MAX_2(0, Kw*(z0-z_sill1)); 
		if (decreasable_width) width1 = widthlaw; else width1 = MAX_2(width1, widthlaw); 

		z_sill1   += -erosrate1*dt/secsperyr;
		erostotal +=  erosrate1*dt/secsperyr;
		//z_sill1 = MAX_2(z_sill1, hypso0_z[0]);
		time += dt;
//		if (timeend<timeini && niters>900 && disch1<disch1ant) {timeend=time+(time-timeini)*1.618033;}
//		if (timeend<timeini && niters>900 && (z_sill1-hypso0_z[0]<1e-3 || erosrate1<1e-6)) timeend=time+vol0/disch1*3; 
//fprintf(stderr, "\ntimeend=%.1e %.1e %.1e %.1e h\n", time/3600, timeend/3600, timeini/3600, dt/3600);
		if (niters>1e6) {fprintf (stdout, "\n\aERROR: TOO MANY ITERATIONS!"); break;}
		disch1ant = disch1;
		niters++;
	} while (time<=timeend || timeend<timeini);

	fprintf (stdout, "\n"CAPTION); fflush(stdout);
	if (verbose_level>=2) {
		fprintf (stderr, "\ntimeend = %.1e at time=%.1e dt= %.1e h; %d iters", timeend/3600, time/3600, dt/3600, niters);
		fprintf (stderr, "\nErosion total = %.1f m\n", erostotal);
		fprintf (stderr, "\n");
	}
	exit(1);
}



int level_and_area_from_volume (float *hypso_z, float *hypso_a, int np, double vol, double *z, double *area) {
	int i;
	double vol_aux=0; 
	if (vol<0) fprintf(stderr, "\nERROR in level_and_area_from_volume: negative volume passed: %.2e m3", vol);
	for (i=1; i<np; i++) {
		double dz, dvol;
		dvol = (hypso_z[i]-hypso_z[i-1])*(hypso_a[i]+hypso_a[i-1])/2;
		//fprintf(stderr, "\n************** %.2f %.2e   %.2f  %e     %d, %d     %.3e  %.3e", hypso_z[i], hypso_a[i], *z, *area, i, np, vol, vol_aux+dvol);
		if (vol_aux+dvol >= vol) {
			double a = (hypso_a[i]-hypso_a[i-1])/(hypso_z[i]-hypso_z[i-1])/2, b=hypso_a[i-1], c=-(vol-vol_aux);
			if (a==0) {
				dz = -c/b;
			}
			else {
				if (b==0)
					dz = sqrt(-c/a);
				else
					dz = (-b+sqrt(b*b-4*a*c))/(2*a);
			}
			*z = hypso_z[i-1] + dz;
			*area = b+dz*a;
			break;
		}
		vol_aux += dvol;
	}
	/*Assumes constant area above last point in hypso_z*/
	if (i==np) {
		*z = hypso_z[np-1] + (vol-vol_aux)/hypso_a[np-1];
		*area = hypso_a[np-1];
	}
//fprintf(stderr, "\n!!!!!!!!!!!!!! %.2f %.2e   %.2f  %e", hypso_z[i], hypso_a[i], *z, *area);
	return (1);
}




int volume_and_area_from_level (float *hypso_z, float *hypso_a, int np, double z, double *vol, double *area) {
	int i;
	*vol=0;
	if (z<hypso_z[0]) fprintf(stderr, "\nERROR in volume_and_area_from_level: passed level below basin floor: %.2f m", z);
	for (i=1; i<np; i++) {
		if (hypso_z[i]<=z) {
			*area = hypso_a[i];
			*vol    += (hypso_z[i]-hypso_z[i-1])*(hypso_a[i]+hypso_a[i-1])/2; 
		}
		else {
			*area = hypso_a[i-1]+((hypso_a[i]-hypso_a[i-1])/(hypso_z[i]-hypso_z[i-1]))*(z-hypso_z[i-1]);
			*vol    += (z-hypso_z[i-1])*(*area+hypso_a[i-1])/2; 
			break;
		}
	}
	return (1);
}


int read_basin_geometry(double z_sill, char *prm, float *hypso_z, float *hypso_a, int *np) {
	int i;
	double area, z_bott;
	char 	*ptr, filename[MAXLENLINE];
	FILE	*file = NULL;

	sscanf(prm, "%s", filename);
	area = atof(strtok(prm, "/"));
	ptr=strtok(NULL, "/");
	if (ptr != NULL) z_bott = atof(ptr);
	ptr=strtok(NULL, "/");
	if (ptr != NULL) ;
	else {
		if ((file = fopen(filename, "rt")) == NULL) {
			if (verbose_level>=3) fprintf(stderr, "\nWarning: Input basin hypsometry file '%s' not found.", filename);
		}
		else {
			if (verbose_level>=3) fprintf(stderr, "\nHypsometry file: '%s'.", filename);
			for (i=0;;i++) {
				TAKE_LINE_2(hypso_z[i], hypso_a[i])
			}
			*np=i;
		}
		return(1);
	}
//fprintf(stderr, "\n>>>>> '%c'  %e %f %s", *ptr, area, z_bott, filename);
	if (*ptr == 'b') {
		/*Default box-like hypsometry*/
		hypso_z[0]=z_bott; 	hypso_a[0]= 0;
		hypso_z[1]=z_bott; 	hypso_a[1]= area;
		hypso_z[2]=z_sill; 	hypso_a[2]= area;
		hypso_z[3]=z_sill+100e3; 	hypso_a[3]= area;
		*np = 4;
		if (verbose_level>=3) fprintf(stderr, "\nAutomatic box-like hypsometry points for basin: %d   type: %c", *np, *ptr);
	}
	else {
		/*Default linear hypsometry*/
		hypso_z[0]=z_bott; 	hypso_a[0]=0;
		hypso_z[1]=z_sill; 	hypso_a[1]= area;
		hypso_z[2]=z_sill+100e3; 	hypso_a[2]= area;
		*np = 3;
		if (verbose_level>=2) fprintf(stderr, "\nAutomatic linear hypsometry points for basin: %d   type: %c", *np, *ptr);
	}
	return(1);
}
	


int syntax (int argc, char **argv) {
	fprintf(stderr, "\n"
			"\nSyntax: %s "
			"\n\t-S[<area0>/<z_bott0>/<b|l>|<hypsometry0>] -O<z_sill1>[/<hl1>/<dist1>] "
			"\n\t-B[<area1>/<z_bott1>/<b|l>|<hypsometry1>] -o<z_sill2>[/<hl2>/<dist2>] "
			"\n\t-b[<area2>/<z_bott2>/<b|l>|<hypsometry2>] "
			"\n\t-e<e0>/<e1>/<e2> [-h] -k<Ke>/<expe> -M<model_eros>[/<model_vel>] -p<p0>/<p1>/<p2> "
			"\n\t-R<roughness> -r<r0>/<r1>/<r2> -t<timeini>/<timeend>/<dt> "
			"\n\t[-V<level>] -W -w<Kw>/<expw>/[d|c] -z<z0>/<z1>/<z2> ", argv[0]);
	fprintf(stderr, "\n"
			"\nThis program calculates the water transfer from a source basin0 overflowing into"
			"a receiving basin1, and if -b is specified also from basin1 to a basin2. Water "
			"flows only in that sense 0->1->2. The erosion produced at sill1 (between basin0 "
			"and basin1) due to water flow is also calculated. "
			"\n\nNotation/legend: "
			"\nz is elevation (growing positive upwards) [m]. "
			"\nbasin0, basin1, and basin2 refer to the source basin, the upper flooded basin, and the lower flooded basin. "
			"\n"
			"\n\t-S indicates that the source basin0 is not infinite (default) but it has either a triangular (l) or box-like (b) hypsometry defined by <area0> and <z_bott0> or follows the hypsometry0 file."
			"\n\t-O provides geometry of outlet from basin0 to basin1. Sill elevation (must be below z0), and, only if model_vel=0 or model_eros!=2, headloss [m] (<0) and slope distance [m]. Headloss is also the maximum effective difference in water level between basins. "
			"\n\t-B provides basin1 geometry. If no hypsometry file: area [m3], bottom z of basin1 [m], and 'b' or 'l' for box or triangular hypsometry. If 2 basins, area1 includes the whole basin above sill2."
			"\n\t-o provides geometry of outlet from basin1 to basin2. Format as -O. So far, this sill is NOT eroded as sill1, and hl2,dist2 have no effect. "
			"\n\t-b to account for a subsidiary basin2. Format as -B. "
			"\nz_sill1, z_sill2 are the initial sill elevations [m]."
			"\nhl1, hl2 are the initial difference in altitude (headloss) between basins 0 and 1 and basins 1 and 2 [m]."
			"\ndist1, dist2 are the distances over which that headloss takes place [m]."
			"\n"
			"\n\t-e to give the surface evaporation at each basin lake [m/y]."
			"\n\t-f to make basins 1 and 2 full since the beggining."
			"\n\t-h for help."
			"\n\t-k to specify the two constants of the erosion law: for model_eros=0 Ke is in units [m/yr/Pa^a]; for model_eros=1 Ke is in [m/yr/(m/s*Pa)^a]"
			"\n\t-M to specify the erosion <model_eros> (0 e=Ke*tau^expe for tau=rho*g*D*S; 1 for stream power per unit area: e=Ke*(tau*V)^expe; 2 as 0 but tau using Darcy-Weisbach equation with Chezy's coeff). <model_vel> is the hydraulic model (0 for Manning's eq.; 1 for critical flow)."
//			"\n\t-P to call the gmt script and produce a postscript with graphics. "
			"\n\t-p to give the precipitation fallen on basins [m/y]."
			"\n\t-R to change the default coefficient of roughness [adimensional]."
			"\n\t-r to give the runoff collected by rivers to basins [m3/s]."
			"\n\t-t sets initial and final run time [hours]. <dt> is the time step. A negative timeend will find an automatic value after discharge values start decreasing. Negative dt gives approximate number of automatically adjusted time steps."
			"\n\t-u sets the sill uplift rate [m/y]. Current sea level rise is 3mm/yr. Last postglacial maximum was 12 mm/yr."
			"\n\t-W to drop the water evaporated in basins 1 and 2 at basin 0 (water conservative)."
			"\n\t-w to specify the two constants of the width law W=Kw*Q^expw [-]. Add 'd' to allow decrease of width during decreasing discharges. Add 'c' for Kw to be the radius of curvature of the outlet section. Add 't' for triangular section with width=Kw*water_depth"
			"\n\t-z to specify the water level in each basin [default is rim of basin0 and bottom for basin1 and basin2]."
			"\n"
			"\nIncision law: e=Ke*tau^expe (for model_eros=0) ; e=Ke*(tau*V)^expe (for model_eros=1), where V is the water velocity. Note that for model_eros=1, expe is equivalent to n (exponent of slope) in the stream power law of river incision. For model=0, expe=3/2*n. So, according to Whipple & Tucker (1999) expe=[2/3-2] for model=1 and expe=[1-3] for model=0."
			"\n"
			"\nReads hypsometry from files <hypsometry0> (source basin), <hypsometry1> (basin1), <hypsometry2> (basin2), in two columns (z, basin_area_above_z). z runs from basin bottom (first line) to the highest z0 [m]. First area [m2] is thus 0. "
			"\nThe source basin0 has an infinite hypsometry by default (ocean)."
			"\nArea of basin1 jumps to the total hypsometry above the sill with basin2. "
			"\nArea of basin2 zeroes when reaching the sill with basin1. "
			"\n"
			"\nWrites in stdout: "
			CAPTION
			"\nThese are: time [hours], sill channel elevation [m], hydraulic radius [m], slope, "
			"\nwater velocity, discharge, channel width, incision rate, transferred volume, basin level");
	fprintf(stderr, "\n-V provides additional verbose information such as the polygon area [1-4].\n");
	fprintf(stderr, "\n"
			"\n#Examples: "
			"\n#Mediterranean, three basins (data from Blanc, 2002):"
			"\nspillover -k8e-6/1.5 -S3.6e14/-5000/l -O-10/-1000/100000 -B1.0e12/-1509/l -o-430/-1000/100000 -b1.51e12/-1718/l -t0/.6/.0001\n"
			"\n"
			"\n#Mediterranean, one basin:"
			"\nspillover -O-10/-1000/100000 -B2.51e12/-1509/l -t0/2.5/.01"
			"\n"
			"\n#Black Sea:"
			"\nspillover -k8e-5/1.5 -O-2/-1000/100000 -B4.4e11/-100/b -t0/200/.1"
			"\n"
			"\nBelow follows a list of default parameters and the corresponding first time steps of an example run:"
			"\n");
	AUTHORSHIP;
	return(0);
}
