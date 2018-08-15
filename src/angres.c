/*--------------------------------------------------------------------
 * $Id: angres.c 3254 2017-02-05 14:47:34Z bernhard.mayer $
 * 
 * This file is part of libRadtran.
 * Copyright (c) 1997-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * ######### Contact info: http://www.libradtran.org #########
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License   
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.        
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   
 * GNU General Public License for more details.                    
 * 
 * You should have received a copy of the GNU General Public License          
 * along with this program; if not, write to the Free Software                
 * Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA.
 *--------------------------------------------------------------------*/

/* 03.04.2002: Fixed problem with finding correct grid cell for direct beam, AKy. */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef _AIX
  #include <getopt.h>
#else
  #include <unistd.h>
#endif

#if HAVE_LIBGSL 
  #include <gsl/gsl_math.h>
  #include <gsl/gsl_monte.h>
  #include <gsl/gsl_monte_plain.h>
  #include <gsl/gsl_monte_miser.h>
  #include <gsl/gsl_monte_vegas.h>
#endif

#include "ascii.h"

#define PROGRAM "angres"
#define VERSION "0.3"
#define DATE    "3 Feb, 2003"

#define DIM_ANG_INT 2  /* Integration to be made over 2 dimensions */

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif

typedef struct {
  int nphi;
  int numu;
  double *umu;
  double *phi;
  double **prod;
} mcinput;


/*************************************************************/
/* print usage information                                   */
/*************************************************************/

/*
  Documentation in latex for the User's Guide. Extracted by sdoc2.awk 
  and input in the doc/tools.tex.
*/
/*
  <lpdoc>

\subsection{Angular response and tilted surfaces - \codeidx{angres}}
The \code{angres} tool takes a precalculated radiance field and integrates 
it over a given angular area using any angular response. Typical usages of
\code{angres} are calculation of radiation on tilted surfaces and estimation 
of effects of imperfect angular response functions.

The \code{angres} tool is invoked as follows:
\begin{Verbatim}[fontsize=\footnotesize]
 angres angres_file raddis_file
\end{Verbatim}
The two required input files will be read by the \code{angres} tool.
\begin{description}
\item[angres\_file] is a two column file with the first column holding
     the angle and the second column the angular response, e.g. a
     measured cosine response. To generate standard angular response 
     function see the \code{make\_angres\_func} tool.

\item[raddis\_file] holds the radiance distribution as output from uvspec
     with the disort solvers for one single wavelength.
\end{description}
After reading the two input files the angular response will be tilted
and rotated if specified with the \code{-t} and \code{-r} options respectively. Finally
the product of the resulting angular response and radiance distribution
field are integrated using Monte Carlo methods to yield the effective
response. The integration is done for the diffuse radiation field only. To
include the direct contribution the -s and -z options must be set to give
the direction of the sun.

Output is 3 numbers:
\begin{enumerate}
\item The integral of the diffuse radiation field times angular response.
\item Estimated absolute error of the above integral.
\item The integral of the diffuse+direct radiation field times angular response
       (requires that \code{-s} and \code{-z} are specified, otherwise same as first number
\end{enumerate}

The angles in the \file{angres\_file} must be in radians if not
the \code{-a} option is used. The \file{raddis\_file} must contain
output from \code{uvspec} run for one single wavelength with one of
the disort solvers and with \code{phi} and \code{umu} set.
Note that the angles in the \file{angres\_file} must follow the same
conventions as for the disort algorithm. This is different from
that typically used when reporting measurements of the angular response.

The \code{angres} tool accepts the following command line options:

\begin{description}
\item[-h]  show this page.
\item[-c]  number of random points used for Monte Carlo integration.
\item[-i]  The diffuse radiation is assumed to be isotropic.
\item[-a]  angular response angle given in degrees and not cosine of angle.
\item[-r]  rotation angle in degrees.
\item[-t]  tilt angle in degrees.
\item[-s]  solar zenith angle in degrees.
\item[-z]  solar azimuth angle in degrees.
\item[-p]  pgm files are made of the angular response before
      and after tilt and rotate.
\end{description}

Sample \code{angres} input and output files are found in the \file{examples} directory.
The following
\begin{Verbatim}[fontsize=\footnotesize]
 angres examples/ANGRES_1_ANG.DAT  \
        examples/ANGRES_RADDIS_1.DAT -a -t -r 0 -s 32 -z 0
\end{Verbatim}
calculates the radiation on a horisontal surface given the angular response 
in \file{examples/ANGRES_1_ANG.DAT}. The input used to calculate the radiance file is
given in the start of \file{examples/ANGRES_RADDIS_1.DAT}.

An example of the use of \code{angres} together with \code{uvspec} is given in 
\citet[section 4.6]{mayer2005}.
  </lpdoc>
*/

static void print_usage (char *filename)
{
  fprintf (stderr, "%s - %s   \n\n", PROGRAM, VERSION);
  fprintf (stderr, "written by  Arve Kylling,\n");
  fprintf (stderr, "            Norwegian Institute for Air Research (NILU)\n");
  fprintf (stderr, "            P.O. Box 100, 2027 Kjeller, Norway\n");
  fprintf (stderr, "Version %s finished %s\n\n", VERSION, DATE);
  fprintf (stderr, "Be aware that this program is a beta version! Please report any\n");
  fprintf (stderr, "kind of error to the author, including the error message and the\n");
  fprintf (stderr, "database being processed when the error occured. Thanks!\n\n");

  fprintf (stderr, "USAGE: %s [options] angres_file raddis_file \n\n",
	   filename);
  fprintf (stderr, "%s will read the two required input files:\n", PROGRAM);
  fprintf (stderr, "  angres_file is a two column file with the first column holding\n");
  fprintf (stderr, "     the angle and the second column the angular response, e.g. a\n");
  fprintf (stderr, "     measured cosine response.\n");
  fprintf (stderr, "  raddis_file holds the radiance distribution as output from uvspec\n");
  fprintf (stderr, "     with the disort solvers for one single wavelength.\n");
  fprintf (stderr, "After reading the two input files the angular response will be tilted\n");
  fprintf (stderr, "and rotated if specified with the -t and -r options respectively. Finally\n");
  fprintf (stderr, "the product of the resulting angular response and radiance distribution\n");
  fprintf (stderr, "field are integrated using Monte Carlo methods to yield the effective\n");
  fprintf (stderr, "response. The integration is done for the diffuse radiation field only. To\n");
  fprintf (stderr, "include the direct contribution the -s and -z options must be set to give\n");
  fprintf (stderr, "the direction of the sun. \n");
  fprintf (stderr, "\n");
  fprintf (stderr, "Output is 3 numbers:\n");
  fprintf (stderr, "   1: The integral of the diffuse radiation field times angular response.\n");
  fprintf (stderr, "   2: Estimated absolute error of the above integral.\n");
  fprintf (stderr, "   3: The integral of the diffuse+direct radiation field times angular response\n");
  fprintf (stderr, "       (requires that -s and -z are specified, otherwise same as 1).\n");
  fprintf (stderr, "\n");
  fprintf (stderr, "The angles in the angres_file must be in radians if not\n");
  fprintf (stderr, "the -a option is used. The raddis_file must contain\n");
  fprintf (stderr, "output from uvspec run for one single wavelength with one of\n");
  fprintf (stderr, "the disort solvers and with phi and umu set.\n");
  fprintf (stderr, "Note that the angles in the angres_file must follow the same\n");
  fprintf (stderr, "conventions as for the disort algorithm. This is different from\n");
  fprintf (stderr, "that typically used when reporting measurements of the angular response.\n");
  fprintf (stderr, "\n");
  fprintf (stderr, " command line options:\n");
  fprintf (stderr, "  -h  show this page\n");
  fprintf (stderr, "  -c  number of random points used for Monte Carlo integration\n");
  fprintf (stderr, "  -i  The diffuse radiation is assumed to be isotropic\n");
  fprintf (stderr, "  -a  angular response angle given in degrees and not cosine of angle\n");
  fprintf (stderr, "  -r  rotation angle in degrees\n");
  fprintf (stderr, "  -t  tilt angle in degrees\n");
  fprintf (stderr, "  -s  solar zenith angle in degrees\n");
  fprintf (stderr, "  -z  solar azimuth angle in degrees\n");
  fprintf (stderr, "  -p  pgm files are made of the angular response before\n");
  fprintf (stderr, "      and after tilt and rotate\n");
}
  
static int get_options (int argc, char **argv,
			char *programname,
			char *angresfilename,
			char *raddisfilename,
			int *pgm,
			int *dev,
			int *ang,
			size_t *calls,
			double *sza,
			double *phi0,
			double *theta_t,
			double *phi_r,
			int *incdir,
			int *isotropic
			)
{
  int c=0;

  /* default settings */

  /* save name of program */
  strncpy (programname, argv[0], FILENAME_MAX);
  
  while ((c=getopt (argc, argv, "ahipc:d:s:t:r:z:")) != EOF)  {
    switch(c)  {
    case 'h':  /* help */
      print_usage (programname);
      return (-1);
      break;
    case 'd': 
      *dev = atoi(optarg);
      break;
    case 'a': 
      *ang = 1;
      break;
    case 'i': 
      *isotropic = 1;
      break;
    case 'p': 
      *pgm = 1;
      break;
    case 'c': 
      *calls = atoi(optarg);
      break;
    case 'r': 
      *phi_r = atof(optarg);
      break;
    case 's': 
      *sza = atof(optarg);
      *incdir=1;
      break;
    case 'z': 
      *phi0 = atof(optarg);
      *incdir=1;
      break;
    case 't': 
      *theta_t = atof(optarg);
      break;
    case '?':
      print_usage (programname);
      return (-1);
      break;
    default:
      print_usage (programname);
      return (-1);
    }
  }
      
  /* check number of remaining command line arguments */
  if (argc - optind != 2)  {
    print_usage (programname);
    return -1;
  }
  
  /* save command line arguments */
  strncpy (angresfilename,    argv[optind+0], FILENAME_MAX);
  strncpy (raddisfilename,    argv[optind+1], FILENAME_MAX);

  return 0;  /* if o.k. */
}

double  **tilt_and_rotate(double *aangres, 
			  int nphi, double *phi, int numu, double *umu, int pgm,
			  double theta_t, double phi_r, int dev) {
  int status = 0, i=0, j=0, k=0, it=0, jt=0, ia=0, ja=0, nt;
  /* double umu_t; */
  double theta=0, x=0, y=0, z=0, xp=0, yp=0, zp=0, r=0;
  double  a11, a12, a13, a21, a22, a23, a31, a32, a33;
  double rn=0, theta_n=0, phi_n=0;
  double **angres=NULL, **angrestilted=NULL;
  double **angresindcnt=NULL;
  double delphirad=0;
  FILE *pgmfile = NULL;

  /* allocate memory for the angular response function and the tilted and rotated one */
  status = ASCII_calloc_double(&angres, numu, nphi);
  if (status!=0) { 
    fprintf (stderr, "tilt_and_rotate: error allocating memory for angular response function\n");
    /*    return status; */
  }
  status = ASCII_calloc_double(&angrestilted, numu, nphi);
  if (status!=0) { 
    fprintf (stderr, "tilt_and_rotate: error allocating memory for angular response function\n");
    /*    return status; */
  }
  
  status = ASCII_calloc_double(&angresindcnt, numu, nphi);
  if (status!=0) { 
    fprintf (stderr, "tilt_and_rotate: error allocating memory for angular response function\n");
    /*    return status; */
  }
  
  /* The read angular response is assumed to have no azimuth dependence, */
  /* hence it is the same for all phi                                    */
  for (i=0;i<numu;i++)
    for (j=0;j<nphi;j++)
      angres[i][j] = aangres[i];

  if (pgm) {
    if ((pgmfile = fopen("angres_orginal.pgm", "w"))==NULL) {
      fprintf (stderr, "Error opening angres_orginal.pgm for writing\n");
      return NULL;
    }
    
    fprintf(pgmfile,"P2\n");
    fprintf(pgmfile,"%d %d\n", nphi, numu);
    fprintf(pgmfile,"256\n");
    for (i=0;i<numu;i++) {
      for (j=0;j<nphi;j++) {
	k = (int) (angres[i][j]*256);
	fprintf(pgmfile,"%d\n", k);
      } 
    }
    fclose(pgmfile); 
  }

  if (dev==1)  
    fprintf (stdout, "Tilt: %f Rotate: %f\n", theta_t, phi_r);
  /* commented out by RB, because was not used but caused compiler warnings */
  /* umu_t = cos(theta_t*PI/180.0); */

  /* Rotation matrix. */
  /* Eq 4-47., p. 147, Goldstein, Classical Mechanics. psi=0 for our purposes */
      
  /*  a11 =  cos(phi_r*PI/180.0);                      */
  /*  a12 = -cos(theta_t*PI/180.0)*sin(phi_r*PI/180.0); */
  /*  a13 =  sin(theta_t*PI/180.0)*sin(phi_r*PI/180.0); */
  /*  a21 =  sin(phi_r*PI/180.0);                      */
  /*  a22 =  cos(theta_t*PI/180.0)*cos(phi_r*PI/180.0); */
  /*  a23 = -sin(theta_t*PI/180.0)*cos(phi_r*PI/180.0); */
  /*  a31 = 0;                                        */
  /*  a32 =  sin(theta_t*PI/180.0);                    */
  /*  a33 =  cos(theta_t*PI/180.0);                    */

  /* First only tilt */
  a11 =  1;
  a12 =  0;
  a13 =  0;
  a21 =  0;
  a22 =  cos(theta_t*PI/180.0);
  a23 = -sin(theta_t*PI/180.0);
  a31 =  0;
  a32 =  sin(theta_t*PI/180.0);
  a33 =  cos(theta_t*PI/180.0);

  if (dev==1) {
    printf("a11: %7.3f, a12: %7.3f a13: %7.3f\n", a11, a12, a13); 
    printf("a21: %7.3f, a22: %7.3f a23: %7.3f\n", a21, a22, a23); 
    printf("a31: %7.3f, a32: %7.3f a33: %7.3f\n", a31, a32, a33); 
  }
  if (dev==1)  
    fprintf(stdout, "Tilt: %f Rotate: %f\n", theta_t, phi_r);

  /* First tilt the angular response */
  for (i=0;i<numu;i++) {
    theta = acos(umu[i]);
    for (j=0;j<nphi;j++) {
      rn = 1;
      xp = rn*cos(phi[j]*PI/180)*sin(theta);
      yp = rn*sin(phi[j]*PI/180)*sin(theta);
      zp = rn*cos(theta);
      if (dev==1) {
	printf("\nxp:  %7.3f, yp:  %7.3f zp:  %7.3f theta: %7.4f phi: %5.1f r: %6.3f umu: %f\n", 
	       xp, yp, zp, theta, phi[j],r, umu[i]); 
      }

      if (dev==1)  
	fprintf (stdout,"Tilt: %f Rotate: %f\n", theta_t, phi_r);
      
      x = a11*xp+a12*yp+a13*zp;
      y = a21*xp+a22*yp+a23*zp;
      z = a31*xp+a32*yp+a33*zp;

      rn = sqrt(x*x+y*y+z*z);
      if (rn > 0) 
	theta_n = acos(z/rn);              else theta_n = 0;

      if (rn > 0) {
	phi_n   = acos(x/(rn*sin(theta_n)));
      }
      else {
	phi_n = 0;
      }

      if (dev==1) {
	printf("x: %7.3f, y: %7.3f z: %7.3f theta: %7.4f phi: %5.1f\n", x, y, z, theta*180/PI, phi[j]); 
	printf("rn: %7.3f, theta_n: %7.4f phi_n: %5.1f phi_r: %f\n", rn, theta_n*180/PI, phi_n*180/PI, phi_r); 
      }

      if (dev==1)  
	fprintf (stdout,"Tilt: %f Rotate: %f\n", theta_t, phi_r);

      it = 0;
      ia = 0;

      while ((theta_n*180/PI-acos(umu[ia++])*180/PI) < 0 && ia < numu); 
      it=ia-1; /* This one is fun:-) */
      
      /* Make sure we get the closest grid point */
      if (fabs(acos(umu[it])*180/PI-theta_n*180/PI) > fabs(acos(umu[it-1])*180/PI-theta_n*180/PI))
	it=it-1;


      jt = 0;
      ja = 0;
      delphirad = fabs(phi[1]-phi[0])*PI/180;
      while (fabs(phi_n-phi[ja++]*PI/180) > delphirad) {
	if (dev==1) {
	  printf("%d %d %d %f %f %f %f\n", j, ja, jt, phi_n, phi[ja-1]*PI/180, fabs(phi_n-phi[ja-1]*PI/180),
	    delphirad);
	}
      } 
      jt=ja-1;
      angrestilted[i][j] = angres[it][jt];      
    }
  }

  /* Now we rotate */
  nt = (phi_r*PI/180)/delphirad;
  if (jt >= nphi) jt-=nphi;
  for (i=0;i<numu;i++) {
    theta = acos(umu[i]);
    for (j=0;j<nphi;j++) {
      jt = j+nt;
      if (jt >= nphi) jt-=nphi;
      if (dev==1 || dev==10) fprintf(stdout,"%3d %3d nt: %d %3d jt: %3d %3d, nphi: %d, phi_r: %5.1f %3.1f\n", 
				      i,j,nt, it, jt, ja-1, nphi, phi_r, angres[it][jt]);

      angres[i][j] = angrestilted[i][jt];
      if (dev==1) {
	printf("angrestilted: %7.4f %7.4f %7.4f\n", acos(umu[it])*180/PI,phi[jt],angrestilted[i][j]); 
      }
    }
  }

  if (pgm) {

    if ((pgmfile = fopen("angres_tilted.pgm", "w"))==NULL) {
      fprintf (stderr, "Error opening angres_tilted.pgm for writing\n");
      return NULL;
    }

    fprintf(pgmfile,"P2\n");
    fprintf(pgmfile,"%d %d\n", nphi, numu);
    fprintf(pgmfile,"256\n");
    for (i=0;i<numu;i++) {
      for (j=0;j<nphi;j++) {
	k = (int) (angres[i][j]*256);
	fprintf(pgmfile,"%d\n", k);
      } 
    }
    fclose(pgmfile); 
  }

  return angres;

}

double linpol (double x1, double x2, double y1, double y2, double x) {
  /* Linearly interpolate between two points to get wanted y value for x. */
  double y;
  double a=0, b=0;

  if (x2 > 0 && x1 > 0 && fabs(x2-x1) < 0.0000001) {
    y = y1;
  }
  else {
    a  = (y2-y1)/(x2-x1);
    b  = y1 - a*x1;
    y  = a*x + b;
  }
  return y;
}


double mc_func (double *k, size_t dim, void *params)
{
  int i1=0, i2=0, j1=0, j2=0;
  double p1, p2, p3;
  double umu=0, phi=0;
  mcinput lmcinp;

  lmcinp =  *(mcinput *)params;

  umu = k[0];
  phi = k[1];
  i2  = 0;
  while (umu>lmcinp.umu[i2++]);
  i1 = i2-1;
  if (i2 ==0) {i2=1; i1=0;}
  if (i2 >= lmcinp.numu) {i2=lmcinp.numu-1; i1=i2-1;}

  j2  = 0;
  while (phi>lmcinp.phi[j2++]);
  j1 = j2-1;
  if (j2 ==0) {j2=1; j1=0;}
  if (j2 >= lmcinp.nphi) {j2=lmcinp.nphi-1; j1=j2-1;}

  p1 = linpol(lmcinp.phi[j1], lmcinp.phi[j2], lmcinp.prod[i1][j1], lmcinp.prod[i1][j2], phi);
  p2 = linpol(lmcinp.phi[j1], lmcinp.phi[j2], lmcinp.prod[i2][j1], lmcinp.prod[i2][j2], phi);

  p3 = linpol(lmcinp.umu[i1], lmcinp.umu[i2], p1, p2, umu);

  return p3;

}

double mc_funciso (double *k, size_t dim, void *params)
{
  int i1=0, i2=0, j1=0, j2=0;
  double p1, p2, p3;
  double umu=0, phi=0;
  mcinput lmcinp;

  lmcinp =  *(mcinput *)params;

  umu = k[0];
  phi = k[1];
  i2  = 0;
  while (umu>lmcinp.umu[i2++]);
  i1 = i2-1;
  if (i2 ==0) {i2=1; i1=0;}
  if (i2 >= lmcinp.numu) {i2=lmcinp.numu-1; i1=i2-1;}

  j2  = 0;
  while (phi>lmcinp.phi[j2++]);
  j1 = j2-1;
  if (j2 ==0) {j2=1; j1=0;}
  if (j2 >= lmcinp.nphi) {j2=lmcinp.nphi-1; j1=j2-1;}

  p1 = linpol(lmcinp.phi[j1], lmcinp.phi[j2], lmcinp.prod[i1][j1], lmcinp.prod[i1][j2], phi);
  p2 = linpol(lmcinp.phi[j1], lmcinp.phi[j2], lmcinp.prod[i2][j1], lmcinp.prod[i2][j2], phi);

  p3 = linpol(lmcinp.umu[i1], lmcinp.umu[i2], p1, p2, umu);

  return p3;

}


int main(int argc, char **argv) 
{
  int i=0, j=0, k=0, status=0, pgm=0, dev=0, ang=0, incdir=0, isotropic=0;
  int ia=0, it=0, ja=0, jt=0;
  double delphirad=0, dir=0;

  char programname     [FILENAME_MAX]="";
  char angresfilename  [FILENAME_MAX]=""; /* File holding angular response      */
  char raddisfilename  [FILENAME_MAX]=""; /* File holding radiance distribution */

  double *phi=NULL, *rphi=NULL, *aumu=NULL, *rumu=NULL, *aangres=NULL;
  double theta_t=0, phi_r=0, sza=0, phi0=0, max=0;
  double dphi=0, dtheta=0, tmpprod=0;
  double uavgso;
  /* commented out by RB, because was not used but caused compiler warnings */
  /* double rfldir, rfldn, flup, uavgso, uavgdn, uavgup; */

  int arows=0, amin_columns=0, amax_columns=0;  /* Angular response file      */
  int rrows=0, rmin_columns=0, rmax_columns=0;  /* Radiance distribution file */
  int nphi=0, numu=0;

  mcinput mcinp;
  /* void *pmcinp; */
  mcinput mcinpiso;
  /* void *pmcinpiso; */

  double **angres=NULL, **avalue=NULL, **rvalue=NULL, **raddis=NULL, **prod=NULL;

  FILE *pgmfile = NULL;

  /* Monte Carlo integration stuff */
  double res=0, err=0;
     
  double xl[DIM_ANG_INT];
  double xu[DIM_ANG_INT];

  /*  size_t calls = 10; */
  size_t calls = 5000000;

#if HAVE_LIBGSL 
  const gsl_rng_type *T;
  const gsl_rng_type *Tiso;
  gsl_rng *r;
  gsl_rng *riso;

  gsl_monte_function G    = { &mc_func, DIM_ANG_INT, &mcinp };
  gsl_monte_function Giso = { &mc_funciso, DIM_ANG_INT, &mcinpiso };
     

  gsl_monte_plain_state *s    = gsl_monte_plain_alloc (DIM_ANG_INT);
  gsl_monte_plain_state *siso = gsl_monte_plain_alloc (DIM_ANG_INT);
#endif
  /* commented out by RB, because was not used but caused compiler warnings */
  /* pmcinp    = &mcinp; */
  /* pmcinpiso = &mcinpiso; */

  /* get command line options  */
  status=get_options (argc, argv, programname, angresfilename, raddisfilename, &pgm,
		      &dev, &ang, &calls, &sza, &phi0, &theta_t, &phi_r, &incdir,
		      &isotropic);

  if (status!=0)
    return status;
  
#if !HAVE_LIBGSL 
  fprintf(stderr,"%s was built without the gsl-library. Thus NO SUPPORT\n", programname);
  fprintf(stderr,"for integration of product of angular response and radiance field.\n");
#endif
  

  /* read angular response file  */
  status = ASCII_file2double (angresfilename,   
			      &arows, &amax_columns, &amin_columns, 
			      &avalue);

  if (status!=0) {
    fprintf (stderr, "error %d reading input file %s\n", status, angresfilename);
    return status;
  }

  if (amin_columns<2) {
    fprintf (stderr, "error, input file must have at least two columns\n");
    return -1;
  }

  if (amin_columns != amax_columns) {
    fprintf (stderr, " WARNING! %s is not a rectangular matrix; only the first\n", angresfilename);
    fprintf (stderr, " %d columns will be used for the analysis\n", amin_columns);
  }
  numu = arows;
  
  /* Read radiance field file  */
  status = ASCII_file2double (raddisfilename,   
			      &rrows, &rmax_columns, &rmin_columns, 
			      &rvalue);

  if (status!=0) {
    fprintf (stderr, "error %d reading input file %s\n", status, raddisfilename);
    return status;
  }

  /* Get umu and phi from radiance field  */
  rphi  = ASCII_row(rvalue, rmax_columns, 1);
  j=0;
  for (i=0;i<rmax_columns;i++) {
    if (rphi[i] >= 0) {
      nphi++;
    }
  }
  phi   = (double *) calloc(nphi, sizeof(double));
  j=0;
  for (i=0;i<nphi;i++) {
    if (rphi[i] != NAN) { /* ??? is that correct? Not even NAN == NAN ??? */
      phi[j++] = rphi[i];
    }
  }
  rumu = ASCII_column(rvalue, rrows, 0);
  for (i=0;i<rrows;i++) { rumu[i] = rumu[i+2]; } /* Shift to skip two first lines */

  /* allocate memory for the radiance field  */
  status = ASCII_calloc_double(&raddis, numu, nphi);
  if (status!=0) { 
    fprintf (stderr, "angres: allocating memory for produce\n");
    return status;
  }
  for (i=0;i<numu;i++) {
    for (j=0;j<nphi;j++) {
      raddis[i][j] = rvalue[i+2][j+2];
    }
  }
  /* commented out by RB, because was not used but caused compiler warnings */
  /* rfldir = rvalue[0][1]; */
  /* rfldn   = rvalue[0][2]; */
  /* flup   = rvalue[0][3]; */
  uavgso = rvalue[0][4];
  /* uavgdn = rvalue[0][5]; */
  /* uavgup = rvalue[0][6]; */

  aumu = ASCII_column(avalue, arows, 0);
  if (ang)
    for (i=0;i<numu;i++)
      aumu[i] = cos(aumu[i]*PI/180);

  aangres = ASCII_column(avalue, arows, 1);

  /* compare umu with similar from angular response file, differences should be small,  */
  /* otherwise interpolate angular response to radiance field.                          */
  for (i=0;i<numu;i++)
    if (fabs(aumu[i]-rumu[i]) > 0.0001)
      fprintf(stderr, "umu %f %f %f\n", aumu[i], rumu[i], aumu[i]-rumu[i]);

  /* Make diffuse radiation field isotropic */
  if (isotropic) {
    /* Integrate (angular response)*radiance over all angles */

    /* Upper and lower integration bounds for each dimension */
    xl[0] = aumu[0];   xu[0] = aumu[numu-1];
    xl[1] = phi[0]*PI/180;    xu[1] = phi[nphi-1]*PI/180;

    mcinpiso.nphi = nphi;
    mcinpiso.numu = numu;
    mcinpiso.phi  = (double *) calloc(nphi, sizeof(double));
    mcinpiso.umu  = (double *) calloc(numu, sizeof(double));
    status = ASCII_calloc_double(&mcinpiso.prod, numu, nphi);
    for (i=0;i<nphi;i++) mcinpiso.phi[i] = phi[i]*PI/180; 
    for (i=0;i<numu;i++) mcinpiso.umu[i] = aumu[i]; 
    for (i=0;i<numu;i++) {
      for (j=0;j<nphi;j++) {
	mcinpiso.prod[i][j] = raddis[i][j];
      }
    }
  
#if HAVE_LIBGSL 
    gsl_rng_env_setup ();
    
    Tiso = gsl_rng_default;
    riso = gsl_rng_alloc (Tiso);
    
    gsl_monte_plain_integrate (&Giso, xl, xu, DIM_ANG_INT, calls, riso, siso,
			       &res, &err);
    gsl_monte_plain_free (siso);
#endif

    dphi    = phi[1]-phi[0];
    dtheta  = fabs(acos(aumu[1])*180/PI-acos(aumu[0])*180/PI);
    tmpprod = dtheta*dphi;

    for (i=0;i<numu;i++) {
      for (j=0;j<nphi;j++) {
	raddis[i][j] = res/tmpprod;
      }
    }
  }



  /* Tilt and rotate angular response */
  angres = tilt_and_rotate(aangres, nphi, phi, numu, aumu, pgm, theta_t, phi_r, dev);
  if (angres==NULL) {
    fprintf (stderr, "Error calling tilt_and_rotate()\n");
    return -1;
  }

  /* Multiply angular response and radiance field */
  /* allocate memory for the product              */

  status = ASCII_calloc_double(&prod, numu, nphi);
  if (status!=0) { 
    fprintf (stderr, "angres: error allocating memory for product\n");
    return status;
  }
  for (i=0;i<numu;i++) {
    for (j=0;j<nphi;j++) {
      prod[i][j] = angres[i][j]*raddis[i][j];
    }
  }

  if (pgm) {

    if ((pgmfile = fopen("prod.pgm","w"))==NULL) {
      fprintf (stderr, "Error opening prod.pgm for writing\n");
      return -1;
    }

    fprintf(pgmfile,"P2\n");
    fprintf(pgmfile,"%d %d\n", nphi, numu);
    fprintf(pgmfile,"256\n");
    max = -9999;
    for (i=0;i<numu;i++)
      for (j=0;j<nphi;j++)
	if (prod[i][j] > max)
	  max = prod[i][j]; 

    for (i=0;i<numu;i++)
      for (j=0;j<nphi;j++) {
	k = (int) (prod[i][j]*256/max);
	fprintf(pgmfile,"%d\n", k);
      } 

    fclose(pgmfile); 
    
    if ((pgmfile = fopen("raddis.pgm","w"))==NULL) {
      fprintf (stderr, "Error opening raddis.pgm for writing\n");
      return -1;
    }

    fprintf(pgmfile,"P2\n");
    fprintf(pgmfile,"%d %d\n", nphi, numu);
    fprintf(pgmfile,"256\n");
    max = -9999;
    for (i=0;i<numu;i++)
      for (j=0;j<nphi;j++)
	if (raddis[i][j] > max) 
	  max = raddis[i][j];

    for (i=0;i<numu;i++)
      for (j=0;j<nphi;j++) {
	k = (int) (raddis[i][j]*256/max);
	fprintf(pgmfile,"%d\n", k);
      } 

    fclose(pgmfile); 
  }

  /* Find angular response for direct beam contribution */
  it = 0;
  ia = 0;

  while (((180-sza)-acos(aumu[ia++])*180/PI) < 0 && ia < numu)
    ;

  it=ia-1; /* This one is fun:-) */

  /* Make sure we get the closest grid point */
  if (fabs(acos(aumu[it])*180/PI-(180-sza)) > fabs(acos(aumu[it-1])*180/PI-(180-sza)))
    it=it-1;

  jt = 0;
  ja = 0;
  delphirad = fabs(phi[1]-phi[0])*PI/180;
  while (fabs(phi0-phi[ja++])*PI/180 > delphirad)
    ;

  jt=ja-1;


  if (incdir) {
    dir = angres[it][jt]*uavgso*4*PI;
  }
  else {
    dir=0;
  }

  ASCII_free_double (avalue, arows);
  ASCII_free_double (rvalue, rrows);
  ASCII_free_double (angres, numu);

  /* Integrate (angular response)*radiance over all angles */

  /* Upper and lower integration bounds for each dimension */
  xl[0] = aumu[0];   xu[0] = aumu[numu-1];
  xl[1] = phi[0]*PI/180;    xu[1] = phi[nphi-1]*PI/180;

  mcinp.nphi = nphi;
  mcinp.numu = numu;
  mcinp.phi  = (double *) calloc(nphi, sizeof(double));
  mcinp.umu  = (double *) calloc(numu, sizeof(double));
  status = ASCII_calloc_double(&mcinp.prod, numu, nphi);
  for (i=0;i<nphi;i++) mcinp.phi[i] = phi[i]*PI/180; 
  for (i=0;i<numu;i++) mcinp.umu[i] = aumu[i]; 
  for (i=0;i<numu;i++) {
    for (j=0;j<nphi;j++) {
      mcinp.prod[i][j] = prod[i][j];
    }
  }
  
#if HAVE_LIBGSL 
  gsl_rng_env_setup ();
  
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  gsl_monte_plain_integrate (&G, xl, xu, DIM_ANG_INT, calls, r, s,
			     &res, &err);
  gsl_monte_plain_free (s);
#endif
  /*  fprintf(stdout,"Result: %f,  error: %f\n", res, err); */
  fprintf(stdout,"%f %f %f\n", res, err, res+dir);

  ASCII_free_double (raddis, numu);

  return 0;   /* if o.k. */
}
