/*--------------------------------------------------------------------
 * $Id: rayleigh.c 2839 2013-01-31 15:13:42Z svn-kylling $
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "mystic.h"
#include "uvspecrandom.h"
#include "f77-uscore.h"
#include "fortran_and_c.h"
#include "miecalc.h"
#include "./ascii.h"

/***********************************************************************************/
/* Function: new_direction_raytracing                                     @62_30i@ */
/* Description:                                                                    */
/*  Calculate new direction at angles (mu, phi) to the original direction.         */
/*  This function calculates the new direction using the geometric optics          */
/*  raytracing approach based on Snells law and Fresnel equations.                 */
/*  The refractive index is retrieved from Warren and Wiscombe's LUT for           */
/*  the respective wavelength with fixed temperature of 273 K.                      */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Linda Forster                                                           */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void F77_FUNC (wrefice, WREFICE) (float *wl, float *temp, float *re, float *im);

double get_crystal_geometry(double **p1, double **p2, double **p3, double **p1_col, double **p2_col, double **p3_col, double **p1_pla, double **p2_pla, double **p3_pla, double fraction_col, int *deg_of_freedom, double fraction_col_orient, double fraction_pla_orient, double dist_col, double dist_pla, double *crystal_axis_tilt, const int um, int *vm);

void raytracing_crystal (double refre, double **p1, double **p2, double **p3, double *dx, double *u, double *v, int deg_of_freedom, double dist, double crystal_axis_tilt);

void random_crystal_rotation(double **p1,double **p2, double **p3, double **p4, double **p5, double **p6, double *sp1, double *sp2, double *sp3, double *xmin, double *xmax, double *ymin, double *ymax, double *zmin, double *zmax, const int um, int *vm, int deg_of_freedom, double dist, double crystal_axis_tilt);

double reflexion_coefficient(double alpha, double n);

double refraction(double *ref1, double *ref2, double *ref3, double *ref6, double theta_i, double *krel, double *k2);

void cross_product(double *a, double *b, double *cross);

void solve_equation_system(double *co11,double *co12, double *co13, double *co21, double *co22, double *co23, double *co31, double *co32, double *co33, double *co14, double *co24, double *co34, double *vf, double *x1, double *x2, double *x3);

void plane_equation(double **p1, double **p2, double **p3, double *nk1, double *nk2, double *nk3, double *fn1, double *fn2, double *fn3, double **kano1, double **kano2, double **kano3, const int um, int *vm);

void PTXA(double *xold, double *yold, double *zold, double *ev1, double *ev2, double *ev3, double *nk1, double *nk2, double *nk3,  double *xs1, double *xs2, double *xs3, double **kano1, double **kano2, double **kano3, double **p1, double **p2, double **p3, double *lq, double *lq_max, int *umtr, int *utr, int *utrold, int *vm, int un, const int um);

void PTXB(double *xold, double *yold, double *zold, double *ev1, double *ev2, double *ev3, double *nk1, double *nk2, double *nk3,  double *xs1, double *xs2, double *xs3, double **kano1, double **kano2, double **kano3, double **p1, double **p2, double **p3, double *lq, double *lq_max, int *umtr, int *utr, int *utrold, int *vm, int un, const int um);


void new_direction_raytracing (atmosphere_struct *atmos, double mu, double phi,
			  direction *dir, double phi_inc, float wvnmlo, float wvnmhi)
{

  int i=0, j=0;
  int l=0;
  double u[3], v[3];
  double u2inv=0;

  float temperature=273.;
  float refre=0;
  float refim=0;
  float wavelength = 0; // in [um]

  static int first_call=1;
  double *first=NULL;
  double *second=NULL;
  double *third=NULL;
  int value=0;
  int um=8;
  int vm[um];
  for (i=0;i<um;i++){
    if (i<2) vm[i]=6;
    else
        vm[i]=4;
  }
  int un = 30;
  int vn = 10;

  double **p1=NULL;
  double **p2=NULL;
  double **p3=NULL;
  p1 = malloc(un*sizeof(double*));
  for(i=0;i<un;i++) p1[i]=malloc(vn*sizeof(double));
  p2 = malloc(un*sizeof(double*));
  for(i=0;i<un;i++) p2[i]=malloc(vn*sizeof(double));
  p3 = malloc(un*sizeof(double*));
  for(i=0;i<un;i++) p3[i]=malloc(vn*sizeof(double));

  static double **p1_col=NULL;
  static double **p2_col=NULL;
  static double **p3_col=NULL;

  static double **p1_pla=NULL;
  static double **p2_pla=NULL;
  static double **p3_pla=NULL;

  double fraction_col = 0.;
  double fraction_col_orient = 0;
  double dist_col = 0.;
//   int dof_col = 0;
//   double fraction_pla = 0;
  double fraction_pla_orient = 0;
  double dist_pla = 0;
//   int dof_pla = 0;
  double dist=0;
  int deg_of_freedom=0;
  double crystal_axis_tilt=0;
  double sum_fraction=0.;
  // start with raytracing code

  wavelength = 1e4*2./(wvnmlo + wvnmhi);
  F77_FUNC (wrefice, WREFICE) (&wavelength, &temperature, &refre, &refim);

  for (i=0; i<atmos->n_raytracing_prop; i++){
      if (strcmp("columns", atmos->raytracing_prop[i].name)==0){
          fraction_col = atmos->raytracing_prop[i].fraction;
          fraction_col_orient = atmos->raytracing_prop[i].oriented_fraction;
          dist_col = atmos->raytracing_prop[i].angdist_width;
//           dof_col = atmos->raytracing_prop[i].orientation_dof;
        }
        if (strcmp(atmos->raytracing_prop[i].name, "plates")==0){
//             fraction_pla = atmos->raytracing_prop[i].fraction;
            fraction_pla_orient = atmos->raytracing_prop[i].oriented_fraction;
            dist_pla = atmos->raytracing_prop[i].angdist_width;
//             dof_pla = atmos->raytracing_prop[i].orientation_dof;
        }
        sum_fraction+=atmos->raytracing_prop[i].fraction;
    }
    if (sum_fraction>1){
        fprintf(stderr, "ERROR crystal fraction sum up to %f, should be 1!\nPlease correct the input parametes in ic_raytracing_file!\n", sum_fraction);
    }

  if (first_call) {
    p1_col = malloc(un*sizeof(double*));
    for(i=0;i<un;i++) p1_col[i]=malloc(vn*sizeof(double));
    p2_col = malloc(un*sizeof(double*));
    for(i=0;i<un;i++) p2_col[i]=malloc(vn*sizeof(double));
    p3_col = malloc(un*sizeof(double*));
    for(i=0;i<un;i++) p3_col[i]=malloc(vn*sizeof(double));

    p1_pla = malloc(un*sizeof(double*));
    for(i=0;i<un;i++) p1_pla[i]=malloc(vn*sizeof(double));
    p2_pla = malloc(un*sizeof(double*));
    for(i=0;i<un;i++) p2_pla[i]=malloc(vn*sizeof(double));
    p3_pla = malloc(un*sizeof(double*));
    for(i=0;i<un;i++) p3_pla[i]=malloc(vn*sizeof(double));

    // read in crystal points
    read_3c_file ("/usr/users/forster/libRadtran/data/ic/crystal_geometry/columns.crystal_points", &first, &second, &third, &value);
    l=0;
    for(i=0;i<um;i++){
        for(j=0;j<vm[i];j++){
            p1_col[i][j]=first[l];
            p2_col[i][j]=second[l];
            p3_col[i][j]=third[l];
            l+=1;
        }
    }
    free(first);
    free(second);
    free(third);
    read_3c_file ("/usr/users/forster/libRadtran/data/ic/crystal_geometry/plates4040.crystal_points", &first, &second, &third, &value);
    l=0;
    for(i=0;i<um;i++){
	for(j=0;j<vm[i];j++){
	    p1_pla[i][j]=first[l];
	    p2_pla[i][j]=second[l];
	    p3_pla[i][j]=third[l];
	    l+=1;
	}
    }
    first_call=0;
    free(first);
    free(second);
    free(third);
  }

  /********************************/
  /* u,v are vectors of length 1; */
  /* u,v,dx form a left-handed    */
  /* orthogonal system            */
  /********************************/

  u[2] = + sqrt( dir->dx[0] * dir->dx[0] + dir->dx[1] * dir->dx[1] );

  if (u[2] == 0.) {
    v[0] = -cosd(phi_inc)*dir->dx[2];
    v[1] = sind(phi_inc)*dir->dx[2];
    v[2] = 0.;
  } else {
    u2inv = 1./u[2];

    v[0] = - dir->dx[1] * u2inv;
    v[1] = + dir->dx[0] * u2inv;
    v[2] =   0.;
  }

  u[0] = - dir->dx[2] * v[1];
  u[1] = + dir->dx[2] * v[0];

  dist = get_crystal_geometry(p1, p2, p3, p1_col, p2_col, p3_col, p1_pla, p2_pla, p3_pla, fraction_col, &deg_of_freedom, fraction_col_orient, fraction_pla_orient, dist_col, dist_pla, &crystal_axis_tilt, um ,vm);

  raytracing_crystal (refre, p1, p2, p3, dir->dx, u, v, deg_of_freedom, dist, crystal_axis_tilt);


  /****************************/
  /* create a new vector from */
  /* u, v, and dx             */
  /****************************/

  /* cosine of solar zenith angle for radiance calculation */
  dir->cotheta = fabs (dir->dx[2]);

  hitflag (dir);
}


double get_crystal_geometry(double **p1, double **p2, double **p3, double **p1_col, double **p2_col, double **p3_col, double **p1_pla, double **p2_pla, double **p3_pla, double fraction_col, int *deg_of_freedom, double fraction_col_orient, double fraction_pla_orient, double dist_col, double dist_pla, double *crystal_axis_tilt, const int um, int *vm){
    int i=0, j=0;
    double rnd=0;
    double dist=0;
    // habit mixture

    rnd = uvspec_random();
    if(rnd<fraction_col){
        // hexagonal columns
        for(i=0;i<um;i++){
            for(j=0;j<vm[i];j++){
            p1[i][j]=p1_col[i][j];  // always start with original coordinates
            p2[i][j]=p2_col[i][j];
            p3[i][j]=p3_col[i][j];
            }
        }

        rnd = uvspec_random();
        if(rnd<fraction_col_orient){
            (*crystal_axis_tilt) = M_PI/2.;
            (*deg_of_freedom) = 2; //oriented columns
            dist = dist_col;
        }
        else{
            (*deg_of_freedom) = 3; // randomly oriented
            dist = 0;
        }
    }
    else{
        // hexagonal plates
        for(i=0;i<um;i++){
            for(j=0;j<vm[i];j++){
            p1[i][j]=p1_pla[i][j];
            p2[i][j]=p2_pla[i][j];
            p3[i][j]=p3_pla[i][j];
            }
        }

        rnd = uvspec_random();
        if(rnd<fraction_pla_orient){
            (*crystal_axis_tilt) = 0;
            (*deg_of_freedom) = 2; //oriented plates
            dist = dist_pla;
        }
        else{
            (*deg_of_freedom) = 3; // randomly oriented
            dist = 0;
        }
    }
    return dist;
}

void raytracing_crystal (double refre, double **p1, double **p2, double **p3, double *dx, double *u, double *v, int deg_of_freedom, double dist, double crystal_axis_tilt){
    int i, j, k;
    int iextref=0, itrans=0, iref=0;

    double ev1[3]={0,0,0};
    double ev2[3]={0,0,0};
    double ev3[3]={0,0,0};
    double vk7[3]={0,0,0};
    double vk9[3]={0,0,0};
    double vk12[3]={0,0,0};
    double vk7_check[3]={0,0,0};
    double cross9[3]={0,0,0};
    double cross12[3]={0,0,0};

    const int un = 30;
    const int vn = 10;
    double sp1, sp2, sp3;
    double **p4=NULL;
    double **p5=NULL;
    double **p6=NULL;
    p4 = malloc(un*sizeof(double*));
    for(i=0;i<un;i++) p4[i]=malloc(vn*sizeof(double));
    p5 = malloc(un*sizeof(double*));
    for(i=0;i<un;i++) p5[i]=malloc(vn*sizeof(double));
    p6 = malloc(un*sizeof(double*));
    for(i=0;i<un;i++) p6[i]=malloc(vn*sizeof(double));

    double *fn1=NULL;
    double *fn2=NULL;
    double *fn3=NULL;
    fn1=malloc(un*sizeof(double));
    fn2=malloc(un*sizeof(double));
    fn3=malloc(un*sizeof(double));
    double *nk1=NULL;
    double *nk2=NULL;
    double *nk3=NULL;
    nk1=malloc(un*sizeof(double));
    nk2=malloc(un*sizeof(double));
    nk3=malloc(un*sizeof(double));

    double **kano1=NULL;
    kano1=malloc(un*sizeof(double*));
    for(i=0;i<un;i++) kano1[i]=malloc(vn*sizeof(double));
    double **kano2=NULL;
    kano2=malloc(un*sizeof(double*));
    for(i=0;i<un;i++) kano2[i]=malloc(vn*sizeof(double));
    double **kano3=NULL;
    kano3=malloc(un*sizeof(double*));
    for(i=0;i<un;i++) kano3[i]=malloc(vn*sizeof(double));

    double *xs1=NULL;
    double *xs2=NULL;
    double *xs3=NULL;
    xs1=malloc(un*sizeof(double));
    xs2=malloc(un*sizeof(double));
    xs3=malloc(un*sizeof(double));

    double *lq=NULL;
    lq=malloc(un*sizeof(double));
    int *umtr=NULL;
    umtr=malloc(un*sizeof(double));

    const int um=8;
    int vm[um];
    for (i=0;i<um;i++){
	if (i<2)
            vm[i]=6;
	else
	    vm[i]=4;
    }
    double xmax=0,xmin=0,ymax=0,ymin=0,zmax=0,zmin=0;
    double lq_max=0;
    double xold=0,yold=0,zold=0;
    double rnd=0;
    double R=0;
    double alpha=0, beta=0, beta_c=0;
    double w=0, nf=0;
    double help1=0, help4=0;
//     double sca_extrefl=0;
    int utr=0, utrold=0;
    double phi=0, mu=0;

    // PEX determine min max
    xmax = p1[0][0];
    ymax = p2[0][0];
    zmax = p3[0][0];
    xmin = p1[0][0];
    ymin = p2[0][0];
    zmin = p3[0][0];
    for (i=0;i<um;i++){
	for(j=0;j<vm[i];j++){
	    if (xmax < p1[i][j]) xmax = p1[i][j];  
	    if (ymax < p2[i][j]) ymax = p2[i][j];
	    if (zmax < p3[i][j]) zmax = p3[i][j];
	    if (xmin > p1[i][j]) xmin = p1[i][j];
	    if (ymin > p2[i][j]) ymin = p2[i][j];
	    if (zmin > p3[i][j]) zmin = p3[i][j];
	}
    }
    lq_max = (xmax - xmin)*(xmax - xmin) + (ymax - ymin)*(ymax - ymin) + (zmax - zmin)*(zmax - zmin);

    random_crystal_rotation(p1, p2, p3, p4, p5, p6, &sp1, &sp2, &sp3, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax, um, vm, deg_of_freedom, dist, crystal_axis_tilt);

    plane_equation(p1,p2,p3,nk1,nk2,nk3,fn1,fn2,fn3,kano1,kano2,kano3, um, vm);

    for(k=0;k<3;k++){
//       ev0[k] = dx[k];
      ev3[k] = dx[k];
      ev1[k] = u[k];
      ev2[k] = v[k];
    }

    utrold = 999;
    utr=-999;
    while(utr==-999){
      phi = uvspec_random()*2.*M_PI;
      mu = sqrt( uvspec_random()); // interval [0,1]

      xold = sp1-dx[0]*sqrt(lq_max)/2.+(cos(phi)*u[0] + sin(phi)*v[0])*sqrt(lq_max)/2.*mu;
      yold = sp2-dx[1]*sqrt(lq_max)/2.+(cos(phi)*u[1] + sin(phi)*v[1])*sqrt(lq_max)/2.*mu;
      zold = sp3-dx[2]*sqrt(lq_max)/2.+(cos(phi)*u[2] + sin(phi)*v[2])*sqrt(lq_max)/2.*mu;

      PTXA(&xold,&yold,&zold,ev1,ev2,ev3,nk1,nk2,nk3,xs1,xs2,xs3,kano1,kano2,kano3,p1,p2,p3,lq,&lq_max,umtr,&utr,&utrold, vm, un, um);
    }

    alpha = acos(-(ev3[0]*fn1[utr]+ev3[1]*fn2[utr]+ev3[2]*fn3[utr]));//*180./pi
    beta = asin(1./(refre) * sin(alpha));
    R = reflexion_coefficient(alpha, refre);

    // PRTA
    w = ev3[0]*fn1[utr] + ev3[1]*fn2[utr] + ev3[2]*fn3[utr];
    if (w < -1)
      w = -1.;
    w = -fabs(w);

    rnd=uvspec_random();
    if(rnd<R){
//       fprintf(stderr, "external reflection\n");
      iextref+=1;
      //********************************************
      // external reflection
      //********************************************
      help1 = -2.*w; // Attention: if w not < 0 in PRTA then @$
      vk9[0] = ev3[0] + help1*fn1[utr]; // direction vector for reflection 
      vk9[1] = ev3[1] + help1*fn2[utr];
      vk9[2] = ev3[2] + help1*fn3[utr];
      nf = sqrt(vk9[0]*vk9[0] + vk9[1]*vk9[1] + vk9[2]*vk9[2]);
      for(k=0;k<3;k++)
	vk9[k] = vk9[k] / nf; 

      if(alpha==0.){
	for(k=0;k<3;k++){
	  vk7[k] = -ev1[k];
	  cross9[k] = -ev2[k];
	}
      }
      else{
	vk7[0] = fn2[utr]*ev3[2] - fn3[utr]*ev3[1]; // vector perp. to plane of incidence
	vk7[1] = fn3[utr]*ev3[0] - fn1[utr]*ev3[2]; // perp = f x prop_in
	vk7[2] = fn1[utr]*ev3[1] - fn2[utr]*ev3[0];
	nf = sqrt(vk7[0]*vk7[0] + vk7[1]*vk7[1] + vk7[2]*vk7[2]);
	help4 = -fabs(w)/w; // sim. f pointing in(out)ward for ex(in)ternal reflection
	for(k=0;k<3;k++)
	  vk7[k] = help4*vk7[k]/nf; 

	cross_product(vk9, vk7, cross9);
      }

      // calculate scattering angle
//       *mu_s = fabs(ev3[0]*vk9[0] +ev3[1]*vk9[1] +ev3[2]*vk9[2]);//*180./M_PI;

      for(k=0;k<3;k++){
	dx[k] = vk9[k];
	u[k] = vk7[k];
	v[k] = cross9[k];	
      }

//       sca_extrefl = acos(*mu_s)*180./M_PI;
//       check = 180.-2.*alpha/M_PI*180.;

    }
    else{
      itrans+=1;
//       fprintf(stderr, "1. transmission\n");
      //********************************************
      // 1. transmission
      //********************************************
      if(alpha==0)
	fprintf(stdout,"Error\n");

      if(alpha!=0){
	vk12[0] = sin(beta)/sin(alpha) * ev3[0] + (sin(beta)/tan(alpha)-cos(beta))* fn1[utr];
	vk12[1] = sin(beta)/sin(alpha) * ev3[1] + (sin(beta)/tan(alpha)-cos(beta))* fn2[utr];
	vk12[2] = sin(beta)/sin(alpha) * ev3[2] + (sin(beta)/tan(alpha)-cos(beta))* fn3[utr];

	nf = sqrt(vk12[0] * vk12[0] + vk12[1] * vk12[1] + vk12[2] * vk12[2]);
	for(k=0;k<3;k++)
	  vk12[k] = vk12[k] / nf; 

	//third solution
	vk7[0] = fn2[utr]*ev3[2] - fn3[utr]*ev3[1]; // vector perp. to plane of incidence
	vk7[1] = fn3[utr]*ev3[0] - fn1[utr]*ev3[2]; // perp = f x prop_in
	vk7[2] = fn1[utr]*ev3[1] - fn2[utr]*ev3[0];
	nf = sqrt(vk7[0]*vk7[0] + vk7[1]*vk7[1] + vk7[2]*vk7[2]);
	help4 = -fabs(w)/w; // sim. f pointing in(out)ward for ex(in)ternal reflection
	for(k=0;k<3;k++)
	  vk7[k] = help4*vk7[k]/nf; 

	// test if every vector is perpendicular
	if((vk7[0]*vk12[0]+vk7[1]*vk12[1]+vk7[2]*vk12[2])>1e-6)
	  fprintf(stdout,"1 Error\n");
	if((vk7[0]*dx[0]+vk7[1]*dx[1]+vk7[2]*dx[2])>1e-6)
	  fprintf(stdout,"2 Error\n");
	if((vk7[0]*fn1[utr]+vk7[1]*fn2[utr]+vk7[2]*fn3[utr])>1e-6)
	  fprintf(stdout,"3 Error\n");
	
	if((vk7[0]*ev1[0]+vk7[1]*ev1[1]+vk7[2]*ev1[2])<0)
	  for(k=0;k<3;k++)
	    vk7[k]=-vk7[k];

	if((vk7[0]*ev1[0]+vk7[1]*ev1[1]+vk7[2]*ev1[2])<0)
	  fprintf(stdout,"error\n");

	cross_product(vk12, vk7, cross12);
	
	if((cross12[0]*ev2[0]+cross12[1]*ev2[1]+cross12[2]*ev2[2])<0)
	  for(k=0;k<3;k++)
	    cross12[k]=-cross12[k];

	if((cross12[0]*vk7[0]+cross12[1]*vk7[1]+cross12[2]*vk7[2])>1e-5)
	  fprintf(stdout,"error\n");
	if((cross12[0]*ev2[0]+cross12[1]*ev2[1]+cross12[2]*ev2[2])<0)
	  fprintf(stdout,"error\n");
      }
//       sca_extrefl = acos(ev3[0]*vk12[0] +ev3[1]*vk12[1] +ev3[2]*vk12[2])*180./M_PI;
//       check = alpha/M_PI*180.-asin(1./(refre)*sin(alpha))/M_PI*180.;
//       check = (alpha-beta)/M_PI*180.;

      for(k=0;k<3;k++){
	ev3[k]=vk12[k];
	ev1[k]=vk7[k];
	ev2[k]=cross12[k];
      }

      rnd=0.;
      while(rnd<R){

        xold = xs1[utr];
        yold = xs2[utr];
        zold = xs3[utr];
	
        utrold = utr;
        // PTXB
        PTXB(&xold,&yold,&zold,ev1,ev2,ev3,nk1,nk2,nk3,xs1,xs2,xs3,kano1,kano2,kano3,p1,p2,p3,lq,&lq_max,umtr,&utr,&utrold, vm, un, um);

	w = ev3[0]*fn1[utr] + ev3[1]*fn2[utr] + ev3[2]*fn3[utr]; 
        if (w>1.)
	  w=1.;
        w = fabs(w);

	beta = acos(ev3[0]*fn1[utr]+ev3[1]*fn2[utr]+ev3[2]*fn3[utr]);//*180./pi

	// Total reflection!!!
	beta_c = asin(1./(refre));
	if ( beta >= beta_c){
	  rnd = 0;
	}
	else{
	  alpha = asin((refre) * sin(beta));
        R = reflexion_coefficient(alpha, refre);
	rnd=uvspec_random();
	}

	if(rnd < R){
//             fprintf(stderr, "internal reflection\n");
	  //********************************************
	  // internal reflection
	  //********************************************
	  iref+=1;
	   //  ---- internal reflection ----  
	  help1 = -2.*w; //! Attention: if w not < 0 in PRTA then @//$
	  vk9[0] = ev3[0] + help1*fn1[utr]; //! direction vector for reflection 
	  vk9[1] = ev3[1] + help1*fn2[utr];
	  vk9[2] = ev3[2] + help1*fn3[utr];   
	  nf = sqrt(vk9[0]*vk9[0] + vk9[1]*vk9[1] + vk9[2]*vk9[2]);
	  vk9[0] = vk9[0] / nf; 
	  vk9[1] = vk9[1] / nf; 
	  vk9[2] = vk9[2] / nf;

	  if(alpha==0.){
	    for(k=0;k<3;k++){
	      vk7[k] = -ev1[k];
	      cross9[k] = -ev2[k];
	    }
	  }
	  else{
	    vk7[0] = fn2[utr]*ev3[2] - fn3[utr]*ev3[1]; // vector perp. to plane of incidence
	    vk7[1] = fn3[utr]*ev3[0] - fn1[utr]*ev3[2]; // perp = f x prop_in
	    vk7[2] = fn1[utr]*ev3[1] - fn2[utr]*ev3[0];
	    nf = sqrt(vk7[0]*vk7[0] + vk7[1]*vk7[1] + vk7[2]*vk7[2]);
	    help4 = -fabs(w)/w; // sim. f pointing in(out)ward for ex(in)ternal reflection
	    for(k=0;k<3;k++)
	      vk7[k] = help4*vk7[k]/nf;
	    cross_product(vk9, vk7, cross9);
	  }
// 	  sca_extrefl = acos(ev3[0]*vk9[0] +ev3[1]*vk9[1] +ev3[2]*vk9[2])*180./M_PI;
// 	  check = 180.-2.*beta/M_PI*180.;
	
	  for(k=0;k<3;k++){
	    ev3[k] = vk9[k];
	    ev1[k] = vk7[k]; // e senkrecht
	    ev2[k] = cross9[k];  // e parallel
	  }
	} // if, internal reflection
      }// while

      //********************************************
      // 2. transmission
      //********************************************
//       fprintf(stderr, "2. transmission\n");
      itrans+=1;

      if(alpha==0){
	for(k=0;k<3;k++)
	  vk12[k]=vk12[k];
      }
      else{
	//  ---- direction of refraction ----  
	vk12[0] = sin(alpha)/sin(beta)*ev3[0] + (cos(alpha)-sin(alpha)/tan(beta))*fn1[utr];
	vk12[1] = sin(alpha)/sin(beta)*ev3[1] + (cos(alpha)-sin(alpha)/tan(beta))*fn2[utr];
	vk12[2] = sin(alpha)/sin(beta)*ev3[2] + (cos(alpha)-sin(alpha)/tan(beta))*fn3[utr];

	nf = sqrt(vk12[0] * vk12[0] + vk12[1] * vk12[1] + vk12[2] * vk12[2]);
	for(k=0;k<3;k++)
	  vk12[k] = vk12[k] / nf; 
      }
      // calculate scattering angle
//       *mu_s = (ev0[0]*vk12[0] +ev0[1]*vk12[1] +ev0[2]*vk12[2]);//*180./M_PI; changed to ev0, since this is the original photon direction
      //third solution
      vk7[0] = fn2[utr]*ev3[2] - fn3[utr]*ev3[1]; // vector perp. to plane of incidence
      vk7[1] = fn3[utr]*ev3[0] - fn1[utr]*ev3[2]; // perp = f x prop_in
      vk7[2] = fn1[utr]*ev3[1] - fn2[utr]*ev3[0];
      nf = sqrt(vk7[0]*vk7[0] + vk7[1]*vk7[1] + vk7[2]*vk7[2]);
      help4 = -fabs(w)/w; // sim. f pointing in(out)ward for ex(in)ternal reflection
      for(k=0;k<3;k++)
	vk7[k] = help4*vk7[k]/nf; 

      if((vk7[0]*vk7_check[0]+vk7[1]*vk7_check[1]+vk7[2]*vk7_check[2])<0)
	for(k=0;k<3;k++)
	  vk7[k]=-vk7[k];

      if((vk7[0]*vk7_check[0]+vk7[1]*vk7_check[1]+vk7[2]*vk7_check[2])<0)
	fprintf(stdout,"error\n");

      cross_product(vk12, vk7, cross12);

      if((cross12[0]*ev2[0]+cross12[1]*ev2[1]+cross12[2]*ev2[2])<0)
	for(k=0;k<3;k++)
	  cross12[k]=-cross12[k];

      if((cross12[0]*vk7[0]+cross12[1]*vk7[1]+cross12[2]*vk7[2])>1e-5)
	fprintf(stdout,"error\n");
      if((cross12[0]*ev2[0]+cross12[1]*ev2[1]+cross12[2]*ev2[2])<0)
	fprintf(stdout,"error\n");

//       sca_extrefl = acos(ev3[0]*vk12[0] +ev3[1]*vk12[1] +ev3[2]*vk12[2])*180./M_PI;
//       check = (alpha-beta)/M_PI*180.;

      for(k=0;k<3;k++){
	dx[k] = vk12[k];
	u[k] = vk7[k];
	v[k] = cross12[k];
      }
    }
    free(fn1);
    free(fn2);
    free(fn3);
    free(nk1);
    free(nk2);
    free(nk3);
    for(i=0;i<un;i++) free(kano1[i]);
    free(kano1);
    for(i=0;i<un;i++) free(kano2[i]);
    free(kano2);
    for(i=0;i<un;i++) free(kano3[i]);
    free(kano3);
    free(xs1);
    free(xs2);
    free(xs3);
    free(lq);
    free(umtr);
    for(i=0;i<un;i++) free(p4[i]);
    free(p4);
    for(i=0;i<un;i++) free(p5[i]);
    free(p5);
    for(i=0;i<un;i++) free(p6[i]);
    free(p6);
    for(i=0;i<un;i++) free(p1[i]);
    free(p1);
    for(i=0;i<un;i++) free(p2[i]);
    free(p2);
    for(i=0;i<un;i++) free(p3[i]);
    free(p3);
}


double reflexion_coefficient(double alpha, double n){
    // Fresnel reflection;
    //alpha = alpha/180.*pi;
//     fprintf(stderr, "reflexion coefficient\n"); 
    double sinalpha=0, sinbeta=0, beta=0;
    double R=0, R1=0, R2=0;
    sinalpha = sin(alpha);
    sinbeta = sinalpha / (n);
    beta = asin(sinbeta);
    if (sinalpha!=0){
      R1 = -sin(alpha-beta)/sin(alpha+beta); // perpendicular
      R2 = tan(alpha-beta)/tan(alpha+beta); // parallel
      R = 0.5*(R1*R1+R2*R2);
    }
    else
	R = (n-1)*(n-1)/(n+1)/(n+1);
    return R;
}

double refraction(double *ref1, double *ref2, double *ref3, double *ref6, double theta_i, double *krel, double *k2){
//     fprintf(stderr, "refraction\n"); 
    double sintiq=0;
    double ref4=0, ref5=0, ref7=0;
    double test1=0;
    //double test2=0;
    double q4=0, q2=0,g=0;
    double rnstar=0, theta_t=0;
    sintiq = sin(theta_i)*sin(theta_i);
    ref4 = 1. - (1. - (*ref2))/ (*ref1) / (*ref3) * sintiq;
    ref5 = 2. * (*krel) / (*ref1) / (*ref3)*sintiq;
    q4 = ref4*ref4 + ref5*ref5;
    q2 = sqrt(q4);
    g = asin(ref5/q2)/2.;
    test1 = acos(ref4/q2)/2.;
    //test2 = atan(ref5/ref4)/2.;
    g = test1;
    ref7 = (cos(g) - (*k2)*sin(g)) * (cos(g) - (*k2)*sin(g));
    rnstar = sqrt(sintiq + (*ref6) * sqrt(ref4*ref4 + ref5*ref5)*ref7);
    theta_t = asin(sin(theta_i)/rnstar);
    return theta_t;
}

void cross_product(double *a, double *b, double *cross){
    double nf=0;
    cross[0] = (a[1])*(b[2]) - (a[2])*(b[1]);
    cross[1] = (a[2])*(b[0]) - (a[0])*(b[2]);
    cross[2] = (a[0])*(b[1]) - (a[1])*(b[0]);
    nf = sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
    cross[0] = cross[0]/nf;
    cross[1] = cross[1]/nf;
    cross[2] = cross[2]/nf;
}
/*
//rotation_matrix(vk71, vk72, vk73, ev11, ev12, e13, ev31, ev32, ev33);
void rotation_matrix(double *en1, double *en2, double *en3, double *eo1, double *eo2, double *eo3, double *ko1, double *ko2, double *ko3){
  myfunc_arg = (*en1)*(*eo1) + (*en2)*(*eo2) + (*en3)*(*eo3); // perp_new * perp_old
  //if (fabs(myfunc_arg) > 1.0) 
  //    myfunc_arg = sign(1.0, myfunc_arg);
  ev21 = ko2*en3 - ko3*en2; 
  ev22 = ko3*en1 - ko1*en3;
  ev23 = ko1*en2 - ko2*en1;
  test1 = myfunc_arg*(eo1*ev21 + eo2*ev22 + eo3*ev23);//  perp_new * par_old 
  help1 = myfunc_arg*myfunc_arg;
  cos_rot = 2.*help1 - 1.;
  help2 = 4.*help1*(1. - help1); // sin^2(2*psi_rot) 
  help3 = sqrt(help2);
  //sin_rot = sign(help3, test1);
  //rot_stokes(2,2) = cos_rot;
  //rot_stokes(2,3) =-sin_rot;
  //rot_stokes(3,2) = sin_rot;
  //rot_stokes(3,3) = cos_rot;
}
*/

void solve_equation_system(double *co11,double *co12, double *co13, double *co21, double *co22, double *co23, double *co31, double *co32, double *co33, double *co14, double *co24, double *co34, double *vf, double *x1, double *x2, double *x3){
  //fprintf(stderr, "solve_equation_system\n"); 
  double d1=0,d2=0,d3=0,d4=0;
  *x1 = 0;
  *x2 = 0;
  *x3 = 0;
  *vf = 99;
  d1 = (*co11) * (*co22) * (*co33)  +  (*co21) * (*co32) * (*co13)  +  (*co31) * (*co12) * (*co23)  -  (*co31) * (*co22) * (*co13)  -  (*co11) * (*co32) * (*co23)  -  (*co21) * (*co12) * (*co33);
  if (d1==0) (*vf) = 2;
  else{
    d2 = (*co14) * (*co22) * (*co33) + (*co24) * (*co32) * (*co13) + (*co34) * (*co12) * (*co23) - (*co34) * (*co22) * (*co13) - (*co14) * (*co32) * (*co23) - (*co24) * (*co12) * (*co33);
    d3 = (*co11) * (*co24) * (*co33) + (*co21) * (*co34) * (*co13) + (*co31) * (*co14) * (*co23) - (*co31) * (*co24) * (*co13) - (*co11) * (*co34) * (*co23) - (*co21) * (*co14) * (*co33);
    d4 = (*co11) * (*co22) * (*co34) + (*co21) * (*co32) * (*co14) + (*co31) * (*co12) * (*co24) - (*co31) * (*co22) * (*co14) - (*co11) * (*co32) * (*co24) - (*co21) * (*co12) * (*co34);
    *x1 = d2 / d1;
    *x2 = d3 / d1;
    *x3 = d4 / d1;
    *vf = 0;
  }
}


void plane_equation(double **p1, double **p2, double **p3, double *nk1, double *nk2, double *nk3, double *fn1, double *fn2, double *fn3, double **kano1, double **kano2, double **kano3, const int um, int *vm){
//     fprintf(stderr, "plane equation\n"); 
// PNK: coefficients of plane equation  nk1(i)*x + nk2(i)*y + nk3(i)*z = 1, i = 1,..,um
    double x1=0,x2=0,x3=0;
    double vf=0;
    double avek1=0,avek2=0,avek3=0,bvek1=0,bvek2=0,bvek3=0,cvek1=0,cvek2=0,cvek3=0;
    double co11=0,co12=0,co13=0,co14=0,co21=0,co22=0,co23=0,co24=0,co31=0,co32=0,co33=0,co34=0;
    double laenge=0;
    double diff_x=0,diff_y=0,diff_z=0;
    int v=0, um1=0, mainu=0, forlim=0,vmu=0, endp=0;

    for (mainu=0; mainu<um; mainu++){
	co11 = p1[mainu][0];
	co12 = p2[mainu][0];
	co13 = p3[mainu][0];
	co14 = 1;
	co21 = p1[mainu][1];
	co22 = p2[mainu][1];
	co23 = p3[mainu][1];
	co24 = 1;
	co31 = p1[mainu][2];
	co32 = p2[mainu][2];
	co33 = p3[mainu][2];
	co34 = 1;

	solve_equation_system( &co11, &co12, &co13, &co21, &co22, &co23, &co31, &co32, &co33, &co14, &co24, &co34, &vf,&x1, &x2, &x3);

	nk1[mainu] = x1;
	nk2[mainu] = x2;
	nk3[mainu] = x3; // Einheitsvektoren, die Ebenengleichung der einzelnen Kristallflaechen aufspannen
	// PFN: outward directed normal vectors of crystal surfaces 
	avek1 = p1[mainu][1] - p1[mainu][0];
	avek2 = p2[mainu][1] - p2[mainu][0];
	avek3 = p3[mainu][1] - p3[mainu][0];
	bvek1 = p1[mainu][2] - p1[mainu][1];
	bvek2 = p2[mainu][2] - p2[mainu][1];
	bvek3 = p3[mainu][2] - p3[mainu][1];
	cvek1 = avek2*bvek3 - avek3*bvek2; //c = a x b 
	cvek2 = avek3*bvek1 - avek1*bvek3; // points inward 
	cvek3 = avek1*bvek2 - avek2*bvek1;
	laenge = sqrt(cvek1*cvek1 + cvek2*cvek2 + cvek3*cvek3);
	fn1[mainu] = -cvek1/laenge; // normalized & points outward 
	fn2[mainu] = -cvek2/laenge;
	fn3[mainu] = -cvek3/laenge;
	// PPD: normal vectors with respect to surface segments 
	um1 = mainu;
	//inside = 'true';
	vmu = vm[um1];
	forlim = vmu;
	for (v=0;v<forlim;v++){
	    endp = v + 1;
	    if (endp == vmu) endp = 0;
	    diff_x = p1[um1][endp] - p1[um1][v]; // vector in the polygon plane 
	    diff_y = p2[um1][endp] - p2[um1][v]; // joining p_n -> p_[n+1] 
	    diff_z = p3[um1][endp] - p3[um1][v];
	    kano1[um1][v] = fn2[um1]*diff_z - fn3[um1]*diff_y; // vector in the plane 
	    kano2[um1][v] = fn3[um1]*diff_x - fn1[um1]*diff_z; // perp. to diff 
	    kano3[um1][v] = fn1[um1]*diff_y - fn2[um1]*diff_x; // pointing outward 
	}
    }
    //return nk1,nk2,nk3,fn1,fn2,fn3,kano1,kano2,kano3,inside;
}

void PTXA(double *xold, double *yold, double *zold, double *ev1, double *ev2, double *ev3, double *nk1, double *nk2, double *nk3,  double *xs1, double *xs2, double *xs3, double **kano1, double **kano2, double **kano3, double **p1, double **p2, double **p3, double *lq, double *lq_max, int *umtr, int *utr, int *utrold, int *vm, int un, const int um){
//     fprintf(stderr, "PTXA\n"); 
    // subroutine PSG: Ray equations 
    // a ray is determined by the intersection of two plane equations, 
    // where plane I is defined by direction of propagation & 1st electric field vector 
    // and plane  II by            "           & 2nd electric field vector
    double co11=0,co12=0,co13=0,co14=0,co21=0,co22=0,co23=0,co24=0,co31=0,co32=0,co33=0,co34=0;	
    double vf=0;
    double snk11=0,snk12=0,snk13=0, snk21=0, snk22=0, snk23=0;
    double lqh=0;
    double direction=0;
    double ac=0,bc=0,cc=0;
    double x1=0, x2=0, x3=0;
    int i=0, j=0,h=0, k=0, ui=0, m=0, vmu=0, v=0;
    double scalpr=0;
    int inside=0;

    co11 = *xold;
    co12 = *yold;
    co13 = *zold;
    co14 = 1.;
    co21 = *xold + ev2[0];
    co22 = *yold + ev2[1];
    co23 = *zold + ev2[2];
    co24 = 1.;
    co31 = *xold + ev3[0];
    co32 = *yold + ev3[1];
    co33 = *zold + ev3[2];
    co34 = 1.;

    solve_equation_system( &co11, &co12, &co13, &co21, &co22, &co23, &co31, &co32, &co33, &co14, &co24, &co34, &vf,&x1, &x2, &x3);
    snk11 = x1;
    snk12 = x2;
    snk13 = x3;

    co11 = *xold + ev1[0];
    co12 = *yold + ev1[1];
    co13 = *zold + ev1[2];
    co14 = 1.;
    co21 = *xold;
    co22 = *yold;
    co23 = *zold;
    co24 = 1.;
    co31 = *xold + ev3[0];
    co32 = *yold + ev3[1];
    co33 = *zold + ev3[2];
    co34 = 1.;

    solve_equation_system( &co11, &co12, &co13, &co21, &co22, &co23, &co31, &co32, &co33, &co14, &co24, &co34, &vf,&x1, &x2, &x3);	
    snk21 = x1;
    snk22 = x2;
    snk23 = x3;

    for (i=0;i<un;i++){
      lq[i]=999;
      umtr[i]=999;
    }
    h=0;
    int u;
    for (i=0;i<um;i++){ // point of intersection with unbounded surface
      u = i;
      if (u!=(*utrold)){
	// points of intersection ray/crystal surface
	co11 = snk11;
	co12 = snk12;
	co13 = snk13;
	co14 = 1.;
	co21 = snk21;
	co22 = snk22;
	co23 = snk23;
	co24 = 1.;
	co31 = nk1[u];
	co32 = nk2[u];
	co33 = nk3[u];
	co34 = 1.;

	solve_equation_system( &co11, &co12, &co13, &co21, &co22, &co23, &co31, &co32, &co33, &co14, &co24, &co34, &vf,&x1, &x2, &x3);
	ac = x1 - (*xold);
	bc = x2 - (*yold);
	cc = x3 - (*zold);
	direction = (ev3[0]*ac + ev3[1]*bc + ev3[2]*cc); 

	if (direction > 0 && vf==0 &&  u!=(*utrold)){ // (*)
	  // !intersection 'in front of' ray & intersection exists & not old plane 
	  xs1[u] = x1; //! sort with respect to distance 
	  xs2[u] = x2;
	  xs3[u] = x3;
	  lqh = ac * ac + bc * bc + cc * cc;
	  lq[h] = lqh;
	  umtr[h] = u;//+1
	  //h = h + 1  // zaehlt h um 1 hoch, wenn bedingung (*) erfuellt
	  // sortiere der groesse von lq nach
	  for (j=0;j<h;j++){
	    if (lq[h] < lq[j]){ // search for shortest distances
	      for (k=h;k>j;k=k-1){
		lq[k] = lq[k-1];
		umtr[k] = umtr[k-1]; // order of hitted surfaces
	      }
	      lq[j] = lqh;
	      umtr[j] = u;
	    }
	  }
	  h = h + 1;  // zaehlt h um 1 hoch, wenn bedingung (*) erfuellt
	}
      }
    }
    m=0;
    *utr=-999;
    // h=6 hier
    while(m<(h-1) && (*utr)!=u && lq[m]<(*lq_max)){ //  determines wether xs is within the boundary of the surface 
      u = umtr[m];
      ui = u;
      inside = 1;
      vmu = vm[ui];
      v = 0;
      while (inside==1 && v!=vmu){
	v=v+1;
	scalpr = kano1[ui][v-1]*(xs1[ui] - p1[ui][v-1]) + kano2[ui][v-1]*(xs2[ui] - p2[ui][v-1]) + kano3[ui][v-1]*(xs3[ui] - p3[ui][v-1]);
	if (scalpr>0.)
	  inside = 0;
      }
      if (inside==1)
	(*utr) = u;
      m=m+1;
    }
}

void PTXB(double *xold, double *yold, double *zold, double *ev1, double *ev2, double *ev3, double *nk1, double *nk2, double *nk3,  double *xs1, double *xs2, double *xs3, double **kano1, double **kano2, double **kano3, double **p1, double **p2, double **p3, double *lq, double *lq_max, int *umtr, int *utr, int *utrold, int *vm, const int un, const int um){
//     fprintf(stderr, "PTXB\n"); 
    // subroutine PSG: Ray equations 
    // a ray is determined by the intersection of two plane equations, 
    // where plane I is defined by direction of propagation & 1st electric field vector 
    // and plane  II by            "           & 2nd electric field vector
    double co11=0,co12=0,co13=0,co14=0,co21=0,co22=0,co23=0,co24=0,co31=0,co32=0,co33=0,co34=0;	
    double vf=0;
    double snk11=0,snk12=0,snk13=0, snk21=0, snk22=0, snk23=0;
    double lqh=0;
    double direction=0;
    double ac=0,bc=0,cc=0;
    double x1=0, x2=0, x3=0;
    int i=0, j=0,h=0, k=0, u=0, ui=0, su=0, m=0, vmu=0, v=0;
    double scalpr=0;
    int inside=0;

    co11 = *xold;
    co12 = *yold;
    co13 = *zold;
    co14 = 1.;
    co21 = *xold + ev2[0];
    co22 = *yold + ev2[1];
    co23 = *zold + ev2[2];
    co24 = 1.;
    co31 = *xold + ev3[0];
    co32 = *yold + ev3[1];
    co33 = *zold + ev3[2];
    co34 = 1.;

    solve_equation_system( &co11, &co12, &co13, &co21, &co22, &co23, &co31, &co32, &co33, &co14, &co24, &co34, &vf,&x1, &x2, &x3);
    snk11 = x1;
    snk12 = x2;
    snk13 = x3;

    co11 = *xold + ev1[0];
    co12 = *yold + ev1[1];
    co13 = *zold + ev1[2];
    co14 = 1.;
    co21 = *xold;
    co22 = *yold;
    co23 = *zold;
    co24 = 1.;
    co31 = *xold + ev3[0];
    co32 = *yold + ev3[1];
    co33 = *zold + ev3[2];
    co34 = 1.;

    solve_equation_system( &co11, &co12, &co13, &co21, &co22, &co23, &co31, &co32, &co33, &co14, &co24, &co34, &vf,&x1, &x2, &x3);	
    snk21 = x1;
    snk22 = x2;
    snk23 = x3;
    //lq = np.zeros(un)
    //umtr = np.zeros(un, dtype=int)
    for (i=0;i<un;i++){
      lq[i]=999;
      umtr[i]=999;
    }
    h=0;
    for (su=0;su<um;su++)
      umtr[su]=0.;

    for (su=0;su<um;su++){ // point of intersection with unbounded surface
      u = su; 
      //if (u!=(*utrold)){
      // points of intersection ray/crystal surface
      co11 = snk11;
      co12 = snk12;
      co13 = snk13;
      co14 = 1.;
      co21 = snk21;
      co22 = snk22;
      co23 = snk23;
      co24 = 1.;
      co31 = nk1[u];
      co32 = nk2[u];
      co33 = nk3[u];
      co34 = 1.;

      solve_equation_system( &co11, &co12, &co13, &co21, &co22, &co23, &co31, &co32, &co33, &co14, &co24, &co34, &vf,&x1, &x2, &x3);

      ac = x1 - (*xold);
      bc = x2 - (*yold);
      cc = x3 - (*zold);
      direction = (ev3[0]*ac + ev3[1]*bc + ev3[2]*cc); 

      if (direction > 0 && vf==0 &&  u!=(*utrold)){ // (*)
	// !intersection 'in front of' ray & intersection exists & not old plane 
	xs1[u] = x1; //! sort with respect to distance 
	xs2[u] = x2;
	xs3[u] = x3;
	lqh = ac * ac + bc * bc + cc * cc;
	lq[h] = lqh;
	umtr[h] = u;//+1
	//h = h + 1  // zaehlt h um 1 hoch, wenn bedingung (*) erfuellt
	// sortiere der groesse von lq nach
	for (j=0;j<h;j++){
	  if (lq[h] < lq[j]){ // search for shortest distances
	    for (k=h;k>j;k=k-1){
	      lq[k] = lq[k-1];
	      umtr[k] = umtr[k-1]; // order of hitted surfaces
	    }
	    lq[j] = lqh;
	    umtr[j] = u;//+1 
	  }
	}
	h = h + 1;  // zaehlt h um 1 hoch, wenn bedingung (*) erfuellt
      }
    }
    m=0;
    *utr=-999;
    // h=6 hier
    while(m<(h-1) && (*utr)!=u){ //  determines wether xs is within the boundary of the surface 
      u = umtr[m];
      ui = u;
      inside = 1;
      vmu = vm[ui];
      v = 0;
      while (inside==1 && v!=vmu){
	v=v+1;
	scalpr = kano1[ui][v-1]*(xs1[ui] - p1[ui][v-1]) + kano2[ui][v-1]*(xs2[ui] - p2[ui][v-1]) + kano3[ui][v-1]*(xs3[ui] - p3[ui][v-1]);
	if (scalpr>0.)
	  inside = 0;
      }
      if (inside==1)
	(*utr) = u;
      m=m+1;
    }
}


void random_crystal_rotation(double **p1,double **p2, double **p3, double **p4, double **p5, double **p6, double *sp1, double *sp2, double *sp3, double *xmin, double *xmax, double *ymin, double *ymax, double *zmin, double *zmax, const int um, int *vm, int deg_of_freedom, double dist, double crystal_axis_tilt){
//   fprintf(stderr, "random crystal rotation\n"); 
  // random rotation of crystal in 
  // Schwerpunkt
  double x_min=0,y_min=0,z_min=0,x_max=0,y_max=0,z_max=0;
  double alpha_euler=0;
  double beta_euler=0;
  double beta=0;
  double gamma_euler=0;
  double c1=0,c2=0,c3=0;
  double r11=0,r12=0,r13=0,r21=0,r22=0,r23=0,r31=0,r32=0,r33=0;
  double sg1=0,sg2=0,sg3=0;
  double s1=0,s2=0,s3=0;
  int n=0;
  int i=0,j=0;
  n=0;
  sg1=0;
  sg2=0;
  sg3=0;
    for (i=0;i<um;i++){
	for (j=0;j<vm[i];j++){
	    n=n+1;
	    sg1= p1[i][j]+sg1;
	    sg2= p2[i][j]+sg2;
	    sg3= p3[i][j]+sg3; 
	}
    }
  sg1=sg1/n;
  sg2=sg2/n;
  sg3=sg3/n;

  // crystal rotation
  if ((deg_of_freedom) == 3){
    //fprintf(stdout,"random\n");
    alpha_euler = 2. * M_PI * uvspec_random();
    beta_euler= acos(1.-2.*sqrt(uvspec_random())); // as in Macke rt.f90
    gamma_euler = 2. * M_PI* uvspec_random();
  }
  //orientation, singly oriented columns
  else if ((deg_of_freedom) == 2){
    //fprintf(stdout,"orient\n");
    alpha_euler = 2. * M_PI * uvspec_random();
    beta = uvspec_random_gauss(dist);
    //     crystal_axis_tilt = M_PI/2.;
    beta_euler = crystal_axis_tilt + beta/180.*M_PI;
    gamma_euler = 2. * M_PI * uvspec_random();
  }
  else
    fprintf(stderr,"Error in deg_of_freedom! %d \n", deg_of_freedom);
  s1=sin(alpha_euler);
  s2=sin(beta_euler);
  s3=sin(gamma_euler);
  c1=cos(alpha_euler);
  c2=cos(beta_euler);
  c3=cos(gamma_euler);

  r11 = -c2*s1*s3 + c1*c3;
  r12 = -c2*s1*c3 - c1*s3;
  r13 =  s2*s1;
  r21 =  c2*c1*s3 + s1*c3;
  r22 =  c2*c1*c3 - s1*s3;
  r23 = -s2*c1;
  r31 =  s2*s3;
  r32 =  s2*c3;
  r33 =  c2;
  for (i=0;i<um;i++){
    for (j=0;j<vm[i];j++){
      // verschiebung in schwerpunkt
      p4[i][j]= p1[i][j]-sg1;
      p5[i][j]= p2[i][j]-sg2;
      p6[i][j]= p3[i][j]-sg3;
    }
  }
  for (i=0;i<um;i++){
    for (j=0;j<vm[i];j++){		
      p1[i][j] =  p4[i][j]*r11 + p5[i][j]*r12 + p6[i][j]*r13;
      p2[i][j] =  p4[i][j]*r21 + p5[i][j]*r22 + p6[i][j]*r23;
      p3[i][j] =  p4[i][j]*r31 + p5[i][j]*r32 + p6[i][j]*r33;
    }
  }
  for (i=0;i<um;i++){
    for (j=0;j<vm[i];j++){
      p4[i][j] = p1[i][j];
      p5[i][j] = p2[i][j];
      p6[i][j] = p3[i][j];
    }
  }
  // min/max
  x_max = p1[0][0];
  y_max = p2[0][0];
  z_max = p3[0][0];
  x_min = p1[0][0];
  y_min = p2[0][0];
  z_min = p3[0][0];
  for (i=0;i<um;i++){
    for (j=0;j<vm[i];j++){
      if (x_max < p1[i][j]) x_max = p1[i][j]; // use maxval, minval 
      if (y_max < p2[i][j]) y_max = p2[i][j];
      if (z_max < p3[i][j]) z_max = p3[i][j];
      if (x_min > p1[i][j]) x_min = p1[i][j];
      if (y_min > p2[i][j]) y_min = p2[i][j];
      if (z_min > p3[i][j]) z_min = p3[i][j];
    }
  }
  for (i=0;i<um;i++){
    for (j=0;j<vm[i];j++){	
      p6[i][j] = p6[i][j] - z_max - 1;
      p3[i][j] = p6[i][j];
    }
  }
  // min/max
  x_max = p1[0][0];
  y_max = p2[0][0];
  z_max = p3[0][0];
  x_min = p1[0][0];
  y_min = p2[0][0];
  z_min = p3[0][0];
  for (i=0;i<um;i++){
    for (j=0;j<vm[i];j++){
      if (x_max < p1[i][j]) x_max = p1[i][j]; // use maxval, minval 
      if (y_max < p2[i][j]) y_max = p2[i][j];
      if (z_max < p3[i][j]) z_max = p3[i][j];
      if (x_min > p1[i][j]) x_min = p1[i][j];
      if (y_min > p2[i][j]) y_min = p2[i][j];
      if (z_min > p3[i][j]) z_min = p3[i][j];
  }
  }
  // Schwerpunkt
  n=0;
  sg1=0;
  sg2=0;
  sg3=0;
  for (i=0;i<um;i++){
    for (j=0;j<vm[i];j++){
      n=n+1;
      sg1= p1[i][j]+sg1;
      sg2= p2[i][j]+sg2;
      sg3= p3[i][j]+sg3;
    }
  }
  (*sp1)=sg1/n;
  (*sp2)=sg2/n;
  (*sp3)=sg3/n;
  (*xmin)=x_min;
  (*ymin)=y_min;
  (*zmin)=z_min;
  (*xmax)=x_max;
  (*ymax)=y_max;
  (*zmax)=z_max;
 }
