/*--------------------------------------------------------------------
 * $Id: ckdfu.c 3307 2017-09-15 14:55:01Z Claudia.Emde $
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

#include "uvspec.h"
#include "fortran_and_c.h"
#include "f77-uscore.h"
#include "ckdfu.h"

/* molecular constants */
#define M_AIR  28.977
#define M_O3   48.000
#define M_H2O  18.015

#define N_MOL_CKDFU 6.0235e23

int ckdfu (float umco2, float umch4, float umn2o, 
	   float umf11, float umf12, float umf22, 
	   float ***temper, float ***press, float *z, 
	   float ****dens, float ****dens_avg, int nlev,
	   int Nx, int Ny,
	   //float *rhoair, float *rhoh2o, float *rhoo3, 
	   //float *pres, float *temp, float *davg, float *z, int nlev,
	   int h2o_cont,
	   ck_profile *ck)
{
  
  #include "fl_radparams.cinc"


  float *tmp_ck   = calloc (nv1x*mbx*mg, sizeof(float));
  float *tmp_delg = calloc (mbx*mg,      sizeof(float));
 
  float ***fck=NULL, **fdelg=NULL;
  int *kg = calloc (mbx, sizeof(int));
  int ig=0, ib=0, lu=0, lc=0, ix=0, iy=0;

  float *rhoair=NULL, *rhoo3=NULL, *rhoh2o=NULL;
  float *pres=NULL, *temp=NULL, *davg=NULL;

  /* allocate memory for density profiles */
  rhoair = calloc (nlev, sizeof(float));
  rhoh2o = calloc (nlev, sizeof(float));
  rhoo3  = calloc (nlev, sizeof(float));
  temp   = calloc (nlev, sizeof(float));
  pres   = calloc (nlev, sizeof(float));
  davg   = calloc (nlev, sizeof(float));

  /* allocate memory for  fck and fdelg */
  ASCII_calloc_float_3D (&fck, mg, nvx, mbx); 
  ASCII_calloc_float (&fdelg, mg, mbx); 
  
  for (ix=0; ix<Nx; ix++){
    for (iy=0; iy<Ny; iy++){
      for (lc=0; lc<nlev; lc++) {
	rhoair[lc] = dens[MOL_AIR][ix][iy][lc] * M_AIR / N_MOL_CKDFU * 1.0e6;
	rhoh2o[lc] = dens[MOL_H2O][ix][iy][lc] * M_H2O / N_MOL_CKDFU * 1.0e6;
	rhoo3 [lc] = dens[MOL_O3][ix][iy][lc] * M_O3  / N_MOL_CKDFU * 1.0e6;
	temp  [lc] = temper[ix][iy][lc];
	pres  [lc] = press[ix][iy][lc];
      }
      
      
      for (lc=0; lc<nlev-1; lc++){
	davg  [lc] = dens_avg[MOL_AIR][ix][iy][lc] * M_AIR  / 
	  N_MOL_CKDFU * 1.0e6 * (z[lc]-z[lc+1]);
      }
      
      
      void F77_FUNC (ckdfuf, CKDFUF) (float *umco2, float *umch4, float *umn2o, 
				      float *umf11, float *umf12, float *umf22, 
				      float *rhoair, float *rhoh2o, float *rhoo3, 
				      float *pres, float *temp, 
				      float *davg, float *z, int *nlev, int *h2o_cont,
				      int *kg, float *tmp_delg, float *tmp_ck);
      
      F77_FUNC  (ckdfuf, CKDFUF) (&umco2, &umch4, &umn2o,
				  &umf11, &umf12, &umf22,
				  rhoair, rhoh2o, rhoo3,
				  pres, temp, davg, z, &nlev, &h2o_cont, kg,
				  tmp_delg, tmp_ck);
  
      fortran2c_3D_float_ary_noalloc (mg, nvx, mbx, tmp_ck, fck);
      fortran2c_2D_float_ary_noalloc (mg, mbx, tmp_delg, fdelg);
    
      /* Loop over bands */
      for (ib=0; ib<mbx; ib++) {
	
	/* copy number of subbands to final destination */
	ck[ib].ngauss = kg[ib];
	
	for (ig=0; ig<kg[ib]; ig++) {
	  
	  /* copy weight of subbands to final destination */
	  ck[ib].weight[ig] = fdelg[ig][ib];
	  
	  /* copy absorption coefficient to final destination */
	  for (lu=0; lu<nlev-1; lu++)
	    ck[ib].crs[ix][iy][lu][ig] = fck[ig][lu][ib];
	}
      }
    }
  }

  free (tmp_ck); free (tmp_delg); 
  free (kg);
  
  free (rhoair);
  free (rhoh2o);
  free (rhoo3);
  free (temp);
  free (pres);
  free (davg);

  ASCII_free_float_3D(fck, mg, nvx);
  ASCII_free_float(fdelg, mg);
  return 0;
}



int ckdfucld (float *z, float *reff, float *lwc, int nlev,
	      float **tau, float **gg, float **ssa)
{
  
  #include "fl_radparams.cinc"

  /* the following statement is meaningless but it prevents compiler warnings */
  /* "warning: variable ‘mg’ set but not used" etc                            */
  int null = nvx + nv1x + mbx + mg - nvx - nv1x - mbx - mg;

  float *tmp_tau   = calloc (nv1x*mbx, sizeof(float));
  float *tmp_ssa   = calloc (nv1x*mbx, sizeof(float));
  float *tmp_ww    = calloc (4*nv1x*mbx, sizeof(float));

  float **tmpc_tau=NULL, **tmpc_ssa=NULL, ***tmpc_ww=NULL;

  int ib=null, lu=0;

  void F77_FUNC (ckdcldf, CKDCLDF) (float *z, float *reff, float *lwc, int *nlev,
			  float *tmp_tau, float *tmp_ww, float *tmp_ssa);
  
  F77_FUNC  (ckdcldf, CKDCLDF) (z, reff, lwc, &nlev, tmp_tau, tmp_ww, tmp_ssa); /* in libsrc_c/ckdcldf.f */

  tmpc_tau  = fortran2c_2D_float_ary (nvx, mbx, tmp_tau);
  tmpc_ssa  = fortran2c_2D_float_ary (nvx, mbx, tmp_ssa);
  tmpc_ww   = fortran2c_3D_float_ary (nvx, mbx, 4, tmp_ww);
  free (tmp_tau); free (tmp_ssa); free (tmp_ww); 

    
  /* Loop over bands */
  for (ib=0; ib<mbx; ib++) {

    /* copy data to final destination */
    for (lu=0; lu<nlev-1; lu++) { 
      tau[ib][lu] = tmpc_tau[lu][ib];
      ssa[ib][lu] = tmpc_ssa[lu][ib];
      gg [ib][lu] = tmpc_ww [lu][ib][0] / 3.0;  /* first moment of the phase function */
    }
  }
  
  /* Free memory - very important! This was a major memory leak with MYSTIC */
  /* where ckdfucld() is called for every grid box                          */      
  ASCII_free_float (tmpc_tau, nvx);
  ASCII_free_float (tmpc_ssa, nvx);
  ASCII_free_float_3D (tmpc_ww, nvx, mbx);

  mg=0;  /* useless statement, but otherwise we get a warning that mg is defined but not used */

  return 0;
}



int ckdfuice (float *z, float *reff, float *lwc, int nlev,
	      float **tau, float **gg, float **ssa, 
	      int unscaled)
{
  
  #include "fl_radparams.cinc"

  /* the following statement is meaningless but it prevents compiler warnings */
  /* "warning: variable ‘mg’ set but not used" etc                            */
  int null = nvx + nv1x + mbx + mg - nvx - nv1x - mbx - mg;

  float *tmp_tau   = calloc (nv1x*mbx, sizeof(float));
  float *tmp_ssa   = calloc (nv1x*mbx, sizeof(float));
  float *tmp_ww    = calloc (4*nv1x*mbx, sizeof(float));

  float **tmpc_tau=NULL, **tmpc_ssa=NULL, ***tmpc_ww=NULL;

  float *deff = calloc (nlev, sizeof(float));

  int ib=null, lu=0;

  void F77_FUNC (ckdicef, CKDICEF) (float *z, float *reff, float *lwc, int *nlev,
			  float *tmp_tau, float *tmp_ww, float *tmp_ssa,
			  int *unscaled);                                 /* in src/ckdicef.f */
  
  /* Fu et al. requires effective diameter */
  for (lu=0; lu<nlev; lu++) {
    deff[lu] = 2.0 * reff[lu];

    /* test if effective diameter within allowed range; unfortunately we need */
    /* to exclude more than necessary: shortwave and longwave properties are  */
    /* valid for slightly different size ranges; however, both are calculated */
    /* at the same time, hence we can only use the region which is available  */
    /* in both parameterizations.                                             */

    if (lwc[lu]>0) {  /* need only test if lwc>0; otherwise optical properties are not required */
      if (deff[lu] < MIN_DEFF_FU96 || deff[lu] > MAX_DEFF_FU96) {
	fprintf (stderr, "Error, effective radius %f um not covered by Fu [1996]\n", deff[lu]/2.0);
	fprintf (stderr, "Allowed range is %7.3f - %7.3f um\n",
		 MIN_DEFF_FU96/2.0, MAX_DEFF_FU96/2.0);
	return -1;
      }
      
      if (deff[lu] < MIN_DEFF_FU98 || deff[lu] > MAX_DEFF_FU98) {
	fprintf (stderr, "Error, effective radius %f um not covered by Fu et al. [1998]\n", deff[lu]/2.0);
	fprintf (stderr, "Allowed range is %7.3f - %7.3f um\n",
		 MIN_DEFF_FU98/2.0, MAX_DEFF_FU98/2.0);
	return -1;
      }
    }
  }


  F77_FUNC  (ckdicef, CKDICEF) (z, deff, lwc, &nlev, tmp_tau, tmp_ww, tmp_ssa, &unscaled);

  tmpc_tau  = fortran2c_2D_float_ary (nvx, mbx, tmp_tau);
  tmpc_ssa  = fortran2c_2D_float_ary (nvx, mbx, tmp_ssa);
  tmpc_ww   = fortran2c_3D_float_ary (nvx, mbx, 4, tmp_ww);
  free (tmp_tau); free (tmp_ssa); free (tmp_ww); 

    
  /* Loop over bands */
  for (ib=0; ib<mbx; ib++) {

    /* copy data to final destination */
    for (lu=0; lu<nlev-1; lu++) { 
      tau[ib][lu] = tmpc_tau[lu][ib];
      ssa[ib][lu] = tmpc_ssa[lu][ib];
      gg [ib][lu] = tmpc_ww [lu][ib][0] / 3.0;  /* first moment of the phase function */
    }
  }

  ASCII_free_float(tmpc_tau, nvx);
  ASCII_free_float(tmpc_ssa, nvx);
  ASCII_free_float_3D(tmpc_ww, nvx, mbx);
  free(deff);

  mg=0;  /* useless statement, but otherwise we get a warning that mg is defined but not used */

  return 0;
}




int ckdfuray (int ib, int ig, float u0, int nlev, float *tau)
{
  
  #include "fl_radparams.cinc"

  /* the following statement is meaningless but it prevents compiler warnings */
  /* "warning: variable ‘mg’ set but not used" etc                            */
  int null = nvx + nv1x + mbx + mg - nvx - nv1x - mbx - mg;

  void F77_FUNC (ckdrayf, CKDRAYF) (int *ib, int *ig, float *u0, int *nlev, float *tau);
  
  F77_FUNC (ckdrayf, CKDRAYF) (&ib, &ig, &u0, &nlev, tau);

  return null;
}




