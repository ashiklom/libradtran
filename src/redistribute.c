/*--------------------------------------------------------------------
 * $Id: redistribute.c 3280 2017-07-11 12:21:41Z Claudia.Emde $
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

#include "uvspec.h"
#include "ascii.h"
#include "errors.h"

#define QSORT_CAST (int (*)(const void *, const void *))
#define EPSILON   1.0e-6


/************************************/
/* prototypes of internal functions */
/************************************/

static int redistribute_moment (float *****moment, int ***nmom, int nlambda, 
                                int nlambda_lower, int nlambda_upper,
				int nlyr_alloc,
                                float *z_old, int nlyr_old, int nphamat,
                                float *z_new, int nlyr_new, int alloc);

int redistribute_caoth ( wl_out_struct     wl,
			 atm_out_struct    atm,
			 int               nphamat,
			 int               alloc_microphys,
			 int               copy_ssprop,
			 int               quiet,
			 /* input / output */
			 caoth_out_struct *caoth );

int redistribute_caoth_microphys ( caoth_microphys_struct *microphys,
				   float                  *zd,
				   int                     nlyr, 
				   float                  *zd_common,
				   int                     nlyr_common );

/* 3DABS moved to redistribute.h */
/* int redistribute_molecular (atm_optprop_struct *optprop, int nlambda,  */
/* 			    int nlambda_lower, int nlambda_upper, */
/* 			    int *nq,  */
/*                             int nipa, */
/* 			    float *zd, int nlyr, int Nx, int Ny,  */
/* 			    float *zd_common, int nlyr_common); */

int redistribute_molecular_crs (crs_out_struct *crs, 
				int nlambda, int nlambda_lower, int nlambda_upper, 
				int nmol, int number_of_ramanshifts,
				float *zd, int nlyr, 
				float *zd_common, int nlyr_common);

int redistribute_amf (atm_microphys_struct *microphys, 
		      crs_out_struct *crs,
		      int nlambda, 
		      int nlambda_lower, int nlambda_upper,
		      float *zd, int nlyr, 
		      float *zd_common, int nlyr_common);

int redistribute_caoth3d ( caoth3d_out_struct *caoth3d, 
			   float              *zd,
			   int                 nlyr,
			   float              *zd_common,
			   int                 nlyr_common,
			   int                 alloc,
			   int                 quiet );

int redistribute_3D (float ****data, int nx, int ny,
		     float *z_old, int nlyr_old, 
		     float *z_new, int nlyr_new,
		     int *threedold, int *threednew,
		     int scale, int alloc);

static int calloc_float_3D_sparse  (float ****value, int nz, int nx, int ny, int *allocz);
static int free_float_3D_sparse (float ***value, int nz, int nx, int*allocz);


/**************************************************************/
/* Setup redistribute irradiance.                         */
/**************************************************************/

int setup_redistribute (input_struct input, output_struct *output)
{
  int    nuser=0, i=0, isp=0;
  float *zuser=NULL;

  int    lc=0, ialt=0, status=0, old_nlev=0, ipa=0, 
    lu=0, li=0, lk=0, lv=0, ntmp_zd=0;
  float *tmp_zd=NULL, *tmp_zd_rgrid=NULL, *tmp_zd_user=NULL;

  int    iv=0;
  double babso_o3=0, babso_co2=0, babso_o2=0, babso_no2=0;
  double babso_bro=0, babso_oclo=0, babso_hcho=0, babso_md=0;
  double babso_amf=0, babso_amftot=0, vertical_column=0, crs_amf=0;
  double babsotot_o3=0, babsotot_co2=0, babsotot_o2=0, babsotot_no2=0;
  double babsotot_bro=0, babsotot_oclo=0, babsotot_hcho=0;
  double vcol_o3=0, vcol_co2=0, vcol_o2=0, vcol_no2=0, vcol_bro=0;
  double vcol_oclo=0, vcol_hcho=0;
 
  double deltaz=0;

  int redist=0;

  if (input.verbose)
    fprintf (stderr, " *** setup_redistribute()\n");

  aer_out_struct *redistributed_aer= calloc(1, sizeof(aer_out_struct));

  
  /********************** add some more vertical levels to common grid ******************/


  /* this cant be done before verbose output as data for equidistant grid would be missing */
  if (input.atm.zout_interpolate == NO_ZOUT_INTERPOLATE) {

    /*  user-defined equidistant altitude grid  */
    if (input.alt.altitude_dz > 0) {
      nuser = (int) ((output->atm.zd[0] - output->atm.zd[output->atm.nlev-1]) / 
                      input.alt.altitude_dz + 0.5) + 1;
      
      zuser = (float *) calloc (nuser, sizeof(float));

      for (i=0; i<nuser; i++)
        zuser[i] = output->atm.zd[output->atm.nlev-1] + input.alt.altitude_dz * (float) i;
  
      while (zuser[nuser-1] > output->atm.zd[0])
        nuser--;

      /* add user-defined grid to atmospheric grid */
      set_common_z_grid (zuser, nuser, &output->atm.zd_common, &output->atm.nlev_common);

      free (zuser);
    }
  }

  
  /************************ Redistribute molecular optical properties *******************/

  if (input.verbose) {
    fprintf (stderr, " ...  redistributing molecular optical properties\n");
    fflush (stderr);
  }
  
  status = redistribute_molecular (&(output->atm.optprop), output->wl.nlambda_r, 
				   output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper,
				   output->atm.nq_r,
				   output->nipa,
				   output->atm.zd, output->atm.nlev-1,
				   output->atm.Nxatm, output->atm.Nyatm, 
				   output->atm.zd_common, output->atm.nlev_common-1);

  if (status!=0) {
    fprintf (stderr, "Error regridding molecular scattering and absorption\n");
    return status;
  }

  /************************ Redistribute molecular cross sections *******************/

  if (input.verbose) {
    fprintf (stderr, " ...  redistributing molecular cross sections\n");
    fflush (stderr);
  }

  status = redistribute_molecular_crs (&(output->crs), output->wl.nlambda_r, 
				       output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper,
				       MOL_NN, output->crs.number_of_ramanshifts,
				       output->atm.zd, output->atm.nlev-1, 
				       output->atm.zd_common, output->atm.nlev_common-1);
  
  if (status!=0) {
    fprintf (stderr, "Error regridding molecular cross sections\n");
    return status;
  }

  /******************** Redistribute aerosol stuff ****************************/

  if (input.verbose) {
    fprintf (stderr, " ...  redistributing aerosol optical properties\n");
    fflush (stderr);
  }
  
  calloc_aer_out (redistributed_aer, output->wl.nlambda_r, output->atm.nlev_common, output->aer.optprop.nphamat);
 
  status = cp_aer_out (&(output->aer), redistributed_aer, output->wl.nlambda_r, 1, input.quiet);
  if (status!=0) {
    fprintf (stderr, "Error %d copying output->aer to redistributed_aer\n", status);
    return status;
  }

  old_nlev = output->aer.nlev;
  tmp_zd = calloc (old_nlev+1, sizeof (float));

  for (lc=0;lc<old_nlev;lc++){
    tmp_zd[lc] = output->aer.zd[lc]; 
  }
  
  free_aer_out(&(output->aer));
  
  status = redistribute_optprop(&(redistributed_aer->optprop),
                                output->wl.nlambda_r,
                                output->wl.nlambda_rte_lower,
                                output->wl.nlambda_rte_upper,
                                tmp_zd, old_nlev-1, redistributed_aer->optprop.nphamat,
                                output->atm.zd_common,
                                output->atm.nlev_common-1, 0);
  if (status!=0) {
    fprintf (stderr, "Error %d during redistribute_optprop (aero) (line %d, function %s in %s)\n", 
	     status, __LINE__, __func__, __FILE__);
    return status;
  }
  
  calloc_aer_out(&(output->aer), output->wl.nlambda_r, output->atm.nlev_common, redistributed_aer->optprop.nphamat);

  status = cp_aer_out(redistributed_aer, &(output->aer), output->wl.nlambda_r, 1, input.quiet);
  if (status!=0) {
    fprintf (stderr, "Error %d copying redistributed_aer to output->aer (line %d, function %s in %s)\n", 
                      status, __LINE__, __func__, __FILE__);
    return status;
  }
  
  if (status!=0) {
    fprintf (stderr, "Error regridding aerosol parameters\n");
    return status;
  }

  /* copy number of levels and re-write altitude grid */
  output->aer.nlev = output->atm.nlev_common;
  free (output->aer.zd);
  output->aer.zd = calloc (output->aer.nlev, sizeof(float));
  for (lc=0; lc<output->aer.nlev; lc++)
    output->aer.zd[lc] = output->atm.zd_common[lc];

  /* free_aer_out(redistributed_aer, output->wl.nlambda_r); */

  free(tmp_zd);
  free_aer_out(redistributed_aer);
  free(redistributed_aer);

  /**************************** Redistribute caoth stuff *************************/

  for (isp=0; isp<input.n_caoth; isp++) {

    if (input.verbose) {
      fprintf (stderr, " ...  redistributing optical properties for %s\n",
	       input.caoth[isp].fullname);
      fflush (stderr);
    }

    /********************************/
    /* loop over independent pixels */
    /********************************/

    /* better would probably be a source switch here SBCA */
    if (input.caoth[isp].ipa)
      for (ipa=0; ipa<output->nipa; ipa++) {
	status = redistribute_caoth ( output->wl,
				      output->atm,
				      output->caoth[isp].optprop.nphamat,
				      output->caoth[isp].microphys.alloc,
				      0,
				      input.quiet,
				      &(output->caoth_ipa[isp][ipa]) );
	if (status)
	  return fct_err_out ( status, "redistribute_caoth", ERROR_POSITION );
      }

    /* also needed in case of 3d files. lwc_layer needs to be
       redistributed. This is kinda overkill, and could be cleaned up
       some day. */
    status = redistribute_caoth ( output->wl,
				  output->atm,
				  output->caoth[isp].optprop.nphamat,
				  output->caoth[isp].microphys.alloc,
				  1,
				  input.quiet,
				  &(output->caoth[isp]) );
    if (status)
      return fct_err_out ( status, "redistribute_caoth", ERROR_POSITION );

    /* redistribute 3D caoths */
    /* XXXXX check if really need to redistribute caoth */
    if (output->caoth3d[isp].nthreed>0) {

      redist=1;  /* only redistribute if layers differ */
      if (output->caoth3d[isp].nlyr==output->atm.nlev_common-1) {
	redist=0;
	for (lc=0;lc<=output->caoth3d[isp].nlyr;lc++) {
	  if ( output->caoth3d[isp].zd[lc] !=
	       output->atm.zd_common[output->caoth3d[isp].nlyr-lc] ) {
	    fprintf (stderr, "LEVELS %f AND %f DIFFER!\n", 
		     output->caoth3d[isp].zd[lc],
		     output->atm.zd_common[output->caoth3d[isp].nlyr-lc]);
	    redist=1;
	    break;
	  }
	}
      }

      if (redist) {
	status = redistribute_caoth3d ( &(output->caoth3d[isp]),
					output->caoth3d[isp].zd,
					output->caoth3d[isp].nlyr, 
					output->atm.zd_common,
					output->atm.nlev_common-1,
					1,
					input.quiet );
	if (status)
	  return fct_err_out ( status, "redistribute_caoth3d", ERROR_POSITION );
      }
    }

  } /* end loop isp */

  /*********************** Interpolate cloud fraction **********************/

  if (output->cf.nlev > 0) {
    /* WHY cf.nlev-1 AND nlev_common-1 ??? UH. 2009-04 */
    status += redistribute_1D ((void *) (&(output->cf.cf)),
			       output->cf.zd, output->cf.nlev-1, 
			       output->atm.zd_common, output->atm.nlev_common-1,
			       0, REDISTRIBUTE_FLOAT, 1);
    if (status != 0) {
      fprintf (stderr, "Error %d during redistribute cloudfraction (line %d, function %s in %s)  \n", 
                       status, __LINE__, __func__, __FILE__);
      return status;
    }
    /* realloc output->cf.zd */
    output->cf.nlev = output->atm.nlev_common-1;
    free(output->cf.zd);
    output->cf.zd = calloc (output->cf.nlev, sizeof (float));
    for (lc=0; lc<output->cf.nlev; lc++)
      output->cf.zd[lc] = output->atm.zd_common[lc];

    if (output->nipa>1) {
      for (ipa=0; ipa<output->nipa; ipa++) {
        status += redistribute_1D ((void *) (&(output->cfipa[ipa].cf)),
                                   output->cfipa[ipa].zd, output->cfipa[ipa].nlev-1, 
                                   output->atm.zd_common, output->atm.nlev_common-1,
                                   0, REDISTRIBUTE_FLOAT, 1);
        if (status != 0) {
          fprintf (stderr, "Error %d during redistribute ipa cloudfraction (line %d, function %s in %s)  \n", 
                            status, __LINE__, __func__, __FILE__);
          return status;
        }
        /* realloc output->cfipa[ipa].zd */
        output->cfipa[ipa].nlev = output->atm.nlev_common-1;
        free(output->cfipa[ipa].zd);
        output->cfipa[ipa].zd = calloc (output->cfipa[ipa].nlev, sizeof (float));
        for (lc=0; lc<output->cfipa[ipa].nlev; lc++)
          output->cfipa[ipa].zd[lc] = output->atm.zd_common[lc];
      }
    }
  }

  /************************ Interpolate atmosphere *************************/

  if (input.verbose) {
    fprintf (stderr, " ...  interpolating atmospheric profiles: p, T, densities \n");
    fflush (stderr);
  }

  status += interpolate_atmosphere (output->atm.zd,
                                    &(output->atm.microphys.press),
                                    &(output->atm.microphys.temper),
                                    &(output->atm.microphys.dens),
                                    &(output->atm.microphys.temper_avg),
                                    &(output->atm.microphys.dens_avg),
                                    output->atm.nlev,
                                    output->atm.zd_common, output->atm.nlev_common,
                                    input.atm.interpol_method_press,
                                    input.atm.interpol_method_temper,
                                    input.atm.interpol_method_gas,
                                    input.quiet);
  
  if (status != 0) {
    fprintf (stderr, "Error %d interpolating atmosphere (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
    return status;
  }
  
  /* additional verbose output */
/*   if (input.verbose) { */
 
/*     fprintf (stderr, "#-------------------------------------------------------------------------------------------------------------------------"); */
/*     if (output->column[MOL_BRO]  !=0) */
/*       fprintf (stderr, "--------------"); */
/*     if (output->column[MOL_OCLO] !=0) */
/*       fprintf (stderr, "--------------"); */
/*     if (output->column[MOL_HCHO] !=0) */
/*       fprintf (stderr, "--------------"); */
/*     fprintf (stderr, "\n"); */
/*     fprintf (stderr, "# lc |  z[km]  |  Pressure  | Temp.  |    Air      |   Ozone     |     O2      | Water vap.  |    CO2      |    NO2      | "); */
/*     if (output->column[MOL_BRO]  !=0) */
/*       fprintf (stderr, "    BRO      |"); */
/*     if (output->column[MOL_OCLO] !=0) */
/*       fprintf (stderr, "    OCLO     |"); */
/*     if (output->column[MOL_HCHO] !=0) */
/*       fprintf (stderr, "    HCHO     |"); */
/*     fprintf (stderr, "\n"); */
/*     fprintf (stderr, "#    |         |   [hPa]    |  [K]   |   [cm-3]    |   [cm-3]    |   [cm-3]    |   [cm-3]    |   [cm-3]    |   [cm-3]    | "); */
/*     if (output->column[MOL_BRO]  !=0) */
/*       fprintf (stderr, "   [cm-3]    |"); */
/*     if (output->column[MOL_OCLO] !=0) */
/*       fprintf (stderr, "   [cm-3]    |"); */
/*     if (output->column[MOL_HCHO] !=0) */
/*       fprintf (stderr, "   [cm-3]    |"); */
/*     fprintf (stderr, "\n"); */
/*     fprintf (stderr, "#-------------------------------------------------------------------------------------------------------------------------"); */
/*     if (output->column[MOL_BRO]  !=0) */
/*       fprintf (stderr, "--------------"); */
/*     if (output->column[MOL_OCLO] !=0) */
/*       fprintf (stderr, "--------------"); */
/*     if (output->column[MOL_HCHO] !=0) */
/*       fprintf (stderr, "--------------"); */
/*     fprintf (stderr, "\n"); */


/*     for (lc=0; lc<output->atm.nlev_common; lc++) { */
/*       if (output->atm.zd[lc] >= output->alt.altitude) { */
/*         fprintf (stderr, "%5d  %7.3f   %10.5f   %6.2f   %.5e   %.5e   %.5e   %.5e   %.5e   %.5e", */
/* 	         lc, output->atm.zd_common[lc], */
/* 	         output->atm.microphys.press[lc], output->atm.microphys.temper[lc],output->atm.microphys.dens[MOL_AIR][lc], */
/* 	         output->atm.microphys.dens[MOL_O3][lc],output->atm.microphys.dens[MOL_O2][lc], output->atm.microphys.dens[MOL_H2O][lc], */
/*                  output->atm.microphys.dens[MOL_CO2][lc],output->atm.microphys.dens[MOL_NO2][lc]); */
/*         if (output->column[MOL_BRO]  !=0) */
/*           fprintf (stderr, "   %.5e",output->atm.microphys.dens[MOL_BRO][lc]); */
/*         if (output->column[MOL_OCLO] !=0) */
/*           fprintf (stderr, "   %.5e",output->atm.microphys.dens[MOL_OCLO][lc]); */
/*         if (output->column[MOL_HCHO] !=0) */
/*           fprintf (stderr, "   %.5e",output->atm.microphys.dens[MOL_HCHO][lc]); */
/*         fprintf (stderr, "\n"); */
/*       } */
/*     } */

/*     fprintf (stderr, "#-------------------------------------------------------------------------------------------------------------------------------------"); */
/*     if (output->column[MOL_BRO]  !=0) */
/*       fprintf (stderr, "--------------"); */
/*     if (output->column[MOL_OCLO] !=0) */
/*       fprintf (stderr, "--------------"); */
/*     if (output->column[MOL_HCHO] !=0) */
/*       fprintf (stderr, "--------------"); */
/*     fprintf (stderr, "\n"); */

/*   } */
 
  /************************ Interpolate index of refraction ************************/

  if (input.verbose) {
    fprintf (stderr, " ...  interpolating index of refraction\n");
    fflush (stderr);
  }

  for (iv=0; iv<output->wl.nlambda_r; iv++)
    status = interpolate_profile (output->atm.zd, &(output->atm.microphys.refind)[iv], output->atm.nlev,
                                  output->atm.zd_common, output->atm.nlev_common, input.atm.interpol_method_refind,
                                  input.quiet);

  if (status!=0) {
    fprintf (stderr, "Error %d regridding refractive index (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
    return status;
  }

  /************************ Redistribute possible amf stuff ************************/

  if (input.verbose) {
    fprintf (stderr, " ...  redistributing airmass stuff\n");
    fflush (stderr);
  }

  status = redistribute_amf (&(output->atm.microphys), 
			     &(output->crs),
			     output->wl.nlambda_r, 
			     output->wl.nlambda_rte_lower, output->wl.nlambda_rte_upper,
			     output->atm.zd, output->atm.nlev-1, 
			     output->atm.zd_common, output->atm.nlev_common-1);

  if (status!=0) {
    fprintf (stderr, "Error %d regridding airmass stuff (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
    return status;
  }

  /* Finally set new number of output level and layers, and set new zd */
  output->atm.nlev = output->atm.nlev_common;
  output->atm.nlyr = output->atm.nlev_common-1;
  free(output->atm.zd);
  output->atm.zd = calloc (output->atm.nlev+1, sizeof (float));
  for (lc=0; lc<output->atm.nlev_common; lc++) 
    output->atm.zd[lc] = output->atm.zd_common[lc];

  /* For Raman scattering we need output at all atmospheric levels */
  if (input.raman) {

    /* Find all unique altitudes */
    ntmp_zd = output->atm.nlev;
    if ((tmp_zd = calloc (ntmp_zd, sizeof (float))) == NULL) {
      fprintf (stderr,"Error, allocating memory for tmp_zd\n");
      fprintf (stderr, "      (line %d, function %s in %s)  \n", __LINE__, __func__, __FILE__);
      return -1;
    }
    if ((tmp_zd_rgrid = calloc (ntmp_zd, sizeof (float))) == NULL) {
      fprintf (stderr,"Error, allocating memory for tmp_zd_rgrid\n");
      fprintf (stderr, "      (line %d, function %s in %s)  \n", __LINE__, __func__, __FILE__);
      return -1;
    }
    lv=0;
    for (lu=0; lu<output->atm.nlev-1; lu++) {
      if ( output->atm.zd[lu] - output->alt.altitude >= 0.0 ) {
	tmp_zd[lv] = output->atm.zd[lu]- output->alt.altitude;
	lv++;
      }
    }
    tmp_zd[lv] = output->atm.zd[output->atm.nlev]- output->alt.altitude;
    if ( output->alt.altitude != 0) {
      ntmp_zd = lv;  /* Number of output levels may be reduced due to altitude option */
    }
    for (lu=0; lu<ntmp_zd; lu++) {
      tmp_zd_rgrid[lu] = tmp_zd[lu];
    }


    if ((tmp_zd_user = calloc (output->atm.nzout, sizeof (float))) == NULL) {
      fprintf (stderr,"Error, allocating memory for tmp_zd_user \n");
      fprintf (stderr, "      (line %d, function %s in %s)  \n", __LINE__, __func__, __FILE__);
      return -1;
    }
    for (lu=0; lu<output->atm.nzout; lu++) 
      tmp_zd_user[lu] = output->atm.zout_sur[lu];

    
    set_common_z_grid (output->atm.zout_sur, output->atm.nzout, &tmp_zd, &ntmp_zd);
    free(output->atm.zout_sur);
    if ((output->atm.zout_sur    = calloc (ntmp_zd, sizeof (float))) == NULL) {
      fprintf (stderr,"Error, allocating memory for zout_sur \n");
      fprintf (stderr, "      (line %d, function %s in %s)  \n", __LINE__, __func__, __FILE__);
      return -1;
    }
    output->atm.nzout_user = output->atm.nzout;
    output->atm.nzout = ntmp_zd;
    
    for (lu=0; lu<output->atm.nzout; lu++)  
      output->atm.zout_sur[lu] = tmp_zd[lu];


    /* Get indices of all user altitudes */
    if ((output->atm.zout_user_index = calloc (output->atm.nzout_user, sizeof (int))) == NULL) {
      fprintf (stderr,"Error, allocating memory for output->atm.zout_user_index (line %d, function %s in %s)\n",
                       __LINE__, __func__, __FILE__);
      return -1;
    }

    if ((output->atm.zout_comp_index = calloc (output->atm.nzout, sizeof (int))) == NULL) {
      fprintf (stderr,"Error, allocating memory for output->atm.zout_comp_index (line %d, function %s in %s)\n",
                       __LINE__, __func__, __FILE__);
      return -1;
    }

    lk = 0;
    for (lu=0;lu<output->atm.nzout_user;lu++) {
      for (li=0;li<output->atm.nzout;li++) {
	if (fabs(output->atm.zout_sur[li]-tmp_zd_user[lu]) < 0.001 ) {
	  output->atm.zout_user_index[lk++] = li;   /* These are the indices for the user altitudes */
	}
      }
    }
    free(tmp_zd);
    free(tmp_zd_user);

    /* Get the indices for the computational altitudes */
    
    lk=0;
    li=0;
    for (lu=0;lu<output->atm.nlev;lu++) {
      for (li=0;li<output->atm.nzout;li++) {
	/*	if (fabs(output->atm.zout_sur[li]-output->atm.zd[lu]) < 0.001 ) { */
	if (fabs(output->atm.zout_sur[li]-tmp_zd_rgrid[lu]) < 0.001 && lk < output->atm.nzout) {
	  output->atm.zout_comp_index[lk++] = li;   /* These are the indices for the computational altitudes */
	}
      }
    }    
    
    free(tmp_zd_rgrid);
  } /* END   if (input.raman) */

  /* last but not least consider altitude */
  if (output->alt.altitude != 0) {

    for (ialt=0; ialt<output->atm.nlev_common; ialt++)
      if (output->alt.altitude == output->atm.zd_common[ialt])
	break;
    
    if (ialt == output->atm.nlev_common) {
      fprintf (stderr, "Error, user-defined altitude %f not found in common grid (line %d, function %s in %s)\n", 
                               output->alt.altitude, __LINE__, __func__, __FILE__);
      return -1;
    }
    
    for (lc=0; lc<output->atm.nlev_common; lc++) {
      output->atm.zd_common[lc] -= output->alt.altitude;
      output->atm.zd[lc]        -= output->alt.altitude;
    }

    output->atm.nlev_common = ialt+1;
    output->atm.nlyr_common = ialt;
    output->atm.nlev        = ialt+1;
    output->atm.nlyr        = ialt;

  }

  /* something special here: we need to set the surface temperature  */
  /* to the atmospheric temperature at altitude; this has to be done */
  /* here rather than in setup_temperature      WHY???? UH           */
  if (output->surface_temperature < 0) {
    output->surface_temperature = output->atm.microphys.temper[0][0][output->atm.nlev-1];
    
    if (input.verbose)
      fprintf (stderr, " ... setup surface temperature\n");

    if (strlen(input.filename[FN_SURFACE_TEMP_MAP]) > 0) {
      /* read surface temperature from an netCDF input map */
      status = get_number_from_netCDF_map (input.latitude, input.longitude, input.UTC, input.atm.time_interpolate,
                                           input.filename[FN_SURFACE_TEMP_MAP], &(output->surface_temperature), TYPE_FLOAT, 
                                           input.netCDF_name_surf_T,
                                           input.verbose, input.quiet);
      if (status!=0) {
        fprintf (stderr, "Error %d reading surface temperature map %s\n", status, input.filename[FN_SURFACE_TEMP_MAP]);
        fprintf (stderr, "      (line %d, function %s in %s)  \n", __LINE__, __func__, __FILE__);
        return status;
      }
    }
    
    if (input.verbose)
      fprintf (stderr, " ... setting surface temperature to %f K\n", output->surface_temperature);
  }


  /* Output vertical columns */
  if (input.verbose) {
    
    /* loop over wavelength */
    for (iv=output->wl.nlambda_rte_lower; iv<=output->wl.nlambda_rte_upper; iv++) {
      
      /* reset data */
      vcol_o3   = 0; 
      vcol_co2  = 0; 
      vcol_o2   = 0; 
      vcol_no2  = 0; 
      vcol_bro  = 0; 
      vcol_oclo = 0; 
      vcol_hcho = 0; 

      babsotot_o3   = 0; 
      babsotot_co2  = 0; 
      babsotot_o2   = 0; 
      babsotot_no2  = 0; 
      babsotot_bro  = 0; 
      babsotot_oclo = 0; 
      babsotot_hcho = 0; 
      babso_amftot  = 0;
      
      vertical_column=0;
      
      fprintf (stderr, "\n*** wavelength: iv = %d, %f nm\n", iv, output->wl.lambda_r[iv]);
      fprintf (stderr, "*** setup_redistribute()\n");
      fprintf (stderr, " --------------------------------------------------------------------------------------------------------------------------------------------------------\n");
      fprintf (stderr, "   lc |      z[km] |\n");
      fprintf (stderr, "      |            |           o3             o2            co2            no2            bro           oclo           hcho        wc.dtau        ic.dtau\n");
      fprintf (stderr, " --------------------------------------------------------------------------------------------------------------------------------------------------------\n");

      /*fprintf(stderr, "3DAbs redistribute atm.nlyr %d \n", output->atm.nlyr); */
      
      for (lc=0; lc<output->atm.nlyr; lc++) {
	
	/* The factor 1.e+5 converts z from km to cm.*/
	deltaz= (output->atm.zd[lc] - output->atm.zd[lc+1])*1e5;
	
	babso_o3 = 0.5*deltaz * (output->atm.microphys.dens[MOL_O3][0][0][lc]*output->crs.crs[0][0][lc][MOL_O3][iv]
				 + output->atm.microphys.dens[MOL_O3][0][0][lc+1]*output->crs.crs[0][0][lc+1][MOL_O3][iv]);
	if (output->atm.microphys.denstab_id != MOL_O3) babso_md += babso_o3;
	vcol_o3    += 0.5*deltaz * (output->atm.microphys.dens[MOL_O3][0][0][lc] + output->atm.microphys.dens[MOL_O3][0][0][lc+1]);
	

	babso_co2 = 0.5*deltaz * (output->atm.microphys.dens[MOL_CO2][0][0][lc]*output->crs.crs[0][0][lc][MOL_CO2][iv]
				 + output->atm.microphys.dens[MOL_CO2][0][0][lc+1]*output->crs.crs[0][0][lc+1][MOL_CO2][iv]);
	if (output->atm.microphys.denstab_id != MOL_CO2) babso_md += babso_co2;
	vcol_co2    += 0.5*deltaz * (output->atm.microphys.dens[MOL_CO2][0][0][lc] + output->atm.microphys.dens[MOL_CO2][0][0][lc+1]);
	

	babso_o2 = 0.5*deltaz * (output->atm.microphys.dens[MOL_O2][0][0][lc]*output->crs.crs[0][0][lc][MOL_O2][iv]
				 + output->atm.microphys.dens[MOL_O2][0][0][lc+1]*output->crs.crs[0][0][lc+1][MOL_O2][iv]);
	if (output->atm.microphys.denstab_id != MOL_O2) babso_md += babso_o2;
	vcol_o2    += 0.5*deltaz * (output->atm.microphys.dens[MOL_O2][0][0][lc] + output->atm.microphys.dens[MOL_O2][0][0][lc+1]);
	

	babso_no2 = 0.5*deltaz * (output->atm.microphys.dens[MOL_NO2][0][0][lc]*output->crs.crs[0][0][lc][MOL_NO2][iv]
				  + output->atm.microphys.dens[MOL_NO2][0][0][lc+1]*output->crs.crs[0][0][lc+1][MOL_NO2][iv]);
	if (output->atm.microphys.denstab_id != MOL_NO2) babso_md += babso_no2;
	vcol_no2   += 0.5*deltaz * (output->atm.microphys.dens[MOL_NO2][0][0][lc] + output->atm.microphys.dens[MOL_NO2][0][0][lc+1]);
	

	babso_bro = 0.5*deltaz * (output->atm.microphys.dens[MOL_BRO][0][0][lc]*output->crs.crs[0][0][lc][MOL_BRO][iv]
				  + output->atm.microphys.dens[MOL_BRO][0][0][lc+1]*output->crs.crs[0][0][lc+1][MOL_BRO][iv]);
	if (output->atm.microphys.denstab_id != MOL_BRO) babso_md += babso_bro;
	vcol_bro   += 0.5*deltaz * (output->atm.microphys.dens[MOL_BRO][0][0][lc] + output->atm.microphys.dens[MOL_BRO][0][0][lc+1]);
	

	babso_oclo = 0.5*deltaz * (output->atm.microphys.dens[MOL_OCLO][0][0][lc]*output->crs.crs[0][0][lc][MOL_OCLO][iv]
				   + output->atm.microphys.dens[MOL_OCLO][0][0][lc+1]*output->crs.crs[0][0][lc+1][MOL_OCLO][iv]);
	if (output->atm.microphys.denstab_id != MOL_OCLO) babso_md += babso_oclo;
	vcol_oclo  += 0.5*deltaz * (output->atm.microphys.dens[MOL_OCLO][0][0][lc] + output->atm.microphys.dens[MOL_OCLO][0][0][lc+1]);
	

	babso_hcho = 0.5*deltaz * (output->atm.microphys.dens[MOL_HCHO][0][0][lc]*output->crs.crs[0][0][lc][MOL_HCHO][iv]
				   + output->atm.microphys.dens[MOL_HCHO][0][0][lc+1]*output->crs.crs[0][0][lc+1][MOL_HCHO][iv]);
	if (output->atm.microphys.denstab_id != MOL_HCHO)  babso_md += babso_hcho;
	vcol_hcho  += 0.5*deltaz * (output->atm.microphys.dens[MOL_HCHO][0][0][lc] + output->atm.microphys.dens[MOL_HCHO][0][0][lc+1]);
	
	babsotot_o3   += babso_o3;
	babsotot_co2  += babso_co2;
	babsotot_o2   += babso_o2;
	babsotot_no2  += babso_no2;
	babsotot_bro  += babso_bro;
	babsotot_oclo += babso_oclo;
	babsotot_hcho += babso_hcho;
	babso_amftot  += babso_amf;
	
        switch(input.ck_scheme) {
        case CK_CRS:
        case CK_RAMAN:
        case CK_REPTRAN:
        case CK_REPTRAN_CHANNEL:
          switch (output->atm.molabs) {
          case MOLABS_CALC:
          case MOLABS_LOOKUP:
	    fprintf (stderr, "%5d | %10.4f | %12.6e | %12.6e | %12.6e | %12.6e | %12.6e | %12.6e | %12.6e | ",
		     lc, output->atm.zd[lc+1]+output->alt.altitude,
		     babso_o3, babso_o2, babso_co2, babso_no2,
		     babso_bro, babso_oclo, babso_hcho);
	    if (input.i_wc!=-1 && output->caoth[input.i_wc].optprop.dtau!=NULL) /* SBCA should be if 1d */
	      fprintf (stderr, "%12.6e | ", output->caoth[input.i_wc].optprop.dtau[iv][lc]);
	    else
	      fprintf (stderr, "%12.6e | ", 0.0);
	    if (input.i_ic!=-1 && output->caoth[input.i_ic].optprop.dtau!=NULL) /* SBCA should be if 1d */
	      fprintf (stderr, "%12.6e\n", output->caoth[input.i_ic].optprop.dtau[iv][lc]);
	    else
	      fprintf (stderr, "%12.6e\n", 0.0);
	    break;
          case MOLABS_FILE_MONO:
          case MOLABS_FILE_SPEC:
	  case MOLABS_NONE:
	    break;
          default:
            fprintf (stderr, "Error, unknown molecular absorption option %d\n", output->atm.molabs);
            fprintf (stderr, "      (line %d, function %s in %s)  \n", __LINE__, __func__, __FILE__);
            return -1;
          }
	  break;
        case CK_KATO:
        case CK_KATO2:
        case CK_KATO2_96:
	case CK_KATO2ANDWANDJI:
        case CK_FU:
        case CK_AVHRR_KRATZ:
        case CK_FILE:
        case CK_LOWTRAN:
	  break;
        default:
          fprintf (stderr, "Error: unimplemented correlated-k scheme\n");
          fprintf (stderr, "      (line %d, function %s in %s)  \n", __LINE__, __func__, __FILE__);
          return -1;
          break;
        }

      }
      fprintf (stderr, " ---------------------------------------------------------------------------------------------------------------------------\n");
      fprintf (stderr, "%6s  %10.4f | %12.6e | %12.6e | %12.6e | %12.6e | %12.6e | %12.6e | %12.6e\n", 
	       "AMFsum", 0.0/0.0,
	       babsotot_o3, babsotot_o2, babsotot_co2, babsotot_no2, babsotot_bro, babsotot_oclo, babsotot_hcho);
      fprintf (stderr, " ---------------------------------------------------------------------------------------------------------------------------\n");
      fprintf (stderr, "AMF column    %4.0f | %12.6e | %12.6e | %12.6e | %12.6e | %12.6e | %12.6e | %12.6e\n", 
	       0.0/0.0, vcol_o3, vcol_o2, vcol_co2, vcol_no2, vcol_bro, vcol_oclo, vcol_hcho);
      fprintf (stderr, " ---------------------------------------------------------------------------------------------------------------------------\n");
      fprintf (stderr, "AMF BOA crs   %4.0f | %12.6e | %12.6e | %12.6e | %12.6e | %12.6e | %12.6e | %12.6e\n", 
	       0.0/0.0,
	       output->crs.crs[0][0][output->atm.nlyr-1][MOL_O3]  [iv], 
	       output->crs.crs[0][0][output->atm.nlyr-1][MOL_O2] [iv], 
	       output->crs.crs[0][0][output->atm.nlyr-1][MOL_CO2]  [iv], 
	       output->crs.crs[0][0][output->atm.nlyr-1][MOL_NO2] [iv], 
	       output->crs.crs[0][0][output->atm.nlyr-1][MOL_BRO] [iv], 
	       output->crs.crs[0][0][output->atm.nlyr-1][MOL_OCLO][iv],
	       output->crs.crs[0][0][output->atm.nlyr-1][MOL_HCHO][iv]);
      fprintf (stderr, " ---------------------------------------------------------------------------------------------------------------------------\n");
      if (output->atm.microphys.denstab_id) {
	fprintf(stderr,"Density matrix absorption column for airmass calculations: %16.10e %16.10e %16.10e\n",
		babso_amftot, vertical_column, crs_amf);
	fprintf (stderr, " ---------------------------------------------------------------------------------------------------------------------------\n");
      }
    }
  }

  return status;
}

int compare_float (void *ap, void *bp) { /* Sort such that altitude is in descending order */
  float *a = (float *) ap;
  float *b = (float *) bp;
  
  if (*a == *b) 
    return 0;

  if (*a > *b) 
    return -1;
  else 
    return 1;
}

/****************************************************************/
/* Regrid data from an input grid to an output grid. All levels */
/* of the input grid must also be available in the output grid; */
/* the scale flag toggles between distributing a quantity over  */
/* the output levels (e.g. an optical thickness, scale==1)      */
/* or simply copying it (e.g. an asymmetry parameter, scale==0) */
/****************************************************************/

int regrid_float (float *xin, float  *yin,  int nin,
		  float *xout, float *yout, int nout,
		  int scale)
{
  int iin=0, iout=0;

  /* first check if both grids are sorted in descending order */
  for (iin=0; iin<nin-1; iin++)
    if (xin[iin] <= xin[iin+1]) {
      fprintf (stderr, "Error, input grid not sorted in descending order! (line %d, function %s in %s)\n",
                        __LINE__, __func__, __FILE__);
      fprintf (stderr, "xin[%d] = %f, xin[%d] = %f\n", iin, xin[iin], iin+1, xin[iin+1]);
      return -1;
    }

  for (iout=0; iout<nout-1; iout++)
    if (xout[iout] <= xout[iout+1]) {
      fprintf (stderr, "Error, output grid not sorted in descending order! (line %d, function %s in %s)\n", 
                       __LINE__, __func__, __FILE__);
      fprintf (stderr, "xout[%d] = %f, xout[%d] = %f\n", iout, xout[iout], iout+1, xout[iout+1]);
      return -1;
    }

  /* first check if all grid points of the input grid */
  /* are also available in the output grid            */

  iout=nout-1;
  for (iin=nin-1; iin>=0; iin--) {
    while (xout[iout]<xin[iin]) 
      if (--iout<0)
	break;
      
    if (xin[iin]!=xout[iout]) {
      fprintf (stderr, "Error, level %f from input not contained in output! (line %d, function %s in %s)\n",
	       xin[iin], __LINE__, __func__, __FILE__);
      return -1;
    }
  }
  
  iin  = nin -1;
  iout = nout-1;

  /* go up to first output level */
  while (xout[iout]<xin[iin]) 
    iout--;

  for (; iout>=0; iout--) {

    while (xout[iout]<xin[iin-1]) {
      if (scale) 
	yout[iout] = yin[iin] * 
	  (xout[iout-1]-xout[iout]) / (xin[iin-1] - xin[iin]);
      else 
	yout[iout] = yin[iin];
      iout--;
    }

    iout++;
    iin--;
    
    if (iin<1)
      break;
  }

  return 0;
}


/*****************************************************************/
/* Regrid data from an input grid to an output grid. All levels  */
/* of the input grid must also be available in the output grid;  */
/* the scale flag toggles between distributing a quantity over   */
/* the output levels (e.g. an optical thickness, scale==1)       */
/* or simply copying it (e.g. an asymmetry parameter, scale==0); */
/*****************************************************************/

int regrid (double *xin, double  *yin,  int nin,
	    double *xout, double *yout, int nout,
	    int scale)
{
  int iin=0, iout=0;

  /* first check if both grids are sorted in descending order */
  for (iin=0; iin<nin-1; iin++)
    if (xin[iin] <= xin[iin+1]) {
      fprintf (stderr, "Error, input grid not sorted in descending order! (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      fprintf (stderr, "xin[%d] = %f, xin[%d] = %f\n", iin, xin[iin], iin+1, xin[iin+1]);
      return -1;
    }

  for (iout=0; iout<nout-1; iout++)
    if (xout[iout] <= xout[iout+1]) {
      fprintf (stderr, "Error, output grid not sorted in descending order! (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      fprintf (stderr, "xout[%d] = %f, xout[%d] = %f\n", iout, xout[iout], iout+1, xout[iout+1]);
      return -1;
    }

  /* second check if all grid points of the input grid */
  /* are also available in the output grid            */

  iout=nout-1;
  for (iin=nin-1; iin>=0; iin--) {
    while (xout[iout]<xin[iin]) 
      if (--iout<0)
	break;
     
    /* this was originally xin[iin] != xout[iout], but this caused rounding error with radiosonde z-grids */
    
    if ( fabs ( xin[iin] - xout[iout] ) > EPSILON  * fabs (xout[iout]))   {
      fprintf (stderr, "Error, level %16.12f from input not contained in output! (line %d, function %s in %s)\n",
               xin[iin], __LINE__, __func__, __FILE__);
      fprintf (stderr, "   input level = %16.12f, nearest output levels = %16.12f, %16.12f\n",
               xin[iin], xout[iout+1], xout[iout]);
      if (xin[iin] == xout[iout])
	fprintf (stderr, "  but both numbers agree when compared with ==\n");
      
      return -1;
    }
    
    /* if floating point rounding variance than set values to equal value */
    /* changes values only in this function, as no pointer are given to the function */
    if ( xin[iin] != xout[iout] && fabs(xin[iin] - xout[iout]) < EPSILON * fabs (xout[iout]) )
      xin[iin] = xout[iout];
    
  }
  
  
  iin  = nin -1;
  iout = nout-1;

  /* go up to first output level */
  while (xout[iout]<xin[iin]) 
    iout--;

  for (; iout>=0; iout--) {

    while (xout[iout]<xin[iin-1]) {
      if (scale) 
	yout[iout] = yin[iin] * 
	  (xout[iout-1]-xout[iout]) / (xin[iin-1] - xin[iin]);
      else 
	yout[iout] = yin[iin];
      iout--;
    }

    iout++;
    iin--;
    
    if (iin<1)
      break;
  }

  return 0;
}


/* copy array source to target and round to 6 digits */
void cp_round (float *source, float *target, int nlev)
{
  int lc=0;
  float pa=0;

  for (lc=0; lc<nlev; lc++) {
    if (source[lc]<=0.0)
      target[lc] = 0.0;
    else {
      pa=pow(10.0,floor(log10(source[lc]))-6.0);
      target[lc] = floor(source[lc]/pa)*pa;
    }
  }
}


void set_common_z_grid (float *zd, int nlev, float **zd_common, int *nlev_common) {

  int i=0, k=0, l=0, included=0;
  float *tmp = NULL;
  float z_min=NOT_DEFINED_FLOAT;
  float z_max=NOT_DEFINED_FLOAT;

  if ( (*nlev_common) != 0 ) {
    if ((*zd_common)[*nlev_common-1] < (*zd_common)[0]) {
      z_min=(*zd_common)[*nlev_common-1];
      z_max=(*zd_common)[0];
    }
    else {
      z_min=(*zd_common)[0];
      z_max=(*zd_common)[*nlev_common-1];    
    }
  }

  if ((tmp = (float *) calloc(nlev+*nlev_common, sizeof(float)))==NULL) {
    fprintf (stderr, "Error: Allocation of tmp in set_common_z_grid() (redistribute.c) \n");
    /* return -1; */
  }

  /* add the existing common grid to tmp */
  for (k=0; k<*nlev_common; k++)
    tmp[k] = (*zd_common)[k];

  /* cp_round (*zd_common, tmp, *nlev_common); */

  for (i=0; i<nlev; i++) {

    /* Is this level already included? */
    included = 0;
    for (l=0; l<k; l++) {
      if ( fabs ( tmp[l] - zd[i] ) <= EPSILON * fabs(zd[i]) ) {
	included=1;
	break;
      }
    }

    /* if not already included then add */
    if (!included) {
      /* if zd is inside the range of zd_common, if not called for the first time              */
      /* e.g. add z-grid of the ozone sonde to atmosphere grid                                 */
      /* for first time (*nlev_common==0) it is OK to add level outside the range of zd_common */
      if ( (*nlev_common)!=0 && z_min < zd[i] && zd[i] < z_max ) {
        tmp[k++] = zd[i];
      }
    }
  }

  *nlev_common = k;

  if (*nlev_common==0) {   /* Special for first time, and no need to sort this */
    (*zd_common) = (float *) calloc (nlev, sizeof(float));

    for (i=0;i<nlev;i++) 
      (*zd_common)[i] = zd[i];      

    *nlev_common = nlev;
  }
  else {
    free(*zd_common);
    (*zd_common) = (float *) calloc (*nlev_common, sizeof(float));
    for (i=0; i<*nlev_common; i++) 
      (*zd_common)[i] = tmp[i];

    qsort ((*zd_common), *nlev_common, sizeof(float), 
	   QSORT_CAST compare_float);

  }

  free(tmp);

}



/* add z levels from a file to an existing grid; */
/* reallocate memory as needed.                  */
int add_file2grid (float **zd, int *nlev, char *filename)
{
  int status=0;
  float *znew=NULL;
  int nnew=0;

  status = read_1c_file_float (filename, &znew, &nnew);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }
  
  /* add levels to aerosol profile */
  set_common_z_grid (znew, nnew, zd, nlev);
  free (znew);

  return 0;
}

/* same as add_file2grid, but from netCDF file */
/* that contains a sigma coordiante            */
int add_sigma_nc_file2grid (float **zd, int *nlev, char *filename, 
                            float latitude, float longitude, struct tm UTC, int time_interpolate,
                            float *z_atm, float *press, int nlev_atm,
                            int verbose, int quiet)
{
  int status=0;

  int ncid=NOT_DEFINED_INTEGER;

  long ilat=NOT_DEFINED_INTEGER,ilon=NOT_DEFINED_INTEGER;                                        
  size_t nlat  = 0;
  size_t nlon  = 0;
  double *lat_grid  = NULL;
  double *lon_grid  = NULL;
  int nt=-1;
  int itime1 = -1, itime2 = -1;
  float dt = NOT_DEFINED_FLOAT;

  size_t nnew=0;              /* number of layers in the netCDF file */
  double *sigma_grid=NULL;
  float *log_p_aero=NULL;
  float *log_p_atm =NULL;

  int lc=NOT_DEFINED_INTEGER;

  float *znew=NULL;
  float tmp;

  char function_name[] = "add_sigma_nc_file2grid";
  char file_name[]     = "redistribute.c";

  status = get_all_netCDF_indices ( filename, latitude, longitude,
                                    &(ncid), &(ilat), &(nlat), &(ilon), &(nlon), 
                                    &(lat_grid), &(lon_grid),
                                    UTC, time_interpolate,
                                    &(nt), &(itime1), &(itime2), &(dt),
                                    verbose, quiet);

  /* read number of layers */
  status = alloc_and_read_netCDF_1D_double(ncid,"ilev", &nnew, "sigma", &(sigma_grid));
  if (status!=0) {
    fprintf (stderr, "Error %d reading 'ilev' from %s in %s (%s)\n", status, filename, function_name, file_name );
    return status;
  }

  if ((log_p_aero = calloc (nnew, sizeof (float))) == NULL) {
    fprintf (stderr,"Error, allocating memory for 'log_p_aero' in %s (%s)\n", function_name, file_name);
    return -10;
  }

  for (lc=0; lc<nnew; lc++) {
    /* convert from sigma(=p/p0) -> log (sigma) */
    log_p_aero[lc] = log(sigma_grid[lc]);
  }

  free(sigma_grid);

  /* alloc one or two timesteps for pressure */
  if ((log_p_atm = calloc (nlev_atm, sizeof (float))) == NULL) {
    fprintf (stderr,"Error, allocating memory for 'log_p_atm' in %s (%s)\n", function_name, file_name);
    return -10;
  }

  for (lc=0; lc<nlev_atm; lc++) {
    /* convert from sigma(=p/p0) -> log (sigma) */
    log_p_atm[lc] = log(press[lc]/press[nlev_atm-1]);
  }

  if ((znew = calloc (nnew, sizeof (float))) == NULL) {
    fprintf (stderr,"Error, allocating memory for 'znew' in %s (%s)\n", function_name, file_name);
    return -10;
  }

  /* linear interpolation of z on the log(p/p0)-grid */
  status += arb_wvn (nlev_atm, log_p_atm,  z_atm, 
                     (int) nnew, log_p_aero, znew, 
                     INTERP_METHOD_LINEAR, 0);
  if (status!=0) {
    fprintf (stderr,"Error, during interpolation of 'z' in %s (%s)\n", function_name, file_name);
    return status;
  }

  free(log_p_aero);
  free(log_p_atm);

  /* reverse order of znew */
  for (lc=0; lc<nnew/2; lc++) {
    tmp             = znew[nnew-1-lc];
    znew[nnew-1-lc] = znew[lc];
    znew[lc]        = tmp;
  }

  /* correct small rounding errors for lowermost level */
  if (fabs(znew[nnew-1]-z_atm[nlev_atm-1]) < 0.00000001 )
    znew[nnew-1]=z_atm[nlev_atm-1];

  /* add levels to aerosol profile */
  set_common_z_grid (znew, nnew, zd, nlev);
  free (znew);

  return status;

}


/****************************************************************/
/* Regrid a 1D float array data from a given input grid (z_old) */
/* to an given output grid. All levels of the input grid must   */
/* also be available in the output grid; the scale flag toggles */
/* between distributing a quantity over the output levels       */
/* (e.g. an optical thickness, scale==1) or simply copying it   */
/* (e.g. an asymmetry parameter, scale==0); memory for data is  */
/* reallocated automatically.                                   */
/****************************************************************/

int redistribute_1D (void **data, 
		     float *z_old, int nlyr_old, 
		     float *z_new, int nlyr_new,
		     int scale, int type, int alloc) 
{
  int lc=0, status=0;
  double *xin =NULL, *yin =NULL;
  double *xout=NULL, *yout=NULL;

  double *tmp=NULL; 

  xin  = (double *) calloc (nlyr_old+1, sizeof(double));
  yin  = (double *) calloc (nlyr_old+1, sizeof(double));

  xout = (double *) calloc (nlyr_new+2, sizeof(double));
  yout = (double *) calloc (nlyr_new+2, sizeof(double));

  /* Allocate memory for temporary array */
  tmp = (double *) calloc (nlyr_old+1, sizeof(double));

  /* copy original stuff to temporary array */
  /* and reallocate with changed resolution */

  switch (type) {
  case REDISTRIBUTE_FLOAT:
    for (lc=0; lc<nlyr_old; lc++)
      tmp[lc] = ((float *)  (*data))[lc];
    
    if(alloc){
    /* free original array and reallocate with changed resolution */
    free (*data);
    *data = calloc (nlyr_new, sizeof(float));
    }
    break;

  case REDISTRIBUTE_DOUBLE:
    for (lc=0; lc<nlyr_old; lc++)
      tmp[lc] = ((double *) (*data))[lc];

    if(alloc){
      free (*data);
      *data = calloc (nlyr_new, sizeof(double));
    }
    break;
    
  default:
    fprintf (stderr, "Error, unknown data type %d in redistribute_1D()\n", type);
    return -1;
  }

  xin[0] = z_old[0];
  yin[0] = tmp[0];
  for (lc=0; lc<nlyr_old; lc++) {
    xin[lc+1]  = z_old[lc+1];
    yin[lc+1]  = tmp[lc];
  }

  for (lc=0; lc<=nlyr_new; lc++)
    xout[lc] = z_new[lc];
  
  status = regrid (xin, yin, nlyr_old+1, xout, yout, nlyr_new+1, scale);
  
  if (status != 0) {
    fprintf(stderr, "redistribute_1D(): Error %d returned by regrid (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
    return status;
  }

  switch (type) {
  case REDISTRIBUTE_FLOAT:
    for (lc=1; lc<=nlyr_new; lc++)
      ((float *) (*data))[lc-1] = yout[lc];      
    break;

  case REDISTRIBUTE_DOUBLE:
    for (lc=1; lc<=nlyr_new; lc++)
      ((double *) (*data))[lc-1] = yout[lc];      
    break;
    
  default:
    fprintf (stderr, "Error, unknown data type %d in redistribute_1D (line %d, function %s in %s)\n", type, __LINE__, __func__, __FILE__);
    return -1;
  }

  free (tmp); 
  free(xin); free(xout);
  free(yin); free(yout);

  return 0;
}



/****************************************************************/
/* Regrid a 2D float array data from a given input grid (z_old) */
/* to an given output grid. All levels of the input grid must   */
/* also be available in the output grid; the scale flag toggles */
/* between distributing a quantity over the output levels       */
/* (e.g. an optical thickness, scale==1) or simply copying it   */
/* (e.g. an asymmetry parameter, scale==0); memory for data is  */
/* reallocated automatically.                                   */
/****************************************************************/

int redistribute_2D (void ***data, int nlambda, 
		     int nlambda_lower, int nlambda_upper,
		     float *z_old, int nlyr_old, 
		     float *z_new, int nlyr_new,
		     int scale, int type, int alloc) 
{
  int lc=0, iv=0, status=0;
  double *xin =NULL, *yin =NULL;
  double *xout=NULL, *yout=NULL;

  double **tmp=NULL; 


  xin  = (double *) calloc (nlyr_old+1, sizeof(double));
  yin  = (double *) calloc (nlyr_old+1, sizeof(double));

  xout = (double *) calloc (nlyr_new+1, sizeof(double));
  yout = (double *) calloc (nlyr_new+1, sizeof(double));



  /* Allocate memory for temporary array */
  if ((status = ASCII_calloc_double (&tmp, nlambda, nlyr_old)) != 0) {
    fprintf (stderr, "Error allocating %d x %d doubles for tmp\n",
	     nlambda, nlyr_old);
    return status;
  }

  /* copy original stuff to temporary array */
  /* and reallocate with changed resolution */
  switch (type) {
  case REDISTRIBUTE_FLOAT:
    for (iv=nlambda_lower; iv<=nlambda_upper; iv++)
      for (lc=0; lc<nlyr_old; lc++) 
	tmp[iv][lc] = ((float **) (*data))[iv][lc];
    
    if(alloc){
      ASCII_free_float ((float **) (*data), nlambda);
      if ((status = ASCII_calloc_float ((float ***) data, nlambda, nlyr_new)) != 0) {
        fprintf (stderr, "Error allocating memory for data\n");
        return status;
      }
    }
    break;

  case REDISTRIBUTE_DOUBLE:
    for (iv=nlambda_lower; iv<=nlambda_upper; iv++)
      for (lc=0; lc<nlyr_old; lc++) 
	tmp[iv][lc] = ((double **) (*data))[iv][lc];
    
    if(alloc){
      ASCII_free_double ((double **) (*data), nlambda);
    if ((status = ASCII_calloc_double ((double ***) data, nlambda, nlyr_new)) != 0) {
      fprintf (stderr, "Error allocating memory for data\n");
      return status;
    }
    }
    break;
    
  default:
    fprintf (stderr, "Error, unknown data type %d (line %d, function %s in %s)\n", type, __LINE__, __func__, __FILE__);
    return -1;
  }
  

  for (iv=nlambda_lower; iv<=nlambda_upper; iv++) {
    xin[0] = z_old[0];
    yin[0] = tmp[iv][0];
    for (lc=0; lc<nlyr_old; lc++) {
      xin[lc+1]  = z_old[lc+1];
      yin[lc+1]  = tmp[iv][lc];
    }
    
    for (lc=0; lc<=nlyr_new; lc++)
      xout[lc] = z_new[lc];
    
    status = regrid (xin, yin, nlyr_old+1,
		     xout, yout, nlyr_new+1, scale);

    if (status != 0) {
      fprintf(stderr, "Error %d returned by regrid (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
      return status;
    }

    switch (type) {
    case REDISTRIBUTE_FLOAT:
      for (lc=1; lc<=nlyr_new; lc++)
	(((float **) (*data))[iv][lc-1]) = yout[lc];      
      break;
      
    case REDISTRIBUTE_DOUBLE:
      for (lc=1; lc<=nlyr_new; lc++)
	(((double **) (*data))[iv][lc-1]) = yout[lc];      
      break;
      
    default:
      fprintf (stderr, "Error, unknown data type %d (line %d, function %s in %s)\n", type, __LINE__, __func__, __FILE__);
      return status;
    }
  }

  ASCII_free_double (tmp, nlambda);

  free(xin); free(xout);
  free(yin); free(yout);
  
  return 0;
}


/***********************************************************************/
/* Regrid a 3D MYSTIC float array [iz, ix, iy] from a given input grid */
/* (z_old) to a given output grid. PLEASE NOTE THAT THESE 3D FIELDS    */
/* ARE COUNTED UPWARDS FROM z=0 TO z=TOA, RATHER THAN DOWNWARDS WHICH  */
/* IS THE CASE FOR ALL OTHER LIBRADTRAN PROFILES.                      */
/* All levels of the input grid must also be available in the output   */ 
/* grid; the scale flag toggles between distributing a quantity over   */
/* the output levels (e.g. an optical thickness, scale==1) or simply   */
/* copying it (e.g. an asymmetry parameter, scale==0); memory data is  */
/* reallocated automatically.                                          */
/***********************************************************************/

int redistribute_3D (float ****data, int nx, int ny,
		     float *z_old, int nlyr_old, 
		     float *z_new, int nlyr_new,
		     int *threedold, int *threednew,
		     int scale, int alloc) 
{
  int lc=0, status=0;
  int ix=0, iy=0;
  double *xin =NULL, *yin =NULL;
  double *xout=NULL, *yout=NULL;

  float ***tmp=NULL; 

  xin  = (double *) calloc (nlyr_old+1, sizeof(double));
  yin  = (double *) calloc (nlyr_old+1, sizeof(double));

  xout = (double *) calloc (nlyr_new+1, sizeof(double));
  yout = (double *) calloc (nlyr_new+1, sizeof(double));


  if ((status = calloc_float_3D_sparse (&tmp, nlyr_old, nx, ny, threedold)) != 0) {
    fprintf (stderr, "Error allocating %d x %d x %d doubles for tmp\n",
	     nlyr_old, nx, ny);
    return status;
  }
 
  /* copy original stuff to temporary array */
  for (lc=0; lc<nlyr_old; lc++)
    if (threedold[lc]) {
      for (ix=0; ix<nx; ix++)
	for (iy=0; iy<ny; iy++)
	  tmp[lc][ix][iy] = (*data)[lc][ix][iy];
    }

  /* and reallocate with changed resolution */
  if(alloc) {
    free_float_3D_sparse (*data, nlyr_old, nx, threedold);
    if ((status = calloc_float_3D_sparse (data, nlyr_new, nx, ny, threednew)) != 0) {
      fprintf (stderr, "Error allocating memory for data\n");
      return status;
    }
  }
  
  for (ix=0; ix<nx; ix++) {
    for (iy=0; iy<ny; iy++) {
      xin[0] = z_old[nlyr_old];
      if (threedold[nlyr_old-1])
	yin[0] = tmp[nlyr_old-1][ix][iy];
      
      for (lc=0; lc<nlyr_old; lc++) {
	xin[lc+1]  = z_old[nlyr_old-lc-1];
	if (threedold[nlyr_old-lc-1])
	  yin[lc+1]  = tmp[nlyr_old-lc-1][ix][iy]; 
      }

      for (lc=0; lc<=nlyr_new; lc++)
	xout[lc] = z_new[lc];
      
      status = regrid (xin, yin, nlyr_old+1,
		       xout, yout, nlyr_new+1, scale);
      
      if (status != 0) {
	fprintf(stderr, "Error %d returned by regrid() in redistribute_2D (redistribute.c)\n", status);
	return status;
      }
      
      for (lc=0; lc<nlyr_new; lc++)
	if (threednew[lc])
	  (*data)[lc][ix][iy] = yout[nlyr_new-lc];      
    }
  }
    
  free_float_3D_sparse (tmp, nlyr_old, nx, threedold);
    
  free(xin); free(xout);
  free(yin); free(yout);
  
  return 0;
}



/****************************************************************/
/* Regrid the threed integer array [iz] from a given input grid */
/* (z_old) to a given output grid. All levels of the input grid */
/* must also be available in the output grid; the scale flag    */
/* toggles between distributing a quantity over the output      */
/* levels (e.g. an optical thickness, scale==1) or simply       */
/* copying it (e.g. an asymmetry parameter, scale==0); memory   */
/* for data is reallocated automatically.                       */
/****************************************************************/

int redistribute_threed (int **threed,
			 float *z_old, int nlyr_old, 
			 float *z_new, int nlyr_new,
			 int alloc) 
{
  int lc=0, status=0;
  double *xin =NULL, *yin =NULL;
  double *xout=NULL, *yout=NULL;

  int *tmp=NULL; 

  xin  = (double *) calloc (nlyr_old+1, sizeof(double));
  yin  = (double *) calloc (nlyr_old+1, sizeof(double));

  xout = (double *) calloc (nlyr_new+1, sizeof(double));
  yout = (double *) calloc (nlyr_new+1, sizeof(double));


  /* Allocate memory for temporary array */
  tmp = calloc (nlyr_old, sizeof(int));
  
  
  /* copy original stuff to temporary array */
  
  for (lc=0; lc<nlyr_old; lc++) 
    tmp[lc] = (*threed)[lc];
  
  /* and reallocate with changed resolution */
  if(alloc) {
    free (*threed);
    *threed = calloc(nlyr_new, sizeof(int));
  }
  
  
  xin[0] = z_old[nlyr_old];
  yin[0] = (double) tmp[nlyr_old-1];
  
  for (lc=0; lc<nlyr_old; lc++) {
    xin[lc+1]  = z_old[nlyr_old-lc-1];
    yin[lc+1]  = tmp[nlyr_old-lc-1]; 
  }
  
  for (lc=0; lc<=nlyr_new; lc++)
    xout[lc] = z_new[lc];
  
  status = regrid (xin, yin, nlyr_old+1,
		   xout, yout, nlyr_new+1, 0);
      
  if (status != 0) {
    fprintf(stderr, "Error %d returned by regrid() in redistribute_2D (redistribute.c)\n", status);
    return status;
  }
  
  for (lc=0; lc<nlyr_new; lc++)
    (*threed)[lc] = (int) (yout[nlyr_new-lc]+0.5);      
  
  free(tmp);
  
  free(xin); free(xout);
  free(yin); free(yout);
  
  return 0;
}
  


/*****************************************************************/
/* Regrid phase function moments from a given input grid (z_old) */
/* to an given output grid. All levels of the input grid must    */
/* also be available in the output grid; the scale flag toggles  */
/* between distributing a quantity over the output levels        */
/* (e.g. an optical thickness, scale==1) or simply copying it    */
/* (e.g. an asymmetry parameter, scale==0); memory for the       */
/* moments is reallocated automatically.                         */
/*****************************************************************/
static int redistribute_moment (float *****moment, int ***nmom, int nlambda, 
                                int nlambda_lower, int nlambda_upper,
				int nlyr_alloc,
                                float *z_old, int nlyr_old, int nphamat, 
                                float *z_new, int nlyr_new, int alloc) 
{
  int lc=0, iv=0, k=0, ip=0, status=0;
  int maxmom=0;
  double *xin =NULL, *yin =NULL;
  double *xout=NULL, *yout=NULL;

  double ***tmp_mom=NULL, ***tmp_newmom=NULL;
  int *tmp_nmom = calloc (nlyr_old, sizeof(int)); 

  xin  = (double *) calloc (nlyr_old+1, sizeof(double));
  yin  = (double *) calloc (nlyr_old+1, sizeof(double));

  xout = (double *) calloc (nlyr_new+1, sizeof(double));
  yout = (double *) calloc (nlyr_new+1, sizeof(double));


  for (iv=nlambda_lower; iv<=nlambda_upper; iv++) {
    
    /* copy original stuff to temporary array */
    for (lc=0; lc<nlyr_old; lc++) 
      tmp_nmom[lc] = (*nmom)[iv][lc];
    
    if(alloc){
      /* free original array and reallocate with changed resolution */
      free ((*nmom)[iv]);
      (*nmom)[iv] = calloc (nlyr_new, sizeof(int));
    }

    /* regrid number of moments */
    xin[0] = z_old[0];
    yin[0] = (float) tmp_nmom[0];
    for (lc=0; lc<nlyr_old; lc++) {
      xin[lc+1]  = z_old[lc+1];
      yin[lc+1]  = tmp_nmom[lc];
    }
    
    for (lc=0; lc<=nlyr_new; lc++)
      xout[lc] = z_new[lc];

    status = regrid (xin, yin, nlyr_old+1,
		     xout, yout, nlyr_new+1, 0);
    
    if (status != 0) {
      fprintf(stderr, "redistribute_moment_noalloc(): Error %d returned by regrid()\n", status);
      return status;
    }
    
    for (lc=1; lc<=nlyr_new; lc++)
      (*nmom)[iv][lc-1] = (int) (yout[lc] + 0.5);

    
    /* determine maximum number of moments */
    maxmom = 0;
    for (lc=0; lc<nlyr_old; lc++) 
      if (tmp_nmom[lc] > maxmom) 
	maxmom = tmp_nmom[lc];
    
    /* Allocate memory for temporary array of moments */
    tmp_mom = calloc (nlyr_old, sizeof (double **));
    for (lc=0; lc<nlyr_old; lc++){
      tmp_mom[lc] = calloc (nphamat, sizeof (double *));
      for(ip=0; ip<nphamat; ip++)
        tmp_mom[lc][ip] = calloc (maxmom+1, sizeof (double));
    }
   
    /* copy original stuff to temporary array */
    for (lc=0; lc<nlyr_old; lc++) 
      for (ip=0; ip<nphamat; ip ++) 
        for (k=0; k<tmp_nmom[lc]; k++) 
          tmp_mom[lc][ip][k] = (*moment)[iv][lc][ip][k];
    
    /* allocate temporary array for regridded moments */
    tmp_newmom = calloc (nlyr_new, sizeof (double **));
    for (lc=0; lc<nlyr_new; lc++){
      tmp_newmom[lc] = calloc (nphamat, sizeof (double *));
      for(ip=0; ip<nphamat; ip++)
        tmp_newmom[lc][ip] = calloc (maxmom+1, sizeof (double));
    }

    /* free original moments array and allocate with changed resolution */
    ASCII_free_float_3D ((*moment)[iv], nlyr_alloc, nphamat);
    (*moment)[iv] = calloc (nlyr_new, sizeof (float **));
    for (lc=0; lc<nlyr_new; lc++){
      (*moment)[iv][lc] = calloc (nphamat, sizeof (float *));
      for(ip=0; ip<nphamat; ip++){
        (*moment)[iv][lc][ip] = calloc ((*nmom)[iv][lc]+1, sizeof (float));
      }
    }
        
    for(ip=0; ip<nphamat; ip++){
      for (k=0; k<maxmom; k++) {
        xin[0] = z_old[0];
        yin[0] = tmp_mom[0][ip][k];
        for (lc=0; lc<nlyr_old; lc++) {
          xin[lc+1]  = z_old[lc+1];
          yin[lc+1]  = tmp_mom[lc][ip][k];
        }
        
        for (lc=0; lc<=nlyr_new; lc++)
          xout[lc] = z_new[lc];
          
        status = regrid (xin, yin, nlyr_old+1,
                         xout, yout, nlyr_new+1, 0);
        if (status != 0) {
          fprintf(stderr, "Error %d returned by regrid()\n", status);
            return status; 
        }
        
        for (lc=1; lc<=nlyr_new; lc++) {
          tmp_newmom[lc-1][ip][k] = yout[lc];
        }
      }
    }

    for (lc=0; lc<nlyr_new; lc++) 
      for (ip=0; ip<nphamat; ip++)
        for (k=0; k<(*nmom)[iv][lc]; k++) 
          (*moment)[iv][lc][ip][k] = (float) tmp_newmom[lc][ip][k];
    
      
    ASCII_free_double_3D (tmp_mom, nlyr_old, nphamat);
    ASCII_free_double_3D (tmp_newmom, nlyr_new, nphamat);
  
  }
  
  free(tmp_nmom);
  
  free(xin); free(xout);
  free(yin); free(yout);

  return 0;
}


/***********************************************************************************/
/* Function: redistribute_phase                                           @30_30i@ */
/* Description: Regrid phase function         from a given input grid (z_old)      */
/*              to an given output grid. All levels of the input grid must         */
/*              also be available in the output grid;                              */
/*              the phase fuznctions are not copied, but pointers are set          */
/*              in a correct way                                                   */
/*                                                                                 */
/* Parameters:                                                                     */
/*   In:    nlambda:                                                               */
/*          nlambda_lower:                                                         */
/*          nlambda_upper:                                                         */
/*          z_old:                                                                 */
/*          nlyr_old:                                                              */
/*          nphamat:                                                               */
/*          z_new:                                                                 */
/*          nlyr_new:                                                              */
/*          alloc:                                                                 */
/*                                                                                 */
/*   Inout: phase:                                                                 */
/*          theta:                                                                 */
/*          mu:                                                                    */
/*          ntheta:                                                                */
/*                                                                                 */
/* Return value: Status: 0 if everything ok                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                   @i30_30@ */
/***********************************************************************************/

static int redistribute_phase (float *****phase, float *****theta, double *****mu,
			       int ****ntheta, int nlambda, 
                               int nlambda_lower, int nlambda_upper,
			       int nlyr_alloc,
			       float *z_old, int nlyr_old, int nphamat, 
			       float *z_new, int nlyr_new, int alloc) 
{
  int status=0, iv=0, lc=0, lo=0, ip=0, in=0, repointed=0;
  float ****phase_old, ****theta_old;
  double ****mu_old;
  int ***ntheta_old;
  int *i_loc;

  /* save pointers to old arrays */
  phase_old  = *phase;
  theta_old  = *theta;
  mu_old     = *mu;
  ntheta_old = *ntheta;

  /* find indices; relate new grid to old grid */
  status = sort_vec (z_old, nlyr_old+1, z_new, nlyr_new+1, &(i_loc));

  if (status!=0) {
    fprintf (stderr, "Error %d sort_vec-ing\n", status);
    return status;
  }

  /* allocate */
  (*ntheta) = calloc(nlambda, sizeof(int    **));
  (*theta)  = calloc(nlambda, sizeof(float  ***));
  (*mu)     = calloc(nlambda, sizeof(double ***));
  (*phase)  = calloc(nlambda, sizeof(float  ***));

  for (iv=0; iv<nlambda; iv++) {
    (*ntheta)[iv] = calloc(nlyr_new, sizeof(int    *));
    (*theta) [iv] = calloc(nlyr_new, sizeof(float  **));
    (*mu)    [iv] = calloc(nlyr_new, sizeof(double **));
    (*phase) [iv] = calloc(nlyr_new, sizeof(float  **));
  }

  for (iv=0; iv<nlambda; iv++)
    for (lc=0; lc<nlyr_new; lc++) {
      (*ntheta)[iv][lc] = calloc(nphamat, sizeof(int    ));
      (*theta) [iv][lc] = calloc(nphamat, sizeof(float  *));
      (*mu)    [iv][lc] = calloc(nphamat, sizeof(double *));
      (*phase) [iv][lc] = calloc(nphamat, sizeof(float  *));

    }

  /* copy grid */
  for (iv=nlambda_lower; iv<=nlambda_upper; iv++) {
    lo=-1;
    repointed=0;
    for (lc=0; lc<nlyr_new; lc++) {

      /* old grid index */
      if (lo != i_loc[lc]) {
	/* new old grid index, set index and make sure we only repoint */
	lo=i_loc[lc];
	repointed=0;
      }

      for(ip=0; ip<nphamat; ip++) {

	/* no theta grid if above highest cloud or below lowest one */
	if (lo != -1)
	  (*ntheta)[iv][lc][ip] = ntheta_old[iv][lo][ip];
	else
	  (*ntheta)[iv][lc][ip] = 0;

	/* exit if no theta grid defined, we do not want to allocate! */
	if (!(*ntheta)[iv][lc][ip])
	  continue;

	if (!repointed) {
	  /* simply point to old arrays */
	  (*theta) [iv][lc][ip] = theta_old[iv][lo][ip];
	  (*mu)    [iv][lc][ip] = mu_old   [iv][lo][ip];
	  (*phase) [iv][lc][ip] = phase_old[iv][lo][ip];
	}
	else {
	  /* old array already pointed to by previous layer, need to allocate and copy now */
	  (*theta) [iv][lc][ip] = calloc((*ntheta)[iv][lc][ip], sizeof(float ));
	  (*mu)    [iv][lc][ip] = calloc((*ntheta)[iv][lc][ip], sizeof(double));
	  (*phase) [iv][lc][ip] = calloc((*ntheta)[iv][lc][ip], sizeof(float ));

	  for (in=0; in<(*ntheta)[iv][lc][ip]; in++) {
	    (*theta) [iv][lc][ip][in] = (*theta) [iv][lc-1][ip][in];
	    (*mu)    [iv][lc][ip][in] = (*mu)    [iv][lc-1][ip][in];
	    (*phase) [iv][lc][ip][in] = (*phase) [iv][lc-1][ip][in];
	  }
	}

      } /* ip loop */
      repointed=1;
    } /* lc loop */
  } /* iv loop */

  /* free old matrices */
  for (iv=0; iv<nlambda; iv++)
    for (lc=0; lc<nlyr_alloc; lc++) {
      free(ntheta_old[iv][lc]);
      free(theta_old [iv][lc]);
      free(mu_old    [iv][lc]);
      free(phase_old [iv][lc]);
    }

  for (iv=0; iv<nlambda; iv++) {
    free(ntheta_old[iv]);
    free(theta_old [iv]);
    free(mu_old    [iv]);
    free(phase_old [iv]);
  }

  free(i_loc);  /* Memory allocated by sort_vec */
  free(ntheta_old);
  free(theta_old);
  free(mu_old);
  free(phase_old);

  return 0;
}


/***********************************************************************************/
/* Function: sort_vec                                                     @62_30i@ */
/* Description:                                                                    */
/*  Locate new grid on old grid and return indices.                                */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int sort_vec (float *xin, int nin, float *xout, int nout, int **i_loc)
{
  int iin=0, iout=0;
  float *xin_c=NULL;
  char function_name[] = "sort_vec";
  char file_name[]     = "redistribute.c";

  xin_c = calloc(nin,sizeof(float));
  for (iin=0; iin<nin; iin++)
    xin_c[iin]=xin[iin];

  /* first check if both grids are sorted in descending order */
  for (iin=0; iin<nin-1; iin++)
    if (xin[iin] <= xin[iin+1]) {
      fprintf (stderr, "Error, input grid not sorted in descending order!\n");
      fprintf (stderr, "xin[%d] = %f, xin[%d] = %f\n", iin, xin[iin], iin+1, xin[iin+1]);
      return -1;
    }

  for (iout=0; iout<nout-1; iout++)
    if (xout[iout] <= xout[iout+1]) {
      fprintf (stderr, "Error, output grid not sorted in descending order!\n");
      fprintf (stderr, "xout[%d] = %f, xout[%d] = %f\n", iout, xout[iout], iout+1, xout[iout+1]);
      return -1;
    }

  /* first check if all grid points of the input grid */
  /* are also available in the output grid            */

  iout=nout-1;
  for (iin=nin-1; iin>=0; iin--) {
    while (xout[iout]<xin[iin]) 
      if (--iout<0)
	break;
     
    /* this was originally xin[iin] != xout[iout], but this caused rounding error with radiosonde z-grids */
    
    if ( fabs ( xin[iin] - xout[iout] ) > EPSILON  * fabs (xout[iout]))   {
      fprintf (stderr, "Error, level %16.12f from input not contained in output, %s (%s)\n",
               xin[iin], function_name, file_name);
      fprintf (stderr, "   input level = %16.12f, nearest output levels = %16.12f, %16.12f\n",
               xin[iin], xout[iout+1], xout[iout]);
      if (xin[iin] == xout[iout])
	fprintf (stderr, "  but both numbers agree when compared with ==\n");
      
      return -1;
    }
    
    /* if floating point rounding variance than set values to equal value */
    /* changes values only in this function, as no pointer are given to the function */
    if ( xin[iin] != xout[iout] && fabs(xin[iin] - xout[iout]) < EPSILON * fabs (xout[iout]) )
      xin_c[iin] = xout[iout];
    
  }
  
  (*i_loc) = calloc(nout-1, sizeof (int));

  /* until first cloud appears, no theta grid defined (-1) */
  iin=-1;
  
  for (iout=0; iout < nout-1; iout++) {
    /* reached next cloud level */
    if (xout[iout] == xin_c[iin+1])
        iin++;
    /* passed lowest cloud level */
    if (iin == nin-1)
      break;
    
    (*i_loc)[iout] = iin;
  }

  /* levels below lowest cloud should have no theta grid defined (-1) */
  for (; iout < nout-1; iout++)
    (*i_loc)[iout] = -1;


  free(xin_c);

  return 0;
}


/**********************/
/* Redistribute caoth */
/**********************/

int redistribute_caoth ( wl_out_struct     wl,
			 atm_out_struct    atm,
			 int               nphamat,
			 int               alloc_microphys,
			 int               copy_ssprop,
			 int               quiet,
			 /* input / output */
			 caoth_out_struct *caoth )
{
  caoth_out_struct *redistributed_caoth = calloc(1, sizeof(caoth_out_struct));
  float *caoth_zd=NULL;
  int lc=0, old_nlev=0, status=0;

  calloc_caoth_out ( redistributed_caoth,
		     (*caoth).name,
		     (*caoth).fullname,
		     wl.nlambda_r,
		     atm.nlev_common,
		     nphamat,
		     1 );

  status = cp_caoth_out ( redistributed_caoth,
			  *caoth,
			  copy_ssprop,
			  1,
			  quiet );
  if (status)
    return mem_err_out ( "output_caoth struct ", ERROR_POSITION );

  old_nlev = (*caoth).nlev;
  caoth_zd = calloc (old_nlev+1, sizeof (float));

  for (lc=0; lc<old_nlev; lc++)
    caoth_zd[lc] = (*caoth).zd[lc];

  free_caoth_out (caoth, 1);

  status = redistribute_optprop ( &(redistributed_caoth->optprop),
				  wl.nlambda_r, 
				  wl.nlambda_rte_lower, 
				  wl.nlambda_rte_upper,
				  caoth_zd,
				  old_nlev-1,
				  nphamat, 
				  atm.zd_common, 
				  atm.nlev_common-1,
				  0 );

  if (status!=0) {
    fprintf (stderr, "Error %d during redistribute_optprop (%s) (line %d, function %s in %s)\n",
	     status, (*caoth).fullname, __LINE__, __func__, __FILE__);
    return status;
  }

  /* also redistribute caoth microphysics (required by MYSTIC) */
  if (alloc_microphys) {
    status = redistribute_caoth_microphys ( &(redistributed_caoth->microphys),
					    caoth_zd,
					    old_nlev-1,
					    atm.zd_common,
					    atm.nlev_common-1 );

    if (status!=0) {
      fprintf (stderr, "Error %d regridding microphysical parameters for %s (line %d, function %s in %s)\n",
	       status, (*caoth).fullname, __LINE__, __func__, __FILE__ );
      return status;
    }
  }

  /* RPB speedup: for some strange reason,
     the following three function calls can NOT be replaced by:
     caoth = redistributed_caoth; */
  calloc_caoth_out ( caoth,
		     (*redistributed_caoth).name,
		     (*redistributed_caoth).fullname,
		     wl.nlambda_r,
		     atm.nlev_common,
		     nphamat,
		     1);

  status = cp_caoth_out ( caoth,
			  *redistributed_caoth,
			  copy_ssprop,
			  1,
			  quiet );

  if (status!=0) {
    fprintf (stderr, "Error %d copying redistributed %s to output structure (line %d, function %s in %s)\n", 
	     status, (*redistributed_caoth).fullname, __LINE__, __func__, __FILE__ );
    return status;
  }

  free_caoth_out(redistributed_caoth, 1);
  free(redistributed_caoth);

  free(caoth_zd);

  /* copy number of levels and re-write altitude grid */
  (*caoth).nlev = atm.nlev_common;
  free ((*caoth).zd);
  (*caoth).zd = calloc ((*caoth).nlev, sizeof(float));
  for (lc=0; lc<(*caoth).nlev; lc++)
    (*caoth).zd[lc] = atm.zd_common[lc];

  return 0;
}


/***********************************/
/* Redistribute optical properties */
/***********************************/

int redistribute_optprop (optprop_struct *optprop, int nlambda, 
			  int nlambda_lower, int nlambda_upper,
			  float *zd, int nlyr, int nphamat,  
			  float *zd_common, int nlyr_common, int alloc)
{
  int status=0;

  /** dtau **/
  status += redistribute_2D ((void *) (&(optprop->dtau)), nlambda, 
			     nlambda_lower, nlambda_upper,
			     zd, nlyr, 
			     zd_common, nlyr_common,
			     1, REDISTRIBUTE_FLOAT, alloc);
  if (status != 0) {
    fprintf(stderr, "Error %d during redistribute dtau (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
    return status;
  }

  /** ssa  **/
  status += redistribute_2D ((void *) (&(optprop->ssa)), nlambda, 
			     nlambda_lower, nlambda_upper,
			     zd, nlyr, 
			     zd_common, nlyr_common,
			     0, REDISTRIBUTE_FLOAT, alloc);
  if (status != 0) {
    fprintf(stderr, "Error %d during redistribute ssa (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
    return status;
  }

  /** f  **/
  status += redistribute_2D ((void *) (&(optprop->f)), nlambda, 
			     nlambda_lower, nlambda_upper,
			     zd, nlyr, 
			     zd_common, nlyr_common,
			     0, REDISTRIBUTE_FLOAT, alloc);
  if (status != 0) {
    fprintf(stderr, "Error %d during redistribute f (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
    return status;
  }

  /**  g1  **/
  status += redistribute_2D ((void *) (&(optprop->g1)), nlambda, 
			     nlambda_lower, nlambda_upper,
			     zd, nlyr, 
			     zd_common, nlyr_common,
			     0, REDISTRIBUTE_FLOAT, alloc);
  if (status != 0) {
    fprintf(stderr, "Error %d during redistribute g1 (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
    return status;
  }

  /**  g2  **/
  status += redistribute_2D ((void *) (&(optprop->g2)), nlambda, 
			     nlambda_lower, nlambda_upper,
			     zd, nlyr, 
			     zd_common, nlyr_common,
			     0, REDISTRIBUTE_FLOAT, alloc);
  if (status != 0) {
    fprintf(stderr, "Error %d during redistribute g2 (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
    return status;
  }

  /**  ff  **/
  status += redistribute_2D ((void *) (&(optprop->ff)), nlambda, 
			     nlambda_lower, nlambda_upper,
			     zd, nlyr, 
			     zd_common, nlyr_common,
			     0, REDISTRIBUTE_FLOAT, alloc);
  if (status != 0) {
    fprintf(stderr, "Error %d during redistribute ff (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
    return status;
  }

  /** dscale **/
  status += redistribute_2D ((void *) (&(optprop->dscale)), nlambda, 
			     nlambda_lower, nlambda_upper,
			     zd, nlyr, 
			     zd_common, nlyr_common,
			     0, REDISTRIBUTE_FLOAT, alloc);
  if (status != 0) {
    fprintf(stderr, "Error %d during redistribute dscale (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
    return status;
  }

  /** moments **/
  /* ??? Wendisch bug
  status += redistribute_moment (&(optprop->moment), &(optprop->nmom), nlambda, 
				 nlambda_lower, nlambda_upper,
				 zd, nlyr, 
				 zd_common, nlyr_common); ??? */
				 
  status += redistribute_moment (&(optprop->moment), &(optprop->nmom), nlambda, 
				 0, nlambda-1, optprop->nlev-1,
				 zd, nlyr, nphamat, 
				 zd_common, nlyr_common, alloc);
  if (status != 0) {
    fprintf(stderr, "Error %d during redistribute moment (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
    return status;
  }
  
  /** phase functions **/
  status += redistribute_phase (&(optprop->phase), &(optprop->theta), &(optprop->mu),
				&(optprop->ntheta), nlambda, 
				0, nlambda-1,
				optprop->nlev-1,
				zd, nlyr, nphamat, 
				zd_common, nlyr_common, alloc);
  if (status != 0) {
    fprintf(stderr, "Error %d during redistribute phase (line %d, function %s in %s)\n", status, __LINE__, __func__, __FILE__);
    return status;
  }
  
  return status;
}



/************************************/
/* Redistribute 3d caoth properties */
/************************************/

int redistribute_caoth3d ( caoth3d_out_struct *caoth3d, 
			   float              *zd,
			   int                 nlyr,
			   float              *zd_common,
			   int                 nlyr_common,
			   int                 alloc,
			   int                 quiet )
{
  int lc=0, status=0;
  int nx=caoth3d->Nx;
  int ny=caoth3d->Ny;
  int *threedold=NULL;

  if (!quiet)
    fprintf (stderr, " ... redistributing %d 3D layers -> %d 3D layers\n", nlyr, nlyr_common);

  /* first redistribute threed, but save a copy in threedold */
  threedold = calloc (nlyr_common, sizeof(int));
  for (lc=0; lc<nlyr; lc++)
    threedold[lc] = caoth3d->threed[lc];
  
  status += redistribute_threed((&(caoth3d->threed)),
				zd, nlyr, 
				zd_common, nlyr_common,
				alloc);

  if (status!=0) {
    fprintf (stderr, "Error %d redistributing threed\n", status);
    return status;
  }

  /** lwc **/
  status += redistribute_3D ((&(caoth3d->lwc)), nx, ny, 
			     zd, nlyr, 
			     zd_common, nlyr_common,
			     threedold, caoth3d->threed,
			     0, alloc);

  /** reff **/
  status += redistribute_3D ((&(caoth3d->reff)), nx, ny, 
			     zd, nlyr, 
			     zd_common, nlyr_common,
			     threedold, caoth3d->threed,
			     0, alloc);

  /** ext **/
  status += redistribute_3D ((&(caoth3d->ext)), nx, ny, 
			     zd, nlyr, 
			     zd_common, nlyr_common,
			     threedold, caoth3d->threed,
			     0, alloc);

  /** ssa **/
  status += redistribute_3D ((&(caoth3d->ssa)), nx, ny, 
			     zd, nlyr, 
			     zd_common, nlyr_common,
			     threedold, caoth3d->threed,
			     0, alloc);

  /** g1 **/
  status += redistribute_3D ((&(caoth3d->g1)), nx, ny, 
			     zd, nlyr, 
			     zd_common, nlyr_common,
			     threedold, caoth3d->threed,
			     0, alloc);

  /** g2 **/
  status += redistribute_3D ((&(caoth3d->g2)), nx, ny, 
			     zd, nlyr, 
			     zd_common, nlyr_common,
			     threedold, caoth3d->threed,
			     0, alloc);

  /** ff **/
  status += redistribute_3D ((&(caoth3d->ff)), nx, ny, 
			     zd, nlyr, 
			     zd_common, nlyr_common,
			     threedold, caoth3d->threed,
			     0, alloc);

  /** f **/
  status += redistribute_3D ((&(caoth3d->f)), nx, ny, 
			     zd, nlyr, 
			     zd_common, nlyr_common,
			     threedold, caoth3d->threed,
			     0, alloc);

  /** dscale **/
  status += redistribute_3D ((&(caoth3d->dscale)), nx, ny, 
			     zd, nlyr, 
			     zd_common, nlyr_common,
			     threedold, caoth3d->threed,
			     0, alloc);


  free(threedold);
  

  /* XXXXXXXXXX and update 3D z-profiles to common grid */
  caoth3d->nlyr = nlyr_common;
  free (caoth3d->zd);
  caoth3d->zd = calloc (nlyr_common+1, sizeof (float));
  for (lc=0; lc<=caoth3d->nlyr; lc++) 
    caoth3d->zd[lc] = zd_common[nlyr_common-lc];

  /* update nthreed - it might have changed */
  caoth3d->nthreed=0;
  for (lc=0; lc<caoth3d->nlyr; lc++)
    if (caoth3d->threed[lc])
      caoth3d->nthreed++;
  
  return status;
}





/***********************************/
/* Redistribute optical properties */
/***********************************/

int redistribute_caoth_microphys ( caoth_microphys_struct *microphys,
				   float                  *zd,
				   int                     nlyr, 
				   float                  *zd_common,
				   int                     nlyr_common )
{
  int status=0;
  
  /** lwc for layers  **/
  status += redistribute_1D((void *) (&(microphys->lwc_layer)),
                            zd, nlyr, 
                            zd_common, nlyr_common,
                            0, REDISTRIBUTE_FLOAT, 0);

  /** reff for layers **/
  status += redistribute_1D((void *) (&(microphys->effr_layer)),
                            zd, nlyr, 
                            zd_common, nlyr_common,
                            0, REDISTRIBUTE_FLOAT, 0);
  

  /** lwc at levels  **/
  status += redistribute_1D((void *) (&(microphys->lwc)),
                            zd, nlyr, 
                            zd_common, nlyr_common,
                            0, REDISTRIBUTE_FLOAT, 0);
  
  /** reff at levels **/
  status += redistribute_1D ((void *) (&(microphys->effr)),
			     zd, nlyr, 
			     zd_common, nlyr_common,
                             0, REDISTRIBUTE_FLOAT, 0);

  /* ????? the regridded lwc and reff are not necessarily correct at the edges ????? */
  /* ????? and should therefore not be used after regridding; a special        ????? */
  /* ????? treatment is required to "conserve" the hard boundaries;            ????? */
  /* ????? use lwc_layer and effr_layer instead!                               ????? */ 

  return status;
}

/**************************************************/
/* Redistribute molecular cross sections          */
/**************************************************/


int redistribute_molecular_crs (crs_out_struct *crs, int nlambda, int nlambda_lower, int nlambda_upper,
				int nmol, int number_of_ramanshifts,
				float *zd, int nlyr, 
				float *zd_common, int nlyr_common)
{
  int lc=0, iv=0, id=0, scale=0, status=0;
  double ***tmpcrs=NULL;

  double *xin=NULL, *yin=NULL, *xout=NULL, *yout=NULL;

  xin  = (double *) calloc (nlyr_common+1, sizeof(double));
  yin  = (double *) calloc (nlyr_common+1, sizeof(double));

  xout = (double *) calloc (nlyr_common+1, sizeof(double));
  yout = (double *) calloc (nlyr_common+1, sizeof(double));

  /* Allocate memory for temporary ary */
  if ((status = ASCII_calloc_double_3D (&tmpcrs, nlyr+1, nmol, nlambda)) != 0) {
    fprintf (stderr, "redistribute_molecular_crs: Error allocating memory for tmpcrs\n");
    return status;
  }

  /* cp original stuff to tmp ary */
  for (lc=0; lc<nlyr; lc++) 
    for (id=0; id<nmol; id++)
      for (iv=nlambda_lower; iv<=nlambda_upper; iv++)
	tmpcrs[lc][id][iv] = crs->crs[0][0][lc][id][iv];
  
  /* free original ary and allocate with changed resolution */
  ASCII_free_float_5D (crs->crs, 1, 1, nlyr+1, nmol);
  if ((status = ASCII_calloc_float_5D(&crs->crs, 1, 1, 
				      nlyr_common+1, nmol, nlambda)) != 0) {
    fprintf (stderr, "redistribute_molecular_crs: Error allocating memory for crs->crs\n");
    return status;
  }

  scale = 0;
  for (iv=nlambda_lower; iv<=nlambda_upper; iv++) {
    for (id=0; id<nmol; id++) {
      xin[0] = zd[0];
      yin[0] = tmpcrs[0][id][iv];

      for (lc=0; lc<nlyr; lc++) {
	xin[lc+1] = zd[lc+1];
	yin[lc+1] = tmpcrs[lc][id][iv];
      }
      
      for (lc=0; lc<=nlyr_common; lc++) 
      	xout[lc] = zd_common[lc];

      status = regrid (xin, yin, nlyr+1,
		       xout, yout, nlyr_common+1, scale);
      if (status != 0) {
	fprintf(stderr, "redistribute_molecular_crs: error calling regrid for crs\n");
	return status;
      }
      for (lc=1; lc<=nlyr_common; lc++) 
	crs->crs[0][0][lc-1][id][iv] = yout[lc];      
    }
  }
  ASCII_free_double_3D (tmpcrs, nlyr+1, nmol);

  /* free memory */
  free (xin); free (xout);
  free (yin); free (yout);


  /* Raman scattering cross section do not need to be redistributed as they are 
     calculated just before the call to qdisort in solve_rte.c. They do however 
     need to be large enough and memomry is thus re-allocated */
  ASCII_free_double (crs->crs_raman_RL, nlyr+1);
  if ((status = ASCII_calloc_double(&crs->crs_raman_RL, nlyr_common+1, number_of_ramanshifts+1)) != 0) {
    fprintf (stderr, "redistribute_molecular_crs: Error allocating memory for crs->crs_raman_RL\n");
    return status;
  }

  ASCII_free_double (crs->crs_raman_RG, nlyr+1);
  if ((status = ASCII_calloc_double(&crs->crs_raman_RG, nlyr_common+1, number_of_ramanshifts+1)) != 0) {
    fprintf (stderr, "redistribute_molecular_crs: Error allocating memory for crs->crs_raman_RG\n");
    return status;
  }


  return status;
}

/**************************************************/
/* Redistribute molecular scattering / absorption */
/**************************************************/

int redistribute_molecular (atm_optprop_struct *optprop, int nlambda, 
			    int nlambda_lower, int nlambda_upper,
			    int *nq, 
                            int nipa, /* ??? CE nipa is not used ! */
			    float *zd, int nlyr, int Nx, int Ny,  
			    float *zd_common, int nlyr_common)
{
  int ix=0, iy=0, lc=0, iv=0, iq=0, scale=0, status=0;
  double *****tmpopt=NULL, ***tmpopt_md=NULL;
  
  double *xin=NULL, *yin=NULL, *xout=NULL, *yout=NULL;

  xin  = (double *) calloc (nlyr_common+1, sizeof(double));
  yin  = (double *) calloc (nlyr_common+1, sizeof(double));

  xout = (double *) calloc (nlyr_common+1, sizeof(double));
  yout = (double *) calloc (nlyr_common+1, sizeof(double));

  /**************** Redistribute molecular optical depth (tau_molabs_r) *********************/

  /* Allocate memory for temporary ary */
  if ((status = ASCII_calloc_double_5D_arylen_restricted (&tmpopt, Nx, Ny, nlyr, 
							  nlambda, nlambda_lower, nlambda_upper, nq)) != 0) {
    fprintf (stderr, "Error allocating memory for tmpopt\n");
    return status;
  }
  
  for (ix=0; ix<Nx; ix++)
    for (iy=0; iy<Ny; iy++)
      /* cp original stuff to tmp ary */
      for (lc=0; lc<nlyr; lc++) 
	for (iv=nlambda_lower; iv<=nlambda_upper; iv++)
	  for (iq=0; iq<nq[iv]; iq++) 
	    tmpopt[ix][iy][lc][iv][iq] = optprop->tau_molabs_r[ix][iy][lc][iv][iq];
  
  
  
  /* free original ary and allocate with changed resolution */
  ASCII_free_double_5D (optprop->tau_molabs_r, Nx, Ny, nlyr, nlambda); 
  if ((status = ASCII_calloc_double_5D_arylen_restricted
       (&optprop->tau_molabs_r, Nx, Ny, nlyr_common,
	nlambda, nlambda_lower, nlambda_upper, nq)) != 0) {
    fprintf (stderr, "Error allocating memory for optprop->tau_molabs_r\n");
    return status;
  }
  
  for (ix=0; ix<Nx; ix++){
    for (iy=0; iy<Ny; iy++){
      scale = 1;
      for (iv=nlambda_lower; iv<=nlambda_upper; iv++) {
	for (iq=0; iq<nq[iv]; iq++) {
	  xin[0] = zd[0];
	  yin[0] = tmpopt[ix][iy][0][iv][iq];
	  for (lc=0; lc<nlyr; lc++) {
	    xin[lc+1] = zd[lc+1];
	    yin[lc+1] = tmpopt[ix][iy][lc][iv][iq];
	  }

	  for (lc=0; lc<=nlyr_common; lc++) 
	    xout[lc] = zd_common[lc];
	  
	  status = regrid (xin, yin, nlyr+1,
			   xout, yout, nlyr_common+1, scale);
	  if (status != 0) {
	    fprintf(stderr, "Redistribute: error calling regrid for tau_molabs_r\n");
	    return status;
	  }
	  for (lc=1; lc<=nlyr_common; lc++) 
	    optprop->tau_molabs_r[ix][iy][lc-1][iv][iq] = yout[lc];      
	}
      }
    }
  }
  
  ASCII_free_double_5D (tmpopt, Nx, Ny, nlyr, nlambda);
  
  /************************ Redistribute molecular optical properties -md *********************/


  /* 3DAbs, this is not used for 3D species, so loop over ix, iy is not included in md part */
  
  /* Allocate memory for temporary ary */
  if ((status = ASCII_calloc_double_3D_arylen_restricted (&tmpopt_md, nlyr, 
							  nlambda, nlambda_lower, nlambda_upper, nq)) != 0) {
    fprintf (stderr, "Error allocating memory for tmpopt\n");
    return status;
  }

  /* cp original stuff to tmp ary */
  for (lc=0; lc<nlyr; lc++) 
    for (iv=nlambda_lower; iv<=nlambda_upper; iv++)
      for (iq=0; iq<nq[iv]; iq++) 
	tmpopt_md[lc][iv][iq] = optprop->tau_molabs_md_r[0][0][lc][iv][iq];
  
  /* free original ary and allocate with changed resolution */
  ASCII_free_double_5D (optprop->tau_molabs_md_r, 1, 1, nlyr, nlambda);
  if ((status = ASCII_calloc_double_5D_arylen_restricted 
       (&optprop->tau_molabs_md_r, 1, 1, nlyr_common,
	nlambda, nlambda_lower, nlambda_upper, nq)) != 0) {
    fprintf (stderr, "Error allocating memory for optprop->tau_molabs_md_r\n");
    return status;
  }

  
  scale = 1;
  for (iv=nlambda_lower; iv<=nlambda_upper; iv++) {
    for (iq=0; iq<nq[iv]; iq++) {
      xin[0] = zd[0];
      yin[0] = tmpopt_md[0][iv][iq];
      for (lc=0; lc<nlyr; lc++) {
	xin[lc+1]  = zd[lc+1];
	yin[lc+1]  = tmpopt_md[lc][iv][iq];
      }
      for (lc=0; lc<=nlyr_common; lc++) 
      	xout[lc] = zd_common[lc];

	status = regrid (xin, yin, nlyr+1,
			 xout, yout, nlyr_common+1, scale);
	if (status != 0) {
	  fprintf(stderr, "Redistribute: error calling regrid for tau_molabs_md\n");
	  return status;
	}
	for (lc=1; lc<=nlyr_common; lc++) 
	  optprop->tau_molabs_md_r[0][0][lc-1][iv][iq] = yout[lc];      
    }
  }
  ASCII_free_double_3D (tmpopt_md, nlyr, nlambda);

  /********************************** Redistribute Rayleigh scattering *********************/

  
  /* Allocate memory for temporary ary */
  if ((status = ASCII_calloc_double_5D_arylen_restricted (&tmpopt, Nx, Ny, nlyr, 
							  nlambda, nlambda_lower, nlambda_upper, nq)) != 0) {
    fprintf (stderr, "Error allocating memory for tmpopt\n");
    return status;
  }

 
  /* cp original stuff to tmp ary */
  for (ix=0; ix<Nx; ix++)
    for (iy=0; iy<Ny; iy++)
      for (lc=0; lc<nlyr; lc++) 
	for (iv=nlambda_lower; iv<=nlambda_upper; iv++)
	  for (iq=0; iq<nq[iv]; iq++) 
	    tmpopt[ix][iy][lc][iv][iq] = optprop->tau_rayleigh_r[ix][iy][lc][iv][iq];
  
  /* free original ary and allocate with changed resolution */
  ASCII_free_double_5D (optprop->tau_rayleigh_r, Nx, Ny, nlyr, nlambda);
  if ((status = ASCII_calloc_double_5D_arylen_restricted 
       (&optprop->tau_rayleigh_r, Nx, Ny, nlyr_common,
	nlambda, nlambda_lower, nlambda_upper, nq)) != 0) {
    fprintf (stderr, "Error allocating memory for tmpopt\n");
    return status;
  }
  
  scale = 1;
  for (ix=0; ix<Nx; ix++){
    for (iy=0; iy<Ny; iy++){
      for (iv=nlambda_lower; iv<=nlambda_upper; iv++) {
	for (iq=0; iq<nq[iv]; iq++) {
	  xin[0] = zd[0];
	  yin[0] = tmpopt[ix][iy][0][iv][iq];
	  for (lc=0; lc<nlyr; lc++) {
	    xin[lc+1]  = zd[lc+1];
	    yin[lc+1]  = tmpopt[ix][iy][lc][iv][iq];
	  }
	  for (lc=0; lc<=nlyr_common; lc++) 
	    xout[lc] = zd_common[lc];
	  
	  status = regrid (xin, yin, nlyr+1,
			   xout, yout, nlyr_common+1, scale);
	  if (status != 0) {
	    fprintf(stderr, "Redistribute: error calling regrid for tau_rayleigh_r\n");
	    return status;
	  }
	  for (lc=1; lc<=nlyr_common; lc++) 
	    optprop->tau_rayleigh_r[ix][iy][lc-1][iv][iq] = yout[lc];      
	}
      }
    }
  }
  ASCII_free_double_5D (tmpopt, Nx, Ny, nlyr, nlambda);
  
  /* free memory */
  free (xin); free (xout);
  free (yin); free (yout);

  return 0;
}



/******************************/
/* Redistribute airmass stuff */
/******************************/

int redistribute_amf (atm_microphys_struct *microphys, 
		      crs_out_struct *crs,
		      int nlambda, 
		      int nlambda_lower, int nlambda_upper,
		      float *zd, int nlyr, 
		      float *zd_common, int nlyr_common)
{
  int is=0, iv=0, lc=0, scale=0, status=0;
  

  float  *xfin=NULL, *yfin=NULL, *xfout=NULL, *yfout=NULL;
  double *xin=NULL, *yin=NULL, *xout=NULL, *yout=NULL;
  double **tmp2=NULL;


  xin  = (double *) calloc (nlyr_common+1, sizeof(double));
  yin  = (double *) calloc (nlyr_common+1, sizeof(double));
  xout = (double *) calloc (nlyr_common+1, sizeof(double));
  yout = (double *) calloc (nlyr_common+1, sizeof(double));

  xfin  = (float *) calloc (nlyr_common+1, sizeof(float));
  yfin  = (float *) calloc (nlyr_common+1, sizeof(float));
  xfout = (float *) calloc (nlyr_common+1, sizeof(float));
  yfout = (float *) calloc (nlyr_common+1, sizeof(float));



  if (microphys->denstab_id > 0) {

    /** crs_amf **/

    /* Allocate memory for temporary ary */
    if ((status = ASCII_calloc_double (&tmp2, nlambda, nlyr+1)) != 0) {
      fprintf (stderr, "Error allocating memory for tmp2\n");
      return status;
    }

    /* cp original stuff to tmp ary */
    for (lc=0; lc<nlyr+1; lc++) 
      for (iv=nlambda_lower; iv<=nlambda_upper; iv++) 
	tmp2[iv][lc] = crs->crs_amf[iv][lc];

    /* allocate needed memory */
    ASCII_free_float (crs->crs_amf, nlambda);
    if ((status = ASCII_calloc_float (&crs->crs_amf, nlambda,
				      nlyr_common+1)) != 0)
      return status;

    scale = 0;
    for (iv=nlambda_lower; iv<=nlambda_upper; iv++) {
      for (lc=0; lc<nlyr+1; lc++) {
	xin[lc]  = zd[lc];
	yin[lc]  = tmp2[iv][lc];
      }
      for (lc=0; lc<=nlyr_common; lc++) 
	xout[lc] = zd_common[lc];
      
      status = regrid (xin, yin, nlyr+1,
		       xout, yout, nlyr_common+1, scale);
      if (status != 0) {
	fprintf(stderr, "Redistribute: error calling regrid for crs.crs_amf\n");
	return status;
      }
      for (lc=0; lc<=nlyr_common; lc++)
	crs->crs_amf[iv][lc] = yout[lc];      
    }
    ASCII_free_double (tmp2, nlambda);

    /** nsza_denstab **/

    /* Allocate memory for temporary ary */
    if ((status = ASCII_calloc_double (&tmp2, 
      				       microphys->nsza_denstab, 
      				       nlyr+1)) != 0)
      return status;

    /* cp original stuff to tmp ary */
    for (lc=0; lc<=nlyr; lc++) 
      for (is=0; is<microphys->nsza_denstab; is++) 
    	tmp2[is][lc] = microphys->denstab_amf[is][lc];

    /* free original ary and allocate with changed resolution */
    ASCII_free_float (microphys->denstab_amf, microphys->nsza_denstab);
    if ((status = ASCII_calloc_float (&microphys->denstab_amf, 
				      microphys->nsza_denstab, 
				      nlyr_common+1)) != 0)
      return status;
    
    for (lc=0; lc<=nlyr; lc++) 
      xfin[nlyr-lc] = zd[lc];
    for (lc=0; lc<=nlyr_common; lc++) 
      xfout[lc]  = zd_common[lc];

    for (is=0; is<microphys->nsza_denstab; is++) {
      for (lc=0; lc<=nlyr; lc++) 
      	yfin[nlyr-lc] = (float) tmp2[is][lc];

      status = arb_wvn (nlyr+1, xfin, yfin,
			nlyr_common+1, xfout, yfout, 1, 0); 
      if (status!=0) {
	fprintf (stderr, "Redistribute: Error %d interpolating denstab_amf\n", status);
	return status;
      }
      for (lc=0; lc<=nlyr_common; lc++) 
	microphys->denstab_amf[is][lc] = yfout[lc];
    }
    
    ASCII_free_double (tmp2, microphys->nsza_denstab);

  }

  /* free memory */
  free (xin); free (xout);
  free (yin); free (yout);

  free (xfin); free (xfout);
  free (yfin); free (yfout);



  return 0;
}



/***********************************************************************************/
/* Function: calloc_float_3D_sparse                                       @30_30i@ */
/* Description: Allocate memory for a three-dimensional array of float, but only   */
/*              if three[iz]==1. Required for allocation of 3D caoth fields.       */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                   @i30_30@ */
/***********************************************************************************/

static int calloc_float_3D_sparse  (float ****value,
				    int nz,
				    int nx,
				    int ny,
				    int *allocz)
{
  int i=0, j=0;

  if ( (*value = (float ***) calloc (nz, sizeof (float **))) == NULL )  
    return -1;

    
  for (i=0; i<nz; i++) {
    if (allocz[i]) {
      if ( ((*value)[i] = (float **) calloc (nx, sizeof (float *))) == NULL )
	return -1;
      for (j=0; j<nx; j++)
	if (((*value)[i][j] = (float *) calloc (ny, sizeof (float))) == NULL) 
	  return -1;
    }
  }

  return 0;
}    


/***********************************************************************************/
/* Function: free_float_3D_sparse                                         @30_30i@ */
/* Description: Free memory, which has been allocated with calloc_float_3D_sparse. */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                   @i30_30@ */
/***********************************************************************************/

static int free_float_3D_sparse (float ***value, 
				 int nz,
				 int nx, int*allocz)
{
  int i=0, j=0;

  for (i=0; i<nz; i++) {
    if (allocz[i]) {
      for (j=0; j<nx; j++)
	free (value[i][j]);
      free (value[i]);
    }
  }

  free (value);
  return 0;
}    
