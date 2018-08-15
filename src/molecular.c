/*--------------------------------------------------------------------
 * $Id: molecular.c 3308 2017-09-18 13:01:41Z Claudia.Emde $
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
#include "numeric.h"
#include "molecular.h"
#include "rayleigh.h"

/************************************/
/* prototypes of internal functions */
/************************************/

static int read_crs_from_files (input_struct input, 
                                output_struct *output,
                                int first,
                                int **molabs_src);

static int crs_o3  (crs_out_struct *crs_out, int first, 
		    float *lambda, int n_lambda,
                    int *molabs_src,
		    float ***temper, int ntemp, int Nxatm, int Nyatm, 
		    char *path,
		    int type);

static int crs_co2 (crs_out_struct *crs_out, int first, 
		    float *lambda, int n_lambda, 
                    int *molabs_src,
		    float ***temper, int ntemp, int Nxatm, int Nyatm, char *path,
		    int type, int quiet);

static int crs_o2 (crs_out_struct *crs_out, int first, 
		   float *lambda, int n_lambda, 
                   int *molabs_src,
		   float ***temper, int ntemp, int Nxatm, int Nyatm, char *path,
		   int type, int quiet);

static int crs_no2 (crs_out_struct *crs_out, int first,
		    float *lambda, int n_lambda, 
                    int *molabs_src,
		    float ***temper, int ntemp, int Nxatm, int Nyatm, char *path,
		    int type);

static int crs_so2 (crs_out_struct *crs_out, int first,
		    float *lambda, int n_lambda, 
                    int *molabs_src,
		    float ***temper, int ntemp, int Nxatm, int Nyatm, char *path,
		    int type);

static int crs_bro (crs_out_struct *crs_out,
		    float *lambda, int n_lambda, 
                    int *molabs_src,
		    float ***temper, int ntemp, int Nxatm, int Nyatm, char *path,
		    int type);

static int crs_oclo (crs_out_struct *crs_out, int first, 
		    float *lambda, int n_lambda, 
                    int *molabs_src,
		    float ***temper, int ntemp, int Nxatm, int Nyatm, char *path,
		    int type);

static int crs_hcho (crs_out_struct *crs_out, int first, 
		    float *lambda, int n_lambda, 
                    int *molabs_src,
		    float ***temper, int ntemp, int Nxatm, int Nyatm, char *path,
		    int type);

static int crs_o4 (crs_out_struct *crs_out, int first, int o4abs,
		   float *lambda, int n_lambda, 
                   int *molabs_src,
		   float ***temper, int ntemp, int Nxatm, int Nyatm, char *path,
		   int type);

static int read_crs (crs_out_struct *crs_out,
		     float *lambda, int n_lambda, 
		     int *molabs_src,
		     float ***temper, int ntemp, int Nxatm, int Nyatm, 
		     char *path, int mol_id);

static double quadratic(double x, double x0, double x1, double x2, 
		       double y0, double y1, double y2);

static double levelsum (atm_out_struct *atm, int mol);

static int setup_crs_from_lookup (float *lambda, 
                                  int n_lambda,
                                  int *iwvl_in_reptran_file,
                                  atm_out_struct *atm,
                                  int quiet,
                                  char *filename, 
                                  int ***molabs_src);

static int calc_crs_from_lookup (float *lambda, 
                                 int n_lambda, 
                                 int *iwvl_in_reptran_file,
                                 int *molabs_src,
                                 atm_out_struct *atm, 
                                 int i_mol,
                                 int quiet,
                                 char *lookup_file,
                                 crs_out_struct *crs_out);


/**************************************************************/
/* Setup absorption by molecules                              */
/**************************************************************/
int setup_crs (input_struct input, output_struct *output)
{
  int status = 0, lc=0, iv=0;
  static int first = 1;
  static int first_3d=1; 
  
  int i_mol,i_wvl;
  int **molabs_src;
  int i_mol_selected;

  char filename[FILENAME_MAX]="";
  
  /* Allocate memory for crs; currently for Rayleigh scattering  */
  /* and the absorption by the MOL_* species defined in uvspec.h */

  
  if (first) {
    status  = ASCII_calloc_float_5D(&output->crs.crs, 
				    output->atm.Nxatm, output->atm.Nyatm, output->atm.nlev, 
				    MOL_NN, output->wl.nlambda_r);
    if (status!=0) { 
      fprintf (stderr, "Error %d allocating memory for output->crs.crs\n", status); 
      return status; 
    }
  }

  /* 3DAbs how to make first only at first call for 3D ??, when is this needed??*/
  if (output->molecular3d)
    first=first_3d; 
  
  if (output->molecular3d){

    if(!input.quiet)
      fprintf(stderr, "3DAbs 3D allocate crs memory Nx %d Ny %d nlev %d MOL_NN %d nlam %d Nzatm %d!\n",
	      output->atm.Nxatm, output->atm.Nyatm, output->atm.nlev, 
	      MOL_NN, output->wl.nlambda_r, output->atm.Nzatm); 
    
    ASCII_free_float_5D(output->crs.crs, 1, 1, output->atm.Nzatm, 
			MOL_NN);

    
    ASCII_calloc_float_5D(&output->crs.crs, 
			  output->atm.Nxatm, output->atm.Nyatm, output->atm.nlev, 
			  MOL_NN, output->wl.nlambda_r);
  }

  
  if (first) {
    status  = ASCII_calloc_float(&output->crs.crs_amf, output->wl.nlambda_r, output->atm.nlev);
    if (status!=0) { 
      fprintf (stderr, "Error %d allocating memory for output->crs.crs_amf\n", status); 
      return status; 
    } 
  }

  switch(input.ck_scheme) {
  case CK_CRS:
  case CK_REPTRAN:
  case CK_REPTRAN_CHANNEL:
  case CK_RAMAN:

    status = ASCII_calloc_int (&molabs_src, MOL_NN, output->wl.nlambda_r);
    if (status!=0) { 
      fprintf (stderr, "Error %d allocating memory for molabs_src.\n", status); 
      return status; 
    }

    switch (output->atm.molabs) {
    case MOLABS_CALC:
      
      for (i_mol=1; i_mol<MOL_NN; i_mol++) 
        for (i_wvl=0; i_wvl<output->wl.nlambda_r; i_wvl++)
          molabs_src[i_mol][i_wvl]=MOLABS_SRC_CALC;
      
      /* Now read the cross sections from the crs_files */
      status = read_crs_from_files (input, output, first, molabs_src);
      if (status!=0) {
        fprintf (stderr, "Error %d reading cross sections from crs_files.\n", status);
	return status;
      }

      break;

    case MOLABS_LOOKUP:

      if (input.verbose)
	fprintf (stderr, " ... setting up cross sections from lookup tables\n");

      /* determine filename with the parameterization */  
      status = reptran_filename (input, -1, filename);

      status = setup_crs_from_lookup (output->wl.lambda_r, 
                                      output->wl.nlambda_r, 
                                      output->wl.iwvl_in_reptran_file_r,
			              &(output->atm),
                                      input.quiet,
                                      filename,
                                      &molabs_src); 
      if (status!=0) {
        fprintf (stderr, "Error %d setting up cross sections from lookup table.\n", status);
	return status;
      }
      
      /* calculate the cross sections from the lookup tables */
      for (i_mol=1; i_mol<MOL_NN; i_mol++) {
     	
	if (levelsum(&(output->atm),i_mol)>0) {
         
          /* determine whether any wavelength is selected for current species */
          i_mol_selected=0;
          for (i_wvl=0; i_wvl<output->wl.nlambda_r; i_wvl++)
            if (molabs_src[i_mol][i_wvl]==MOLABS_SRC_LOOKUP)
              i_mol_selected=1;
         
          if (i_mol_selected){
          
            /* determine filename with the lookup table */  
            status = reptran_filename (input, i_mol, filename);
	    
	    if (input.verbose)
	      fprintf (stderr, " ... calculating absorption cross sections from lookup table %s\n", filename);

            status = calc_crs_from_lookup(output->wl.lambda_r, 
                                          output->wl.nlambda_r, 
                                          output->wl.iwvl_in_reptran_file_r, 
                                          molabs_src[i_mol], 
                                          &(output->atm), 
                                          i_mol, 
                                          input.quiet,
                                          filename,
                                          &output->crs);
            if (status != 0) {
              fprintf (stderr, "Error while calculating absorption cross sections from lookup table %s.\n",filename);
              return -1;
            }
          }
          else
	    if (input.verbose)
	      fprintf (stderr, " ... no absorption cross sections from lookup tables required for %s.\n", gas_number2string(i_mol));
        }
      }

      /* read the cross sections from the crs_files */
      status = read_crs_from_files (input, output, first, molabs_src);
      if (status!=0) {
        fprintf (stderr, "Error %d reading cross sections from crs_files.\n", status);
	return status;
      }

      break;
      
    case MOLABS_FILE_MONO:
    case MOLABS_FILE_SPEC:
    case MOLABS_NONE:
      if (input.verbose)
	fprintf (stderr, " ... not reading cross sections\n");

      /* no need to read anything */
      break;

    default:
      fprintf (stderr, "Error, unknown molecular absorption option %d\n", output->atm.molabs);
      return -1;
    }      

    status = ASCII_free_int (molabs_src, MOL_NN);

    /* Set crs for amf calculations, 3DAbs not included*/
    if (output->atm.microphys.denstab_id > 0) {
      for (lc=0; lc<output->atm.nlev; lc++) 
        for (iv=0; iv<output->wl.nlambda_r; iv++)
          output->crs.crs_amf[iv][lc] = output->crs.crs[0][0][lc][output->atm.microphys.denstab_id][iv];
    }

    break;

  case CK_FU:
  case CK_FILE:          /* ??? need to check if that is appropriate here ??? */
  case CK_AVHRR_KRATZ:   /* ??? need to check if this is appropriate here ??? */
  case CK_LOWTRAN:       /* ??? need to check if this is appropriate here ??? */
  case CK_KATO:
  case CK_KATO2:
  case CK_KATO2_96:
  case CK_KATO2ANDWANDJI:

    if (input.atm.interpol_method_press != INTERP_METHOD_LOG) {
      fprintf (stderr, "pressure interpolation is assumed to be logarithmic\n");
      fprintf (stderr, "please change programm code for kato in crs_ck (in ck.c) \n");
      return -1;
    }

    /* setup correlated k-tables */
    status = crs_ck (&(output->ck),     
		     output->wl.lambda_r, output->wl.nlambda_r, 
		     output->atm.microphys.temper, output->atm.microphys.temper_avg, 
                     output->atm.microphys.press, output->atm.zd,
		     output->atm.microphys.dens, output->atm.microphys.dens_avg,
                     output->atm.nlev, output->atm.Nxatm, output->atm.Nyatm,
		     output->mixing_ratio, input.ck_h2ocont, 
		     input.ck_abs[CK_ABS_O4], input.ck_abs[CK_ABS_N2], 
		     input.ck_abs[CK_ABS_CO], input.ck_abs[CK_ABS_SO2], input.ck_abs[CK_ABS_NH3], input.ck_abs[CK_ABS_NO], input.ck_abs[CK_ABS_HNO3],
		     output->atm.sza_r,
		     input.filename[FN_PATH], input.ck_scheme_filename, first, input.verbose,
                     output->nipa,
		     &output->crs_ck,
                     input, output);   /* in ck.c */
    
    if (status!=0) {
      fprintf (stderr, "Error %d setting up cross sections for mol_abs_param.\n", status);
      return status;
    }

    break;

  default:
    fprintf (stderr, "Error: unsupported mol_abs_param scheme\n");
    return -1;
    break;
  }

  first=0;
  if(output->molecular3d)
    first_3d=0;
  
  return 0;
}

/**************************************************************/
/* Setup molecular (Rayleigh) scattering                      */
/**************************************************************/

int setup_rayleigh (input_struct input, output_struct *output)
{
  int status=0, iv=0;
  float *crs=NULL;

  /* allocate memory for depolarization and initialize */
  output->rayleigh_depol = calloc (output->wl.nlambda_r, sizeof(float));
  crs = calloc (output->wl.nlambda_r, sizeof(float));
  
  for (iv=0; iv<output->wl.nlambda_r; iv++)
    
    if ((input.crs_model[CRS_MOL_RAYLEIGH] == CRS_MODEL_NICOLET && input.rayleigh_depol == (float) -999.0)
	|| ( input.crs_model[CRS_MOL_RAYLEIGH] == CRS_MODEL_BODHAINE29 && input.rayleigh_depol == (float) -999.0) )
      output->rayleigh_depol[iv] = 0;
    else
      output->rayleigh_depol[iv] = input.rayleigh_depol;

  switch (input.crs_model[CRS_MOL_RAYLEIGH]) {
  case CRS_MODEL_BODHAINE:
    if (input.verbose)
      fprintf (stderr, "     using Bodhaine et al. Rayleigh scattering cross section\n");
    status = crs_rayleigh_bodhaine (output->wl.lambda_r, output->wl.nlambda_r, 
				    output->mixing_ratio[MX_CO2], 
                                    input.rayleigh_depol,
				    &(output->rayleigh_depol), &crs);
    if (status!=0) {
      fprintf (stderr, "Error %d setting up Rayleigh cross section from Bodhaine\n", status);
    }
    break;
  case CRS_MODEL_BODHAINE29:
    if (input.verbose)
      fprintf (stderr, "     using Bodhaine et al. Eq. 29 Rayleigh scattering cross section\n");
    status = crs_rayleigh_bodhaine29 (output->wl.lambda_r, output->wl.nlambda_r, 
				      &crs);
    if (status!=0) {
      fprintf (stderr, "Error %d setting up Rayleigh cross section from Bodhaine, Eq. 29.\n", status);
    }
    break;
    
  case CRS_MODEL_NICOLET:
    if (input.verbose)
      fprintf (stderr, "     using Nicolet Rayleigh scattering cross section\n");
    status = crs_rayleigh_nicolet (output->wl.lambda_r, output->wl.nlambda_r, &crs);
    if (status!=0) {
      fprintf (stderr, "Error %d setting up Rayleigh cross section from Nicolet\n", status);
    }
    break;

  case CRS_MODEL_PENNDORF:
    if (input.verbose)
      fprintf (stderr, "     using Penndorf Rayleigh scattering cross section\n");
    status = crs_rayleigh_penndorf (output->wl.lambda_r,  output->wl.nlambda_r, 
				    output->mixing_ratio[MX_CO2], 
                                    input.rayleigh_depol,
				    &(output->rayleigh_depol), &crs);

    if (status!=0) {
      fprintf (stderr, "Error %d setting up Rayleigh cross section from Penndorf\n", status);
    }
    break;

  default:
    fprintf (stderr, "Error, unknown Rayleigh scattering cross section %d\n", status);
  }
  
  for (iv=0; iv<output->wl.nlambda_r; iv++){
    output->crs.crs[0][0][0][MOL_AIR][iv]=crs[iv]; 
    
    if (input.crs_model[CRS_MOL_RAYLEIGH] == CRS_MODEL_NICOLET){
      if (input.rayleigh_depol == (float) -999.0)
        output->rayleigh_depol[iv] = 0;
      else
        output->rayleigh_depol[iv] = input.rayleigh_depol;
    }      
  }
  free(crs);
  
  return status;
}

/**************************************************************/
/* Setup inelastic molecular (Raman) scattering                      */
/**************************************************************/

int setup_raman (input_struct input, output_struct *output)
{
  /* Raman wavelength shifts and cross section are calculated just before the call */
  /* to qdisort, see solve_rte.c                                                   */

  int status=0;

  if ((output->crs.wvl_of_ramanshifts = calloc(output->crs.number_of_ramanwavelengths, sizeof(double)))==NULL) {
    status=-1;
    fprintf (stderr, "Error %d allocating memory for output->crs.wvl_of_ramanshifts\n", status); 
    fprintf (stderr,"       (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return status; 
  }  

  if ((status = ASCII_calloc_double (&output->crs.crs_raman_RL , output->atm.nlev, 
					output->crs.number_of_ramanshifts+1)) !=0) {
    fprintf (stderr, "Error %d allocating memory for output->crs.crs_raman \n", status); 
    fprintf (stderr,"       (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return status; 
  }

  if ((status = ASCII_calloc_double (&output->crs.crs_raman_RG , output->atm.nlev, 
					output->crs.number_of_ramanshifts+1)) !=0) {
    fprintf (stderr, "Error %d allocating memory for output->crs.crs_raman \n", status); 
    fprintf (stderr,"       (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__ );
    return status; 
  }
  
  return status;

}

/*****************************************************/
/* Read the cross sections from the crs_files.       */
/* Cross sections are only read if the concentration */
/* is larger than 0 in at least one of the levels    */
/*****************************************************/
int read_crs_from_files (input_struct input, 
                         output_struct *output, 
                         int first, 
                         int **molabs_src)
{

  //3DAbs loop over ix iy 
  int status=0, i=0;

  if (input.verbose)
    fprintf (stderr, " ... reading cross sections from crs_files\n");

  /* O3 absorption */
  if (levelsum (&(output->atm), MOL_O3) > 0) {
    status = crs_o3 (&output->crs, first, 
                     output->wl.lambda_r, output->wl.nlambda_r,
                     molabs_src[MOL_O3],
                     output->atm.microphys.temper, output->atm.nlev, 
		     output->atm.Nxatm, output->atm.Nyatm,
                     input.filename[FN_PATH], input.crs_model[CRS_MOL_O3]);
	
    if (status!=0) {
      fprintf (stderr, "Error %d setting up O3 cross section.\n", status);
      return status;
    }
  }
  else {
    if (input.verbose)
	  fprintf (stderr, " ... not reading O3 cross section because concentration is zero\n");
  }      

  /* CO2 absorption */
  if (levelsum (&(output->atm), MOL_CO2) > 0) {
    status = crs_co2 (&output->crs, first, 
                      output->wl.lambda_r, output->wl.nlambda_r, 
                      molabs_src[MOL_CO2],
                      output->atm.microphys.temper, output->atm.nlev,
		      output->atm.Nxatm, output->atm.Nyatm,
                      input.filename[FN_PATH], input.co2_spectrum, input.quiet);
      
    if (status!=0) {
      fprintf (stderr, "Error %d setting up CO2 cross section.\n", status);
      return status;
    }
  }
  else {
    if (input.verbose)
      fprintf (stderr, " ... not reading CO2 cross section because concentration is zero\n");
  }      

  /* O2 absorption */
  if (levelsum (&(output->atm), MOL_O2) > 0) {
    status = crs_o2 (&output->crs, first,
                     output->wl.lambda_r, output->wl.nlambda_r, 
                     molabs_src[MOL_O2],
                     output->atm.microphys.temper, output->atm.nlev,
		     output->atm.Nxatm, output->atm.Nyatm,
                     input.filename[FN_PATH], input.o2_spectrum, input.quiet);
	
    if (status!=0) {
      fprintf (stderr, "Error %d setting up O2 cross section.\n", status);
      return status;
    }
  }
  else {
    if (input.verbose)
      fprintf (stderr, " ... not reading O2 cross section because concentration is zero\n");
  }      

  /* NO2 absorption */
  if (levelsum (&(output->atm), MOL_NO2) > 0) {
    status = crs_no2 (&output->crs, first,
                      output->wl.lambda_r, output->wl.nlambda_r, 
                      molabs_src[MOL_NO2],
                      output->atm.microphys.temper, output->atm.nlev,
		      output->atm.Nxatm, output->atm.Nyatm,
                      input.filename[FN_PATH], input.crs_model[CRS_MOL_NO2]);

    if (status!=0) {
      fprintf (stderr, "Error %d setting up NO2 cross section.\n", status);
      return status;
    }
  }
  else {
    if (input.verbose)
      fprintf (stderr, " ... not reading NO2 cross section because concentration is zero\n");
  }      

  /* BRO absorption */
  if (levelsum (&(output->atm), MOL_BRO) > 0) {
    status = crs_bro (&output->crs,
                      output->wl.lambda_r, output->wl.nlambda_r, 
                      molabs_src[MOL_BRO],
                      output->atm.microphys.temper, output->atm.nlev,
		      output->atm.Nxatm, output->atm.Nyatm,
                      input.filename[FN_PATH], input.bro_spectrum);

    if (status!=0) {
      fprintf (stderr, "Error %d setting up BrO cross section.\n", status);
      return status;
    }
  }
  else {
    if (input.verbose)
      fprintf (stderr, " ... not reading BRO cross section because concentration is zero\n");
  }      

  /* OCLO absorption */
  if (levelsum (&(output->atm), MOL_OCLO) > 0) {
    status = crs_oclo (&output->crs, first, 
		       output->wl.lambda_r, output->wl.nlambda_r, 
		       molabs_src[MOL_OCLO],
		       output->atm.microphys.temper, output->atm.nlev,
		       output->atm.Nxatm, output->atm.Nyatm,
		       input.filename[FN_PATH], input.oclo_spectrum);

    if (status!=0) {
      fprintf (stderr, "Error %d setting up OClO cross section.\n", status);
      return status;
    }
  }
  else {
    if (input.verbose)
      fprintf (stderr, " ... not reading OCLO cross section because concentration is zero\n");
  }      
	
  /* HCHO absorption */
  if (levelsum (&(output->atm), MOL_HCHO) > 0) {
    status = crs_hcho (&output->crs, first, 
                       output->wl.lambda_r, output->wl.nlambda_r,
                       molabs_src[MOL_HCHO],
                       output->atm.microphys.temper, output->atm.nlev,
		       output->atm.Nxatm, output->atm.Nyatm,
                       input.filename[FN_PATH], input.hcho_spectrum);
	
    if (status!=0) {
      fprintf (stderr, "Error %d setting up HCHO cross section.\n", status);
      return status;
    }
  }
  else {
    if (input.verbose)
      fprintf (stderr, " ... not reading HCHO cross section because concentration is zero\n");
  }      

  /* O4 absorption */
  if (levelsum (&(output->atm), MOL_O4) > 0) {
    status = crs_o4 (&output->crs, first, input.ck_abs[CK_ABS_O4],
                     output->wl.lambda_r, output->wl.nlambda_r, 
                     molabs_src[MOL_O4],
                     output->atm.microphys.temper, output->atm.nlev,
		     output->atm.Nxatm, output->atm.Nyatm,
                     input.filename[FN_PATH], input.o4_spectrum);
	
    if (status!=0) {
      fprintf (stderr, "Error %d setting up O4 cross section.\n", status);
      return status;
    }
  }
  else {
    if (input.verbose)
      fprintf (stderr, " ... not reading O4 cross section because concentration is zero\n");
  }      

  /* SO2 absorption */
  if (levelsum (&(output->atm), MOL_SO2) > 0) {
    status = crs_so2 (&output->crs, first,
                      output->wl.lambda_r, output->wl.nlambda_r, 
                      molabs_src[MOL_SO2],
                      output->atm.microphys.temper, output->atm.nlev,
		      output->atm.Nxatm, output->atm.Nyatm,
                      input.filename[FN_PATH], input.so2_spectrum);
    if (status!=0) {
      fprintf (stderr, "Error %d setting up SO2 cross section.\n", status);
      return status;
    }
  }
  else {
    if (input.verbose)
      fprintf (stderr, " ... not reading SO2 cross section because concentration is zero\n");
  }      

  /* Next see if any cross section are given with crs_file, if so, overwrite */
      
  for (i=0; i<MOL_NN; i++) {
    if (strlen(output->crs.filename[i]) > 0) {
      status = read_crs (&output->crs, 
                         output->wl.lambda_r, output->wl.nlambda_r,
                         molabs_src[i],
                         output->atm.microphys.temper, output->atm.nlev,
			 output->atm.Nxatm, output->atm.Nyatm,
                         input.filename[FN_PATH], i);
      if (status!=0) {
        fprintf (stderr, "Error %d reading cross section file %s.\n", 
                 status, output->crs.filename[i]);
        return status;
      }
    }
  }

  return 0;

}

/***********************************************************************************/
/* Function: crs_o3                                                                */
/* Description:                                                                    */
/*  Calculates ozone cross sections as a function of wavelength and temperature.   */ 
/*                                                                                 */
/* Parameters:                                                                     */
/*  crs_out_struct *crs_out: structure holding the cross sections                  */
/*  int first: If 1, temperature coefficients for ozone will be calculated at      */
/*             wavelengths lambda; if 0, the temperature coefficients from an      */
/*             earlier call will be used. first should obviously be 1 the first    */
/*             time crs_o3 is called. It should also be 1 if lambda is changed.    */
/*             It does not need to be 1 if temper is changed.                      */
/*  float *lambda: Array of wavelengths [nm].                                      */
/*  int n_lambda:  Number of wavelengths.                                          */
/*  int *molabs_src: If !=MOLABS_SRC_CALC for a wavelength, cross section for      */
/*                   this wavelength is not written to crs_out->crs                */
/*  float *temper: Temperatures [K] at the atmospheric levels.                     */
/*  int n_temp:    Number of temperatures.                                         */ 
/*  char *path:    Location of the libRadtran data files.                          */
/*  int type:      The ozone cross section data to be used,                        */
/*                 1 (Molina and Molina, 1986);                                    */
/*                 2 (Daumont et al., 1992; Malicet et al., 1995)                  */
/*                 3 (Paur and Bass, 1985)                                         */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                                 */
/***********************************************************************************/

static int crs_o3 (crs_out_struct *crs_out, int first, 
		   float *lambda, int n_lambda, 
                   int *molabs_src,
		   float ***temper, int ntemp, int Nxatm, int Nyatm,
		   char *path,
		   int type)
{
  int ia=0, ib=0, it=0, iv=0, ix=0, iy=0;
  int n_crs=0, status=0;
  static float fact=0;
  char filenm[FILENAME_MAX+200]="", fullnm[FILENAME_MAX+200]="";
  static float *lambda_raw=NULL, *coeff_crs[3]={NULL,NULL,NULL};
  float a=0.0, b=0.0, gradient=0.0;

  if (first) {
    if (n_lambda <= 0 && ntemp != 0)
      return 2;
    else { 
      if (n_lambda != 0 &&  ntemp <= 0)
	return 3;
    }

    switch (type) {
    case CRS_MODEL_MOLINA:
      strcpy(filenm,"crs/crs_o3_mol_cf.dat");
      fact     = 1.0E-20;
      break;
    case CRS_MODEL_DAUMONT:
      strcpy(filenm,"crs/crs_o3_dau_cf.dat");
      fact     = 1.0E-20;
      break;
    case CRS_MODEL_BASS_AND_PAUR:
      strcpy(filenm,"crs/crs_o3_pab_cf.dat");
      fact     = 1.0E-20;
      break;
    case CRS_MODEL_BOGUMIL:
      strcpy(filenm,"crs/crs_O3_UBremen_cf.dat");
      fact     = 1.0E-20;
      break;
    default:
      fprintf (stderr, "Error, unknown ozone absorption cross section %d\n", status);
      return status;
    }
    
    strcpy(fullnm,path);
    strcat(fullnm,filenm);

    status = read_4c_file_float (fullnm, 
				 &lambda_raw, 
				 &(coeff_crs[0]), &(coeff_crs[1]), &(coeff_crs[2]), 
				 &n_crs);

    if (status!=0) { 
      fprintf (stderr, "Error %d opening file %s\n", status, fullnm); 
      return status; 
    } 

  } /* if (first) */
  
  for (iv=0; iv<n_lambda; iv++) {
   
    if (molabs_src[iv]==MOLABS_SRC_CALC && lambda[iv] < lambda_raw[n_crs-1]){

      ia=closest_above(lambda[iv], lambda_raw, n_crs);
      ib=closest_below(lambda[iv], lambda_raw, n_crs);

      if (ib<0) {
        fprintf (stderr, "Error (crs_o3), wavelength %.1f nm out of range!\n", lambda[iv]);
        fprintf (stderr, "lambda: %f n_crs: %d ia: %d ib: %d\n",lambda[iv], n_crs, ia, ib);
        return 9;
      }
      
      for (ix=0; ix<Nxatm; ix++){
	for (iy=0; iy<Nyatm; iy++){
	  for (it=0; it<ntemp; it++) {
	    if (ia > (n_crs-1) || ib<0)  
	      crs_out->crs[ix][iy][it][MOL_O3][iv] = 0;
	    else {
	      a = coeff_crs[0][ia] + 
		coeff_crs[1][ia] * (temper[ix][iy][it]-273.15) + 
		coeff_crs[2][ia]*(temper[ix][iy][it]-273.15)*(temper[ix][iy][it]-273.15);
	      if (ia == ib) {
		crs_out->crs[ix][iy][it][MOL_O3][iv] = a*fact;
	      }
	      else {
		b=coeff_crs[0][ib] + 
		  coeff_crs[1][ib] * (temper[ix][iy][it]-273.15) + 
		  coeff_crs[2][ib]*(temper[ix][iy][it]-273.15)*(temper[ix][iy][it]-273.15);
		gradient=(a-b)/(lambda_raw[ia]-lambda_raw[ib]);
		crs_out->crs[ix][iy][it][MOL_O3][iv] = (gradient*(lambda[iv]-lambda_raw[ib])+b) * fact;
		/*printf("%d %d %f %e\n",it,iv,lambda[iv],crs_out->crs[it][MOL_O3][iv]);*/
	      }
	    }
	  }
	}
      }
    }
  }

  return 0;
}



/***********************************************************************************/
/* Function: crs_co2                                                               */
/* Description:                                                                    */
/*  Calculates CO2 cross sections as a function of wavelength and temperature.     */ 
/*                                                                                 */
/* Parameters:                                                                     */
/*  crs_out_struct *crs_out: structure holding the cross sections                  */
/*  int first: If 1, temperature coefficients for ozone will be calculated at      */
/*             wavelengths lambda; if 0, the temperature coefficients from an      */
/*             earlier call will be used. first should obviously be 1 the first    */
/*             time crs_co2 is called. It should also be 1 if lambda is changed.   */
/*             It does not need to be 1 if temper is changed.                      */
/*  float *lambda: Array of wavelengths [nm].                                      */
/*  int n_lambda:  Number of wavelengths.                                          */
/*  int *molabs_src: If !=MOLABS_SRC_CALC for a wavelength, cross section for      */
/*                   this wavelength is not written to crs_out->crs                */
/*  float *temper: Temperatures [K] at the atmospheric levels.                     */
/*  int n_temp:    Number of temperatures.                                         */ 
/*  char *path:    Location of the libRadtran data files.                          */
/*  int type:      The CO2 cross section data to be used,                          */
/*                 1 (Lewis and Carver, 1983);                                     */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                                 */
/***********************************************************************************/

static int crs_co2 (crs_out_struct *crs_out, int first, 
		    float *lambda, int n_lambda, 
                    int *molabs_src,
		    float ***temper, int ntemp, int Nxatm, int Nyatm, char *path,
		    int type, int quiet)
{
  int ia=0, ib=0, it=0, iv=0, ix=0, iy=0;
  int n_crs=0, status=0;
  static float fact=0;
  char filenm[FILENAME_MAX+200]="", fullnm[FILENAME_MAX+200]="";
  static float *lambda_raw=NULL, *coeff_crs=NULL;
  float a=0, b=0, gradient=0;

  if (first) {
    if (n_lambda <= 0 && ntemp != 0)
      return 2;
    else { 
      if (n_lambda != 0 &&  ntemp <= 0)
	return 3;
    }
    
    switch (type) {
    case CO2_YOSHINO:
      strcpy(filenm,"crs/crs_co2_yoshino96_295K.dat");
      fact     = 1.0E-21;
      break;
    default:
      fprintf (stderr, "Error, unknown CO2 absorption cross section %d\n", type);
      return status;
    }

    
    strcpy(fullnm,path);
    strcat(fullnm,filenm);
    
    /* CO2 absorption cross section for 295K from Yoshino et al. [1996]             */
    /* ??? currently temperature dependence not considered; should change that! ??? */
    if (lambda[0] <= 200.0 && lambda[n_lambda-1] >= 118.0) { 
      if (!quiet)
	fprintf (stderr, " ... reading CO2 absorption cross section for 295K from Yoshino et al. [1996]\n");
      
      status = read_2c_file_float (fullnm, 
				   &lambda_raw, 
				   &coeff_crs, 
				   &n_crs);
      
      if (status!=0) { 
	fprintf (stderr, "Error %d reading file %s\n", status, fullnm); 
	return status; 
      } 
    }
  } /* if (first) */
  
  /* CO2 absorption cross section for 295K from Yoshino et al. [1996]             */
  /* ??? currently temperature dependence not considered; should change that! ??? */
  if (lambda[0] <= 200.0 && lambda[n_lambda-1] >= 118.0) { 
    for (iv=0; iv<n_lambda; iv++) {
    
      if (molabs_src[iv]==MOLABS_SRC_CALC && lambda[iv] < lambda_raw[n_crs-1]){

        ia=closest_above(lambda[iv], lambda_raw, n_crs);
        ib=closest_below(lambda[iv], lambda_raw, n_crs);
      
        if (ib<0) {
          fprintf (stderr, "Error (crs_co2), wavelength %.1f nm out of range!\n", lambda[iv]);
          fprintf (stderr, "lambda: %f n_crs: %d ia: %d ib: %d\n",lambda[iv], n_crs, ia, ib);
          return 9;
        }
      
	for (ix=0; ix<Nxatm; ix++){
	  for (iy=0; iy<Nyatm; iy++){
	    for (it=0; it<ntemp; it++) {
	      a = coeff_crs[ia];
	      if (ia == ib) {
		crs_out->crs[ix][iy][it][MOL_CO2][iv] = a*fact;
	      }
	      else {
		b = coeff_crs[ib];
		gradient = (a-b) / (lambda_raw[ia]-lambda_raw[ib]);
		crs_out->crs[ix][iy][it][MOL_CO2][iv] = (gradient*(lambda[iv]-lambda_raw[ib])+b) * fact;
	      }
	    }
	  }
	}
      }
    }
  }

  return 0;
}



/***********************************************************************************/
/* Function: crs_o2                                                                */
/* Description:                                                                    */
/*  Calculates O2 cross sections as a function of wavelength and temperature.      */ 
/*                                                                                 */
/* Parameters:                                                                     */
/*  crs_out_struct *crs_out: structure holding the cross sections                  */
/*  int first: If 1, temperature coefficients for O2 will be calculated at         */
/*             wavelengths lambda; if 0, the temperature coefficients from an      */
/*             earlier call will be used. first should obviously be 1 the first    */
/*             time crs_o2 is called. It should also be 1 if lambda is changed.    */
/*             It does not need to be 1 if temper is changed.                      */
/*  float *lambda: Array of wavelengths [nm].                                      */
/*  int n_lambda:  Number of wavelengths.                                          */
/*  int *molabs_src: If !=MOLABS_SRC_CALC for a wavelength, cross section for      */
/*                   this wavelength is not written to crs_out->crs                */
/*  float *temper: Temperatures [K] at the atmospheric levels.                     */
/*  int n_temp:    Number of temperatures.                                         */ 
/*  char *path:    Location of the libRadtran data files.                          */
/*  int type:      The O2 cross section data to be used,                           */
/*                 1 (default);                                                    */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                                 */
/***********************************************************************************/

static int crs_o2 (crs_out_struct *crs_out, int first, 
		   float *lambda, int n_lambda, 
                   int *molabs_src,
		   float ***temper, int ntemp, int Nxatm, int Nyatm, char *path,
		   int type, int quiet)
{
  int ia=0, ib=0, it=0, iv=0, i=0, j=0, ix=0, iy=0;
  int n_crs_SR=0, n_crs_Yo88=0, n_crs_Yo2005=0, n_crs_Og75=0, n_crs_Le83=0, status=0, n_tmp=0;
  static float fact_SR=0;
  char filenm[3][FILENAME_MAX+200]={"","",""}, fullnm[FILENAME_MAX+200]="";
  static float *lambda_SR=NULL, *lambda_Yo88=NULL, *lambda_Yo2005=NULL, *lambda_Og75=NULL, *lambda_Le83=NULL;
  static float *crs_SR[3][3], *crs_Yo88=NULL, *crs_Yo2005=NULL, *crs_Og75=NULL, *crs_Le83_288=NULL, *crs_Le83_alpha=NULL;
  float *tmp=NULL;
  float *tmpw=NULL, *tmp0=NULL, *tmp1=NULL, *tmp2=NULL;
  float a=0.0, b=0.0, d=0.0, gradient=0.0;
  
  if (first) {
    if (n_lambda <= 0 && ntemp != 0)
      return 2;
    else { 
      if (n_lambda != 0 &&  ntemp <= 0)
	return 3;
    }

    /* read Schumann-Runge from Minschwaner et al [1992] */
    if (lambda[0] <= (float) (1e7/49000.5) && lambda[n_lambda-1] >= (float) (1e7/57000.0)) {

      if (!quiet) {
	fprintf (stderr, " ... reading temperature-dependent O2 Schumann-Runge cross sections\n");
	fprintf (stderr, " ... by Minschwaner et al. [1992]\n");
      }

      switch (type) {
      case O2_DEFAULT:
	strcpy(filenm[0],"crs/Minschwaner_O2_fitcoef_cold.txt");
	strcpy(filenm[1],"crs/Minschwaner_O2_fitcoef_mid.txt");
	strcpy(filenm[2],"crs/Minschwaner_O2_fitcoef_hot.txt");
	fact_SR = 1.0E-20;
	break;
      default:
	fprintf (stderr, "Error, unknown oxygen absorption cross section %d\n", status);
	return status;
      }
      
      for (i=0;i<3;i++) {
	strcpy(fullnm,path);
	strcat(fullnm,filenm[i]);
	status = read_4c_file_float (fullnm, 
				     &tmp, 
				     &(crs_SR[i][0]), &(crs_SR[i][1]), &(crs_SR[i][2]), 
				     &n_tmp);
	
	if (status!=0) { 
	  fprintf (stderr, "Error %d reading file %s\n", status, fullnm); 
	  return status; 
	} 

	/* sort by wavelength, assuming that Minschwaner was sorted by wavenumber */
	tmpw = calloc (n_tmp, sizeof(float));
	tmp0 = calloc (n_tmp, sizeof(float));
	tmp1 = calloc (n_tmp, sizeof(float));
	tmp2 = calloc (n_tmp, sizeof(float));
	
	for (j=0; j<n_tmp; j++) {
	  tmpw[j] = tmp[n_tmp-1-j];
	  tmp0[j] = crs_SR[i][0][n_tmp-1-j];
	  tmp1[j] = crs_SR[i][1][n_tmp-1-j];
	  tmp2[j] = crs_SR[i][2][n_tmp-1-j];
	}
	
	for (j=0; j<n_tmp; j++) {
	  tmp[n_tmp-1-j]       = tmpw[j];
	  crs_SR[i][0][j] = tmp0[j];
	  crs_SR[i][1][j] = tmp1[j];
	  crs_SR[i][2][j] = tmp2[j];
	}
	
	free(tmpw); free(tmp0); free(tmp1); free(tmp2); 
	

	if (i==0) {
	  n_crs_SR = n_tmp;
	  lambda_SR = calloc (n_crs_SR, sizeof(float));
	  for (iv=0; iv<n_crs_SR; iv++)
	    lambda_SR[iv] = 1.0E7/tmp[n_crs_SR-1-iv];
	  
	  /* check if the data is sorted by wavelength */
	  for (iv=1;iv<n_crs_SR;iv++)
	    if (lambda_SR[iv] <= lambda_SR[iv-1]) {
	      fprintf (stderr, "Error, %s not sorted by wavenumber\n", fullnm);
	      return -1;
	    }
	}
	else {
	  if (n_crs_SR!=n_tmp) {
	    fprintf (stderr, "Error, Minschwaner et al. [1992] wavelength grids\n");
	    fprintf (stderr, "not equal: different number of grid points.\n");
	    return -1;
	  }
	  for (iv=0;iv<n_crs_SR;iv++) {
	    if (fabs((lambda_SR[iv]-1.0E7/tmp[n_crs_SR-1-iv])/lambda_SR[iv]) > 1e-6)  {
	      fprintf (stderr, "Error, Minschwaner et al. [1992] wavelength grids\n");
	      fprintf (stderr, "not equal: %f vs. %f\n", 
		       lambda_SR[iv], (float) (1.0E7/tmp[n_crs_SR-1-iv]));
	      return -1;
	    }
	  }
	}
	free (tmp);
      }
    }

    /* Yoshino et al. [1988] */
    if (lambda[0] <= 240.0 && lambda[n_lambda-1] >= 205.0) { 

      if (!quiet)
	fprintf (stderr, " ... reading O2 Herzberg continuum, 205-240nm, from Yoshino et al. [1988]\n");

      switch (type) {
      case O2_DEFAULT:
	strcpy(filenm[0],"crs/crs_o2_yoshino1988.dat");
	break;
      default:
	fprintf (stderr, "Error, unknown oxygen absorption cross section %d\n", status);
	return status;
      }

      strcpy(fullnm,path);
      strcat(fullnm,filenm[0]);

      status = read_2c_file_float (fullnm, 
				   &lambda_Yo88, 
				   &crs_Yo88, &n_crs_Yo88);

      if (status!=0) { 
	fprintf (stderr, "Error %d reading file %s\n", status, fullnm); 
	return status; 
      } 
    }

    /* Yoshino et al. [2005], 160nm - 175.5 nm */
    if (lambda[0] <= (float) (1e7/57000.0) && lambda[n_lambda-1] >= 160.0) { 

      if (!quiet)
	fprintf (stderr, " ... reading O2 Schumann-Runge continuum, 160-175nm, from Yoshino et al. [2005]\n");

      switch (type) {
      case O2_DEFAULT:
	strcpy(filenm[0],"crs/crs_o2_yoshino2005.dat");
	break;
      default:
	fprintf (stderr, "Error, unknown oxygen absorption cross section %d\n", status);
	return status;
      }

      strcpy(fullnm,path);
      strcat(fullnm,filenm[0]);

      status = read_2c_file_float (fullnm, 
				   &lambda_Yo2005, 
				   &crs_Yo2005, &n_crs_Yo2005);

      if (status!=0) { 
	fprintf (stderr, "Error %d reading file %s\n", status, fullnm); 
	return status; 
      } 
    }      

    /* Ogawa and Ogawa [1975], 108nm - 160 nm */
    if (lambda[0] <= (160.0) && lambda[n_lambda-1] >= 108.0) { 

      if (!quiet)
	fprintf (stderr, " ... reading Schumann-Runge continuum, 108-160nm, from Ogawa and Ogawa [1975]\n");

      switch (type) {
      case O2_DEFAULT:
	strcpy(filenm[0],"crs/crs_o2_ogawa1975.dat");
	break;
      default:
	fprintf (stderr, "Error, unknown oxygen absorption cross section %d\n", status);
	return status;
      }

      strcpy(fullnm,path);
      strcat(fullnm,filenm[0]);

      status = read_2c_file_float (fullnm, 
				   &lambda_Og75, 
				   &crs_Og75, &n_crs_Og75);

      if (status!=0) { 
	fprintf (stderr, "Error %d reading file %s\n", status, fullnm); 
	return status; 
      } 
    }      

    /* Lewis [1983], 121.4nm - 121.9 nm (Lyman alpha region) */
    if (lambda[0] <= (121.9) && lambda[n_lambda-1] >= 121.4) { 

      if (!quiet)
	fprintf (stderr, " ... reading Schumann-Runge continuum, 121.4-121.9nm, from Lewis [1983]\n");

      switch (type) {
      case O2_DEFAULT:
	strcpy(filenm[0],"crs/crs_o2_lewis1983.dat");
	break;
      default:
	fprintf (stderr, "Error, unknown oxygen absorption cross section %d\n", status);
	return status;
      }

      strcpy(fullnm,path);
      strcat(fullnm,filenm[0]);

      status = read_3c_file_float (fullnm, 
				   &lambda_Le83, 
				   &crs_Le83_288, &crs_Le83_alpha, &n_crs_Le83);

      if (status!=0) { 
	fprintf (stderr, "Error %d reading file %s\n", status, fullnm); 
	return status; 
      } 
    }
  } /* if (first) */
  

  /* calculate absorption cross section */
  for (iv=0; iv<n_lambda; iv++) {
    
    if (molabs_src[iv]==MOLABS_SRC_CALC){

      /* Schumann-Runge range */
      if (lambda[iv] <= (float) (1e7/49000.5) && lambda[iv] >= (float) (1e7/57000.0)) {

        ia=closest_above(lambda[iv], lambda_SR, n_crs_SR);
        ib=closest_below(lambda[iv], lambda_SR, n_crs_SR);
      
	for (ix=0; ix<Nxatm; ix++){
	  for (iy=0; iy<Nyatm; iy++){
	    for (it=0; it<ntemp; it++) {
	      if (temper[ix][iy][it]>=130&&temper[ix][iy][it]<=190)
		i=0;
	      
	      if (temper[ix][iy][it]>190&&temper[ix][iy][it]<=280)
		i=1;
	      
	      if (temper[ix][iy][it]>280&&temper[ix][iy][it]<=500)
		i=2;
	      
	      if (temper[ix][iy][it]<130||temper[ix][iy][it]>500) {
		fprintf (stderr, "Error, temperature %f out of range of the\n", temper[ix][iy][it]);
		fprintf (stderr, "Minschaner et al. [1997] Schumann-Runge parameterization\n");
		return -1;
	      }
	      
	      /* d is the transformed temperature variable defined in eq. (6): */
	      /* d = [(T-100)/10]**2.                                          */
	      
	      d = ((temper[ix][iy][it]-100.0)/10.0)*((temper[ix][iy][it]-100.0)/10.0);
	
	      /* assume that cross section outside the wavelength region is zero */
	      if (ia > (n_crs_SR-1) || ib<0)
		crs_out->crs[ix][iy][it][MOL_O2][iv]=0;
	      else {
		a = crs_SR[i][0][ia] * d*d + crs_SR[i][1][ia] * d + crs_SR[i][2][ia];
		if (ia == ib) {
		  crs_out->crs[ix][iy][it][MOL_O2][iv] = a * fact_SR;
		}
		else {
		  b=crs_SR[i][0][ib] * d*d + crs_SR[i][1][ib] * d + crs_SR[i][2][ib];
		  gradient=(a-b)/(lambda_SR[ia]-lambda_SR[ib]);
		  crs_out->crs[ix][iy][it][MOL_O2][iv] = (gradient*(lambda[iv]-lambda_SR[ib])+b) * fact_SR;
		}
	      }
	    }
	  }
	}
      }
      /* add Herzberg continuum to Schumann-Runge range; eq (3) from Yoshino et al [1988] */
      if (lambda[iv] >= 194.0 && lambda[iv] < 205.0)
	for (ix=0; ix<Nxatm; ix++){
	  for (iy=0; iy<Nyatm; iy++){
	    for (it=0; it<ntemp; it++)
	      crs_out->crs[ix][iy][it][MOL_O2][iv] += 6.884*(204.87/lambda[iv]) * exp(-69.7374*pow(log(204.87/lambda[iv]),2))*1e-24;
	  }
	}
      
      /* Yoshino et al. [1988] */
      if (lambda[iv] >= 205 && lambda[iv] <= 240) { 

        ia = closest_above(lambda[iv], lambda_Yo88, n_crs_Yo88);
        ib = closest_below(lambda[iv], lambda_Yo88, n_crs_Yo88);
	
	for (ix=0; ix<Nxatm; ix++){
	  for (iy=0; iy<Nyatm; iy++){
	    for (it=0; it<ntemp; it++) {
	      if (ia > (n_crs_Yo88-1) || ib<0)
		crs_out->crs[ix][iy][it][MOL_O2][iv]=0;
	      else {
		if (ia == ib)
		  crs_out->crs[ix][iy][it][MOL_O2][iv] = crs_Yo88[ia];
		else {
		  gradient=(crs_Yo88[ia] - crs_Yo88[ib]) / (lambda_Yo88[ia] - lambda_Yo88[ib]);
		  crs_out->crs[ix][iy][it][MOL_O2][iv] = (gradient * (lambda[iv]-lambda_Yo88[ib]) + crs_Yo88[ib]);
		}
	      }
	    }
	  }
	}
      }

      /* Yoshino et al. [2005] */
      if (lambda[iv] <= (float) (1e7/57000.0) && lambda[iv] >= 160.0) { 

        ia = closest_above(lambda[iv], lambda_Yo2005, n_crs_Yo2005);
        ib = closest_below(lambda[iv], lambda_Yo2005, n_crs_Yo2005);

	for (ix=0; ix<Nxatm; ix++){
	  for (iy=0; iy<Nyatm; iy++){
	    
	    for (it=0; it<ntemp; it++) {
	      if (ia > (n_crs_Yo2005-1) || ib<0)
		crs_out->crs[ix][iy][it][MOL_O2][iv]=0;
	      else {
		if (ia == ib)
		  crs_out->crs[ix][iy][it][MOL_O2][iv] = crs_Yo2005[ia];
		else {
		  gradient=(crs_Yo2005[ia] - crs_Yo2005[ib]) / (lambda_Yo2005[ia] - lambda_Yo2005[ib]);
		  crs_out->crs[ix][iy][it][MOL_O2][iv] = (gradient * (lambda[iv]-lambda_Yo2005[ib]) + crs_Yo2005[ib]);
		}
	      }
	    }
	  }
	}
      }

      /* Ogawa and Ogawa [1975], 108nm - 160 nm */
      if (lambda[iv] <= 160.0 && lambda[iv] >= 108.0) { 

        ia = closest_above(lambda[iv], lambda_Og75, n_crs_Og75);
        ib = closest_below(lambda[iv], lambda_Og75, n_crs_Og75);
	
	for (ix=0; ix<Nxatm; ix++){
	  for (iy=0; iy<Nyatm; iy++){
	    
	    for (it=0; it<ntemp; it++) {
	      if (ia > (n_crs_Og75-1) || ib<0)
		crs_out->crs[ix][iy][it][MOL_O2][iv]=0;
	      else {
		if (ia == ib)
		  crs_out->crs[ix][iy][it][MOL_O2][iv] = crs_Og75[ia];
		else {
		  gradient=(crs_Og75[ia] - crs_Og75[ib]) / (lambda_Og75[ia] - lambda_Og75[ib]);
		  crs_out->crs[ix][iy][it][MOL_O2][iv] = (gradient * (lambda[iv]-lambda_Og75[ib]) + crs_Og75[ib]);
		}
	      }
	    }
	  }
	}
      }
      /* Lewis [1983], 121.4nm - 121.9 nm (Lyman alpha region) */
      if (lambda[iv] <= (121.9) && lambda[iv] >= 121.4) { 

        ia = closest_above(lambda[iv], lambda_Le83, n_crs_Le83);
        ib = closest_below(lambda[iv], lambda_Le83, n_crs_Le83);

        /* temperature dependence in the Lyman alpha region */
        for (ix=0; ix<Nxatm; ix++){
	  for (iy=0; iy<Nyatm; iy++){
	    for (it=0; it<ntemp; it++) {

	      a = crs_Le83_288[ia] + crs_Le83_alpha[ia] * (temper[ix][iy][it]-288.0);
	      if (ia > (n_crs_Le83-1) || ib<0)
		crs_out->crs[ix][iy][it][MOL_O2][iv]=0;
	      else {
		if (ia == ib)
		  crs_out->crs[ix][iy][it][MOL_O2][iv] = a;
		else {
		  b = crs_Le83_288[ib] + crs_Le83_alpha[ib] * (temper[ix][iy][it]-288.0);
		  gradient=(a-b)/(lambda_Le83[ia]-lambda_Le83[ib]);
		  crs_out->crs[ix][iy][it][MOL_O2][iv] = (gradient*(lambda[iv]-lambda_Le83[ib])+b);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  return 0;
}



/***********************************************************************************/
/* Function: crs_no2                                                               */
/* Description:                                                                    */
/*  Calculates NO2 cross sections as a function of wavelength and temperature.     */ 
/*                                                                                 */
/* Parameters:                                                                     */
/*  crs_out_struct *crs_out: structure holding the cross sections                  */
/*  float *lambda: Array of wavelengths [nm].                                      */
/*  int n_lambda:  Number of wavelengths.                                          */
/*  int *molabs_src: If !=MOLABS_SRC_CALC for a wavelength, cross section for      */
/*                   this wavelength is not written to crs_out->crs                */
/*  float *temper: Temperatures [K] at the atmospheric levels.                     */
/*  int n_temp:    Number of temperatures.                                         */ 
/*  char *path:    Location of the libRadtran data files.                          */
/*  int type:      The  cross section data to be used,                             */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                                 */
/***********************************************************************************/

static int crs_no2 (crs_out_struct *crs_out, int first,
		    float *lambda, int n_lambda, 
                    int *molabs_src,
		    float ***temper, int ntemp, int Nxatm, int Nyatm, char *path,
		    int type)
{
  int ia=0, ib=0, it=0, iv=0, ix=0, iy=0;
  int n_raw=0, n_crs=0, status=0;
  char filenm[FILENAME_MAX+200]="", fullnm[FILENAME_MAX+200]="";
  
  static double fact=0;
  double a=0.0, b=0.0, gradient=0.0;
  double crs=0;
  static double *lambda_raw=NULL, *coeff_crs[3]={NULL,NULL,NULL};

  double *coeffc_raw=NULL;
  double *a0=NULL, *a1=NULL, *a2=NULL, *a3=NULL;

  if (n_lambda <= 0 && ntemp != 0)
    return 2;
  else { 
    if (n_lambda != 0 &&  ntemp <= 0)
      return 3;
  }

  switch (type) {
  case CRS_MODEL_BOGUMIL:
    if (first) {
      strcpy(filenm,"crs/crs_NO2_UBremen_cf.dat");
      strcpy(fullnm,path);
      fact     = 1.0E-20;
      strcat(fullnm,filenm);
      status = read_4c_file (fullnm, 
			     &lambda_raw, 
			     &(coeff_crs[0]), &(coeff_crs[1]), &(coeff_crs[2]), 
			     &n_crs);
      
      if (status!=0) { 
	fprintf (stderr, "Error %d opening file %s\n", status, fullnm); 
	return status; 
      } 
    }
    for (iv=0; iv<n_lambda; iv++) {
      if (molabs_src[iv]==MOLABS_SRC_CALC){
        ia=closest_above_double(lambda[iv], lambda_raw, n_crs);
        ib=closest_below_double(lambda[iv], lambda_raw, n_crs);
        if (ia > (n_crs-1) || ib<0) {
	  fprintf (stderr, "Error (crs_no2), wavelength %.1f nm out of range!\n", lambda[iv]);
	  fprintf (stderr, "Cannot do a monochromatic calculation at this wavelength;\n");
	  fprintf (stderr, "for infrared, have a look at the 'mol_abs_param' option!\n");
	  fprintf (stderr, "lambda: %f n_crs: %d ia: %d ib: %d\n",lambda[iv], n_crs, ia, ib);
	  return 9;
        }
	for (ix=0; ix<Nxatm; ix++){
	  for (iy=0; iy<Nyatm; iy++){
	    for (it=0; it<ntemp; it++) {
	      if (ia > (n_crs-1) || ib<0)
		crs_out->crs[ix][iy][it][MOL_NO2][iv] = 0;
	      else {
		a = coeff_crs[0][ia] + 
		  coeff_crs[1][ia] * (temper[ix][iy][it]-273.15) + 
		  coeff_crs[2][ia]*(temper[ix][iy][it]-273.15)*(temper[ix][iy][it]-273.15);
		if (ia == ib) {
		  crs_out->crs[ix][iy][it][MOL_NO2][iv] = a*fact;
		}
		else {
		  b=coeff_crs[0][ib] + 
		    coeff_crs[1][ib] * (temper[ix][iy][it]-273.15) + 
		    coeff_crs[2][ib]*(temper[ix][iy][it]-273.15)*(temper[ix][iy][it]-273.15);
		  gradient=(a-b)/(lambda_raw[ia]-lambda_raw[ib]);
		  crs_out->crs[ix][iy][it][MOL_NO2][iv] = (gradient*(lambda[iv]-lambda_raw[ib])+b) * fact;
		}
	      }
	    }
	  }
	}
      }
    }
    break;
  case CRS_MODEL_BURROWS:
    strcpy(filenm,"crs/crs_no2_gom.dat");
    strcpy(fullnm,path);
    strcat(fullnm,filenm);
    /* read cross section data */
    status = read_2c_file (fullnm,
			   &lambda_raw, &coeffc_raw, &n_raw);    
    if (status!=0) { 
      fprintf (stderr, "Error %d reading %s\n", status, fullnm); 
      return status; 
    } 
    /* linear interpolation to user-defined wavelength grid */
    status = linear_coeffc (lambda_raw, coeffc_raw, n_raw, 
			    &a0, &a1, &a2, &a3);
    if (status!=0) { 
      fprintf (stderr, "Error %d calculating interpolation coefficients for NO2 cross section\n", 
	       status); 
      return status; 
    } 
    for (iv=0; iv<n_lambda; iv++) {
      if (molabs_src[iv]==MOLABS_SRC_CALC){
        if (lambda[iv] < lambda_raw[0] || lambda[iv] > lambda_raw[n_raw-1])
	  crs=0;
        else
	  status = calc_splined_value (lambda[iv], &crs, lambda_raw, n_raw,
				       a0, a1, a2, a3);
        if (status!=0) {
	  fprintf (stderr, "Error %d interpolating NO2 cross section to %f nm\n", 
		   status, lambda[iv]);
	  return status;
        }
        /* ??? No temperature dependence for GOME cross sections ; ??? */
        /* ??? use constant cross section for all levels  ??? */
	for (ix=0; ix<Nxatm; ix++){
	  for (iy=0; iy<Nyatm; iy++){
	    for (it=0; it<ntemp; it++) {
	      crs_out->crs[ix][iy][it][MOL_NO2][iv] = crs;
	    }
	  }
	}
      }
    }
    free(a0); free(a1); free(a2); free(a3);
    free(lambda_raw); free(coeffc_raw);
    break;
  default:
    fprintf (stderr, "Error, unknown NO2 absorption cross section %d\n", type);
    return status;
  }
    
  return 0;
}


/***********************************************************************************/
/* Function: crs_so2                                                               */
/* Description:                                                                    */
/*  Calculates SO2 cross sections as a function of wavelength and temperature.     */ 
/*                                                                                 */
/* Parameters:                                                                     */
/*  crs_out_struct *crs_out: structure holding the cross sections                  */
/*  float *lambda: Array of wavelengths [nm].                                      */
/*  int n_lambda:  Number of wavelengths.                                          */
/*  int *molabs_src: If !=MOLABS_SRC_CALC for a wavelength, cross section for      */
/*                   this wavelength is not written to crs_out->crs                */
/*  float *temper: Temperatures [K] at the atmospheric levels.                     */
/*  int n_temp:    Number of temperatures.                                         */ 
/*  char *path:    Location of the libRadtran data files.                          */
/*  int type:      The  cross section data to be used,                             */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                                 */
/***********************************************************************************/

static int crs_so2 (crs_out_struct *crs_out, int first,
		    float *lambda, int n_lambda, 
                    int *molabs_src,
		    float ***temper, int ntemp, int Nxatm, int Nyatm, char *path,
		    int type)
{
  int ia=0, ib=0, it=0, iv=0, ix=0, iy=0;
  int n_crs=0, status=0;
  char filenm[FILENAME_MAX+200]="", fullnm[FILENAME_MAX+200]="";
  
  static double fact=0;
  double a=0.0, b=0.0, gradient=0.0;
  static double *lambda_raw=NULL, *coeff_crs[3]={NULL,NULL,NULL};

  if (n_lambda <= 0 && ntemp != 0)
    return 2;
  else { 
    if (n_lambda != 0 &&  ntemp <= 0)
      return 3;
  }

  switch (type) {
  case SO2_BREMEN:
    if (first) {
      strcpy(filenm,"crs/crs_SO2_UBremen_cf.dat");
      strcpy(fullnm,path);
      fact     = 1.0E-20;
      strcat(fullnm,filenm);
      status = read_4c_file (fullnm, 
			     &lambda_raw, 
			     &(coeff_crs[0]), &(coeff_crs[1]), &(coeff_crs[2]), 
			     &n_crs);
      
      if (status!=0) { 
	fprintf (stderr, "Error %d opening file %s\n", status, fullnm); 
	return status; 
      } 
    }
    for (iv=0; iv<n_lambda; iv++) {
      if (molabs_src[iv]==MOLABS_SRC_CALC){
        ia=closest_above_double(lambda[iv], lambda_raw, n_crs);
        ib=closest_below_double(lambda[iv], lambda_raw, n_crs);
        if (ia > (n_crs-1) || ib<0) {
	  fprintf (stderr, "Error (crs_so2), wavelength %.1f nm out of range!\n", lambda[iv]);
	  fprintf (stderr, "Cannot do a monochromatic calculation at this wavelength;\n");
	  fprintf (stderr, "for infrared, have a look at the 'mol_abs_param' option!\n");
	  fprintf (stderr, "lambda: %f n_crs: %d ia: %d ib: %d\n",lambda[iv], n_crs, ia, ib);
	  return 9;
        }
	for (ix=0; ix<Nxatm; ix++){
	  for (iy=0; iy<Nyatm; iy++){
	    for (it=0; it<ntemp; it++) {
	      if (ia > (n_crs-1) || ib<0)
		crs_out->crs[ix][iy][it][MOL_SO2][iv] = 0;
	      else {
		a = coeff_crs[0][ia] + 
		  coeff_crs[1][ia] * (temper[ix][iy][it]-273.15) + 
		  coeff_crs[2][ia]*(temper[ix][iy][it]-273.15)*(temper[ix][iy][it]-273.15);
		if (ia == ib) {
		  crs_out->crs[ix][iy][it][MOL_SO2][iv] = a*fact;
		}
		else {
		  b=coeff_crs[0][ib] + 
		    coeff_crs[1][ib] * (temper[ix][iy][it]-273.15) + 
		    coeff_crs[2][ib]*(temper[ix][iy][it]-273.15)*(temper[ix][iy][it]-273.15);
		  gradient=(a-b)/(lambda_raw[ia]-lambda_raw[ib]);
		  crs_out->crs[ix][iy][it][MOL_SO2][iv] = (gradient*(lambda[iv]-lambda_raw[ib])+b) * fact;
		}
	      }
	    }
	  }
	}
      }
    }
    break;
  default:
    fprintf (stderr, "Error, unknown SO2 absorption cross section %d\n", type);
    return status;
  }
    
  return 0;
}

/***********************************************************************************/
/* Function: crs_bro                                                               */
/* Description:                                                                    */
/*  Calculates BrO cross sections as a function of wavelength and temperature.     */ 
/*  ??? Temperature dependence is not yet considered ???                           */
/*                                                                                 */
/* Parameters:                                                                     */
/*  crs_out_struct *crs_out: structure holding the cross sections                  */
/*  float *lambda: Array of wavelengths [nm].                                      */
/*  int n_lambda:  Number of wavelengths.                                          */
/*  int *molabs_src: If !=MOLABS_SRC_CALC for a wavelength, cross section for      */
/*                   this wavelength is not written to crs_out->crs                */
/*  float *temper: Temperatures [K] at the atmospheric levels.                     */
/*  int n_temp:    Number of temperatures.                                         */ 
/*  char *path:    Location of the libRadtran data files.                          */
/*  int type:      The  cross section data to be used,                             */
/*                 no choices yet                                                  */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                                 */
/***********************************************************************************/

static int crs_bro (crs_out_struct *crs_out,
		    float *lambda, int n_lambda, 
                    int *molabs_src,
		    float ***temper, int ntemp, int Nxatm, int Nyatm, char *path,
		    int type)
{
  int it=0, iv=0, ix=0, iy=0;
  int n_raw=0, status=0;
  char filenm[FILENAME_MAX+200]="", fullnm[FILENAME_MAX+200]="";
  
  double crs=0;

  double *lambda_raw=NULL, *coeffc_raw=NULL;
  double *a0=NULL, *a1=NULL, *a2=NULL, *a3=NULL;

  if (n_lambda <= 0 && ntemp != 0)
    return 2;
  else { 
    if (n_lambda != 0 &&  ntemp <= 0)
      return 3;
  }

  switch (type) {
  case BRO_IASB:
    strcpy(filenm,"crs/crs_bro_wahner_228.dat");
    break;
  default:
    fprintf (stderr, "Error, unknown BrO absorption cross section %d\n", type);
    return status;
  }
    
  strcpy(fullnm,path);
  strcat(fullnm,filenm);

  /* read cross section data */
  status = read_2c_file (fullnm,
			 &lambda_raw, &coeffc_raw, &n_raw);
  
  if (status!=0) { 
    fprintf (stderr, "Error %d reading %s\n", status, fullnm); 
    return status; 
  } 
  
  /* linear interpolation to user-defined wavelength grid */
  status = linear_coeffc (lambda_raw, coeffc_raw, n_raw, 
			  &a0, &a1, &a2, &a3);
  if (status!=0) { 
    fprintf (stderr, "Error %d calculating interpolation coefficients for BrO cross section\n", 
	     status); 
    return status; 
  } 

  for (iv=0; iv<n_lambda; iv++) {

    if (molabs_src[iv]==MOLABS_SRC_CALC){
    
      if (lambda[iv] < lambda_raw[0] || lambda[iv] > lambda_raw[n_raw-1])
        crs=0;
      else
        status = calc_splined_value (lambda[iv], &crs, lambda_raw, n_raw,
				     a0, a1, a2, a3);
    

      if (status!=0) {
        fprintf (stderr, "Error %d interpolating BrO cross section to %f nm\n", 
		 status, lambda[iv]);
        return status;
      }

      /* ??? temperature dependence not yet considered; ??? */
      /* ??? use constant cross section for all levels  ??? */
      for (ix=0; ix<Nxatm; ix++){
	for (iy=0; iy<Nyatm; iy++){
	  for (it=0; it<ntemp; it++) 
	    crs_out->crs[ix][iy][it][MOL_BRO][iv] = crs;
	}
      }
    }
  }

  free(a0); free(a1); free(a2); free(a3);
  free(lambda_raw); free(coeffc_raw);

  return 0;
}



/***********************************************************************************/
/* Function: crs_oclo                                                              */
/* Description:                                                                    */
/*  Calculates oclo  cross sections as a function of wavelength and temperature.   */ 
/*                                                                                 */
/* Parameters:                                                                     */
/*  crs_out_struct *crs_out: structure holding the cross sections                  */
/*  int first: If 1, temperature coefficients for ozone will be calculated at      */
/*             wavelengths lambda; if 0, the temperature coefficients from an      */
/*             earlier call will be used. first should obviously be 1 the first    */
/*             time crs_oclo is called. It should also be 1 if lambda is changed.  */
/*             It does not need to be 1 if temper is changed.                      */
/*  float *lambda: Array of wavelengths [nm].                                      */
/*  int n_lambda:  Number of wavelengths.                                          */
/*  int *molabs_src: If !=MOLABS_SRC_CALC for a wavelength, cross section for      */
/*                   this wavelength is not written to crs_out->crs                */
/*  float *temper: Temperatures [K] at the atmospheric levels.                     */
/*  int n_temp:    Number of temperatures.                                         */ 
/*  char *path:    Location of the libRadtran data files.                          */
/*  int type:      The cross section data to be used,                              */
/*                 1 (Molina and Molina, 1986);                                    */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                                 */
/***********************************************************************************/

static int crs_oclo (crs_out_struct *crs_out, int first, 
		     float *lambda, int n_lambda, 
		     int *molabs_src,
		     float ***temper, int ntemp, int Nxatm, int Nyatm, char *path,
		     int type)
{
  int ia=0, ib=0, it=0, iv=0, ix=0, iy=0;
  int n_crs=0, status=0;
  static float fact=1.e-20;
  char filenm[FILENAME_MAX+200]="", fullnm[FILENAME_MAX+200]="";
  static float *lambda_raw=NULL, *crs_at_T[3]={NULL,NULL,NULL};
  float a=0.0, b=0.0, gradient=0.0;
  double T0=204, T1=296, T2=378;
  double s0=0, s1=0, s2=0;

  if (first) {
    if (n_lambda <= 0 && ntemp != 0)
      return 2;
    else { 
      if (n_lambda != 0 &&  ntemp <= 0)
	return 3;
    }

    switch (type) {
    case OCLO_WAHNER:
      strcpy(filenm,"crs/crs_oclo.dat");
      break;
    default:
      fprintf (stderr, "Error, unknown oclo absorption cross section %d\n", status);
      return status;
    }
    
    strcpy(fullnm,path);
    strcat(fullnm,filenm);

    status = read_4c_file_float (fullnm, 
				 &lambda_raw, 
				 &(crs_at_T[0]), &(crs_at_T[1]), &(crs_at_T[2]), 
				 &n_crs);

    if (status!=0) { 
      fprintf (stderr, "Error %d opening file %s\n", status, fullnm); 
      return status; 
    } 

  } /* if (first) */
  
  for (iv=0; iv<n_lambda; iv++) {

    if (molabs_src[iv]==MOLABS_SRC_CALC){

      ia=closest_above(lambda[iv], lambda_raw, n_crs);
      ib=closest_below(lambda[iv], lambda_raw, n_crs);

      if (ia > (n_crs-1) || ib<0) {
        fprintf (stderr, "Error (crs_oclo), wavelength %.1f nm out of range!\n", lambda[iv]);
        fprintf (stderr, "Cannot do a monochromatic calculation at this wavelength;\n");
        fprintf (stderr, "for infrared, have a look at the 'mol_abs_param' option!\n");
        fprintf (stderr, "lambda: %f n_crs: %d ia: %d ib: %d\n",lambda[iv], n_crs, ia, ib);
        return 9;
      }

      for (ix=0; ix<Nxatm; ix++){
	for (iy=0; iy<Nyatm; iy++){
	  for (it=0; it<ntemp; it++) {
	    if (fabs(temper[ix][iy][it]-T0) < 0.001) {
	      a = crs_at_T[0][ia];
	      b = crs_at_T[0][ib];
	    }
	    else if (fabs(temper[ix][iy][it]-T1) < 0.001) {
	      a = crs_at_T[1][ia];
	      b = crs_at_T[1][ib];
	    }
	    else if (fabs(temper[ix][iy][it]-T2) < 0.001) {
	      a = crs_at_T[2][ia];
	      b = crs_at_T[2][ib];
	    }
	    else {
	      s0 = crs_at_T[0][ia]*1e+20;
	      s1 = crs_at_T[1][ia]*1e+20;
	      s2 = crs_at_T[2][ia]*1e+20;
	      a = quadratic(temper[ix][iy][it], T0, T1, T2, s0, s1, s2);
	      s0 = crs_at_T[0][ib]*1e+20;
	      s1 = crs_at_T[1][ib]*1e+20;
	      s2 = crs_at_T[2][ib]*1e+20;
	      b = quadratic(temper[ix][iy][it], T0, T1, T2, s0, s1, s2);
	      a *= fact;
	      b *= fact;
	    }
	    
	    if (ia == ib) {
	      crs_out->crs[ix][iy][it][MOL_OCLO][iv] = a;
	    }
	    else {
	      gradient=(a-b)/(lambda_raw[ia]-lambda_raw[ib]);
	      crs_out->crs[ix][iy][it][MOL_OCLO][iv] = (gradient*(lambda[iv]-lambda_raw[ib])+b);
	    }
	  }
	}
      }
    }
  }
  return 0;
}

/***********************************************************************************/
/* Function: crs_hcho                                                              */
/* Description:                                                                    */
/*  Calculates hcho  cross sections as a function of wavelength and temperature.   */ 
/*                                                                                 */
/* Parameters:                                                                     */
/*  crs_out_struct *crs_out: structure holding the cross sections                  */
/*  int first: If 1, temperature coefficients for ozone will be calculated at      */
/*             wavelengths lambda; if 0, the temperature coefficients from an      */
/*             earlier call will be used. first should obviously be 1 the first    */
/*             time crs_hcho is called. It should also be 1 if lambda is changed.  */
/*             It does not need to be 1 if temper is changed.                      */
/*  float *lambda: Array of wavelengths [nm].                                      */
/*  int n_lambda:  Number of wavelengths.                                          */
/*  int *molabs_src: If !=MOLABS_SRC_CALC for a wavelength, cross section for      */
/*                   this wavelength is not written to crs_out->crs                */
/*  float *temper: Temperatures [K] at the atmospheric levels.                     */
/*  int n_temp:    Number of temperatures.                                         */ 
/*  char *path:    Location of the libRadtran data files.                          */
/*  int type:      The cross section data to be used,                              */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                                 */
/***********************************************************************************/

static int crs_hcho (crs_out_struct *crs_out, int first, 
		     float *lambda, int n_lambda, 
                     int *molabs_src,
		     float ***temper, int ntemp, int Nxatm, int Nyatm, char *path,
		     int type)
{
  int it=0, iv=0, ix=0, iy=0;
  int n_raw=0, status=0;
  char filenm[FILENAME_MAX+200]="", fullnm[FILENAME_MAX+200]="";
  
  double crs=0;

  double *lambda_raw=NULL, *coeffc_raw=NULL;
  double *a0=NULL, *a1=NULL, *a2=NULL, *a3=NULL;

  if (n_lambda <= 0 && ntemp != 0)
    return 2;
  else { 
    if (n_lambda != 0 &&  ntemp <= 0)
      return 3;
  }

  switch (type) {
  case HCHO_IASB:
    strcpy(filenm,"crs/crs_hcho_cantrell.dat");
    break;
  default:
    fprintf (stderr, "Error, unknown HCHO absorption cross section %d\n", status);
    return status;
  }
    
  strcpy(fullnm,path);
  strcat(fullnm,filenm);
  
  /* read cross section data */
  status = read_2c_file (fullnm,
			 &lambda_raw, &coeffc_raw, &n_raw);
  
  if (status!=0) { 
    fprintf (stderr, "Error %d reading %s\n", status, fullnm); 
    return status; 
  } 
  
  /* linear interpolation to user-defined wavelength grid */
  status = linear_coeffc (lambda_raw, coeffc_raw, n_raw, 
			  &a0, &a1, &a2, &a3);
  if (status!=0) { 
    fprintf (stderr, "Error %d calculating interpolation coefficients for HCHO cross section\n", 
	     status); 
    return status; 
  } 
  
  for (iv=0; iv<n_lambda; iv++) {
    
    if (molabs_src[iv]==MOLABS_SRC_CALC){

      if (lambda[iv] < lambda_raw[0] || lambda[iv] > lambda_raw[n_raw-1])
        crs=0;
      else
        status = calc_splined_value (lambda[iv], &crs, lambda_raw, n_raw,
				     a0, a1, a2, a3);
    
    
      if (status!=0) {
        fprintf (stderr, "Error %d interpolating HCHO cross section to %f nm\n", 
	         status, lambda[iv]);
        return status;
      }
    
      /* ??? temperature dependence not yet considered; ??? */
      /* ??? use constant cross section for all levels  ??? */
      for (ix=0; ix<Nxatm; ix++){
	for (iy=0; iy<Nyatm; iy++){
	  for (it=0; it<ntemp; it++) 
	    crs_out->crs[ix][iy][it][MOL_HCHO][iv] = crs;
	}
      }
    }
  }
  
  return 0;
}

/***********************************************************************************/
/* Function: crs_o4                                                                */
/* Description:                                                                    */
/*  Calculates o4  cross sections as a function of wavelength.                     */ 
/*                                                                                 */
/* Parameters:                                                                     */
/*  crs_out_struct *crs_out: structure holding the cross sections                  */
/*  int first: If 1, temperature coefficients for ozone will be calculated at      */
/*             wavelengths lambda; if 0, the temperature coefficients from an      */
/*             earlier call will be used. first should obviously be 1 the first    */
/*             time crs_hcho is called. It should also be 1 if lambda is changed.  */
/*             It does not need to be 1 if temper is changed.                      */
/*  float *lambda: Array of wavelengths [nm].                                      */
/*  int n_lambda:  Number of wavelengths.                                          */
/*  int *molabs_src: If !=MOLABS_SRC_CALC for a wavelength, cross section for      */
/*                   this wavelength is not written to crs_out->crs                */
/*  float *temper: Temperatures [K] at the atmospheric levels.                     */
/*  int n_temp:    Number of temperatures.                                         */ 
/*  char *path:    Location of the libRadtran data files.                          */
/*  int type:      The cross section data to be used,                              */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                                 */
/***********************************************************************************/

static int crs_o4 (crs_out_struct *crs_out, int first, int o4abs,
		   float *lambda, int n_lambda, 
                   int *molabs_src,
		   float ***temper, int ntemp, int Nxatm, int Nyatm, char *path,
		   int type)
{
  int it=0, iv=0, ix=0, iy=0;
  int n_raw=0, status=0;
  char filenm[FILENAME_MAX+200]="", fullnm[FILENAME_MAX+200]="";
  
  double crs=0;

  double *lambda_raw=NULL, *coeffc_raw=NULL;
  double *a0=NULL, *a1=NULL, *a2=NULL, *a3=NULL;

  if (n_lambda <= 0 && ntemp != 0)
    return 2;
  else { 
    if (n_lambda != 0 &&  ntemp <= 0)
      return 3;
  }

  switch (type) {
  case O4_GREENBLATT:
    strcpy(filenm,"crs/crs_o4_greenblatt.dat");
    break;
  default:
    fprintf (stderr, "Error, unknown O4 absorption cross section %d\n", status);
    return status;
  }
    
  strcpy(fullnm,path);
  strcat(fullnm,filenm);
  
  /* read cross section data */
  status = read_2c_file (fullnm,
			 &lambda_raw, &coeffc_raw, &n_raw);
  
  if (status!=0) { 
    fprintf (stderr, "Error %d reading %s\n", status, fullnm); 
    return status; 
  } 
  
  /* linear interpolation to user-defined wavelength grid */
  status = linear_coeffc (lambda_raw, coeffc_raw, n_raw, 
			  &a0, &a1, &a2, &a3);
  if (status!=0) { 
    fprintf (stderr, "Error %d calculating interpolation coefficients for O4 cross section\n", 
	     status); 
    return status; 
  } 
  
  for (iv=0; iv<n_lambda; iv++) {
    
    if (molabs_src[iv]==MOLABS_SRC_CALC){

      if (lambda[iv] < lambda_raw[0] || lambda[iv] > lambda_raw[n_raw-1])
        crs=0;
      else
        status = calc_splined_value (lambda[iv], &crs, lambda_raw, n_raw,
				     a0, a1, a2, a3);
    
    
      if (status!=0) {
        fprintf (stderr, "Error %d interpolating O4 cross section to %f nm\n", 
		 status, lambda[iv]);
        return status;
      }
    
      /* ??? temperature dependence not yet considered; ??? */
      /* ??? use constant cross section for all levels  ??? */
      for (ix=0; ix<Nxatm; ix++){
	for (iy=0; iy<Nyatm; iy++){
	  for (it=0; it<ntemp; it++) 
	    if (o4abs)
	      crs_out->crs[ix][iy][it][MOL_O4][iv] = crs;
	    else
	      crs_out->crs[ix][iy][it][MOL_O4][iv] = 0.0;
	}
      }
    }
  }

  free(lambda_raw);
  free(coeffc_raw);
  free(a0);
  free(a1);
  free(a2);
  free(a3);

  return 0;
}

/***********************************************************************************/
/* Function: read_crs                                                              */
/* Description:                                                                    */
/*  Read a cross section versus wavelength file. No temperature dependence may     */
/*  be specified.                                                                  */
/*                                                                                 */
/* Parameters:                                                                     */
/*  crs_out_struct *crs_out: structure holding the cross sections                  */
/*  float *lambda: Array of wavelengths [nm].                                      */
/*  int n_lambda:  Number of wavelengths.                                          */
/*  int *molabs_src: If !=MOLABS_SRC_CALC for a wavelength, cross section for      */
/*                   this wavelength is not written to crs_out->crs                */
/*  float *temper: Temperatures [K] at the atmospheric levels.                     */
/*  int n_temp:    Number of temperatures.                                         */ 
/*  char *path:    Location of the libRadtran data files.                          */
/*  int mol_id:    Molecule identifier                                             */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                                 */
/***********************************************************************************/

static int read_crs (crs_out_struct *crs_out,
		     float *lambda, int n_lambda, 
                     int *molabs_src,
		     float ***temper, int ntemp, int Nxatm, int Nyatm,
		     char *path, int mol_id)
{
  int it=0, iv=0, ix=0, iy=0;
  int n_raw=0, status=0;
  char filenm[FILENAME_MAX+200]="", fullnm[FILENAME_MAX+200]="";
  
  double crs=0;

  double *lambda_raw=NULL, *coeffc_raw=NULL;
  double *a0=NULL, *a1=NULL, *a2=NULL, *a3=NULL;

  if (n_lambda <= 0 && ntemp != 0)
    return 2;
  else { 
    if (n_lambda != 0 &&  ntemp <= 0)
      return 3;
  }

  strcpy(filenm,crs_out->filename[mol_id]);
    
  strcat(fullnm,filenm);

  /* read cross section data */
  status = read_2c_file (fullnm, &lambda_raw, &coeffc_raw, &n_raw);
  
  if (status!=0) { 
    fprintf (stderr, "Error %d reading %s\n", status, fullnm); 
    return status; 
  } 
  
  /* linear interpolation to user-defined wavelength grid */
  status = linear_coeffc (lambda_raw, coeffc_raw, n_raw, 
			  &a0, &a1, &a2, &a3);
  if (status!=0) { 
    fprintf (stderr, "Error %d calculating interpolation coefficients for cross section file %s\n", 
	     status, crs_out->filename[mol_id]); 
    return status; 
  } 

  for (iv=0; iv<n_lambda; iv++) {

    if (molabs_src[iv]==MOLABS_SRC_CALC){
    
      if (lambda[iv] < lambda_raw[0] || lambda[iv] > lambda_raw[n_raw-1])
        crs=0;
      else
        status = calc_splined_value (lambda[iv], &crs, lambda_raw, n_raw,
				   a0, a1, a2, a3);
    
      if (status!=0) {
        fprintf (stderr, "Error %d interpolating cross section file %s to %f nm\n", 
	       status, crs_out->filename[mol_id], lambda[iv]);
        return status;
      }

      /* ??? temperature dependence not yet considered; ??? */
      /* ??? use constant cross section for all levels  ??? */
      for (ix=0; ix<Nxatm; ix++){
	for (iy=0; iy<Nyatm; iy++){
	  for (it=0; it<ntemp; it++) 
	    crs_out->crs[ix][iy][it][mol_id][iv] = crs;
	}
      }
    }
  }

  free(lambda_raw);
  free(coeffc_raw);
  free(a0);
  free(a1);
  free(a2);
  free(a3);

  return 0;
}

int closest_above (float lambda, float *lambda_raw, int n_crs) 
{
  int iv=0;
  
  for (iv=n_crs-1;iv>=0;iv--)  {
    if ((lambda_raw[iv]-lambda) == 0.0) 
      return iv;
    else {
      if ((lambda_raw[iv]-lambda) < 0.0) 
	return iv+1;
    }
  }
  return 0;
}

int closest_above_double (double lambda, double *lambda_raw, int n_crs) 
{
  int iv=0;
  
  for (iv=n_crs-1;iv>=0;iv--)  {
    if ((lambda_raw[iv]-lambda) == 0.0) 
      return iv;
    else {
      if ((lambda_raw[iv]-lambda) < 0.0) 
	return iv+1;
    }
  }
  return 0;
}

int closest_below (float lambda, float *lambda_raw, int n_crs) 
{
  int iv=0;

  for (iv=0; iv<n_crs; iv++) {
    if ((lambda_raw[iv]-lambda) == 0.0) 
      return iv;
    else  {
      if ((lambda_raw[iv]-lambda) > 0.0) 
	return iv-1;
    }
  }
  return 0;
}

int closest_below_double (double lambda, double *lambda_raw, int n_crs) 
{
  int iv=0;

  for (iv=0; iv<n_crs; iv++) {
    if ((lambda_raw[iv]-lambda) == 0.0) 
      return iv;
    else  {
      if ((lambda_raw[iv]-lambda) > 0.0) 
	return iv-1;
    }
  }
  return 0;
}

static double quadratic(double x, double x0, double x1, double x2, double y0, double y1, double y2)
{
  double y=0.0;
  double a=0, b=0, c=0;
  
  c = (1/(x2-x1))*(((y2-y0)/(x2-x0))-((y1-y0)/(x1-x0)));
  b = (y1-y0)/(x1-x0)-c*(x1+x0);
  a = y0 - b*x0 -c*x0*x0;

  y = a + b*x + c*x*x;

  return y;
}



static double levelsum (atm_out_struct *atm, int mol)
{
  int lc=0, ix=0, iy=0;
  double sum=0;
 
  for (ix=0; ix<atm->Nxatm; ix++){
    for (iy=0; iy<atm->Nyatm; iy++){
      for (lc=0;lc<atm->nlev;lc++){
	sum+=atm->microphys.dens[mol][ix][iy][lc];
      }
    }
  }
  return sum;
}


/******************************************************************************************/
/* This function checks whether all required species are considered in the representative */
/*   wavelength file and sets up the variable molabs_src that determine whether and from  */
/*   which source the spectral absorption cross sections for the species are obtained.    */
/******************************************************************************************/
int setup_crs_from_lookup (float *lambda, 
                           int n_lambda, 
                           int *iwvl_in_reptran_file,
                           atm_out_struct *atm,
                           int quiet,
                           char *filename,
                           int ***molabs_src)
{
#if HAVE_LIBNETCDF
 
  int i,j,i_mol;

  int status, ncid;
  int varid;

  size_t start[5] = {0,0,0,0,0};
  size_t count[5] = {0,0,0,0,0};

  size_t nspecies_in_reptran_file;
  size_t max_len_species_name;
  size_t nwvl;
  double *wvl_in_reptran_file;
  char **species_name_in_reptran_file;
  int **cross_section_source;

  int ispecies_in_reptran_file[MOL_NN];
  
  char *gas = NULL;


  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status) {
    fprintf (stderr, "Netcdf error %d opening %s\n", status, filename);
    return -1;
  }
  else {
    
    /* read the data from the netcdf file */
    status =0;

    status += nc_inq_dimid (ncid, "nspecies", &varid);
    status += nc_inq_dimlen (ncid, varid, &nspecies_in_reptran_file);

    status += nc_inq_dimid (ncid, "max_len_species_name", &varid);
    status += nc_inq_dimlen (ncid, varid, &max_len_species_name);

    status += nc_inq_dimid (ncid, "nwvl", &varid);
    status += nc_inq_dimlen (ncid, varid, &nwvl);

    status += nc_inq_varid (ncid, "species_name", &varid);
    status += ASCII_calloc_char(&species_name_in_reptran_file,nspecies_in_reptran_file,max_len_species_name);
    for (i=0; i<nspecies_in_reptran_file; i++) {
      start[0]=i; start[1]=0;
      count[0]=1; count[1]=max_len_species_name;
      status += nc_get_vara_text (ncid, varid, start, count, species_name_in_reptran_file[i]);
    }

    status += nc_inq_varid (ncid, "cross_section_source", &varid);
    status += ASCII_calloc_int(&cross_section_source,nwvl,nspecies_in_reptran_file);
    for (i=0; i<nwvl; i++) {
      start[0]=i; start[1]=0;
      count[0]=1; count[1]=nspecies_in_reptran_file;
      status += nc_get_vara_int (ncid, varid, start, count, cross_section_source[i]);
    }
    
    status += nc_inq_varid (ncid, "wvl", &varid);
    wvl_in_reptran_file = (double *) calloc (nwvl, sizeof(double));
    status += nc_get_var_double (ncid, varid, wvl_in_reptran_file);

    status += nc_close(ncid);
  
    if (status) {
      fprintf (stderr, "Error %d while reading %s.\n", status, filename);
      return -1;
    }

  }

  /* find out which species are considered in the reptran_file */
  for (i_mol=1; i_mol<MOL_NN; i_mol++) {

    gas=gas_number2string(i_mol);

    for (j=0; j<n_lambda; j++)  /* by default absorption from crs_files */
      (*molabs_src)[i_mol][j]=MOLABS_SRC_CALC;

    ispecies_in_reptran_file[i_mol]=-1;

    for (j=0; j<nspecies_in_reptran_file; j++)
      if (strncmp(gas,species_name_in_reptran_file[j],4) == 0) 
        ispecies_in_reptran_file[i_mol]=j;

    if (ispecies_in_reptran_file[i_mol]>-1) 
      for (j=0; j<n_lambda; j++) 
        if (cross_section_source[iwvl_in_reptran_file[j]][ispecies_in_reptran_file[i_mol]]==1) 
          (*molabs_src)[i_mol][j]=MOLABS_SRC_LOOKUP;

    free(gas);

  }
 
  status = ASCII_free_int(cross_section_source,nwvl);
  status = ASCII_free_char(species_name_in_reptran_file,nspecies_in_reptran_file);
  free(wvl_in_reptran_file);

  return 0;
#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif

}


/****************************************************************************/
/* Interpolate the spectral absorption cross sections of a given species    */
/* from lookup table file. Variable 'molabs_src' determines the wavelengths */
/* to be considered in this function.                                       */
/****************************************************************************/
//3DAbs!!
int calc_crs_from_lookup (float *lambda, 
                          int n_lambda,
                          int *iwvl_in_reptran_file,
                          int *molabs_src,
		          atm_out_struct *atm, 
                          int i_mol,
                          int quiet,
                          char *lookup_file,
                          crs_out_struct *crs_out)
{
  
#if HAVE_LIBNETCDF

  double *wvl;
  int *wvl_index;
  double *pressure;
  double *vmrs;
  double *t_ref;
  double *t_pert;
  float ****xsec;

  double **xsec_pre_interpolated;
  double *xsec_pre_p;
  double *xsec_pre_vmrs;
  double *xsec_pre_t;

  size_t start[5] = {0,0,0,0,0};
  size_t count[5] = {0,0,0,0,0};

  size_t nwvl, n_vmrs, n_pressure, n_t_pert;

  double press_tmp;

  int *iwvl_in_lookup_file;
  int *iwvl_active; // internal counter for active wavelengths
  int n_active_wvl;

  int status, ncid, varid;
  int nwvl_dimid, n_vmrs_dimid, n_pressure_dimid, n_t_pert_dimid;

  int i_lev=0, ix=0, iy=0;

  double dt,vmr=0;

  int i, i_wvl, i_t_pert, i_lower_p, i_lower_vmr, i_lower_t, i_p, i_t, i_vmr;


  /* initialize arrays required for linking the wavelengths */
  iwvl_in_lookup_file = (int *) calloc (n_lambda, sizeof(int));
  iwvl_active = (int *) calloc (n_lambda, sizeof(int));
  for (i_wvl=0; i_wvl<n_lambda; i_wvl++) {
    iwvl_in_lookup_file[i_wvl]=-1;
    iwvl_active[i_wvl]=-1;
  }

  /* open the netcdf file */
  status = nc_open (lookup_file, NC_NOWRITE, &ncid);
  if (status!=0) {
    fprintf (stderr, "Netcdf error %d opening %s\n", status, lookup_file);
    return -1;
  }
  else{

    /* read dimensions */
    status = 0;

    status += nc_inq_dimid (ncid, "nwvl", &nwvl_dimid);
    status += nc_inq_dimid (ncid, "n_vmrs", &n_vmrs_dimid);
    status += nc_inq_dimid (ncid, "n_pressure", &n_pressure_dimid);
    status += nc_inq_dimid (ncid, "n_t_pert", &n_t_pert_dimid);
    
    status += nc_inq_dimlen (ncid, nwvl_dimid, &nwvl);
    status += nc_inq_dimlen (ncid, n_vmrs_dimid, &n_vmrs);
    status += nc_inq_dimlen (ncid, n_pressure_dimid, &n_pressure);
    status += nc_inq_dimlen (ncid, n_t_pert_dimid, &n_t_pert);
  
    if (status!=0) {
      fprintf (stderr, "Error %d while reading dimensions in %s.\n", status, lookup_file);
      return -1;
    }

    /* allocate arrays and read data from file */
    status = 0;

    status += nc_inq_varid (ncid, "wvl", &varid);
    wvl = (double *) calloc (nwvl, sizeof(double));
    status += nc_get_var_double (ncid, varid, wvl);

    status += nc_inq_varid (ncid, "wvl_index", &varid); /* wvl_index links the wavelengths in the lookup table file to the wavelengths in the representative wavelengths file */
    wvl_index = (int *) calloc (nwvl, sizeof(int));
    status += nc_get_var_int (ncid, varid, wvl_index);

    /* find the wavelengths required from the lookup table */
    n_active_wvl=0;
    for (i_wvl=0; i_wvl<n_lambda; i_wvl++)
      for (i=0; i<nwvl; i++)
        if (molabs_src[i_wvl]==MOLABS_SRC_LOOKUP && iwvl_in_reptran_file[i_wvl]==(wvl_index[i]-1)) { 
          iwvl_in_lookup_file[i_wvl]=i;
          iwvl_active[i_wvl]=n_active_wvl;
          n_active_wvl++;
        }
  
    if (n_active_wvl==0) /* this should not happen because calc_crs_from_lookup is called only if at least one wavelength was selected */
      err_out("Error in calc_crs_from_lookup: No active wavelength found.\n", -1);

    status += nc_inq_varid (ncid, "pressure", &varid);
    pressure = (double *) calloc (n_pressure, sizeof(double));
    status += nc_get_var_double (ncid, varid, pressure);

    status += nc_inq_varid (ncid, "vmrs", &varid);
    vmrs = (double *) calloc (n_vmrs, sizeof(double));
    status += nc_get_var_double (ncid, varid, vmrs);

    status += nc_inq_varid (ncid, "t_ref", &varid);
    t_ref = (double *) calloc (n_pressure, sizeof(double));
    status += nc_get_var_double (ncid, varid, t_ref);

    status += nc_inq_varid (ncid, "t_pert", &varid);
    t_pert = (double *) calloc (n_t_pert, sizeof(double));
    status += nc_get_var_double (ncid, varid, t_pert);

    status += nc_inq_varid (ncid, "xsec", &varid);
    status += ASCII_calloc_float_4D(&xsec,n_t_pert,n_vmrs,n_active_wvl,n_pressure);

    for (i_t_pert=0; i_t_pert<n_t_pert; i_t_pert++)
      for (i_wvl=0; i_wvl<n_lambda; i_wvl++)
        if (iwvl_active[i_wvl]>-1)                // read the cross sections only at the wavelengths where they are needed
          for (i_vmr=0; i_vmr<n_vmrs; i_vmr++) {
            start[0]=i_t_pert; start[1]=i_vmr; start[2]=iwvl_in_lookup_file[i_wvl]; start[3]=0;
            count[0]=1; count[1]=1; count[2]=1; count[3]=n_pressure;
            status += nc_get_vara_float (ncid, varid, start, count, xsec[i_t_pert][i_vmr][iwvl_active[i_wvl]]);
          }
    
    /* close netcdf file */
    nc_close(ncid);

    if (status!=0) {
      fprintf (stderr, "Error %d while reading data from %s\n", status, lookup_file);
      return -1;
    }


    /* Initialize some additional arrays */
    status = ASCII_calloc_double(&xsec_pre_interpolated,n_active_wvl,2);
    xsec_pre_t = (double *) calloc (2, sizeof(double));
    xsec_pre_vmrs = (double *) calloc (2, sizeof(double));
    xsec_pre_p = (double *) calloc (2, sizeof(double));

    
    /* calculate the absorption cross sections for each level */
    for (ix=0; ix<atm->Nxatm; ix++){
      for(iy=0; iy<atm->Nyatm; iy++){
	for (i_lev=0; i_lev<atm->nlev; i_lev++) {
	  
	  press_tmp=atm->microphys.press[ix][iy][i_lev]*100; /* factor 100 comes from conversion of hPa to Pa */
   
	  /* check if pressure is covered by lookup table */
	  if (press_tmp<pressure[n_pressure-1]){ 
	    if (!quiet) { 
	      fprintf(stderr,"     Warning: Pressure at lev %d is %e Pa, thus lower than %e Pa, which is the lower limit in the absorption lookup table.\n",i_lev, press_tmp, pressure[n_pressure-1]);
	      fprintf(stderr,"       Using %e Pa for calculating the absorption cross sections at this level.\n", pressure[n_pressure-1]);
	    }
	    press_tmp = pressure[n_pressure-1];
	  }
	  else if (press_tmp>pressure[0]){
	    if (!quiet) { 
	      fprintf(stderr,"     Warning: Pressure at lev %d is %e Pa, thus higher than %e Pa, which is the upper limit in the absorption lookup table.\n",i_lev,  press_tmp, pressure[0]);
	      fprintf(stderr,"       Using %e Pa for calculating the absorption cross sections at this level.\n", pressure[0]);
	    }
	    press_tmp = pressure[0];
	  }
	  
	  /* Find the pressure in the lookup table */
	  i_lower_p = locate(pressure,n_pressure,press_tmp);
	  
	  /* n_vmrs>1 means that absorption depends on volume mixing ratio */
	  if (n_vmrs>1) {
	    
	    if (i_mol!=MOL_H2O) {
	      fprintf(stderr,"Error: The absorption lookup table %s considers mixing ratio dependence for a gas species other than H2O.\n", lookup_file);
	      return -1;
	    }
	    
	    vmr=atm->microphys.dens[i_mol][ix][iy][i_lev]/atm->microphys.dens[MOL_AIR][ix][iy][i_lev];

	    /* check if vmr is covered by lookup table */
	    if (vmr<vmrs[0]){
	      if (!quiet) {
		fprintf(stderr,"     Warning: VMR at lev %d is lower than %e, which is the lower limit in the absorption lookup table.\n",i_lev, vmrs[0]);
		fprintf(stderr,"       Using %e for calculating the absorption cross sections at this level.\n", vmrs[0]);
	      }
	      vmr = vmrs[0];
	    }
	    else if (vmr>vmrs[n_vmrs-1]){
	      if (!quiet) {
		fprintf(stderr,"     Warning: VMR at lev %d is higher than %e, which is the upper limit in the absorption lookup table.\n",i_lev, vmrs[n_vmrs-1]);
		fprintf(stderr,"       Using %e for calculating the absorption cross sections at this level.\n", vmrs[n_vmrs-1]);
	      }
	      vmr = vmrs[n_vmrs-1];
	    }           
	    
	    /* find the gas mixing ratio in the lookup table */
	    i_lower_vmr = locate(vmrs,n_vmrs,vmr);
	    
	  }
	  else
	    i_lower_vmr=0;
	  
	  for (i_p=0; i_p<2; i_p++) {
	    dt=atm->microphys.temper[ix][iy][i_lev]-t_ref[i_lower_p+i_p]; // temperature pertubation
	    
	    /* check if temperature pertubation is covered by lookup table */
	    if (dt<t_pert[0]){
	      if (!quiet) {
		fprintf(stderr,"     Warning: Temperature pertubation at lev %d is %e K, thus lower than %e K, which is the lower limit in the absorption lookup table.\n",i_lev, dt, t_pert[0]);
		fprintf(stderr,"       Using %e K for calculating the absorption cross sections at this level.\n", t_pert[0]);
	      }
	      dt = t_pert[0];
	    }
	    else if (dt>t_pert[n_t_pert-1]){
	      if (!quiet) {
		fprintf(stderr,"     Warning: Temperature pertubation at lev %d is %e K, thus higher than %e K, which is the upper limit in the absorption lookup table.\n",i_lev, dt, t_pert[n_t_pert-1]);
		fprintf(stderr,"       Using %e K for calculating the absorption cross sections at this level.\n", t_pert[n_t_pert-1]);
	      }
	      dt = t_pert[n_t_pert-1];
	    }
	    
	    /* find the temperature perturbation in the lookup table */
	    i_lower_t = locate(t_pert,n_t_pert,dt);
	    
	    /* now do the interpolations */
	    for (i_wvl=0; i_wvl<n_lambda; i_wvl++) {
	      if (iwvl_active[i_wvl]>-1) {
		for (i_t=0; i_t<2; i_t++) {
		  
		  if (n_vmrs>1) {
		    for (i_vmr=0; i_vmr<2; i_vmr++){
		      xsec_pre_vmrs[i_vmr]=(double) xsec[i_lower_t+i_t][i_lower_vmr+i_vmr][iwvl_active[i_wvl]][i_lower_p+i_p];
		    }
		    xsec_pre_t[i_t] = xsec_pre_vmrs[0] + (xsec_pre_vmrs[1] - xsec_pre_vmrs[0]) * (vmr - vmrs[i_lower_vmr]) / (vmrs[i_lower_vmr+1] - vmrs[i_lower_vmr]);
		  }
		  else{
		    xsec_pre_t[i_t]=(double) xsec[i_lower_t+i_t][0][iwvl_active[i_wvl]][i_lower_p+i_p];
		  }
		  
		}
		
		xsec_pre_interpolated[iwvl_active[i_wvl]][i_p] = xsec_pre_t[0] + (xsec_pre_t[1] - xsec_pre_t[0]) * (dt - t_pert[i_lower_t]) / (t_pert[i_lower_t+1] - t_pert[i_lower_t]);
	      }
	    }
	  }
	  
	  for (i_wvl=0; i_wvl<n_lambda; i_wvl++) {
	    if (iwvl_active[i_wvl]>-1) {
	      
	      for (i_p=0; i_p<2; i_p++)
		xsec_pre_p[i_p]=xsec_pre_interpolated[iwvl_active[i_wvl]][i_p];
	      
	      //crs_out->crs[i_lev][i_mol][i_wvl] = xsec_pre_p[0] + (xsec_pre_p[1] - xsec_pre_p[0]) * (log(press_tmp) - log(pressure[i_lower_p])) / (log(pressure[i_lower_p+1]) - log(pressure[i_lower_p]));
	      crs_out->crs[ix][iy][i_lev][i_mol][i_wvl] = xsec_pre_p[0] + (xsec_pre_p[1] - xsec_pre_p[0]) * (press_tmp - pressure[i_lower_p]) / (pressure[i_lower_p+1] - pressure[i_lower_p]);
	      crs_out->crs[ix][iy][i_lev][i_mol][i_wvl] *= 1e-16; /* The lookup table file contains the cross sections in units of 10^(-20)m^2; here we need cm^2, thus we multiply with 10^(-16). */
	      
	      /* fprintf(stderr, " 3DAbs calc_crs_from_lookup %d %d %d %d %d  crs_out->crs %g \n",  ix,iy,i_lev,i_mol,i_wvl,  crs_out->crs[ix][iy][i_lev][i_mol][i_wvl]);  */
	      if (crs_out->crs[ix][iy][i_lev][i_mol][i_wvl]<0){  /* Negative absorption cross sections are set to zero */
		crs_out->crs[ix][iy][i_lev][i_mol][i_wvl]=0;
		if(!quiet) fprintf(stderr,"Warning in calc_crs_from_lookup: The absorption cross section of %s in grid cell %d,%d,%d at %e nm wavelength is negative; set to zero!\n", gas_number2string(i_mol), ix, iy, i_lev, lambda[i_wvl]);
	      }
	    }
	  }
	}
      }
    }
    
    /* free the arrays */
    free(wvl);
    free(wvl_index);
    free(pressure);
    free(vmrs);
    free(t_ref);
    free(t_pert);
    status = ASCII_free_float_4D(xsec,n_t_pert,n_vmrs,n_active_wvl);
    status = ASCII_free_double(xsec_pre_interpolated,n_active_wvl);
    free(xsec_pre_t);
    free(xsec_pre_vmrs);
    free(xsec_pre_p);
    free(iwvl_in_lookup_file);
    free(iwvl_active);
    
    return 0;

  }
#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif

}
