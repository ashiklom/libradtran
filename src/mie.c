/*--------------------------------------------------------------------
 * $Id: mie.c 3247 2016-08-31 23:38:54Z bernhard.mayer $
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

#include "ascii.h"
#include "mie.h"
#include "numeric.h"
#include "getopt.h"

#if HAVE_LIBNETCDF
#include <netcdf.h>
#endif


#define PI 3.1415926535897932
#define MAXLINE 255



int aerosol_input_data(/* Input */
                       char* data_files_path, int aerosol_type, 
                       int verbose, 
                       /* Output */
                       aerdist_struct *aerdistr);

int aerosol_refractive_index(/* Input */
                             double* lambda,
                             int nlambda, 
                             char* data_files_path, 
                             int aerosol_type,
                             int n_hum,
                             int verbose,
                             /* Output */
                             double*** ref_real,
                             double*** ref_imag);

int aerosol_size_distribution(/* Input */
                              double rmin,
                              double rmax,
                              double rmod,   
                              double sigma,
                              double lambda,
                              int verbose, 
                              /* Output */
                              int* ndens,
                              double** radius,
                              double** number_dens);


int calc_distribution_radius_grid(/* Input */
                                  int distribution, 
                                  char* sd_filename,
                                  double r_eff_min,
                                  double r_eff_max,
                                  double lambda,      /*wavelength in micrometer*/
                                  double dx_max,
                                  double n_r_max,
                                  int verbose,
                                  /* Output */
                                  int* n_dens,
                                  double** radius, 
                                  double** number_dens);

int calc_distribution_number_dens (/* Input */
                                   int distribution,
                                   double alpha,
                                   double r_eff,
                                   double rho_medium,
                                   int n_dens,
                                   double* radius,
                                   int verbose,
                                   /* Output */
                                   double** number_dens
                                   );

int calc_phase_function(/* Input */
                        int  nmom,
                        double lambda, 
                        double reff,
                        int nstokes, 
                        int nthetamax,
                        double accuracy,
                        /* Output */
                        float*** pmom,
                        int* nmommax,
                        float*** phase, 
                        float*** theta,
                        int** ntheta,
                        int* warnings_nmom,
                        int* warnings_ntheta);

double double_min(double x, double y);

double double_max(double x, double y);

double gamma_ln(float xx);

int integrate_sizedist(/* Input */
                       float* qext_array, 
                       float* qsca_array,
                       float* gsca_array,
                       float* qback_array,
                       float*** pmom_array,
                       int nstokes, int nmom,
                       double *x_size, double *y_size, int n_size, 
                       int ir, int verbose, 
                       /* Output */
                       double *beta,
                       double *omega, 
                       double *g,
                       double *back,
                       float*** pmom);

int mie_calc_all_radii(/* Input */
                       int program, int medium, mie_complex crefin, 
                       double wavelength, float temperature, int nstokes,
                       double *x_size, int n_size, 
                       mie_inp_struct input, 
                       int verbose,
                       /* Output */
                       mie_out_struct *output, mie_complex *ref,
                       float **qext_array, float **qsca_array, float **gsca_array,
                       float **qback_array, float ****pmom_array);


int read_refractive_index(/* Input */
                          char* filename,
                          int verbose, 
                          /* Output */
                          int* nlambda, double **lambda,
                          double ***ri_real, double ***ri_imag);

int setup_reff_grid(/* Input */
                    int distribution,
                    double r_eff_step,
                    double r_eff_min, 
                    double r_eff_max, 
		    int r_eff_log,
                    int n_hum,           /* nhum is used for aerosol*/
                    double* hum,
                    int verbose,
                    /* Output */
                    int* n_r_eff,  
                    double** r_eff);

int setup_wavelength_grid(/* Input */
                         double wl_start, double wl_end, double wl_step,
                         char* data_files_path, char* filename_ck_generic,
                         int ck_scheme, int lambda_i_start, int lambda_i_end,
                         int verbose, 
                         /* Output */
                         int* nlambda, double** lambda);


int write_output (input_struct input, output_struct *output, int ir, int iv);


#if HAVE_LIBNETCDF

int create_netcdf_file(/* Input */
                       char *basename,
                       int medium, 
                       int aerosol_type,
                       int distribution, 
                       int alpha,
                       int nlam,
                       double* lambda,
                       int nreff, 
                       double* reff, 
                       int nmommax, 
                       int nstokes,
                       int nthetamax, 
                       int nrho, 
                       int complete_moments,
                       double n_r_max,
                       double dx_max,
                       char *sd_filename,
                       /* Output */
                       netcdf_ids* ncid,
                       char** outputfile);

int write_output_netcdf(netcdf_ids ncid, input_struct input, output_struct *output, int iv, int ir);

#endif

/**
 * mie
 * 
 * Main function of mie tool.
 *
 *
 * @author Bernhard Mayer, Ulrich Hamann, Claudia Emde
 * @date  ????? 
 *        2010-02-11 Restructured and cleaned (CE)        
 */
int mie (input_struct input, output_struct *output)
{
 
  int ip=0, im=0, ir=0, iv=0, status=0;
  int first=1;
  
  double wavelength=0;

  float *qext_array;
  float *qsca_array;
  float *gsca_array;
  float *qback_array;
  float ***pmom_array;  /* Phase function moments */
  double norm=0.0; 
  
  double **ref_real;
  double **ref_imag;

  double ext, vol;
  int lambda_i_start=0;
  int lambda_i_end=0;
  
  mie_complex crefin;
  crefin.re=0.0;
  crefin.im=0.0;
    
  output->ref.re = 0.0;
  output->ref.im = 0.0; 

#if HAVE_LIBNETCDF
  netcdf_ids ncid;
  int nmommax;
  int complete_moments=0; 
#endif


  /***********************************************************/
  /* Setup                                                   */
  /***********************************************************/
  
  /* Read input data and specify wavelengths grid*/ 
  if (strlen(input.data_files_path) == 0)
    input.data_files_path="../data/"; 


  /* Check whether refractive index file is defined and if yes get
     wavelength grid from refractive index file */
  if(strlen (input.ri_filename)>0) {
    status = read_refractive_index(input.ri_filename, input.verbose, 
                                   &(output->wl.nlambda), &(output->wl.lambda),
                                   &ref_real, &ref_imag );
    if (status!=0) {
      fprintf (stderr, "Error %d returned by read_refractive_index()\n", status);
      return status;
    }
  }
  /* Else set up the wavelenth grid from input file specification */
  else {
    status = setup_wavelength_grid(input.wl.wl_start, input.wl.wl_end, input.wl.wl_step, 
                                   input.data_files_path, input.filename_ck_generic,
                                   input.ck_scheme, input.lambda_i_start, input.lambda_i_end,
                                   input.verbose, 
                                   &(output->wl.nlambda), &(output->wl.lambda));
     if (status!=0) {
      fprintf (stderr, "Error %d returned by setup_wavelengths_grid()\n", status);
      return status;
     }
  }
  if (input.lambda_i_start != -999){
    lambda_i_start = input.lambda_i_start;
    lambda_i_end = input.lambda_i_end;
  }
  else{
    lambda_i_start = 0;
    lambda_i_end = output->wl.nlambda-1;                                
  }
  
  /* Read aerosol distribution parameters from OPAC*/
  if (input.distribution==DIST_AER){
    status = aerosol_input_data(input.data_files_path, input.aerosol_type, input.verbose, 
                                &(output->aerdistr));                               
    if (status!=0) {
      fprintf (stderr, "Error %d returned by aerosol_input_data()\n", status);
      return status;
    }
    status = aerosol_refractive_index(output->wl.lambda, output->wl.nlambda, 
                                      input.data_files_path, 
                                      input.aerosol_type, output->aerdistr.n_hum,
                                      input.verbose, 
                                      &ref_real, &ref_imag);
    if (status!=0) {
      fprintf (stderr, "Error %d returned by aerosol_refractive_index()\n",
               status);
      return status;
    }
  }
  
  /* Set up the effective radius grid */
  status=setup_reff_grid(input.distribution, input.r_eff_step, input.r_eff_min, 
                         input.r_eff_max, input.r_eff_log,
                         output->aerdistr.n_hum, output->aerdistr.rel_hum,
                         input.verbose, 
                         &(output->n_r_eff),
                         &(output->r_eff));
  if(status!=0){
    fprintf (stderr, "Error %d returned by setup_reff_grid()\n", status);
    return status;
  }
  
  /**********************************************************/
  /* Memory allocation                                      */
  /**********************************************************/
  
  /* Legendre polynomials of all phase matrix elements */
  if(input.mie.nmom > 0){
    /* Output of Mie program for one size */
    output->mie.pmom = calloc (input.nstokes, sizeof(float *));
    for (ip=0; ip<input.nstokes; ip++)
      output->mie.pmom[ip] = calloc (input.mie.nmom+1, sizeof(float));
    
    /* Integrated, normalized and converted moments */
    output->pmom = calloc (input.nstokes, sizeof(double *));
    for (ip=0; ip<input.nstokes; ip++) 
      output->pmom [ip] = calloc (input.mie.nmom+1, sizeof(double));
  } 
  
  /* Warnings */
  output->warnings_nmom = calloc (output->wl.nlambda, sizeof(int *));
  output->warnings_ntheta = calloc (output->wl.nlambda, sizeof(int *));
  for (iv=0; iv<output->wl.nlambda; iv++){
    output->warnings_nmom[iv] = calloc (output->n_r_eff, sizeof(int));
    output->warnings_ntheta[iv] = calloc (output->n_r_eff, sizeof(int));
  }
  
  if(input.output_format == OUTPUT_FORMAT_NETCDF){
    
    output->outputfile=calloc(FILENAME_MAX, sizeof(char)); 
    
    /* allocate memory for phase function*/
    output->ntheta=calloc(input.nstokes, sizeof(int));
    output->theta=calloc(input.nstokes, sizeof(float *));
    output->phase=calloc(input.nstokes, sizeof(float *));
    for  (ip=0; ip<input.nstokes; ip++){
      output->theta[ip]=calloc(input.nthetamax, sizeof(float));
      output->phase[ip]=calloc(input.nthetamax, sizeof(float));
    }
  }
  
  /**************************************************************/
  /* Start calculation                                          */
  /**************************************************************/     
  
  
   /* loop over wavelengths */
  for (iv=lambda_i_start; iv<=lambda_i_end; iv++) {
    
    if(input.verbose)
      fprintf (stderr, " ... start Mie-calculation for lambda = %10.3f nm  \n", 
               output->wl.lambda[iv]);
    
    
    /* Convert from nm to micro meter */
    wavelength = (output->wl.lambda[iv])*1.e-03;  
    
    /* Aerosol treatment different from clouds for consistency reasons
       with OPAC, size distribution calculated below in loop over
       relative humudities. For each relative humidity a new Mie
       calculation is required, because the refractive index changes. */
    if (input.distribution != DIST_AER){
      /* generate the radius grid of the size distribution */
      status = calc_distribution_radius_grid (input.distribution, input.sd_filename,
                                              input.r_eff_min, input.r_eff_max,
                                              wavelength, input.dx_max, input.n_r_max,
                                              input.verbose, 
                                              &output->dist.n_dens, &output->dist.radius, 
                                              &output->dist.number_dens); 
      if (status!=0) {
        fprintf (stderr, "Error %d returned by calc_distribution_radius_grid()\n", status);
        return status;
      }
    }
    
       
    /* Loop over effective radii, relative humidities (aerosol)*/
    for (ir=0; ir<output->n_r_eff; ir++) {
    
      /* refractive index */ 
      if (strlen (input.ri_filename)>0) {
        crefin.re = ref_real[0][iv];
        crefin.im = ref_imag[0][iv];
      }
      else if (input.distribution == DIST_AER){
        crefin.re = ref_real[ir][iv];
        crefin.im = ref_imag[ir][iv];
        fprintf(stderr, " ... aerosol refractive index %g %g \n", 
                crefin.re, crefin.im);
      }
      else{
        crefin.re = input.crefin.re;
        crefin.im = input.crefin.im;
      }
        
      
      if (input.distribution != DIST_AER){
        /* number density distribution per radius intervall*/ 
        status = calc_distribution_number_dens (input.distribution,
                                                input.alpha, output->r_eff[ir], 
                                                output->rho_medium, output->dist.n_dens,
                                                output->dist.radius, input.verbose,  
                                                &(output->dist.number_dens));
        if (status!=0) {
          fprintf (stderr, "Error %d returned by calc_distribution_number_dens()\n", status);
          return status;
        }
      }
      else{
        /* Set up aerosol size distribution radius grid and calculate
           number densities.*/ 
        status = aerosol_size_distribution(output->aerdistr.r_min[ir], 
                                           output->aerdistr.r_max[ir],
                                           output->aerdistr.rmod[ir],
                                           output->aerdistr.sigma,
                                           output->wl.lambda[iv],
                                           input.verbose, 
                                           &output->dist.n_dens, 
                                           &(output->dist.radius),
                                           &(output->dist.number_dens)
                                           );
        if (status!=0) {
          fprintf (stderr, "Error %d returned by aerosol_size_distribution()\n",status);
          return status;
        }
        /* refractive index for the required relative humidity */
        output->rho_medium=output->aerdistr.rho[ir];
        
        
      } 
  
      /* Allocate memory for optical properties for each size in the
         size distribution */
      if (ir==0 || input.distribution == DIST_AER){
        qext_array = (float *) calloc (output->dist.n_dens, sizeof(float));
        qsca_array = (float *) calloc (output->dist.n_dens, sizeof(float));
        gsca_array = (float *) calloc (output->dist.n_dens, sizeof(float));
        qback_array = (float *) calloc (output->dist.n_dens, sizeof(float));
        
        pmom_array = calloc (input.nstokes, sizeof(float **));
        for (ip=0; ip<input.nstokes; ip++) {
          pmom_array[ip] = calloc (input.mie.nmom+1, sizeof(float *));
          for (im=0; im<input.mie.nmom+1; im++)
            pmom_array[ip][im] = calloc (output->dist.n_dens, sizeof(float));
        }
      }

      /***************************************************************/
      /* call Mie programs                                           */
      /***************************************************************/
      
      switch(input.distribution) {
      case DIST_NONE:
          /* one single radius */
        status = mie_calc (input.mie, &(output->mie),
                           input.program, input.medium, crefin,
                           wavelength, output->r_eff[ir], input.temperature, 
                           input.nstokes,
                           &(output->ref));
        if (status!=0) {
          fprintf (stderr, "Error %d returned by mie_calc()\n", status);
          return status;
        }
        /* calculate beta and omega */
        ext = PI * output->r_eff[ir] * output->r_eff[ir]* output->mie.qext *1e-12;
	//20120816ak sca is not used anywhere, commented
	//        sca = PI * output->r_eff[ir]* output->r_eff[ir] * output->mie.qsca *1e-12;
        vol = (4.0/3.0) * PI * pow((output->r_eff[ir]*1E-6),3.0); 
        output->beta  = ext / vol / 1000.0;
        output->omega = output->mie.qsca / output->mie.qext ;
        output->back = PI * output->r_eff[ir] * output->r_eff[ir]* output->mie.qback *1e-12 /vol / 1000.0; 
        
        /* Moments of phase matrix */
        if (input.mie.nmom>0)
          {
            if(input.nstokes==4){
              /* Copy moments */
              for (ip=0; ip<input.nstokes; ip++) 
                for (im=0; im<=input.mie.nmom; im++)
                  output->pmom[ip][im] = (output->mie.pmom)[ip][im];
              
              /*  Convert from Wiscombe definition to Evans definition of Mueller matrix */
              norm=output->mie.pmom[0][0]+output->mie.pmom[1][0];
              for (im=0; im<=input.mie.nmom; im++){ 
                output->pmom[0][im]=(output->mie.pmom[0][im]+output->mie.pmom[1][im])/norm ;        
                output->pmom[1][im]=(output->mie.pmom[1][im]-output->mie.pmom[0][im])/norm;  
                output->pmom[2][im]= 2.*output->mie.pmom[2][im]/norm;        
                output->pmom[3][im]= 2.*output->mie.pmom[3][im]/norm;      
              }
            }
            else if(input.nstokes==1){
              norm=output->mie.pmom[0][0];
              for (im=0; im<=input.mie.nmom; im++)
                output->pmom[0][im]=output->mie.pmom[0][im]/norm;
            }
            else{
              fprintf(stderr, "Error, nstokes must be either 1 or 4\n");
              return -1;
            }
          }
        
        
        break;
      case DIST_FILE:
      case DIST_GAMMA:
      case DIST_LOGNORMAL:
      case DIST_AER:  
        
        /* droplet radius distribution */
        if(ir==0 || input.distribution==DIST_AER){
          status = mie_calc_all_radii (input.program, input.medium, crefin, 
                                       wavelength, input.temperature, input.nstokes,
                                       output->dist.radius, 
                                       output->dist.n_dens,
                                       input.mie, 
                                       input.verbose,
                                       &(output->mie), 
                                       &(output->ref),
                                       &qext_array, &qsca_array, 
                                       &gsca_array, &qback_array, &pmom_array);
          
          if (status!=0) {
            fprintf (stderr, "Error %d returned by mie_calc_all_radii()\n", status);
            return status;
          }
        }
        
        status = integrate_sizedist (qext_array, qsca_array, 
                                     gsca_array, qback_array, pmom_array,
                                     input.nstokes, input.mie.nmom, 
                                     output->dist.radius, 
                                     output->dist.number_dens, output->dist.n_dens,
                                     ir, input.verbose, 
                                     &(output->beta), &(output->omega),
                                     &(output->g), &(output->back), &(output->pmom)); 
        
        if (status!=0) {
          fprintf (stderr, "Error %d returned by integrate_sizedist()\n", status);
          return status;
        }
        break;
      default:
        fprintf (stderr, "Error, distribution number in mie.c %d\n", input.distribution);
        return -1;
        
        break;
      }

      /***************************************************************/
      /* write output                                                */
      /***************************************************************/
      
      
      if(input.output_format != OUTPUT_FORMAT_NETCDF){
        /* write output */
        status = write_output (input, output, ir, iv);
        if (status!=0) {
          fprintf (stderr, "Error %d returned by write_output()\n", status);
          return status;
        }
      }
#if HAVE_LIBNETCDF
      else{
        
        status=calc_phase_function(input.mie.nmom, 
                                   wavelength, output->r_eff[ir], input.nstokes,
                                   input.nthetamax, input.accuracy,
                                   &(output->pmom), &(output->nmommax), 
                                   &(output->phase), &(output->theta), 
                                   &(output->ntheta),
                                   &(output->warnings_nmom)[iv][ir],
                                   &(output->warnings_ntheta)[iv][ir]);
        if (status!=0) {
          fprintf (stderr, "Error %d returned by calc_phase_function()\n", status);
          return status;
        }
        
       
        
        if (iv == 0  && ir ==0){
          if (input.nmom_netcdf==0){
            nmommax=input.mie.nmom;
            complete_moments=1;
          }
          else{
            nmommax=input.nmom_netcdf;
            complete_moments=0;
          }
          
          status=create_netcdf_file(input.basename, input.medium, input.aerosol_type, 
                                    input.distribution, input.alpha, 
                                    output->wl.nlambda, output->wl.lambda,
                                    output->n_r_eff, output->r_eff,
                                    nmommax, 
                                    input.nstokes, 
                                    input.nthetamax,
                                    output->aerdistr.n_hum,
                                    complete_moments,
                                    input.n_r_max,
                                    input.dx_max,
                                    input.sd_filename,
                                    &ncid, &(output->outputfile));
          if (status!=0) {
            fprintf (stderr, "Error %d returned by create_netcdf_file()\n", status);
            return status;
          }
        }
        
        write_output_netcdf(ncid, input, output, iv, ir);
        
      }
#else 
      fprintf (stderr, " ***********************************************************************\n");
      fprintf (stderr, " * You have built mie without libnetcdf and hence cannot               *\n");
      fprintf (stderr, " * use the option *output_user netcdf*. Please get netcdf and rebuild. *\n");
      fprintf (stderr, " ***********************************************************************\n");
      return -1;    
#endif
      
      /* Free variables, for aerosol the number of points in the size
         distribution is different for each relative humidity (ir),
         therefore these variables are allocated for each ir.*/
      if (ir==output->n_r_eff-1 || input.distribution==DIST_AER){
        for (ip=0; ip<input.nstokes; ip++) {
          for (im=0; im<input.mie.nmom+1; im++)
            free(pmom_array[ip][im]);
          free(pmom_array[ip]);
        }
        
        free(pmom_array);
        
        free(qext_array);
        free(qsca_array);
        free(gsca_array);
        
        free(output->dist.radius);
        free(output->dist.number_dens);
      }
      
    } /* for (ir=0; ir<output->n_r_eff; ir++) { */ 
    
  }/* loop over iv */
  
  /* write warnings*/
  if (input.output_format == OUTPUT_FORMAT_CLOUDPROP ||
      input.output_format == OUTPUT_FORMAT_AERPROP ||
      input.output_format == OUTPUT_FORMAT_NETCDF) {
    
    for (iv=0; iv<output->wl.nlambda; iv++) {
      for (ir=0; ir<output->n_r_eff; ir++) {
        if (output->warnings_nmom[iv][ir]==1) {
          if (first == 1) { /* true */
	      fprintf (stderr, "*** Warning: maxmimum accuracy of the phase function not reached\n");
              fprintf (stderr, "*** with %4d Legendre coefficients.\n",input.mie.nmom);
              fprintf (stderr, "*** In order to get more accurate results\n"); 
              fprintf (stderr, "*** please increase 'nmom' in the mie input file!\n"); 
              first = 0;
          }
          if(input.medium !=AEROSOL)
            fprintf (stderr, "          lambda = %10.3f nm, r_eff = %8.3f micro meter\n", output->wl.lambda[iv], output->r_eff[ir]);
          else
            fprintf (stderr, "          lambda = %10.3f nm, rel_hum = %2f percent\n", output->wl.lambda[iv], output->aerdistr.rel_hum[ir]);
        }
      } 
    } 
    
    first=1; 
    for (iv=0; iv<output->wl.nlambda; iv++) {
      for (ir=0; ir<output->n_r_eff; ir++) {
        if (output->warnings_ntheta[iv][ir]==1) {
          if (first == 1) {  
            fprintf (stderr, "*** Warning: maxmimum accuracy of the phase function not reached\n");
            fprintf (stderr, "*** with %4d scattering angles.\n",input.nthetamax);
            fprintf (stderr, "*** In order to get more accurate results\n"); 
            fprintf (stderr, "*** please increase 'nthetamax' in the mie input file!\n"); 
            first = 0;
          }
          if(input.medium !=AEROSOL)
            fprintf (stderr, "          lambda = %10.3f nm, r_eff = %8.3f micro meter\n", output->wl.lambda[iv], output->r_eff[ir]);
          else
            fprintf (stderr, "          lambda = %10.3f nm, rel_hum = %2f percent\n", output->wl.lambda[iv], output->aerdistr.rel_hum[ir]);
        }
      } 
    } 
    
  }
  
  /*free memory */
  for (iv=0; iv<output->wl.nlambda; iv++){
    free(output->warnings_nmom[iv]); 
    free(output->warnings_ntheta[iv]); 
  }    
  free(output->warnings_nmom);
  free(output->warnings_ntheta);
  
  if(input.output_format == OUTPUT_FORMAT_NETCDF){
   
    free(output->outputfile);
    for (ip=0; ip<input.nstokes; ip++){
      free(output->phase[ip]); 
      free(output->theta[ip]);
    }
    free(output->phase); 
    free(output->theta);
    free(output->ntheta); 
  }
  
  return 0;
}


/**
 * setup_reff_grid
 * 
 * Set up the effective radius grid for the Mie calculation from user
 * input parameters. 
 * 
 * @author      ..., Claudia Emde
 * @date        2010-02-11 put into subroutine by C. Emde
 */
int setup_reff_grid(/* Input */
                  int distribution,
                  double r_eff_step,
                  double r_eff_min, 
                  double r_eff_max, 
		  int r_eff_log,
                  int n_hum,           /* nhum is used for aerosol*/
                  double* hum, 
                  int verbose, 
                  /* Output */
                  int* n_r_eff,  
                  double** r_eff)
{
  int status=0; 

  int ir=0; 
  
  /* determine the number of r_eff to be simulated */
  switch(distribution) {
    
  case DIST_NONE:
  case DIST_GAMMA:
  case DIST_LOGNORMAL:
    /* determine number of r_eff to be calculated */
    if (r_eff_step != 0.0)  {
      if (!r_eff_log)
	/* 1e-7 avoids rouding errors */
	*n_r_eff = (int) floor((r_eff_max - r_eff_min+1e-7) / r_eff_step)+1; 
      else {
	if (r_eff_step<=1) {
	  fprintf (stderr, "Error, effective radius increment needs to be larger\n");
	  fprintf (stderr, "than 1 for log steps.\n");
	  return -1;
	}
	*n_r_eff = (int) (ceil(log (r_eff_max/r_eff_min)/log(r_eff_step)))+1;
      }
    }
    else
      *n_r_eff = 1;

    /* setup r_eff grid */  
    *r_eff = (double *) calloc (*n_r_eff, sizeof(double));
    if(verbose)
      fprintf(stderr, "... effective radii: \n");

    (*r_eff)[0] = r_eff_min;
    for (ir=1; ir<*n_r_eff; ir++) {
      if (!r_eff_log)
	(*r_eff)[ir] = r_eff_min + ir * r_eff_step;
      else 
	(*r_eff)[ir] = (*r_eff)[ir-1] * r_eff_step;

      if (verbose)
	fprintf(stderr, "%g ", (*r_eff)[ir]);
    }

    if (verbose)
      fprintf(stderr, "\n");
        
    break;
    
  case DIST_FILE:
    /* size distribution file means only one r_eff */
    *n_r_eff = 1;
    *r_eff = (double *) calloc (*n_r_eff, sizeof(double));
    (*r_eff)[0] = -999.0; /* effective radius value is not required */
    break;
    
  case DIST_AER:
    /* Effective radius is "replaced" by relative humidity for aerosols */
    *n_r_eff=n_hum;
    *r_eff = (double *) calloc (*n_r_eff, sizeof(double));
    for (ir=0; ir<n_hum; ir++)
      (*r_eff)[ir]=hum[ir];
    fprintf(stderr, "aerosol nhum %d\n", *n_r_eff);
    
    break; 
  default:
    fprintf (stderr, "Error, distribution number in mie.c %d\n", distribution);
    return -3;
    break;
  }
  return status;
}

/**
 * read_refractive_index
 *
 * Read the refractive index from a file
 *
 * @author ...., Claudia Emde
 * @date   2010-02-11 restructured by C. Emde 
 */
int read_refractive_index(/* Input */
                          char* filename, 
                          int verbose, 
                          /* Output */
                          int* nlambda, double **lambda,
                          double ***ref_real, double ***ref_imag)
{
  int status = 0;
  
  int iv=0; 
  
  /* Here we have only one refractive index for each wavelength, this
     is different for aerosols. */
  *ref_real = calloc(1, sizeof(double *)); 
  *ref_imag = calloc(1, sizeof(double *)); 
  status = read_3c_file (filename, &(*lambda), &(*ref_real)[0], &(*ref_imag)[0], 
                         &(*nlambda));
  
  if (status != 0) {
    fprintf (stderr, "error %d reading %s\n", status, filename);
    return status;
  }
  
  fprintf (stderr, " ... read %d data points from %s\n", 
           *nlambda, filename);
  
  /* check if data points are sorted in ascending order */
  for (iv=0; iv < *nlambda-1; iv++){
    if ( (*lambda)[iv+1] < (*lambda)[iv]) {
      fprintf (stderr, "Error: %s not sorted in ascending order\n", 
               filename);
      return -1;
    }
  }
  
  return 0;
}


/**
 * setup_wavelength_grid
 *
 * Setup the wavelength grid for the Mie calculation from user
 * specified input parameters, read wavelengths from file, or specify
 * wavelengths for k-distributions.
 *
 * @author ...., Claudia Emde
 * @date   2010-02-11 restructured by C. Emde
 * 
 */

int setup_wavelength_grid(/* Input */
                         double wl_start, double wl_end, double wl_step,
                         char* data_files_path, char* filename_ck_generic,
                         int ck_scheme, int lambda_i_start, int lambda_i_end,   
                         int verbose, 
                         /* Output */
                         int* nlambda, double** lambda)
{  
  int status = 0;
  int iv=0;
  char *wvl_filename=NULL; 

  int min_columns=0;
  int max_columns=0;
  int max_length=0;

  char ***string=NULL;
  char *dummy=NULL;

  int unit_factor=1;    /* specifies if wvl contain micro meter or nm */

  char function_name[]="setup_wavelength_grid";
  char file_name[]="mie.c";
  
  if (ck_scheme == CK_CRS) {
    if (wl_step != 0.0)
      *nlambda = (int) (wl_end-wl_start) / wl_step + 1;
    else
      *nlambda = 1;
    
    *lambda = (double *) calloc (*nlambda, sizeof(double));
    
    if (verbose)
      fprintf(stderr,"... wavelengths:\n");
    for (iv=0; iv<*nlambda; iv++){
      (*lambda)[iv] = wl_start + wl_step*iv;
      if(verbose)
        fprintf(stderr, "%g ", (*lambda)[iv]);
    }
    if (verbose)
      fprintf(stderr,"\n"); 
  }
  
  else {
    
    wvl_filename = calloc (FILENAME_MAX,sizeof(char));
    
    strcpy(wvl_filename, data_files_path);

    /* Note: For the correlated-k schemes the optical properties are
       calculated for the center wavelengths. Evans' cloudprop does a
       weighting over the Planck functions (Sun, Earth). */ 
    switch(ck_scheme) {
    case CK_KATO:
      strcat(wvl_filename,"correlated_k/kato/wvl.dat");
      fprintf (stderr, " ... Kato wavelength bands\n");
      unit_factor=1;  /* Kato wvl.dat contains nano meter */
      break;
    case CK_KATO2:
      strcat(wvl_filename,"correlated_k/kato2/wvl.dat");
      fprintf (stderr, " ... Kato2 wavelength bands\n");
      unit_factor=1;  /* Kato2 wvl.dat contains nano meter */
      break;
    case CK_KATO2ANDWANDJI:
      strcat(wvl_filename,"correlated_k/kato2andwandji/wvl.dat");
      fprintf (stderr, " ... Kato2andwandji wavelength bands\n");
      unit_factor=1;  /* Kato2 wvl.dat contains nano meter */
      break;
    case CK_KATO2_96:
      strcat(wvl_filename,"correlated_k/kato2.hitran96/wvl.dat");
      fprintf (stderr, " ... Kato2_96 wavelength bands\n");
      unit_factor=1;  /* Kato2_96 wvl.dat contains nano meter */
      break;
    case CK_FU:
      strcat(wvl_filename,"correlated_k/fu/wvl.dat");
      fprintf (stderr, " ... Fu wavelength bands\n");
      unit_factor=1000;  /* Fu wvl.dat contains micro meter */
      break;
    case CK_AVHRR_KRATZ:
      strcat(wvl_filename,"correlated_k/kratz/avhrr_wvl.dat");
      fprintf (stderr, " ... AVHRR-Kratz wavelength bands\n");
      unit_factor=1000;  /* AVHRR-Kratz wvl.dat contains micro meter */
      break;
    case CK_GENERIC:
      strcpy (wvl_filename, filename_ck_generic);
      fprintf (stderr, " ... generic correlated_k wavelength bands\n");
      break;
    default:
      fprintf (stderr, "Error: unsupported correlated-k scheme, nr %2d\n", ck_scheme);
      return -1;
    }
    
    status = ASCII_checkfile(wvl_filename, &(*nlambda), &min_columns, &max_columns, &max_length);
    if (status!=0) {
      fprintf(stderr,"Error %d reading wavelength_file %s in %s (%s) \n", status, wvl_filename, function_name, file_name);
      return status;
    }

    if ( min_columns < 2 ) {
      fprintf(stderr,"Error, wavelength_file %s must have at least two columns! \n", wvl_filename);
      return -1;
    }
    
    if ((status = ASCII_calloc_string(&string, *nlambda, max_columns, FILENAME_MAX+1)) != 0) {
      fprintf(stderr,"Error %d allocating string in %s (%s)\n", status, function_name, file_name);
      return status;
      }
    
    if ((status = ASCII_readfile(wvl_filename, string)) != 0) {
      fprintf(stderr,"Error %d reading wavelength-file: %s in %s (%s) \n", status, wvl_filename, function_name, file_name);
      return status;
    }
    
    if ((*lambda = (double *) calloc (*nlambda, sizeof(double))) == NULL) {
      fprintf(stderr,"Error %d allocating output->wl.lambda in %s (%s) \n", status, function_name, file_name);
      return status;
    }
    
    /* second column is center wavelength */
    for (iv=0; iv< *nlambda; iv++) {
      (*lambda)[iv] = (double) strtod (string[iv][1], &dummy) * unit_factor;
    }
    
    if(verbose)
      for (iv=0; iv < *nlambda; iv++)
        fprintf (stderr, "   lambda[%3d] = %10.3f nm\n", iv+1, (*lambda)[iv]);
    
    
    /* Check validity of specified wavelength indices */
    if (lambda_i_start != -999 && lambda_i_start < 1) {
      fprintf(stderr,"  wavelength_index_start = %d has to be 1 or larger \n", lambda_i_start);
      status=-1;
    }
    if (*nlambda < lambda_i_start-1) {
      fprintf(stderr,"  wavelength_index_start = %d, but there are only %d wavelength bands \n",
              lambda_i_start, *nlambda);
      status=-2;
    }
    
    
    if (lambda_i_end != -999  && lambda_i_end < lambda_i_start) {
      fprintf(stderr,"  wavelength_index_end = %d has to be larger than wavelength_index_start = %d\n",
              lambda_i_end, lambda_i_start);
      status=-3;
    }
    if (*nlambda < lambda_i_end) {
      fprintf(stderr,"  wavelength_index_end = %d, but there are only %d wavelength bands \n",
                        lambda_i_end, *nlambda);
      status=-4;
    }
  }
  return status;
}

/**
 * calc_distribution_radius_grid
 * 
 * Calculate a radius grid for the size distribution over which the
 * optical properties are intgrated. This grid must be sufficiently
 * fine, otherwise "ripples" in the phase function are not averaged
 * out. Here the sampling is done in the same way as in cloudprp by
 * F. Evans.
 *
 * @author ???
 * @date   ???
 */
int calc_distribution_radius_grid(/* Input */
                                  int distribution, 
                                  char* sd_filename,
                                  double r_eff_min,
                                  double r_eff_max,
                                  double lambda,      /*wavelength in micrometer*/
                                  double dx_max,
                                  double n_r_max,
                                  int verbose,
                                  /* Output */
                                  int* n_dens,
                                  double** radius, 
                                  double** number_dens 
                                  )
{
  double radmin=0.0, radmax=0.0, rad=0.0, delta_rad=0.0;
  double x=0.0, delta_x=0.0, dx=0.0;
  int is=0;
  int status=0;
  
  /* fprintf(stderr,"calculate droplet radius distribution\n"); */

  switch(distribution) {
  
  case DIST_NONE:
    /* single radius calculation, nothing to do */
    *n_dens = 1;
    break;
    
  case DIST_FILE:
    /* read size distribution from file */
    status = read_2c_file (sd_filename, &(*radius), &(*number_dens), &(*n_dens));
    
    if (status != 0) {
      fprintf (stderr, "error %d reading %s\n", status, sd_filename);
      return status;
    }

    if (verbose)
      fprintf (stderr, " ... read %d data points from %s\n", 
               *n_dens, sd_filename);
    
    break;
  case DIST_GAMMA:
  case DIST_LOGNORMAL:
    /* ????? are these reasonable values? Evans used 0.02 reff_min - 5 reff_max; ????? */
    /* ????? BM found that 5 is not enough for small particles - but somebody    ????? */
    /* ????? increased the upper boundary to 8 reff_max which is VERY EXPENSIVE  ????? */
    
    radmin = 0.02 * r_eff_min; /* cover the whole droplet spectrum   */
    radmax = n_r_max * r_eff_max; /* in order to save mie-calculations  */
    
    /* CE used n_r_max = 5 to calculate new mie data */
    /* radmax = 5.00 * input.r_eff_max;*/
    
    /* CE ????? */
    /* dx is the maximum bin size (size parameter) in size distribution grid */
    /* Original value is 0.03, this seems to be unsufficient, especially */
    /* for backward scattering region, also in rainbow region lots of wiggles. */
    /* May be we should think of another way to sample the size distribution */
    dx=dx_max;
    
    /* count number of radii */
    rad = radmin;
    is = 0;
    while (rad <= radmax) {
      x = 2*PI*rad / lambda; 
      delta_x = double_max(dx,dx*sqrt(x));
      delta_rad = delta_x * lambda / (2*PI);
      rad = rad + delta_rad;
      is = is + 1;
    }
    *n_dens = is;
    
    if (verbose)
      fprintf(stderr, "size distributions: r_min %g um, r_max %g um, number of grid points %d \n", radmin, radmax, *n_dens); 
    
    if (((*radius) = calloc(*n_dens, sizeof(double)))==NULL) {
      fprintf (stderr, "Error: Allocation of output->dist.radius (make_distribution)\n");
      return -2;
    }
    if (((*number_dens) = calloc(*n_dens, sizeof(double)))==NULL) {
      fprintf (stderr, "Error: Allocation of output->dist.numb_dens (make_distribution)\n");
      return -3;
    } 

    /* calculate radius grid */
    rad = radmin;
    is = 0;
    while (rad <= radmax) {
      (*radius)[is]      = rad;
      x = 2*PI*rad / lambda; 
      delta_x = double_max(dx,dx*sqrt(x));
      delta_rad = delta_x * lambda / (2*PI);
      rad = rad + delta_rad;
      is = is + 1;
    }
    break; 
    
  case DIST_AER: 
    /* Do nothing here. Radius grid and number densities are calculated in
       aerosol_size_distribution. */
    break;
  default:
    fprintf (stderr, "Error in calc_distribution_radius_grid: unrecognized distribution %d\n", distribution);
    return -1;
    break;
  } 
  
  /* check if data points are sorted in ascending order */
  switch(distribution) {
  case DIST_NONE:
    break;
  case DIST_FILE:
  case DIST_GAMMA:
  case DIST_LOGNORMAL:
    /* for gamma and lognormal this is not nessesary, but it costs almost no time */
    
    for (is=0; is<*n_dens-1; is++)
      if ( (*radius)[is+1]< (*radius)[is]) {
        if (distribution == DIST_FILE)
          fprintf (stderr, "Error: %s not sorted in ascending order\n", 
                   sd_filename);
        else
          fprintf (stderr, "Error: radius in distribution not sorted in ascending order, internal program bug!!!\n");
        return -5;
      }
    
    if ((*radius)[0]==0) {
      (*radius)[0] = 0.1 * (*radius)[1];
      fprintf (stderr, " ... cannot handle r=0; setting first data point to r=%g\n",
               (*radius)[0]);
    }
    break;
  default:
    fprintf (stderr, "Error in calc_distribution_radius_grid: unrecognized distribution %d\n", distribution);
    break;
  }
  return status;
}

/**
 * calc_distribution_number_dens
 *
 * Calculate particle number densities for each radius in the size
 * distribution.
 *
 * @author ???
 * @date   ???
 */
int calc_distribution_number_dens (/* Input */
                                   int distribution,
                                   double alpha,
                                   double r_eff,
                                   double rho_medium,
                                   int n_dens,
                                   double* radius,
                                   int verbose, 
                                   /* Output */
                                   double** number_dens
                                   )
{
  double A=0.0, B=0.0;
  float LWC=1.0;
  int is=0;
  int status=0;
  double lnA=0;

  /* fprintf(stderr,"calculate droplet radius distribution\n"); */

  if (rho_medium <= 0.0) {
    fprintf (stderr, "Error, rho_medium = %f is 0 or negative!\n", rho_medium);
    return -1;
  }
  
  
  switch(distribution) {
  case DIST_NONE:
  case DIST_FILE:
    /* single radius calculation, nothing to do */
    break;
    
  case DIST_GAMMA:
    
    /* if (verbose) */
    /*       fprintf(stderr, "Gamma size distribution, alpha=%f \n", alpha); */
  
    B = (alpha+3.0)/r_eff;
    /* original version by Claudia Emde */
    /* A = (0.75/PI)*pow(B,(alpha+4))*1000./exp(gamma_ln(alpha+4.)); */
    /* A = A*LWC/rho_medium;                                         */
    
    /* new version by Bernhard Mayer, 31.8.2016 */
    lnA = log(0.75/PI*1000.0*LWC/rho_medium) + (alpha+4)*log(B) - gamma_ln(alpha+4.);

    /* calculate number density */
    for (is=0;is<n_dens;is++) {
      /* replaced gamma function by exp (ln(gamma)) since otherwise we often get a numerical overflow or NAN.         */
      /* instead of multiplying a very large number by a very small number, we calculate the logarithm and then exp() */
      /* original version by Claudia Emde */
      // (*number_dens)[is] = A * pow(radius[is], alpha) * exp(-B*radius[is]);

      /* new version by Bernhard Mayer, 31.8.2016 */
      (*number_dens)[is] = exp (lnA + alpha * log(radius[is]) - B*radius[is]);

      /*  if(verbose) */
      /*         fprintf(stderr," %3d %g %g \n", is, radius[is], (*number_dens)[is]);    */
    }
    break;
    
  case DIST_LOGNORMAL:
    
    if (verbose)
      fprintf(stderr, "log-normal distribution, alpha=%f, r_mod=%f \n", alpha, r_eff);
    
    /* specified r_eff is taken as r_mod for log normal distribution */
    for (is=0;is<n_dens;is++) {
      A=(log((radius)[is])-log(r_eff))/log(alpha); 
      (*number_dens)[is] = 1.0 /
      (sqrt(2.0*PI)*log(alpha) * radius[is]) * exp(-0.5*A*A); 
    }
    
    /* old implementation:  */
    /*   CE: this was not consistent with the commonly used log-normal size distribution as for  */
    /*       instance in OPAC (Hess et al. 1998). Now made implementation consistent with OPAC. */
    /* B = r_eff*exp(-2.5*alpha*alpha); */
    /*     A = 1000./((4./3.*PI)* sqrt(2*PI)*alpha * B*B*B *exp(4.5*alpha*alpha)); */
    
    /*     A = A*LWC/rho_medium; */
    
    /* calculate number density */
    /* for (is=0;is<n_dens;is++){ */
    /*       (*number_dens)[is] = A/radius[is] *  */
    /*         exp(-0.5* pow(log(radius[is]/B),2) / pow(alpha,2)); */
    /*       /\* if(verbose ) *\/ */
    /*         fprintf(stderr," %3d %g %g \n", is, radius[is], (*number_dens)[is]);    */
    /*     } */
    break;
  default:
    fprintf (stderr, "Error in make_distribution : unrecognized distribution %d\n", distribution);
    status=-1;
    break;
  } /* switch(distribution) */
  
  return status;
}


double gamma_ln(float xx)
{
  int j=0;
  double *cof;
  double stp = 2.50662827465;
  double half = 0.5, one=1.0, fpf = 5.5, x, tmp, ser;
  double gammaln = 0.0;

  cof  = (double *) calloc (6, sizeof(double));
  cof[0] =  76.18009173;
  cof[1] = -86.50532033;
  cof[2] =  24.01409822;
  cof[3] =  -1.231739516;
  cof[4] =   0.120858003;
  cof[5] =  -0.53638200;

  x = xx - one;
  tmp = x + fpf;
  tmp = (x+half)*log(tmp)-tmp;
  ser=one;
  for (j=0; j<=5; j++) {
    x = x + one;
    ser = ser + cof[j]/x; 
  }
  gammaln=tmp+log(stp*ser);
 
  free(cof); 
  
  return gammaln;
}


double double_min(double x, double y)
{
  if (x <= y)
    return x;
  else
    return y;
}

double double_max(double x, double y)
{
  if (x >= y)
    return x;
  else
    return y;
}


/**
 * mie_calc_all_radii
 * 
 * Perform Mie calculations for all radii of all size
 * distributions. The resulting optical properties are afterwards
 * integrated over the size distributions. This function needs to be
 * called only once for each wavelength and optical properties for all
 * specified effective radii can be obtained by integrating the
 * resulting optical properties arrays over the various size
 * distributions. 
 *
 * @author ...
 * @date   2010-02-11 restructured by C. Emde
 */
int mie_calc_all_radii(/* Input */
                       int program, int medium, mie_complex crefin, 
                       double wavelength, float temperature, int nstokes,
                       double *x_size, int n_size, 
                       mie_inp_struct input, 
                       int verbose,
                       /* Output */
                       mie_out_struct *output, mie_complex *ref,
                       float **qext_array, float **qsca_array, float **gsca_array,
                       float **qback_array, float ****pmom_array)
{
  int status=0;
  
  int is=0, ip=0, im=0; 

  for (is=0; is<n_size; is++) {
    
      if (is % 100 == 0 && verbose)
        fprintf (stderr, "%d ", is);
      
      status = mie_calc (input, output, program, medium, crefin, 
	  	         wavelength, x_size[is], temperature, nstokes, ref);
      if (status!=0) {
        fprintf (stderr, "Error %d returned by mie_calc()\n", status);
        return status;
      }
      
      /* save results for further calculations */
      (*qext_array)[is] = output->qext;
      (*qsca_array)[is] = output->qsca;
      (*gsca_array)[is] = output->gsca;
      (*qback_array)[is] = output->qback;

      if (!(output->qext > 0) && verbose){
        fprintf(stderr, "Warning, something strange happened in Mie calculation: \n");
        fprintf(stderr, "is %d wavelength %g x_size %g qext %g \n", is,  wavelength,  x_size[is], output->qext);
      }
      
      if (input.nmom>0)
        for (ip=0; ip<nstokes; ip++)
          for (im=0; im<=input.nmom; im++)
	    (*pmom_array)[ip][im][is] = output->pmom[ip][im];
  }
  if(verbose)
    fprintf (stderr, "\n");
  
  return status; 
} 
  
/**
 * integrate_sizedist
 * 
 * Integrate optical properties arrays (output of mie_calc_all_radii)
 * over a size distribution. 
 * 
 * @author ...
 * @date   2010-02-11 restructured by C.Emde
 */
int integrate_sizedist(/* Input */
                       float* qext_array, 
                       float* qsca_array,
                       float* gsca_array,
                       float* qback_array,
                       float*** pmom_array,
                       int nstokes, int nmom,
                       double *x_size, double *y_size, int n_size, 
                       int ir, int verbose, 
                       /* Output */
                       double *beta,
                       double *omega, 
                       double *g,
                       double *back,
                       float*** pmom
                       )
{
                     
  int is=0, ip=0, im=0;
  
  double xsquared=0, xcubed=0, norm=0;
  double ext=0, sca=0, asy=0, vol=0, backsca;
  
  /* allocate memory for fields */
  double *extinction = calloc (n_size, sizeof (double));
  double *scattering = calloc (n_size, sizeof (double));
  double *asymmetry  = calloc (n_size, sizeof (double));
  double *volume     = calloc (n_size, sizeof (double));
  double *area       = calloc (n_size, sizeof (double));
  double *backscattering = calloc (n_size, sizeof (double));

  double ***pmom_tmp=NULL;
  
  if (nmom>0) {
    
    pmom_tmp = calloc (nstokes, sizeof(double **));
    for (ip=0; ip<nstokes; ip++) {
      pmom_tmp[ip] = calloc (nmom+1, sizeof(double *));
      for (im=0; im<=nmom; im++)
        pmom_tmp[ip][im] = calloc (n_size, sizeof(double));
    }
  }
  
  if(verbose)
    fprintf(stderr, "Integration over size distribution %d\n", ir);

  for (is=0; is<n_size; is++) {
    xsquared = x_size[is] * x_size[is] * 1e-12;
    xcubed   = xsquared   * x_size[is] * 1e-6;
    
    extinction[is] = qext_array[is] * xsquared * y_size[is]; 
    scattering[is] = qsca_array[is] * xsquared * y_size[is];
    asymmetry [is] = gsca_array[is] * scattering[is];
    backscattering[is] = qback_array[is] * xsquared * y_size[is]; 
    volume    [is] = xcubed * y_size[is];
    area      [is] = xsquared * y_size[is];

    if (nmom>0)
      for (ip=0; ip<nstokes; ip++)
        for (im=0; im<=nmom; im++)
	  pmom_tmp[ip][im][is] = pmom_array[ip][im][is] * y_size[is];
  }
  
  norm = integrate (x_size, y_size,     n_size);
  
  ext  = integrate (x_size, extinction, n_size) / norm;
  sca  = integrate (x_size, scattering, n_size) / norm;
  asy  = integrate (x_size, asymmetry,  n_size) / norm;
  backsca = integrate (x_size, backscattering, n_size) / norm;
  vol  = integrate (x_size, volume,     n_size) / norm;
  //20120816ak: are is not in use at the moment, thus commented.
  //  are  = integrate (x_size, area,       n_size) / norm;

  if (nmom>0)
    {
      for (ip=0; ip<nstokes; ip++)
        for (im=0; im<=nmom; im++)
          (*pmom)[ip][im] = integrate (x_size, pmom_tmp[ip][im], n_size) / norm;
      
      if(nstokes==4){
        /* Copy moments */
        for (ip=0; ip<nstokes; ip++) 
          for (im=0; im<=nmom; im++)
            pmom_tmp[ip][im][0] = (*pmom)[ip][im];
        
        /*  Convert from Wiscombe definition to Evans definition of Mueller matrix */
        norm=pmom_tmp[0][0][0]+pmom_tmp[1][0][0];
        for (im=0; im<=nmom; im++){ 
          (*pmom)[0][im]=(pmom_tmp[0][im][0]+pmom_tmp[1][im][0])/norm ;        
          (*pmom)[1][im]=(pmom_tmp[1][im][0]-pmom_tmp[0][im][0])/norm;  
          (*pmom)[2][im]= 2.*pmom_tmp[2][im][0]/norm;        
          (*pmom)[3][im]= 2.*pmom_tmp[3][im][0]/norm;      
        }
      }
      else if(nstokes==1){
        norm=(*pmom)[0][0];
        for (im=0; im<=nmom; im++)
          (*pmom)[0][im]=(*pmom)[0][im]/norm;
      }
      else{
        fprintf(stderr, "Error, nstokes must be either 1 or 4\n");
        return -1;
      }
    }

  /* averaged extinction efficiency qext, should be included as alternative option: */
  /* *beta  = ext/are; */

  /* extinction coefficient */
  *beta  = 3.0*ext/(4.0*vol)/1000.0;
  /* Single scattering albedo */
  *omega = sca/ext;
  /* Asymmetry parameter */
  *g     = asy/sca;
  /* Backscattering coefficient */
  *back  =  3.0*backsca/(4.0*vol)/1000.0;

  /* free memory */
  free (extinction);
  free (scattering);
  free (asymmetry);
  free (volume);

  if (nmom>0) {
    for (ip=0; ip<nstokes; ip++) {
      for (im=0; im<=nmom; im++)
	free(pmom_tmp[ip][im]);
      free(pmom_tmp[ip]);
    }    
    free(pmom_tmp);
  }
  
  return 0;
}
    

/** 
 * aerosol_input_data
 *
 * Reads aeosol input data from OPAC datafile 
 * data/aerosol/OPAC/size_distr.cfg. 
 * The data includes size distribution parameters 
 * and the density of the medium. 
 * 
 * @author Claudia Emde
 */
int aerosol_input_data(/* Input */
                       char* data_files_path, int aerosol_type, 
                       int verbose, 
                       /* Output */
                       aerdist_struct* aerdistr)
{

  int i=0, ih=0;
  int status; 
  int n_dist, n_hum; 
  double *aerosol=NULL;
  double *rel_hum=NULL; 
  double *r_min=NULL;
  double *r_max=NULL; 
  double *rmod=NULL; 
  double *rho=NULL; 
  double *sigma=NULL; 
  char *aerdistr_filename=NULL; 

  aerdistr_filename= calloc (FILENAME_MAX,sizeof(char));
  
  strcpy(aerdistr_filename, data_files_path);
  strcat(aerdistr_filename, "aerosol/OPAC/size_distr.cfg");
  
  /* allocate memory for size distribution data */
  n_dist = 39; /* number of distributions in datafile */
   
  n_hum=0;
  status = read_7c_file(aerdistr_filename, 
                        &aerosol, &rel_hum, &r_min, &r_max, &rmod, &rho, &sigma, &n_dist);
  if (status != 0) {
    fprintf (stderr, "error %d reading %s\n", status, aerdistr_filename );
    return status;
  }
  
  /* Check whether aerosol type is water soluble or dry. In the first case size distributions
     for 8 relative humidities are needed; in the second case only 1. */
  if(aerosol_type == TYPE_INSO || aerosol_type == TYPE_MIAM || 
     aerosol_type == TYPE_MICM || aerosol_type == TYPE_MINM || 
     aerosol_type == TYPE_MITR || aerosol_type == TYPE_SOOT){
    n_hum= 1;
  }
  else{
    n_hum = 8; 
  }
  
  /* allocate aerdistr_struct */
  (*aerdistr).n_hum = n_hum;
  (*aerdistr).rel_hum = calloc(n_hum, sizeof(double));
  (*aerdistr).r_min = calloc(n_hum, sizeof(double));
  (*aerdistr).r_max = calloc(n_hum, sizeof(double));
  (*aerdistr).rmod = calloc(n_hum, sizeof(double));
  (*aerdistr).rho = calloc(n_hum, sizeof(double));
  
  (*aerdistr).type = aerosol_type;
  while(i<n_dist){
    if((int)aerosol[i] == aerosol_type) {
      (*aerdistr).rel_hum[ih]=rel_hum[i];
      (*aerdistr).r_min[ih]=r_min[i]; 
      (*aerdistr).r_max[ih]=r_max[i];
      (*aerdistr).rmod[ih]=rmod[i];
      (*aerdistr).rho[ih]=rho[i];
    if(ih==0)
      (*aerdistr).sigma=sigma[i];
    ih++;
    }
    i++;
  }
  
  /* Check, whether all relative humidities have been  read. */
  if (ih != (*aerdistr).n_hum ){
    fprintf(stderr, "Something went wrong while reading aerosol input.ih %d nhum %d \n", ih, (*aerdistr).n_hum); 
    return -1;
  }
  
  free(aerdistr_filename);
  free(aerosol); free(rel_hum); free(r_min); free(r_max);
  free(rmod); free(rho); free(sigma);
  return 0; 
}


/** 
 * aerosol_size_distribution
 *
 * Aerosol size distributions are generated using the OPAC parameters. The 
 * distributions are log-normal distributions. The radius grid is generated 
 * using the method proposed by Evans (cloudprp.f), which seems to be an efficient 
 * and accurate method.
 *
 * @author Claudia Emde
 */
int aerosol_size_distribution(/* Input */
                              double rmin,
                              double rmax,
                              double rmod,
                              double sigma,
                              double lambda,
                              int verbose, 
                              /* Output */
                              int* ndens,
                              double** radius,
                              double** number_dens)
{

  int i=1; 
  double rstep=0.0, a=0.0, xstep=0.0, x=0.0; 
  double rad=0.0, nd=0.0;
  
  /* Calculate number of radius grid points */
  rad = rmin;
  /* The next points are determined iteratively */
  while (rad <= rmax &&
         /* Cutoff for small particle number densities. */
         /* Assure that calculation is not stopped too early.*/ 
         (nd > 1e-12 || rad < rmod ) ) { 
    
    /* Size parameter */
    x = 2*PI* rad /(lambda*1.e-03);  /* 1E-3 = nm -> micro m */ 
    xstep = double_max(0.03,0.03*sqrt(x)); 
    /* This is the same rule as proposed by Evans in cloudprp.f */
    rstep = xstep * lambda*1.e-03 / (2.0*PI);
    rad += rstep;
    a=(log(rad)-log(rmod))/log(sigma); 
    nd = 1.0 / (sqrt(2.0*PI)*log(sigma)*log(10.0)* rad) * exp(-0.5*a*a); 
    i++;
  } 
  *ndens=i;
  
  /* Allocate memory */
  *radius = calloc(*ndens, sizeof(double));      
  *number_dens = calloc(*ndens, sizeof(double)); 
  
  /* Calculate size distribution, same as above, but now we may fill
     the arrays */
  (*radius)[0]=rmin;
  a=(log((*radius)[0])-log(rmod))/log(sigma);
  (*number_dens)[0]=1.0/(sqrt(2.0*PI)*log(sigma)/log(10))*exp(-0.5*a*a); 
 
  if(verbose) 
  /*     fprintf(stderr, "aerosol size distribution: \n");  */
    fprintf(stderr, " ... %d radius grid points in aerosol size distribution \n", *ndens);

  i=1;  
  /* The next points are determined iteratively */
  while (i < *ndens){
    /* Size parameter */
    x = 2*PI* (*radius)[i-1] / (lambda*1.e-03);  /* 1E-3 = nm -> micro m */ 
    xstep = double_max(0.03,0.03*sqrt(x)); 
    /* This is the same rule as proposed by Evans in cloudprp.f */
    rstep = xstep * lambda*1.e-03 / (2.0*PI);
    (*radius)[i]= (*radius)[i-1] + rstep;
    a=(log((*radius)[i])-log(rmod))/log(sigma); 
    (*number_dens)[i] = 1.0 /
      (sqrt(2.0*PI)*log(sigma)*log(10.0)* (*radius)[i]) * exp(-0.5*a*a); 
    
    /* if(verbose) */
    /*       fprintf(stderr, "i %d %g   %g \n", i, (*radius)[i], (*number_dens)[i]); */

    i++;
  } 
  
  
  /* An alternative would be to use a logarithmic radius. This is done in OPAC.f. */
  /* The solution above should be more accurate */
  /*  int n;  */
  /*   double rlog; */
  /*   n=1000;  */
  
  /*   rstep=(log(rmax)-log(rmin))/(double)n;  */
  
  /*   (*radius)[0]=rmin;   */
  /*   rlog=log(rmin);   */
  /*   i=1;    */
  /*   while (i<n+1){   */
  /*     rlog+=rstep;   */
  /*     /\* radius grid *\/  */
  /*     (*radius)[i]=exp(rlog);   */
  /*     /\* particle number densities, total number density is set to 1.*\/   */
  /*     a=(log((*radius[i])-log(rmod)))/  */
  /*       log(sigma);   */
  /*     (*number_dens)[i]=1/(sqrt(2*PI)*log(sigma)*log(10))*  */
  /*       exp(-0.5*a*a);    */
  /*     i++;   */
  /*   }   */
  /*   n_dens=n+1;   */
  
  return 0; 
}

/** 
 * aerosol_refractive_index
 *
 * This function reads the refractive index from the OPAC data 
 * in data/aerosol/OPAC for the requested aerosol type and 
 * interpolates it on the required wavelength.
 * 
 * 
 * @author Claudia Emde
 */
int aerosol_refractive_index(/* Input */
                             double* lambda,
                             int nlambda, 
                             char* data_files_path, 
                             int aerosol_type,
                             int n_hum,
                             int verbose, 
                             /* Output */
                             double*** ref_real,
                             double*** ref_imag)
  
{
  int status=0;
  double* lambda_data=NULL; 
  double* ref_real_data=NULL; 
  double* ref_imag_data=NULL;
  int nlambda_data=0, ih=0, iv=0;
  char* type;
  char* rel_hum;
  char* filename;
  
  /* Coefficents for interpolation */
  double *a0=NULL, *a1=NULL, *a2=NULL, *a3=NULL;
  
  if (strlen(data_files_path) == 0)
    data_files_path="../data/"; 
  
  filename= calloc (FILENAME_MAX,sizeof(char));
  
  if(aerosol_type == TYPE_INSO)
    type="inso"; 
  else if (aerosol_type == TYPE_WASO)
    type="waso";
  else if (aerosol_type == TYPE_SOOT)
    type="soot";
  else if (aerosol_type == TYPE_SSAM)
    type="ssam";
  else if (aerosol_type == TYPE_SSCM)
    type="sscm";
  else if (aerosol_type == TYPE_MINM)
    type="minm";
  else if (aerosol_type == TYPE_MIAM)
    type="miam";
  else if (aerosol_type == TYPE_MICM)
    type="micm";
  else if (aerosol_type == TYPE_MITR)
    type="mitr";
  else if (aerosol_type == TYPE_SUSO)
    type="suso";
  else{
    fprintf(stderr, "Error: Unknown aerosol type \n");
    return -1; 
  }
  
  
  /* Allocate refractive index array*/
  (*ref_real) = calloc(n_hum, sizeof(double *)); 
  (*ref_imag) = calloc(n_hum, sizeof(double *)); 
  for (ih=0; ih<n_hum; ih++){
    (*ref_real)[ih] = calloc(nlambda, sizeof(double)); 
    (*ref_imag)[ih] = calloc(nlambda, sizeof(double)); 
  }
  
  for (ih=0; ih<n_hum; ih++){
    
    if (ih == 0) rel_hum="00";
    else if (ih == 1) rel_hum="50";
    else if (ih == 2) rel_hum="70";
    else if (ih == 3) rel_hum="80";
    else if (ih == 4) rel_hum="90";
    else if (ih == 5) rel_hum="95";
    else if (ih == 6) rel_hum="98";
    else if (ih == 7) rel_hum="99";
    else{
      fprintf(stderr, "Error: Humidity index ih=%d is wrong. \n", ih);
      return -1; 
    } 
    
    /* Create datafile name */
    strcpy(filename, data_files_path);
    strcat(filename,"aerosol/OPAC/refractive_indices/");
    strcat(filename, type);
    strcat(filename, rel_hum);
    strcat(filename, "_refr.dat");


    status = read_3c_file (filename, &lambda_data, &(ref_real_data), &(ref_imag_data), &nlambda_data);
    if (status != 0) {
      fprintf (stderr, "error %d reading %s\n", status, filename);
      return status;
    }
    
    /* for (iv=0; iv<nlambda_data; iv++) */
    /*       fprintf(stderr, "%g %g %g \n", lambda_data[iv], ref_real_data[iv], ref_imag_data[iv]); */
 
    /* Interpolate refractive index on the required wavelength. */
    /* real part*/
    status = linear_coeffc (lambda_data, ref_real_data, nlambda_data, &a0, &a1, &a2, &a3);
    for (iv=0; iv<nlambda; iv++){
      status = calc_splined_value (lambda[iv]*1e-3, &(*ref_real)[ih][iv], lambda_data, nlambda_data,
                                   a0, a1, a2, a3);
      if (status!=0)  {
        fprintf (stderr, "Sorry cannot interpolate refractive index (real)\n");
        return status;
      }
    }
    free(a0); free(a1); free(a2); free(a3); 
    
    
    /* imaginary part */
    status = linear_coeffc (lambda_data, ref_imag_data, nlambda_data, &a0, &a1, &a2, &a3);
    for (iv=0; iv<nlambda; iv++){
      status = calc_splined_value (lambda[iv]*1e-3, &(*ref_imag)[ih][iv], lambda_data, nlambda_data,
                                   a0, a1, a2, a3);
      if (status!=0)  {
        fprintf (stderr, "Sorry cannot interpolate refractive index (imaginary)\n");
        return status;
      }
      (*ref_imag)[ih][iv]*=-1; /* Mie tool requires positive values*/
    }
    free(a0); free(a1); free(a2); free(a3); 

  } 
  
  /* Free variables */
  
  free(filename);
  free(lambda_data); 
  
  free(ref_real_data);
  free(ref_imag_data);
    
  return 0;
}

/**
 * write_output
 * 
 * Write output to the screen or into output files. The output file
 * format corresponds to the output of the cloudprp program by
 * F. Evans. The screen output can be specified by the user. 
 * 
 * @author ???, C. Emde (phase matrix output) 
 * @date ???
 */

int write_output (input_struct input, output_struct *output, int ir, int iv)
{

  int  io=0;
  int  nmom=-999, nstokes=1;
  int  im=0;
  char outputfile[255];
  char ck_band_str[4];
  int  status=0;
  int  ip=0, nphamat=0;
  int default_output=0;

  if (input.output_format == OUTPUT_FORMAT_NORMAL) {
    
    if (input.n_output_user == 0) {  
      /* default output */
      input.n_output_user = 8;
      input.output_user = calloc (input.n_output_user, sizeof (int));
      input.output_user[0]=OUTPUT_USER_WAVE;
      input.output_user[1]=OUTPUT_USER_REFRAC_REAL;
      input.output_user[2]=OUTPUT_USER_REFRAC_IMAG;
      input.output_user[3]=OUTPUT_USER_QEXT;
      input.output_user[4]=OUTPUT_USER_OMEGA;
      input.output_user[5]=OUTPUT_USER_GG;
      input.output_user[6]=OUTPUT_USER_SPIKE;
      input.output_user[7]=OUTPUT_USER_PMOM;
      default_output=1;
    }

    for (io=0; io<input.n_output_user; io++) {
      switch (input.output_user[io]) {
      case OUTPUT_USER_NONE:
        fprintf (stdout, "%9s ", " ");
        break;
      case OUTPUT_USER_WAVE:
        fprintf (stdout, "%9.3f ", output->wl.lambda[iv]);
        break;
      case OUTPUT_USER_WAVENUMBER:
        fprintf (stdout, "%9.3f ", 1E7/output->wl.lambda[iv]);
        break;
      case OUTPUT_USER_R_EFF:
        fprintf (stdout, "%12.6f ", output->r_eff[ir]);
	  break;
      case OUTPUT_USER_REFRAC_REAL:
	  fprintf (stdout, "%12.6e ", output->ref.re);
	  break;
      case OUTPUT_USER_REFRAC_IMAG:
        fprintf (stdout, "%12.6e ", output->ref.im);
        break;
      case OUTPUT_USER_QEXT:
        if (input.distribution == DIST_NONE)
          fprintf (stdout, "%12.6e ", output->mie.qext);
        else
          fprintf (stdout, "%12.6e ", output->beta);
        break;
      case OUTPUT_USER_QSCA:
        if (input.distribution == DIST_NONE)
          fprintf (stdout, "%12.6e ", output->mie.qsca);
        else
          fprintf (stdout, "%12.6e ", output->omega*output->beta);
        break;
      case OUTPUT_USER_OMEGA:
        if (input.distribution == DIST_NONE)
          fprintf (stdout, "%12.6e ", output->mie.qsca/output->mie.qext);
        else
          fprintf (stdout, "%12.6e ", output->omega);
        break;
      case OUTPUT_USER_GG:
        if (input.distribution == DIST_NONE)
          fprintf (stdout, "%12.6e ", output->mie.gsca);
        else
	    fprintf (stdout, "%12.6e ", output->g);
        break;
      case OUTPUT_USER_SPIKE:
        if (input.distribution == DIST_NONE )
          fprintf (stdout, "%12.6e ", output->mie.spike);
        else
          fprintf (stdout, "%12.6e ", -1.0);
        break;
      case OUTPUT_USER_SFORW:
	if (input.distribution == DIST_NONE)
          fprintf (stdout, "%12.6e %12.6e ", output->mie.sforw.re, output->mie.sforw.im);
        else
          fprintf (stdout, "%12.6e %12.6e ", NAN, NAN);  /* not yet implemented*/
        break;
      case OUTPUT_USER_SBACK:
	if (input.distribution == DIST_NONE)
          fprintf (stdout, "%12.6e %12.6e ", output->mie.sback.re, output->mie.sback.im);
        else
          fprintf (stdout, "%12.6e %12.6e ", NAN, NAN);  /* not yet implemented*/
        break;
      case OUTPUT_USER_QBACK:
        if (input.program == BHMIE){
          if (input.distribution == DIST_NONE)
            fprintf (stdout, "%12.6e ", output->mie.qback);
          else
            fprintf (stdout, "%12.6e ", output->back);
        }
        else
          fprintf (stdout, "%12.6e ", NAN); /* not yet implemented */
        break;
      case OUTPUT_USER_PMOM:
        if (input.mie.nmom>0){
          if (input.nstokes>1 || default_output==0){
            fprintf(stdout , "\n");
            for (im=0; im<=input.mie.nmom; im++){
              for (ip=0; ip < input.nstokes; ip++)
                fprintf (stdout, " %14.6e", output->pmom[ip][im]);
              fprintf(stdout , "\n");
            }
          }
          else{
            for (im=0; im<=input.mie.nmom; im++)
              fprintf (stdout, " %14.6e", output->pmom[0][im]);
            fprintf(stdout , "\n");
          }
        }
        break;
      default:
        fprintf (stderr, "Error, unknown user output %d\n", input.output_user[io]);
        break;
      } /* switch (input.output_user[io]) { */
    } /* for (io=0; io<input.n_output_user; io++) { */
    fprintf (stdout, "\n");
    fflush (stdout);
  } /* if (input.output_user[0] != OUTPUT_USER_CLOUDPROP) { */

  else if (input.output_format == OUTPUT_FORMAT_CLOUDPROP || 
        input.output_format == OUTPUT_FORMAT_AERPROP   ) {
    
    if ( input.distribution == DIST_NONE ){
      fprintf(stderr, "The option output_user cloudprop works only for size distributions, not for mono-disperse particles. \n");
      return -1;
    }
      
    /* For polarized calculations, all four elements are calculated for mie files */
    if (input.nstokes > 1)
      nstokes =4; 
    
    if (ir == 0) {
      /* generate file name */
      switch (input.medium) {
      case WATER:
        strcpy (outputfile, "wc.");          
	break;
      case ICE:
        strcpy (outputfile, "ic.");   
        break;
      case USER:
        strcpy (outputfile, "uc.");   
	break;
      case AEROSOL:
        if(input.aerosol_type == TYPE_INSO)
          strcpy (outputfile, "inso."); 
        else if (input.aerosol_type == TYPE_WASO)
          strcpy (outputfile, "waso."); 
        else if (input.aerosol_type == TYPE_SOOT)
          strcpy (outputfile, "soot."); 
        else if (input.aerosol_type == TYPE_SSAM)
          strcpy (outputfile, "ssam."); 
        else if (input.aerosol_type == TYPE_SSCM)
          strcpy (outputfile, "sscm.");
        else if (input.aerosol_type == TYPE_MINM)
          strcpy (outputfile, "minm.");
        else if (input.aerosol_type == TYPE_MIAM)
          strcpy (outputfile, "miam."); 
        else if (input.aerosol_type == TYPE_MICM)
          strcpy (outputfile, "micm."); 
        else if (input.aerosol_type == TYPE_MITR)
          strcpy (outputfile, "mitr.");
        else if (input.aerosol_type == TYPE_SUSO)
          strcpy (outputfile, "suso.");
        else{
          fprintf(stderr, "Error: Unknown aerosol type \n");
          return -1; 
        }
        break;
      }
      
      sprintf(ck_band_str, "%03d",iv+1);
      strcat (outputfile, ck_band_str);
      strcat (outputfile, ".mie");
      fprintf (stderr, " ... generating file: %s\n", outputfile);
          
      if ((output->ptr_output = fopen (outputfile, "w")) == NULL) {
        fprintf (stderr, "Error, cannot open file: %s in write_output! \n", outputfile);   
        return -1;         
      }
    }
    
    if(input.medium != AEROSOL){
      fprintf(stderr,"   lambda = %10.3f nm (%3d of %3d), r_eff = %8.3f micro meter (%3d of %3d) \n", 
              output->wl.lambda[iv], iv+1,output->wl.nlambda+1, output->r_eff[ir], ir+1,output->n_r_eff+1);
    }
    else{
      fprintf(stderr,"   lambda = %10.3f nm (%3d of %3d), rel_hum = %4.0f percent (%3d of %3d) \n", 
               output->wl.lambda[iv], iv+1,output->wl.nlambda+1, output->aerdistr.rel_hum[ir], ir+1, output->aerdistr.n_hum);
    }
      
    
    if (ir == 0) {
      switch (input.medium) {
      case WATER:
        fprintf (output->ptr_output, "! Mie table vs. effective radius (LWC=1 g/m^3)\n");
	break;
      case ICE:
        fprintf (output->ptr_output, "! Mie table vs. effective radius (IWC=1 g/m^3)\n");
	break;
      case USER:
        fprintf (output->ptr_output, "! Mie table vs. effective radius (ATTENTION!!! assume rho_medium = %7.4f g/cm3)\n", 
                                      output->rho_medium);
        break;
      case AEROSOL: 
        fprintf (output->ptr_output, "! Mie table vs. relative humidity (N=1, rho_medium = %7.4f g/cm3)\n", 
                 output->rho_medium);
	break;
      }
      
      fprintf (output->ptr_output," %12.6e  wavelength (micron) \n",output->wl.lambda[iv]/1000.);
      fprintf (output->ptr_output," %12.6e  %12.6e  index of refraction \n", output->ref.re, output->ref.im);
      if (input.distribution == DIST_NONE) {
        fprintf (output->ptr_output,"S  distribution type ! one Single radius  \n");
        fprintf   (output->ptr_output,"-999.99  distribution shape parameter !(dummy to fit the format)\n");
      }
      else if (input.distribution == DIST_GAMMA) {
        fprintf (output->ptr_output,"G  distribution type\n");
        fprintf   (output->ptr_output,"%7.4f  distribution shape parameter\n", input.alpha);
      }
      else if (input.distribution == DIST_LOGNORMAL) {
        fprintf (output->ptr_output,"L  distribution type\n");
        fprintf   (output->ptr_output,"%7.4f  distribution shape parameter\n", input.alpha);
      }
      else if (input.distribution == DIST_FILE) {
        fprintf (output->ptr_output,"U  distribution type ! User defined, from file \n");
        fprintf   (output->ptr_output,"-999.99  distribution shape parameter !(dummy to fit the format)\n");
      }
      else if (input.distribution == DIST_AER) {
        fprintf (output->ptr_output,"A  distribution type \n");
        fprintf (output->ptr_output,"Size distribution parameters from OPAC database \n");
      }

      if (input.medium != AEROSOL){
        fprintf (output->ptr_output,"%4d %8.4f %8.4f  number, starting, ending effective radius \n", 
                 output->n_r_eff+1, output->r_eff[0], output->r_eff[output->n_r_eff]); 
      }
      else{
        fprintf (output->ptr_output,"%4d %8.4f %8.4f  number, starting, ending relative humidity \n", 
                 output->aerdistr.n_hum, output->aerdistr.rel_hum[0], output->aerdistr.rel_hum[output->aerdistr.n_hum-1]); 
      }
      
      
      } /* if (ir == 0) { */

    
    nmom=0;
    for (ip=0; ip<nstokes; ip++){
      for (im=0;im<=input.mie.nmom;im++){
        /* if (fabs((2*im+1)*output->pmom[ip][im]/output->pmom[ip][0]) < 0.000005) { */
        if (fabs((2*im+1)*output->pmom[ip][im]) < 1E-5){
          if ( im > nmom )
            nmom=im;
          break;
        }
      }
    }
    
    /* remember to write a warning, if there could be more Legendre coefficients */
    if (nmom == input.mie.nmom) 
      output->warnings_nmom[iv][ir]=1;
    
    if(nstokes==1) nphamat = 1;
    else nphamat = 4;
    
    if (input.medium != AEROSOL){
    fprintf (output->ptr_output,"%9.4f  %12.5e  %12.5e %6d  %4d Reff  Ext  Alb  Nleg Nphamat\n", output->r_eff[ir], 
             output->beta/output->rho_medium, 
             output->omega, nmom, nphamat);
    }
    else{
      fprintf (output->ptr_output,"%9.4f  %12.5e  %12.5e %6d  %4d RelHum  Ext  Alb  Nleg Nphamat\n", output->aerdistr.rel_hum[ir], 
               output->beta/output->rho_medium, 
               output->omega, nmom, nphamat);
    }
    
    if (input.mie.nmom>0) {
      fprintf (output->ptr_output, "  ");
      for (ip=0; ip < input.nstokes; ip++){
        for (im=0; im<=nmom; im++)
          fprintf (output->ptr_output, " %9.5f", (2*im+1) * output->pmom[ip][im]);
        fprintf( output->ptr_output, "\n");
      }
    }
    fprintf (output->ptr_output, "\n");
    fflush (output->ptr_output);  
   
    if ((ir == output->n_r_eff && input.output_format == OUTPUT_FORMAT_CLOUDPROP) ||
        (ir == output->aerdistr.n_hum && input.output_format == OUTPUT_FORMAT_AERPROP) )
      fclose(output->ptr_output);
    
  } /* else if (input.output_format == OUTPUT_FORMAT_CLOUDPROP || OUTPUT_FORMAT_AERPROP) ... */ 
  
  return status;
}
 
 
#if HAVE_LIBNETCDF

/** 
 * create_netcdf_file 
 *
 * Create netcdf file.
 * This can directly be used as input for uvspec. 
 * 
 *
 * @author  Claudia Emde
 * @date    2010-01-09 Created
 */         
int create_netcdf_file(/* Input */
                       char *basename,
                       int medium, 
                       int aerosol_type,
                       int distribution, 
                       int alpha,
                       int nlam,
                       double* lambda,
                       int nreff, 
                       double* reff, 
                       int nmommax, 
                       int nstokes,
                       int nthetamax,
                       int nrho,
                       int complete_moments,
                       double n_r_max,
                       double dx_max,
                       char *sd_filename,
                       /* Output */
                       netcdf_ids* ncid,
                       char** outputfile
                       )
{
  /* Indices*/
  int retval=0, status=0, i=0;
  
  /* Netcdf dimension IDs */
  int dimids2D[2], dimids3D[3], dimids4D[4];
  
  float fill_value[1];
  char name[MAXLINE], unit[23], varname[10];
  
  char medium_str[MAXLINE];
  char distr[MAXLINE], param[MAXLINE]; 
  double distr_param[1]; 
  long version[1];
  int comp_moments[1]; 
  
  int nphamat=nstokes; 
  double lambda_tmp[nlam];
  
  /*-------------------------------------------------------------------*/
  /* Generate filename                                                 */
  /*-------------------------------------------------------------------*/

  if(strlen(basename) != 0){
    strcpy (*outputfile, basename);
  }
  else{
    switch (medium) {
    case WATER:
      strcpy (*outputfile, "wc.");
      strcpy (medium_str, "water");
      break;
    case ICE:
    strcpy (*outputfile, "ic.");   
    strcpy (medium_str, "ice");
    break;
    case USER:
      strcpy (*outputfile, "uc.");
      strcpy (medium_str, "user defined refractive index");
      break;
    case AEROSOL:
      
      if(aerosol_type == TYPE_INSO){
        strcpy (*outputfile, "inso."); 
        strcpy (medium_str, "Aerosol: insoluble"); 
      }
      else if (aerosol_type == TYPE_WASO){
        strcpy (*outputfile, "waso."); 
        strcpy (medium_str, "Aerosol: water soluble");
      }
      else if (aerosol_type == TYPE_SOOT){
        strcpy (*outputfile, "soot."); 
        strcpy (medium_str, "Aerosol: soot");
      }
      else if (aerosol_type == TYPE_SSAM){
        strcpy (*outputfile, "ssam."); 
        strcpy (medium_str, "Aerosol: sea salt accumulated mode");
      }    
      else if (aerosol_type == TYPE_SSCM){
        strcpy (*outputfile, "sscm.");
        strcpy (medium_str, "Aerosol: sea salt coarse mode");
      }  
      else if (aerosol_type == TYPE_MINM){
        strcpy (*outputfile, "minm.");
        strcpy (medium_str, "Aerosol: mineral nucleation mode");
      }
      else if (aerosol_type == TYPE_MIAM){
        strcpy (*outputfile, "miam.");
        strcpy (medium_str, "Aerosol: mineral accumulated mode");
      }
      else if (aerosol_type == TYPE_MICM){
        strcpy (*outputfile, "micm."); 
        strcpy (medium_str, "Aerosol: mineral coarse mode");
      }
      else if (aerosol_type == TYPE_MITR){
        strcpy (*outputfile, "mitr.");
        strcpy (medium_str, "Aerosol: mineral transported");
      } 
      else if (aerosol_type == TYPE_SUSO){
        strcpy (*outputfile, "suso.");
        strcpy (medium_str, "Aerosol: sulfate");
      }
      else{
        fprintf(stderr, "Error: Unknown aerosol type \n");
        return -1; 
      }
      break;
    default:
      fprintf (stderr, "Error in write_output_netcdf : unrecognized medium %d\n", medium); 
      return -1; 
    }
  }
  
  switch(distribution){
  case DIST_NONE:
    strcpy (distr, "Calculation for single radius.");
    distr_param[0]=-999.;
    break;
  case DIST_GAMMA:
    strcpy (distr, "Gamma distribution.");
    distr_param[0]=alpha;
    break; 
  case DIST_LOGNORMAL:
    strcpy (distr, "Lognormal distribution.");
    distr_param[0]=alpha;
    break;
  case DIST_FILE:
    strcpy (distr, "User defined size distribution.");
    distr_param[0]=-999.;
    break; 
  case DIST_AER:
    strcpy (distr, 
            "Aerosol size distribution,\nlog-normal with parameters from OPAC");
    distr_param[0]=-999.;
    break;
  default: 
    fprintf (stderr, "Error in write_output_netcdf : unrecognized distribution %d\n", distribution); 
    return -1; 
    break; 
  }
  
  strcat (*outputfile, "mie.cdf");
  fprintf (stderr, " ... generating file: %s\n", *outputfile);
  
  if ((retval = nc_create(*outputfile, NC_CLOBBER, &(*ncid).ncid)))
    {
      fprintf(stderr, "Error: creating netcdf file");
      return -1;
    }

  /* ----------------------------------------------------------------*/
  /* Set global attributes                                           */
  /* ----------------------------------------------------------------*/

  version[0]=20110916;
  strcpy(param, "mie"); 
  comp_moments[0] = complete_moments; 

  strcpy(name, "Netcdf file created using libRadtran-Mie-tool");
  nc_put_att_text((*ncid).ncid, NC_GLOBAL, "file_info", strlen(name), name);  
  
  nc_put_att_long((*ncid).ncid, NC_GLOBAL, "version", NC_LONG, 1, version);  

  nc_put_att_int((*ncid).ncid, NC_GLOBAL, "complete_moments", NC_INT, 1, comp_moments); 

  nc_put_att_text((*ncid).ncid, NC_GLOBAL, "parameterization", strlen(param), param);  
  
  nc_put_att_text((*ncid).ncid, NC_GLOBAL, "size_distr", strlen(distr), distr);
  
  nc_put_att_double((*ncid).ncid, NC_GLOBAL, "param_alpha", NC_DOUBLE, 1, distr_param);
 
  nc_put_att_double((*ncid).ncid, NC_GLOBAL, "n_r_max", NC_DOUBLE, 1, &n_r_max);
  
  nc_put_att_double((*ncid).ncid, NC_GLOBAL, "dx_max", NC_DOUBLE, 1, &dx_max);

  if (strcmp(sd_filename,"") == 1)
    nc_put_att_text((*ncid).ncid, NC_GLOBAL, "size distribution file", strlen(sd_filename), sd_filename);
  
  /* ----------------------------------------------------------------*/
  /* Define dimensions.                                              */
  /* ----------------------------------------------------------------*/
  
  if ((retval = nc_def_dim((*ncid).ncid, "nlam", nlam, &(*ncid).id_nlam))){
    fprintf(stderr, "Error creating dimension nlam.\n" );
    return -1;
  }
  
  if ((retval = nc_def_dim((*ncid).ncid, "nmommax", nmommax, &(*ncid).id_nmommax))){
    fprintf(stderr, "Error creating dimension nmommax.\n" );
    return -1;
  }
  
  if ((retval = nc_def_dim((*ncid).ncid, "nphamat", nphamat, &(*ncid).id_nphamat))){
    fprintf(stderr, "Error creating dimension nphase.\n" );
  }
  
  if (medium==AEROSOL){
    if ((retval = nc_def_dim((*ncid).ncid, "nhum", nreff, &(*ncid).id_nreff))){
      fprintf(stderr, "Error creating dimension nreff.\n" );
      return -1;
    }
  }
  else{
    if((retval = nc_def_dim((*ncid).ncid, "nreff", nreff, &(*ncid).id_nreff))){
      fprintf(stderr, "Error creating dimension nreff.\n" );
      return -1;
    }
  }
  
  if ((retval = nc_def_dim((*ncid).ncid, "nthetamax", nthetamax, &(*ncid).id_nthetamax))){
    fprintf(stderr, "Error creating dimension nthetamax.\n" );
    return -1;
  }
  
  if(medium!=AEROSOL)
    nrho=1;
  
  if ((retval = nc_def_dim((*ncid).ncid, "nrho", nrho, &(*ncid).id_nrho))){
    fprintf(stderr, "Error creating dimension nrho.\n" );
    return -1;
  }
   
  
  /* ---------------------------------------------------------------*/
  /* Define variables                                               */ 
  /* ---------------------------------------------------------------*/
  
  
  /* Wavelength*/
  
  if ((retval = nc_def_var((*ncid).ncid, "wavelen", NC_DOUBLE, 1, 
                           &(*ncid).id_nlam, &(*ncid).id_lam))){
    fprintf(stderr, "Error: defining variable wavelen. \n");
    return -1;
  } 
  strcpy(name,"wavelength");
  strcpy(unit,"micrometer");
  nc_put_att_text((*ncid).ncid, (*ncid).id_lam, "long_name", strlen(name), name);  
  nc_put_att_text((*ncid).ncid, (*ncid).id_lam, "units", strlen(unit), unit);
  
  
  /* effective radius or relative humidity */
  
  if (medium == AEROSOL){
    strcpy(varname,"hum");
    strcpy(name,"relative humidity");
    strcpy(unit,"per cent");
  }
  else{
    strcpy(varname,"reff");
    strcpy(name,"effective radius");
    strcpy(unit,"micrometer");
  }
  if ((retval = nc_def_var((*ncid).ncid, varname, NC_DOUBLE, 1, 
                           &(*ncid).id_nreff, &(*ncid).id_reff))){
    fprintf(stderr, "Error defining variable reff. \n");
    return -1;
  } 
  nc_put_att_text((*ncid).ncid, (*ncid).id_reff, "long_name", strlen(name), name);  
  nc_put_att_text((*ncid).ncid, (*ncid).id_reff, "units", strlen(unit), unit);
  
  /* Number of scattering angles */
  
  dimids3D[0] = (*ncid).id_nlam;
  dimids3D[1] = (*ncid).id_nreff;
  dimids3D[2] = (*ncid).id_nphamat;
  
  if ((retval = nc_def_var((*ncid).ncid, "ntheta", NC_INT, 3, 
                           dimids3D, &(*ncid).id_ntheta))){
    fprintf(stderr, "Error defining variable ntheta. \n");
    return -1;
  } 
  strcpy(name,"number of scattering angles");
  strcpy(unit,"-");
  nc_put_att_text((*ncid).ncid, (*ncid).id_ntheta, "long_name", strlen(name), name);  
  nc_put_att_text((*ncid).ncid, (*ncid).id_ntheta, "units", strlen(unit), unit);
  
  /* Scattering angle grid */
  
  fill_value[0]=-999.0;
  
  /* The dimids array is used to pass the IDs of the dimensions of
   * the variable. */
  dimids4D[0] = (*ncid).id_nlam;
  dimids4D[1] = (*ncid).id_nreff;
  dimids4D[2] = (*ncid).id_nphamat;
  dimids4D[3] = (*ncid).id_nthetamax;
  
  if ((retval = nc_def_var((*ncid).ncid, "theta", NC_FLOAT, 4, 
                           dimids4D, &(*ncid).id_theta))){
    fprintf(stderr, "Error defining variable theta. \n");
    return -1;
  } 
  strcpy(name,"theta");
  strcpy(unit,"degrees");
  nc_put_att_text((*ncid).ncid, (*ncid).id_theta, "long_name", strlen(name), name);  
  nc_put_att_text((*ncid).ncid, (*ncid).id_theta, "units", strlen(unit), unit);
  nc_put_att_float((*ncid).ncid, (*ncid).id_theta, "_FillValue", NC_FLOAT, 1, fill_value);
  
  /* Phase matrix */

  if ((retval = nc_def_var((*ncid).ncid, "phase", NC_FLOAT, 4, 
                           dimids4D, &(*ncid).id_phase))){
    fprintf(stderr, "Error defining variable phase. \n");
    return -1;
  } 
  strcpy(name,"phase");
  strcpy(unit,"-");
  nc_put_att_text((*ncid).ncid, (*ncid).id_phase, "long_name", strlen(name), name);  
  nc_put_att_text((*ncid).ncid, (*ncid).id_phase, "units", strlen(unit), unit);
  nc_put_att_float((*ncid).ncid, (*ncid).id_phase, "_FillValue", NC_FLOAT, 1, fill_value);
  
  /* Number of Legendre polynomials */
  
  dimids3D[0] = (*ncid).id_nlam;
  dimids3D[1] = (*ncid).id_nreff;
  dimids3D[2] = (*ncid).id_nphamat;
  
  if ((retval = nc_def_var((*ncid).ncid, "nmom", NC_INT, 3, 
                           dimids3D, &(*ncid).id_nmom))){
    fprintf(stderr, "Error defining variable nmom. \n");
    return -1;
  } 
  strcpy(name,"number of Legendre polynomials");
  strcpy(unit,"-");
  nc_put_att_text((*ncid).ncid, (*ncid).id_nmom, "long_name", strlen(name), name);  
  nc_put_att_text((*ncid).ncid, (*ncid).id_nmom, "units", strlen(unit), unit);

  /* Legendre moments */
  fill_value[0]=0.0;
  
  dimids4D[0] = (*ncid).id_nlam;
  dimids4D[1] = (*ncid).id_nreff;
  dimids4D[2] = (*ncid).id_nphamat;
  dimids4D[3] = (*ncid).id_nmommax;

  if ((retval = nc_def_var((*ncid).ncid, "pmom", NC_FLOAT, 4, 
                           dimids4D, &(*ncid).id_pmom))){
    fprintf(stderr, "Error defining variable pmom. \n");
    return -1;
  } 
  strcpy(name,"Legendre polynomials");
  strcpy(unit,"including factor 2*l+1");
  nc_put_att_text((*ncid).ncid, (*ncid).id_pmom, "long_name", strlen(name), name);  
  nc_put_att_text((*ncid).ncid, (*ncid).id_pmom, "units", strlen(unit), unit);
  nc_put_att_float((*ncid).ncid, (*ncid).id_pmom, "_FillValue", NC_FLOAT, 1, fill_value);
  
  /* Extinction coefficient */
  
  dimids2D[0] = (*ncid).id_nlam;
  dimids2D[1] = (*ncid).id_nreff; 
  
  if ((retval = nc_def_var((*ncid).ncid, "ext", NC_DOUBLE, 2, 
                           dimids2D, &(*ncid).id_ext))){
    fprintf(stderr, "Error defining variable ext. \n");
    return -1;
  } 
  strcpy(name,"extinction coefficient");
  strcpy(unit,"km^-1/(g/m^3)");
  nc_put_att_text((*ncid).ncid, (*ncid).id_ext, "long_name", strlen(name), name);  
  nc_put_att_text((*ncid).ncid, (*ncid).id_ext, "units", strlen(unit), unit);
   
  
  /* Single scattering albedo */
  
  if ((retval = nc_def_var((*ncid).ncid, "ssa", NC_DOUBLE, 2, 
                           dimids2D, &(*ncid).id_ssa))){
    fprintf(stderr, "Error defining variable ssa. \n");
    return -1;
  } 
  strcpy(name,"single scattering albedo");
  strcpy(unit,"-");
  nc_put_att_text((*ncid).ncid, (*ncid).id_ssa, "long_name", strlen(name), name);  
  nc_put_att_text((*ncid).ncid, (*ncid).id_ssa, "units", strlen(unit), unit);

  /* Asymmetry factor */
  
  if ((retval = nc_def_var((*ncid).ncid, "gg", NC_DOUBLE, 2, 
                           dimids2D, &(*ncid).id_gg))){
    fprintf(stderr, "Error defining variable gg. \n");
    return -1;
  } 
  strcpy(name,"Asymmetry factor");
  strcpy(unit,"-");
  nc_put_att_text((*ncid).ncid, (*ncid).id_gg, "long_name", strlen(name), name);  
  nc_put_att_text((*ncid).ncid, (*ncid).id_gg, "units", strlen(unit), unit);
  
  /* Refractive index (real) */
  
  if(medium != AEROSOL){
    if ((retval = nc_def_var((*ncid).ncid, "refre", NC_DOUBLE, 1, 
                             &(*ncid).id_nlam, &(*ncid).id_refre))){
      fprintf(stderr, "Error defining variable refre. \n");
      return -1;
    }
  }
  else{
    if ((retval = nc_def_var((*ncid).ncid, "refre", NC_DOUBLE, 2, 
                             dimids2D, &(*ncid).id_refre))){
      fprintf(stderr, "Error defining variable refre. \n");
      return -1;
    }
  }
      
  strcpy(name,"refractive index (real)");
  strcpy(unit," ");
  nc_put_att_text((*ncid).ncid, (*ncid).id_refre, "long_name", strlen(name), name);  
  nc_put_att_text((*ncid).ncid, (*ncid).id_refre, "units", strlen(unit), unit);
  
  
  /* Refractive index (imaginary) */
  if(medium != AEROSOL){
    if ((retval = nc_def_var((*ncid).ncid, "refim", NC_DOUBLE, 1, 
                             &(*ncid).id_nlam, &(*ncid).id_refim))){
      fprintf(stderr, "Error defining variable refre. \n");
      return -1;
    }
  }
  else{
    if ((retval = nc_def_var((*ncid).ncid, "refim", NC_DOUBLE, 2, 
                             dimids2D, &(*ncid).id_refim))){
      fprintf(stderr, "Error defining variable refre. \n");
      return -1;
    }
  } 
  
  strcpy(name,"refractive index (imaginary)");
  strcpy(unit," ");
  nc_put_att_text((*ncid).ncid, (*ncid).id_refim, "long_name", strlen(name), name);  
  nc_put_att_text((*ncid).ncid, (*ncid).id_refim, "units", strlen(unit), unit);
  
  if ((retval = nc_def_var((*ncid).ncid, "rho", NC_DOUBLE, 2, 
                           dimids2D, &(*ncid).id_rho))){
    fprintf(stderr, "Error defining variable refre. \n");
    return -1;
  }
  strcpy(name,"density of medium");
  strcpy(unit,"g/cm^3");
  nc_put_att_text((*ncid).ncid, (*ncid).id_rho, "long_name", strlen(name), name);  
  nc_put_att_text((*ncid).ncid, (*ncid).id_rho, "units", strlen(unit), unit);
  
  
  /* End definitions */
  if ((retval = nc_enddef((*ncid).ncid))){
    fprintf(stderr, "Error netcdf definition end. \n");
    return -1;
  } 

  /* Variables that are not depedent on ir, iv  */
  /* wavelengths */
  for (i=0; i<nlam; i++)
    lambda_tmp[i]=lambda[i]*1e-3; 
  if ((retval = nc_put_var_double((*ncid).ncid, (*ncid).id_lam, &lambda_tmp[0]))){
    fprintf(stderr, "Error writing variable lambda.\n");
    return -1;
  } 
  
  /* effective radius */
  if ((retval = nc_put_var_double((*ncid).ncid, (*ncid).id_reff, &reff[0]))){
    fprintf(stderr, "Error writing variable reff.\n");
    return -1;
  } 

  /* Close netcdf file. This frees up any internal netCDF resources
   * associated with the file, and flushes any buffers. */
  if ((status = nc_close((*ncid).ncid))){
    fprintf(stderr, "Error closing netcdf file.\n");
    return -1;
  } 
  
  return 0; 
}


/** 
 * write_output_netcdf 
 *
 * Write the output of the Mie calculation in netcdf format. 
 * This can directly be used as input for uvspec. 
 * 
 *
 * @author Claudia Emde
 * @date   ????
 *         2009-06-29 (modifications to adapt to new format of uvspec files)
 *         2010-02-10 (rewrite)
 */
int write_output_netcdf(netcdf_ids ncid, 
                        input_struct input, output_struct *output,
                        int iv, int ir)
{
  /* Indices*/
  int ip=0, i=0;
  int retval=0, status=0;
 
  int nphamat=input.nstokes; 
   
  /* Temporary 1D arrays for netcdf output */
  int int_array_tmp[1]; int id; 
  double double_array_tmp[1]; 
  float *pmom_tmp; 
  int nmoms=0, nmommax=0; 
  int comp_moments[1]; 
  
  size_t start1D[1], start2D[2], start3D[3], start4D[4];
  ptrdiff_t stride1D[1], stride2D[2], stride3D[3],stride4D[4] ;
  size_t count1D[1], count2D[2], count3D[3], count4D[4];
  
  /* Initialization*/
  start1D[0]=iv; 
  count1D[0]=1; 
  stride1D[0]=1; 

  start2D[0]=iv;       
  start2D[1]=ir; 
  for (i=0; i<2; i++){
    count2D[i]=1;     /* number of indices to be written*/
    stride2D[i]=1;    /* index interval */
  }

  start3D[0]=iv;       
  start3D[1]=ir; 
  start3D[2]=0;   
  for (i=0; i<3; i++){
    count3D[i]=1;     
    stride3D[i]=1; 
  }

  start4D[0]=iv;       
  start4D[1]=ir; 
  start4D[2]=0;   
  start4D[3]=0; 
  for (i=0; i<4; i++){
    count4D[i]=1;       
    stride4D[i]=1;      
  }
  if(input.nmom_netcdf!=0){
    nmommax=input.nmom_netcdf;
    nmoms=nmommax;
  }
  else{
    nmommax=input.mie.nmom;
  }
  
  pmom_tmp=calloc(nmommax, sizeof(float));
    
  /* Open netcdf file. */
  if ((status = nc_open(output->outputfile, NC_WRITE, &id))){
    fprintf(stderr, "Error open netcdf file.\n");
    return -1;
  } 

  /* ---------------------------------------------------------------*/
  /* Write variables                                                */
  /* ---------------------------------------------------------------*/

  if(input.verbose)
    fprintf(stderr, "... write data to netcdf file for iv %d, ir %d \n", iv, ir);  

  for (ip=0; ip<nphamat; ip++){
    
    start3D[2]=ip;
    start4D[2]=ip; 

    /* ntheta */
    int_array_tmp[0]=output->ntheta[ip];
    if ((retval = nc_put_vars_int(id, ncid.id_ntheta, 
                                  start3D, count3D, stride3D, int_array_tmp))){
      fprintf(stderr, "Error writing variable ntheta.\n");
      return -1;
    }
    
    /* theta */
    count4D[3]=output->ntheta[ip];   
    if ((retval = nc_put_vars_float(id, ncid.id_theta,
                                    start4D, count4D, stride4D, output->theta[ip]))){
      fprintf(stderr, "Error writing variable theta.\n");
      return -1;
    }
    
    /* phase functions */
    if ((retval = nc_put_vars_float(id, ncid.id_phase,
                                    start4D, count4D, stride4D, output->phase[ip]))){
      fprintf(stderr, "Error writing variable phase.\n");
      return -1;
    }
    
    /* Number of Legendre polynomials in netcdf file, if option nmom_netcdf is not specified. */
    if( input.nmom_netcdf==0 ){
      nmoms=0; 
      for (i=0;i<=output->nmommax;i++){
        if (fabs((2*i+1)*output->pmom[ip][i]) < 1E-5){
          if ( i > nmoms ){
            nmoms=i;
            break;
          }
        }
      }
    }
    
    if(nmoms==0){
      nmoms=output->nmommax;
      comp_moments[0]=0; 
      nc_put_att_int(id, NC_GLOBAL, "complete_moments", NC_INT, 1, comp_moments);
    }
      
    /* Number of Legendre polynomials  */
    int_array_tmp[0]=nmoms;
    if ((retval = nc_put_vars_int(id, ncid.id_nmom, 
                                  start3D, count3D, stride3D, int_array_tmp))){
      fprintf(stderr, "Error writing variable nmom.\n");
      return -1;
    }
    
    /* Legendre polynomials */
    count4D[3]=nmommax;
    for (i=0; i<nmommax; i++)
      pmom_tmp[i]=output->pmom[ip][i]*(2.0*i+1.0);
    if ((retval = nc_put_vars_float(id, ncid.id_pmom,
                                    start4D, count4D, stride4D, pmom_tmp))){
      fprintf(stderr, "Error writing variable pmom.\n");
      return -1;
    }
    
  }
  
  /* extinction coeffidient */
  double_array_tmp[0]=output->beta/output->rho_medium;
  if ((retval = nc_put_vars_double(id, ncid.id_ext,
                                   start2D, count2D, stride2D, double_array_tmp))){
    fprintf(stderr, "Error writing variable ext.\n");
    return -1;
  }
  
  /* single scattering albedo */
  double_array_tmp[0]=output->omega;
  if ((retval = nc_put_vars_double(id, ncid.id_ssa,
                                   start2D, count2D, stride2D, double_array_tmp))){
    fprintf(stderr, "Error writing variable ssa.\n");
    return -1;
  }

  /* asymmetry factor */
  if (input.distribution == DIST_NONE)
    double_array_tmp[0]=output->mie.gsca;
  else
    double_array_tmp[0]=output->g;
  if ((retval = nc_put_vars_double(id, ncid.id_gg,
                                   start2D, count2D, stride2D, double_array_tmp))){
    fprintf(stderr, "Error writing variable gg.\n");
    return -1;
  }
  
  /* refractive index (real) */ 
  double_array_tmp[0]=output->ref.re;
  if(input.medium != AEROSOL){
    if ((retval = nc_put_vars_double(id, ncid.id_refre,
                                     start1D, count1D, stride1D, double_array_tmp))){
      fprintf(stderr, "Error writing variable refre.\n");
      return -1;
    }
  }
  else{
    if ((retval = nc_put_vars_double(id, ncid.id_refre,
                                     start2D, count2D, stride2D, double_array_tmp))){
      fprintf(stderr, "Error writing variable refre.\n");
      return -1;
    }
  }
  
  /* refractive index (imag) */ 
  double_array_tmp[0]=output->ref.im;
  if(input.medium != AEROSOL){
    if ((retval = nc_put_vars_double(id, ncid.id_refim,
                                     start1D, count1D, stride1D, double_array_tmp))){
      fprintf(stderr, "Error writing variable refim.\n");
      return -1;
    }
  }
  else{
    if ((retval = nc_put_vars_double(id, ncid.id_refim,
                                     start2D, count2D, stride2D, double_array_tmp))){
      fprintf(stderr, "Error writing variable refim.\n");
      return -1;
    }
  }
  
  /* density */ 
  double_array_tmp[0]=output->rho_medium;
  if ((retval = nc_put_vars_double(id, ncid.id_rho,
                                   start2D, count2D, stride2D, double_array_tmp))){
    fprintf(stderr, "Error writing variable rho.\n");
    return -1;
  }

  if ((status = nc_close(id))){
    fprintf(stderr, "Error closing netcdf file.\n");
    return -1;
  } 
  
  return 0; 
}
 
#endif



/**
 * calc_phase_functiom
 *
 * Calculate the phase function from Legendre moments. The phase
 * function is stored on an optimized scattering angle grid. If the
 * accuracy of the phase function is less than 1% (when a specified
 * maximum number of grid points - nthetamax - is used) a warning is
 * printed to the screen. 
 *
 * @author Claudia Emde
 * @date   2010-02-11 restructured 
 */
int calc_phase_function(/* Input */
                        int  nmom,
                        double lambda, 
                        double reff,
                        int nstokes, 
                        int nthetamax,
                        double accuracy,
                        /* Output */
                        float*** pmom,
                        int* nmommax,
                        float*** phase, 
                        float*** theta,
                        int** ntheta,
                        int* warnings_nmom,
                        int* warnings_ntheta)
{
  int ip=0, it=0, im=0; 
  double theta_i=0.0, dtheta=0.0;/* accuracy=0.0; */
  double *pmom_tmp=NULL, *mu_tmp=NULL, **phase_tmp=NULL;
  double *mu_opt=NULL, **phase_opt=NULL; 
  int Nopt=0, Ntheta=0;
  
  int status=0;
  
  
  /* Find number of polynomials that are non-zero, this is done for P11.*/
  /*   For P12, P34 the */
  /* polynomials might oscillate around 0 so it is here not a good idea to cut the Legendre series */
  /* when the values are small */
  *nmommax=0;
  for (im=0; im<=nmom; im++){
    if (fabs( (2*im+1) * (*pmom)[0][im] ) < 1E-5){
      if ( im > (*nmommax) )
        (*nmommax)=im;
      break; 
    }
  }
  
  if(*nmommax !=0 )
    for (ip=0; ip<nstokes; ip++){
      /* set values above nmommax to 0 */
      for (im= *nmommax+1; im <=nmom; im++)
        (*pmom)[ip][im]=0.0;
    }
  
  /* remember to write a warning, if there could be more
     Legendre coefficients, only phase function is checked, the
     other phase matrix elements not */
  if (*nmommax == 0){
    *nmommax=nmom;
    *warnings_nmom=1;
  }

  /* check size parameter, if larger than 250, 2 digits in angular 
     grid are not sufficient */
  if( 2*PI*reff/lambda > 250.0){
    Ntheta=180001;
    dtheta=0.001;
    fprintf(stderr, "3 digits optimization \n");
  }
  else{
    Ntheta=18001; 
    dtheta=0.01; 
    fprintf(stderr, "2 digits optimization \n");
  }
  
  /* Calculate phase functions */
  pmom_tmp = calloc(nmom, sizeof(double)); 
  mu_tmp = calloc(Ntheta, sizeof(double)); 
  phase_tmp = calloc(nstokes, sizeof(double *)); 
  phase_opt= calloc(nstokes, sizeof(double *));
  /* Allocate maximum allowed number of grid points */
  mu_opt= calloc(nthetamax, sizeof(double)); 
  for (ip=0; ip<nstokes; ip++){
    phase_tmp[ip] = calloc(Ntheta,sizeof(double)); 
    phase_opt[ip] = calloc(nthetamax, sizeof(double)); 
  }
  

  for (ip=0; ip<nstokes; ip++){
    
    for (im=0; im<nmom; im++)
      pmom_tmp[im] = (2*im+1)* (*pmom)[ip][im];
      
    for (it=0; it<Ntheta; it++){
      theta_i=180.0-it*dtheta;
      mu_tmp[it]=cos(theta_i/180.0*PI);
      phase_tmp[ip][it]=mom2phase(cos(theta_i/180.0*PI),
                                  pmom_tmp, nmom);
    }
  }
    
  /* Accuracy of phase function: 1% */
  /*accuracy=0.01;*/
  /* Optimize angular grid for phase function */ 
  status=optimize_theta_grid_phamat(&mu_opt, &phase_opt, &Nopt,
                                    mu_tmp, phase_tmp, Ntheta, accuracy,
                                    nthetamax, nstokes, &(*warnings_ntheta));
  if (status != 0) {
    fprintf (stderr, "Error %d returned by optimize_theta_grid()\n", status);
    return status;
  }
  for (ip=0; ip<nstokes; ip++){
    (*ntheta)[ip]=Nopt; 
    for (it=0; it<Nopt; it++){
      /* At the moment the maximum accuracy is 3 digits. */
      /* Print a warning if this accuracy is not sufficient */
      (*theta)[ip][it]=acos(mu_opt[it])*180.0/PI;
      (*phase)[ip][it]=phase_opt[ip][it];
      if((*theta)[ip][it]==0.001 && ip==0)
        fprintf(stderr, "Warning: The accuracy of the angular grid is not sufficient !!!!\n)");
    }
  }  

  for (ip=0; ip<nstokes; ip++){
    free(phase_tmp[ip]); 
    free(phase_opt[ip]); 
  }
  
  free(mu_tmp); free(phase_tmp); free(pmom_tmp); 
  free(mu_opt); free(phase_opt);
  
  return 0; 
}
  
