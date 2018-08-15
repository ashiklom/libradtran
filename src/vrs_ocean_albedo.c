/*--------------------------------------------------------------------
 * $Id: vrs_ocean_albedo.c 2838 2012-12-19 11:06:09Z robert.buras $
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
  
/*
<lpdoc>
\subsection{Calculate albedo of snow - \codeidx{Gen\_snow\_tab}, \codeidx{snowalbedo}}

The \code{vrs_ocean_albedo} program may be used to calculate the ocean
reflectance without and with vibrational Raman scattering as  
formulated by \citet{Vountas2003}.

Note that the \citet{Bricaud1995} regression coefficients for particulate
matter are scaled to agree with those from \citet{Vasilkov2005}. This is done
to avoid unphysical jumps in the calculated albedo. Obviously more data is required
for optical properties of particulate matter.

Run  \code{vrs_ocean_albedo} with \code{-h} to see its usage.

</lpdoc>

*/
/*
  @c17@ 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#include "ascii.h"
#include "ancillary.h"
#include "numeric.h"
#include "table.h"
#include "f77-uscore.h"

#define PROGRAM "vrs_ocean_albedo"

#define NSUM 4

double frequency_redistribution(double nu, double nup)
   {
     double f_R=0;
     double k_R = 0.005152;  /* Constant k_R, Haltrinn and Kattawar, 1993, page 5365.*/
     double arg=0.0;
     int i=0,nsum = NSUM;
     static double alpha[NSUM]; /* alpha_i in Table 1, Haltrinn and Kattawar, 1993, page 5365.*/
     static double sigma[NSUM]; /* sigma_i in Table 1, Haltrinn and Kattawar, 1993, page 5365.*/
     static double dnu[NSUM];   /* dnu_i in Table 1, Haltrinn and Kattawar, 1993, page 5365.*/
     alpha[0] = 0.41;
     alpha[1] = 0.39;
     alpha[2] = 0.10;
     alpha[3] = 0.10;
     dnu[0]   = 3250;
     dnu[1]   = 3425;
     dnu[2]   = 3530;
     dnu[3]   = 3525;
     sigma[0] = 89.179;
     sigma[1] = 74.317;
     sigma[2] = 59.453;
     sigma[3] = 59.453;

     /* Eq. A2, Haltrinn and Kattawar, 1993, page 5365.*/
     f_R=0.0;
     for(i=0;i<nsum;i++){
       arg =((nup-nu-dnu[i])*(nup-nu-dnu[i]))/(2*sigma[i]*sigma[i]);
       f_R += alpha[i]*exp(-arg);
     }
     f_R *= k_R;
     
     return f_R;
   }

static void print_usage (char *name)
{
  fprintf (stderr, "\nusage:   %s  [options]\n\n", name);
  fprintf (stderr, "%s  calculates the ocean reflectance according to\n", name);
  fprintf (stderr, "Vountas et al., Atmos. Chem. Phys., 3, 1365-1375, 2003.\n");
  fprintf (stderr, "Vountas et al., Atmos. Chem. Phys., 3, 1365-1375, 2003.\n");
  fprintf (stderr, "Required arguments:\n");
  fprintf (stderr, "  -f            uvspec output filename. The corresponding input file\n");
  fprintf (stderr, "                must include the following output specification:\n");
  fprintf (stderr, "                output_user lambda eglo\n");
  fprintf (stderr, "                Also, normally you would want to include (for obvious reasons):\n");
  fprintf (stderr, "                zout boa\n");
  fprintf (stderr, "  -s            The first wavelength for which the VRS ocean reflectance\n");
  fprintf (stderr, "                will be calculated (nm). The uvspec output file must\n");
  fprintf (stderr, "                include this wavelength.\n");
  fprintf (stderr, "Optional arguments:\n");
  fprintf (stderr, "  -C            The chlorophyll content in mg/m^3. Default value 1.e-20.\n");
  fprintf (stderr, "                May not be zero.\n");
  fprintf (stderr, "\n");
  fprintf (stderr, "The output has five columns\n");
  fprintf (stderr, "\n");
  fprintf (stderr, "1. column: wavelength (nm)\n");
  fprintf (stderr, "2. column: Total (column 3+4) ocean reflectance.\n");
  fprintf (stderr, "3. column: Elastic scattering ocean reflectance.\n");
  fprintf (stderr, "4. column: vibrational Raman scattering (VRS) ocean reflectance.\n");
  fprintf (stderr, "5. column: Global radiation incident on ocean surface, from uvspec output file.\n");
}

int main (int argc, char **argv)
{
  char *s=NULL;
  char libradtran_data_path[FILENAME_MAX+200] = "";
  char tmp[FILENAME_MAX+200] = "";
  char uvspec_output_file[FILENAME_MAX+200] = "";
  char Bricaud1995_filename[FILENAME_MAX+200] = "albedo/Bricaud1995Table2.dat";
  char Lu2006_filename[FILENAME_MAX+200] = "albedo/Lu2006_TableV1.dat";
  char SB1981_filename[FILENAME_MAX+200] = "albedo/SmithBaker1981Table1.dat";
  char Vasilkov2005_filename[FILENAME_MAX+200] = "albedo/Vasilkov2005Table1.dat";
  double **data=NULL,*R_E=NULL,*R_R=NULL,R_R_int=0.0;
  double fact=0.0;
  double f_R=0,nu=0,nup=0;
  double *eglo=NULL,*lambda=NULL,*aw=NULL,aw_440=0,*bw=NULL,bw_550=0,lambdap=0;
  double *ac=NULL,ac_440=0;
  double adom=0;
  double *lambda_lu=NULL,*aw_lu=NULL, lambda_tmp[1];
  double *lambda_SB=NULL,*bw_SB=NULL;
  double *lambda_Bricaud=NULL,*A_Bricaud=NULL,*B_Bricaud=NULL;
  double *lambda_Vasilkov=NULL,*A_Vasilkov=NULL,*B_Vasilkov=NULL;
  double *lambda_VB=NULL,*ac_VB=NULL;
  double S_Vasilkov=0.014; /* nm-1 */
  double nmtocm_1=1.e-07; 
  double *K=0,*k_E=0,*a=0,*b_b=0,mu_d=0,mu_E_u=0,s_E=0;
  double E_d=0,E_dp=0,b_b_R=0,k_R=0,mu_R_u=0,bracket=0,R_tot=0;
  double b_b_Ep=0,Kp=0,k_Ep=0;
  double C=1.e-20; /* Chlorophyll content in mg/m**3 */
  double lambda_start=-999;
  double sum=0,dlambda=0,nu1=0,nu2=0,dnu=0;
  int c=0;
  int iv=0,ivp=0,ivs=0;
  /* commented out by RB, because was not used but caused compiler warnings */
  /* int print =1; */
  int rows=0, min_columns=0, max_columns=0;
  int nlambda=0,nlambda_lu=0,nlambda_SB=0,nlambda_Bricaud=0,nlambda_Vasilkov=0,nlambda_VB=0;
  int status=0;

  s = NULL; 
  s = getenv("LIBRADTRAN_DATA_FILES");
  if ( s == NULL ) { s = "../data/"; }
  strcpy (libradtran_data_path, s);  
  if (libradtran_data_path[strlen(libradtran_data_path)-1]!='/') {
    strcat(libradtran_data_path,"/");
  }

  strcpy (tmp,libradtran_data_path);
  strcat (tmp,Lu2006_filename);
  strcpy (Lu2006_filename, tmp);

  strcpy (tmp,libradtran_data_path);
  strcat (tmp,SB1981_filename);
  strcpy (SB1981_filename, tmp);

  strcpy (tmp,libradtran_data_path);
  strcat (tmp,Vasilkov2005_filename);
  strcpy (Vasilkov2005_filename, tmp);

  strcpy (tmp,libradtran_data_path);
  strcat (tmp,Bricaud1995_filename);
  strcpy (Bricaud1995_filename, tmp);

  while ((c=getopt (argc, argv, "C:f:s:")) != EOF)  {
    switch(c)  {
    case 'C': 
      C = atof(optarg);
      break;
    case 's': 
      lambda_start = atof(optarg);
      break;
    case 'f': 
      strcpy (uvspec_output_file, optarg);
      break;
    case 'p':  /* no print messages */  
      /* commented out by RB, because was not used but caused compiler warnings */
      /* print = 0; */
      break;   
    case 'h':
      print_usage (argv[0]);
      return (-1);
      break;
    case '?':
      print_usage (argv[0]);
      return (-1);
      break;
    default:
      print_usage (argv[0]);
      return (-1);
    }
  }

  /* Read uvspec output file containing wavelength and global radiation at bottom */
  /* of the atmosphere */
  if ((status = ASCII_file2double (uvspec_output_file, 
				 &rows, &max_columns, &min_columns, &data)) != 0) {
    fprintf (stderr, "Error %d reading %s\n", status, uvspec_output_file);
    return status;
  }
  nlambda=rows;
  if ((lambda  = calloc (nlambda, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'lambda'\n");
    return -1;
  }
  if ((eglo  = calloc (nlambda, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'eglo'\n");
    return -1;
  }
  for (iv=0;iv<nlambda;iv++) {
    lambda[iv] = data[iv][0];
    eglo[iv]   = data[iv][1];
  }
  ASCII_free_double (data, rows);

  if ( lambda_start<0.0 || lambda_start < lambda[0] ) {
    fprintf (stderr,"\nError, 'lambda_start' (option -s) too small.\n\n");
    print_usage (argv[0]);
  }

  /* Read Lu 2006 water absorption coefficients */
  if ((status = ASCII_file2double (Lu2006_filename, 
				 &rows, &max_columns, &min_columns, &data)) != 0) {
    fprintf (stderr, "Error %d reading %s\n", status, Lu2006_filename);
    return status;
  }
  nlambda_lu=rows;
  if ((lambda_lu  = calloc (nlambda_lu, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'lambda'\n");
    return -1;
  }
  if ((aw_lu  = calloc (nlambda_lu, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'aw_lu'\n");
    return -1;
  }
  for (iv=0;iv<nlambda_lu;iv++) {
    lambda_lu[iv] = data[iv][0];
    aw_lu[iv]   = data[iv][1];
  }
  ASCII_free_double (data, rows);
  /* Interpolate aw_lu to wavelengths from uvspec output file */
  if ((aw  = calloc (nlambda, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'aw'\n");
    return -1;
  }
  lambda_tmp[0]=440.0;
  
  status = arb_wvn_double(nlambda_lu, lambda_lu, aw_lu, 1, lambda_tmp, aw, 1, 0);
  aw_440 = aw[0];
  status = arb_wvn_double(nlambda_lu, lambda_lu, aw_lu, nlambda, lambda, aw, 1, 0);

  /* Read Smith and Baker 1981 seawater back scattering coefficients */
  if ((status = ASCII_file2double (SB1981_filename, 
				 &rows, &max_columns, &min_columns, &data)) != 0) {
    fprintf (stderr, "Error %d reading %s\n", status, SB1981_filename);
    return status;
  }
  nlambda_SB=rows;
  if ((lambda_SB  = calloc (nlambda_SB, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'lambda_SB'\n");
    return -1;
  }
  if ((bw_SB  = calloc (nlambda_SB, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'bw_SB'\n");
    return -1;
  }
  for (iv=0;iv<nlambda_SB;iv++) {
    lambda_SB[iv] = data[iv][0];
    bw_SB[iv]     = data[iv][3];
  }
  ASCII_free_double (data, rows);
  /* Interpolate bw_SB to wavelengths from uvspec output file */
  if ((bw  = calloc (nlambda, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'bw'\n");
    return -1;
  }
  /* Get bw at 550 nm */
  lambda_tmp[0]=550.0;
  status = arb_wvn_double(nlambda_SB, lambda_SB, bw_SB, 1, lambda_tmp, bw, 1, 0);
  bw_550 = bw[0]; 
  status = arb_wvn_double(nlambda_SB, lambda_SB, bw_SB, nlambda, lambda, bw, 1, 0);

  /* Read Vasilkov et al. 2005 regression coefficients for particulate matter. */
  if ((status = ASCII_file2double (Vasilkov2005_filename, 
				 &rows, &max_columns, &min_columns, &data)) != 0) {
    fprintf (stderr, "Error %d reading %s\n", status, Vasilkov2005_filename);
    return status;
  }
  nlambda_Vasilkov=rows;
  if ((lambda_Vasilkov  = calloc (nlambda_Vasilkov, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'lambda_Vasilkov'\n");
    return -1;
  }
  if ((A_Vasilkov  = calloc (nlambda_Vasilkov, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'A_Vasilkov'\n");
    return -1;
  }
  if ((B_Vasilkov  = calloc (nlambda_Vasilkov, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'B_Vasilkov'\n");
    return -1;
  }
  for (iv=0;iv<nlambda_Vasilkov;iv++) {
    lambda_Vasilkov[iv] = data[iv][0];
    A_Vasilkov[iv]      = data[iv][1];
    B_Vasilkov[iv]      = data[iv][2];
  }
  ASCII_free_double (data, rows);

  /* Read Bricaud et al. 1995 regression coefficients for particulate matter. */
  if ((status = ASCII_file2double (Bricaud1995_filename, 
				 &rows, &max_columns, &min_columns, &data)) != 0) {
    fprintf (stderr, "Error %d reading %s\n", status, Bricaud1995_filename);
    return status;
  }
  nlambda_Bricaud=rows;
  if ((lambda_Bricaud  = calloc (nlambda_Bricaud, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'lambda_Bricaud'\n");
    return -1;
  }
  if ((A_Bricaud  = calloc (nlambda_Bricaud, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'A_Bricaud'\n");
    return -1;
  }
  if ((B_Bricaud  = calloc (nlambda_Bricaud, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'B_Bricaud'\n");
    return -1;
  }
  for (iv=0;iv<nlambda_Bricaud;iv++) {
    lambda_Bricaud[iv] = data[iv][0];
    A_Bricaud[iv]      = data[iv][1];
    B_Bricaud[iv]      = data[iv][2];
  }
  ASCII_free_double (data, rows);

  /* Combine Bricaud and Vasilkov data. Take 400 nm data from Vasilkov, hence ignore */
  /* first element in Bricaud */
  nlambda_VB=nlambda_Vasilkov+nlambda_Bricaud-1;
  if ((ac_VB  = calloc (nlambda_VB, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'ac_VB'\n");
    return -1;
  }
  if ((lambda_VB  = calloc (nlambda_VB, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'lambda_VB'\n");
    return -1;
  }
  for (iv=0;iv<nlambda_Vasilkov;iv++) {
    lambda_VB[iv] = lambda_Vasilkov[iv];
    /* Vasilkov et al. 2005, page 2864 */
    ac_VB[iv]     = C*A_Vasilkov[iv]*pow(C,-B_Vasilkov[iv]);
  }
  fact = (C*A_Vasilkov[nlambda_Vasilkov-1]*pow(C,-B_Vasilkov[nlambda_Vasilkov-1]))/
    (C*A_Bricaud[0]*pow(C,-B_Bricaud[0]));
  for (iv=1;iv<nlambda_Bricaud;iv++) {
    lambda_VB[nlambda_Vasilkov-1+iv] = lambda_Bricaud[iv];
    /* Vasilkov et al. 2005, page 2864 */
    ac_VB[nlambda_Vasilkov-1+iv]     = fact*C*A_Bricaud[iv]*pow(C,-B_Bricaud[iv]);
  }

  /* Get ac at uvspec wavelengths */
  if ((ac  = calloc (nlambda, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'ac'\n");
    return -1;
  }
  lambda_tmp[0]=440.0;
  status = arb_wvn_double(nlambda_VB, lambda_VB, ac_VB, 1, lambda_tmp, ac, 1, 0);
  ac_440 = ac[0]; 
  status = arb_wvn_double(nlambda_VB, lambda_VB, ac_VB, nlambda, lambda, ac, 1, 0);

  /* Allocate memory for elastic and raman reflectances and more */
  if ((a  = calloc (nlambda, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'a'\n");
    return -1;
  }
  if ((K  = calloc (nlambda, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'K'\n");
    return -1;
  }
  if ((k_E  = calloc (nlambda, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'k_E'\n");
    return -1;
  }
  if ((b_b  = calloc (nlambda, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'b_b'\n");
    return -1;
  }
  if ((R_E  = calloc (nlambda, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'R_E'\n");
    return -1;
  }
  if ((R_R  = calloc (nlambda, sizeof (double)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'R_R'\n");
    return -1;
  }

  s_E = 1.0;  /* According to Vountas et al. 2003, p. 1368 and references therein. */
  mu_d= 0.75; /* According to Vountas et al. 2003, p. 1368 and references therein. */
  mu_E_u= 0.5;  /* According to Vountas et al. 2003, p. 1368 and references therein. */
  mu_R_u= 0.5;  /* According to Vountas et al. 2003, p. 1368 and references therein. */
  
  /* Calculate needed quantities for all wavelengths in uvspec output file*/
  for (iv=0;iv<nlambda;iv++) {
    adom = 0.2*(aw_440+ac_440)*exp(S_Vasilkov*(lambda[iv]-440.0)); /* Vountas  et al. 2003, p. 1368 */
    a[iv]   = aw[iv]+ac[iv]+adom;
    /* Eq 5, Vountas et al. 2003, p. 1369 */
    b_b[iv] = 0.5*bw[iv]+(0.002+0.02*(0.5-0.25*log10(C))*550.0/lambda[iv])*(0.3*pow(C,0.62)-bw_550);

    K[iv]   = (a[iv]+b_b[iv])/mu_d;       /* Vountas et al. 2003, p. 1368 */
    k_E[iv] = (a[iv]+b_b[iv])/mu_E_u;   /* Vountas et al. 2003, p. 1368 */
    R_E[iv] = (s_E * b_b[iv]/mu_d)* (1.0/(K[iv] + k_E[iv]));  /* Eq. 2 Sathyendranath and Platt, 1998. */

    /* Eq 4, Vountas et al. 2003, p. 1368 */
    /* Sum over all primed (p) wavelengths to get vibrational Raman scattering contribution */
  }

  ivs = 0;
  while ( lambda[ivs] < lambda_start) {
    ivs++;
  }
  for (iv=ivs;iv<nlambda;iv++) {
    lambdap=lambda_start-55; /* To go 55 nm shorter should cover most in the UV/visible */
    if ( lambdap < lambda[0] ) {
      fprintf (stderr,"\nError, the wavelength region covered by uvspec file is to small.\n");
      fprintf (stderr,"\nIt starts at %f nm. It needs to start at %f nm..\n", lambda[0], lambdap );
      exit(0);
    }
    sum =0;
    dlambda = lambda[1]-lambda[0];
    nu1 = 1/(lambdap*nmtocm_1);
    nu2 = 1/((lambdap+dlambda)*nmtocm_1);
    dnu = nu1-nu2;
    nu  = 1./ (lambda_start*nmtocm_1) ;
    E_dp  = 0.0;
    b_b_R = 0.0;
    b_b_Ep= 0.0;
    Kp    = 0.0;
    k_Ep  = 0.0;
    ivp   = 0;
    R_R_int = 0.0;
    while ( lambda[ivp] < lambdap) {
	ivp++;
      }
    while (lambdap < lambda[iv]){ /* Only consider shorter wavelengths, Mobley (1994) */ 
      nup = 1. / (lambdap*nmtocm_1); 
      f_R= frequency_redistribution(nu, nup); 
      sum+=dnu*f_R;

      E_dp = eglo[ivp]*dnu;//*f_R;
      Kp   = ((a[ivp]+b_b[ivp])/mu_d);//*dnu*f_R;
      k_Ep = ((a[ivp]+b_b[ivp])/mu_E_u);//*dnu*f_R;
      b_b_Ep= (0.5*bw[ivp]+(0.002+0.02*(0.5-0.25*log10(C))*550.0/lambdap)*(0.3*pow(C,0.62)-bw_550));//*dnu*f_R;

      b_b_R = 2.610e-4*pow(lambdap/488.,5.3); /* m-1. Eq 6, Vountas et al. 2003, p. 1369 */
      E_d = eglo[iv];
      k_R  = (a[iv]+b_b[iv])/mu_R_u;
      bracket = (1+b_b[iv]/k_E[iv]+(2*b_b_Ep/(Kp+k_Ep)));
      R_R_int += (E_dp/E_d)*(b_b_R/mu_d)*(1./(K[iv]+k_R))*bracket*dnu*f_R; 

      lambdap+=dlambda; 
      nu1 = 1/(lambdap*nmtocm_1);
      nu2 = 1/((lambdap+dlambda)*nmtocm_1);
      dnu = nu1-nu2;
      ivp++;
    } 
    
    R_R[iv] = R_R_int;
    R_tot  = R_E[iv]+R_R[iv];//*bracket;

    printf("%12.6f %12.6e %12.6e %12.6e %12.6e\n",
	   lambda[iv],R_tot,eglo[iv],R_E[iv],R_R[iv]);
  }

  free(eglo);
  free(lambda);
  free(lambda_lu);
  free(lambda_SB);
  free(lambda_VB);
  free(lambda_Vasilkov);
  free(a);
  free(aw);
  free(aw_lu);
  free(bw);
  free(bw_SB);
  free(b_b);
  free(ac);
  free(ac_VB);
  free(A_Vasilkov);
  free(B_Vasilkov);
  free(K);
  free(k_E);
  free(R_E);
  free(R_R);
 
  exit(0);

}
