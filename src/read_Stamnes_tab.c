/*--------------------------------------------------------------------
 * $Id: read_Stamnes_tab.c 2623 2011-12-23 10:52:38Z robert.buras $
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

/*  @27c@                

<lpdoc>

\subsection{Stamnes tables for ozone and cloud optical depth \index{Stamnes tables}}

\citet{Stamnes1991a} devised a method to derive the total ozone column 
and cloud optical depth from global irradiance measurements. For ozone 
column  retrieval this method requires a table of irradiance ratios as 
a function of solar zenith angle and ozone column. The irradiance ratio 
is taken as the ratio of irradiances at non-absorbing and ozone-absorbing  
wavelengths. The cloud optical depth is retrieved from tables of cloud/cloudless
irradiance ratios as a function of solar zenith angle and water cloud optical
depth. 

The {\sl libRadtran} package comes with three tools for calculation and reading
of these so-called Stamnes tables.
The Perl script \codeidx{Gen\_o3\_tab.pl} is used to generate a matrix of
ozone values for solar zenith angle versus a chosen ratio of global
irradiances at different wavelengths. For cloud optical depths the 
Perl script \codeidx{Gen\_wc\_tab.pl} may be used to generate a matrix of
cloud optical depth for solar zenith angle versus a chosen global
irradiance at a selected wavelength. Both tables may be read by the C program
\code{read\_Stamnes\_tab} which, for a solar zenith angle and a measured ratio,
returns the overhead ozone column or cloud optical depth. The Perl scripts 
\code{Gen\_o3\_tab.pl} and \code{Gen\_wc\_tab.pl} and the C program are 
briefly described below. For example of their use please see 
\citet{Mayer1998b,Kylling2005,mayer2005}.

\subsubsection{Generation of the Stamnes ozone column table- \codeidx{Gen\_o3\_tab}}

The Perl script \code{Gen\_o3\_tab.pl} is used to generate a matrix of
ozone values for solar zenith angle versus a chosen ratio of global
irradiance at different wavelengths. The table is read by the C program
\code{read\_Stamnes\_tab} which, for a solar zenith angle and a measured 
irradiance ratio, returns the overhead ozone column. The following options
are understood by \code{Gen\_o3\_tab.pl}:
\begin{description}
    \item[--absolute]         The wavelengths in the bandpass files are in
                               absolute units. Default is relative units.
    \item[--albedo $<$value$>$] Lambertian surface albedo. Default is 0.0.
    \item[--alpha $<$value$>$]  Angstrom alpha coefficient. Default is 0.0.
    \item[--beta $<$value$>$]   Angstrom beta coefficient. Default is 0.0.
    \item[--altitude $<$value$>$] Altitude above sea level [km]. Default is 0.0.
    \item[--atmmod $<$name$>$] Name of atmosphere file. Default atmmod/afglus.dat.
    \item[--help]              Prints help message.
    \item[--o3\_crs $<$name$>$] Name of o3 cross section to use. Default is Molina.
                               See \code{uvspec} documentation for other options.
    \item[--slitfunction $<$name$>$] Name of slitfunction file.
    \item[--bandpasslower $<$name$>$] Name of file holding bandpass for lower wavelength.
    \item[--bandpassupper $<$name$>$]  Name of file holding bandpass for upper wavelength.
    \item[--file $<$name$>$]    Name of file where the table will be stored.
    \item[--lower\_lambda $<$value$>$] Value for lower wavelength, in nm.
    \item[--upper\_lambda $<$value$>$] Value for upper wavelength, in nm.
    \item[--zenith]              Calculate zenith sky radiance table.
\end{description}

Two different types of tables may be generated depending on the measurement type
and the preferred analysis method.

\paragraph{Simple wavelength ratios with \code{Gen\_o3\_tab}}
                  
The simplest type of table is made of ratios of the global irradiance 
at two single wavelengths. This is the type of table described by 
\citet{Stamnes1991a} and it is typically used to analyse measurements 
of the global irradiance from spectroradiometers.
It is generated by the following command (\\ is line continuation character) 

\begin{Verbatim}[fontsize=\footnotesize] 
  perl Gen_o3_tab.pl --slitfunction slitfncfile --lower_lambda 305. \
  --upper_lambda 340. --file table.dat
\end{Verbatim}

Here \file{slitfncfile} is the name of the slit function file. It is
a two column file where the first column is the wavelength (nm, in relative
units) and the second column holds the slit function. The slit function must
be normalized to unity at the center wavelength.

The generated table \file{table.dat} is read by \code{read\_Stamnes\_tab} for a
measured ratio, \code{-r 10.0}, and solar zenith angle, \code{-s 30.0}, 
corresponding to the modelled ratio in the table

\begin{Verbatim}[fontsize=\footnotesize] 
read_Stamnes_tab -r 10.0 -s 30.0 table.dat
\end{Verbatim}


\paragraph{Bandpassed wavelength ratios with \code{Gen\_o3\_tab}}

Instead of using single wavelengths it may be of advantage to use ratios
of irradiances covering a certain wavelength range and weighted with
a bandpass function. 
This approach may reduce problems due to changes
in cloud cover and experimental uncertainties. This approach is also suitable
to calculate ozone columns from multichannel, moderate bandwidth filter
instruments  \citep{Dahlback1996}.
Such tables are generated by

\begin{Verbatim}[fontsize=\footnotesize] 
perl Gen_o3_tab.pl --slitfunction slitfncfile  --lower_lambda 305.0 \ 
                   --upper_lambda 320.0 --file table.dat \ 
                   --bandpasslower bplow.dat --bandpassupper bpupp.dat
\end{Verbatim}

Here \file{bplow.dat} and \file{bpupp.dat} are the bandpass function of the 
lower and upper wavelength region respectively. The bandpass files have two
columns. The first column is the wavelength in nm and relative units to 
\code{--lower\_lambda} and \code{--upper\_lambda}. If absolute units are specified as for
filter instruments, use the \code{--absolute option}. The second column is the
bandpass function.

The tables are read in the same way as the simple wavelength ratio tables.


\subsubsection{Generation of the Stamnes cloud optical thickness table - \codeidx{Gen\_wc\_tab}}

The Perl script \code{Gen\_wc\_tab.pl} is used to generate a matrix of
cloud optical depth for solar zenith angle versus a chosen global
irradiance at a selected wavelength. The wavelength should be chosen 
such that it is not affected by ozone, e.g. 380 nm. 
The table is read by the C program
\code{read\_Stamnes\_tab} which, for a solar zenith angle and a measured irradiance,
returns the overhead cloud optical depth. The available options are 
\begin{description}
    \item[--absolute]         The wavelengths in the bandpass file are in
                               absolute units. Default is relative units.
    \item[--albedo $<$value$>$] Lambertian surface albedo. Default is 0.0.
    \item[--alpha $<$value$>$]  Angstrom alpha coefficient. Default is 0.0.
    \item[--beta $<$value$>$]   Angstrom beta coefficient. Default is 0.0.
    \item[--altitude $<$value$>$] Altitude above sea level [km]. Default is 0.0.
    \item[--atmmod $<$name$>$] Name of atmosphere file. Default atmmod/afglus.dat.
    \item[--help]              Prints help message.
    \item[--o3\_crs $<$name$>$] Name of o3 cross section to use. Default is Molina.
                               See \code{uvspec} documentation for other options.
    \item[--slitfunction $<$name$>$] Name of slitfunction file.
    \item[--bandpass $<$name$>$] Name of file holding bandpass for chosen wavelength.
    \item[--file $<$name$>$]    Name of file where the table will be stored.
    \item[--lambda $<$value$>$] Value of chosen wavelength, in nm.
    \item[--wc\_file $<$name$>$] Name of water cloud file. Default none. Must be specified.

\end{description}

The following different types of tables may be generated.
     
\paragraph{Simple wavelength ratios with \code{Gen\_wc\_tab}}

                  
The simplest type of table is made of the global irradiance 
at a single wavelength. This is the type of table described by 
\citet{Stamnes1991a}. This type
of table is typically used to analyse measurements of the global 
irradiance from spectroradiometers.
It is generated by the following command (\ is line continuation character) 

\begin{Verbatim}[fontsize=\footnotesize] 
perl Gen_wc_tab.pl --slitfunction slitfncfile --lambda 380. \
                   --file table.dat --wc_file ../examples/WC.DAT
\end{Verbatim}
 
Here \file{slitfncfile} is the name of the slit function file. It is
a two column file where the first column is the wavelength (nm, in relative
units) and the second column holds the slit function. The slit function must
be normalized to unity at the center wavelength.

The generated table \file{table.dat} is read by \code{read\_Stamnes\_tab} for a
measured global irradiance, \code{-r 10.0}, and solar zenith angle, \code{-s 30.0}, 
corresponding to the modelled ratio in the table. The table must be corrected 
for the Earth--Sun distance for the day of the measurement. This is achieved
by specifying \code{-d 170}, where \code{170} is the day number. The table is
generated for day 1.

\begin{Verbatim}[fontsize=\footnotesize] 
read_o3_tab -r 10.0 -s 30.0 -d 170 table.dat
\end{Verbatim}     
 
\paragraph{Bandpassed wavelength ratios with Gen\_wc\_tab}
                  
Instead of using a single wavelength it may be of advantage to use
irradiances covering a certain wavelength range and weighted with
a bandpass function. 
This approach may reduce problems due to changes
in cloud cover and experimental uncertainties. This approach is also suitable
to calculate cloud optical depth from multichannel, moderate bandwidth filter
instruments \citep{Dahlback1996}.
Such tables are generated by

\begin{Verbatim}[fontsize=\footnotesize] 
perl Gen_wc_tab.pl --slitfunction slitfncfile  --lambda 380.0 \ 
                   --file table.dat --bandpass bp.dat                   
\end{Verbatim}

Here \file{bp.dat} is the bandpass function of the wavelength region. 
The bandpass file have two columns. The first column is the wavelength 
in nm and relative units to --lambda. If absolute units are specified as for
filter instruments, use the --absolute option. The second column is the
bandpass function.

The tables are read in the same way as the simple wavelength irradiance tables.


 </lpdoc>

*/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include "ascii.h"
#include "numeric.h"
#include "table.h"
#include "sun.h"

typedef struct {
  int      n_zenith;
  int      n_ratio;
  double*   zenith;
  double*   ratio;
  double**  ozone;
} STAMNES_TABLE;

/* Error codes */
#define NOT_A_VALID_STAMNES_TABLE -30
#define RATIO_OUT_OF_RANGE        -31
#define ERROR_SPLINING_ZENITH     -32

#define PROGRAM "read_Stamnes_tab"
#define VERSION "1.0"

static void print_usage (char *name)
{
  fprintf (stderr, "\nusage:   %s -r ratio -s sza [-ph] [-d day_of_year] Stamne_table_file\n\n", name);
  fprintf (stderr, "%s calculates the ozone column from a measured irradiance\n", PROGRAM);
  fprintf (stderr, "ratio for a specific solar zenith angle. A Stamnes ozone table file must\n");
  fprintf (stderr, "be present. See Gen_o3_tab.pl for information on how to generate\n");
  fprintf (stderr, "such a table.\n\n");
  fprintf (stderr, "%s may also be used to calculate the cloud optical depth from a \n", PROGRAM);
  fprintf (stderr, "measured irradiance for a specific solar zenith angle. A Stamnes \n");
  fprintf (stderr, "cloud table file must be present. See Gen_wc_tab.pl for information \n");
  fprintf (stderr, "on how to generate such a table. Use the -d option when inferring cloud\n\n");
  fprintf (stderr, "optical depht.\n\n");
  fprintf (stderr, "Required arguments:\n");
  fprintf (stderr, "  -r  <value>    The measured irradiance ratio corresponding to\n");
  fprintf (stderr, "                 the modelled ratio in the Stamnes table.\n");
  fprintf (stderr, "                 Or the measured irradiance when used to derive.\n");
  fprintf (stderr, "                 water cloud optical depths.\n");
  fprintf (stderr, "  -s  <value>    The solar zenith angle.\n");
  fprintf (stderr, "  -c             Return cloud/clear sky factor  .\n");
  fprintf (stderr, "Optional arguments:\n");
  fprintf (stderr, "  -d             Day of year, correct for Earth-Sun distance.\n");
  fprintf (stderr, "  -h             This help message.\n");
  fprintf (stderr, "  -p             Turn of printing of messages.\n\n");
}

/**************************************************************/
/* Read Stamnes table from file filename.                     */
/* zenith [0..n_zenith-1] is an array of solar zenith angles, */ 
/* ratio [0..n_ratio-1] an array of ratios found in the       */
/* file. The two-dimensional integer array ozone holds the    */
/* corresponding ozone columns.                               */
/* Memory for arrays is allocated automatically.              */
/**************************************************************/

int read_stamnes_table (char *filename, 
			STAMNES_TABLE **table)
{
  int i=0, j=0, status=0;
  int rows=0, columns=0, max_columns=0;
  double **value=NULL;


  /* read file to two-dimensional double array value */
  status = ASCII_file2double (filename,
			      &rows, &max_columns, &columns,
			      &value);

  if (status != 0) 
    return status;


  /* check if a rectangular matrix */  
  if (columns != max_columns)  { 
    ASCII_free_double (value, rows);
    fprintf (stderr, " ... read_stamnes_table(): %s is not a rectangular matrix!\n", filename);
    fprintf (stderr, " ... minimum number of columns: %d, maximum number of columns: %d\n", 
	     columns, max_columns);

    return NOT_A_VALID_STAMNES_TABLE;
  }
    

  /* allocate memory for Stamnes table */
  *table = calloc (1, sizeof(STAMNES_TABLE));

  
  /* first row/column contains ratio/zenith */
  (*table)->n_zenith = columns-1;
  (*table)->n_ratio  = rows-1;
 
  /* allocate memory for arrays */
  (*table)->zenith = (double *)  calloc ((*table)->n_zenith, sizeof (double));
  (*table)->ratio  = (double *)  calloc ((*table)->n_ratio,  sizeof (double));
  (*table)->ozone  = (double **) calloc ((*table)->n_ratio,  sizeof (double *));

  for (i=0; i<(*table)->n_ratio; i++)
    (*table)->ozone[i] = (double *) calloc ((*table)->n_zenith,  sizeof (double));
  
  /* copy first row to array zenith */
  for (i=1; i<=(*table)->n_zenith; i++)  
    (*table)->zenith[i-1] = value[0][i];

  /* copy first column to array ratio */
  for (j=1; j<=(*table)->n_ratio; j++)  
    (*table)->ratio[j-1] = value[j][0];



  /* copy the remainder to 2dim double array ozone */
  for (i=1; i<=(*table)->n_ratio; i++) 
    for (j=1; j<=(*table)->n_zenith; j++)
      (*table)->ozone[i-1][j-1] = value[i][j];
    

  /* free memory of double array */
  ASCII_free_double (value, rows);
 
  return 0;  /* if o.k. */
}

/**************************************************************/
/* free memory of data read from Stamnes table.               */
/**************************************************************/

void free_stamnes_table (STAMNES_TABLE *table)
{
  int i=0;

  for (i=0; i<table->n_ratio; i++) 
    free (table->ozone[i]);

  free(table->ozone);

  free(table->zenith);
  free(table->ratio);

  free(table);
}



/***************************************************************/
/* stamnes_calculate uses the data read from a Stamnes table   */
/* to calculate the ozone column for given zenith angle and    */
/* irradiance ratio.                                           */
/*                                                             */
/* zenith[0..n_zenith-1]  zenith angles found in Stamnes table */
/* ratio[0..n_ratio-1]    ratio's found in Stamnes table       */
/* ozone[0..n_ratio-1][0..n_zenith-1]   ozone columns found in */
/*   Stamnes table.                                            */
/* meas_zenith   zenith angle to be used for calculation       */
/* meas_ratio    ratio to be used for calculation              */
/* column                                                      */
/***************************************************************/

int stamnes_calculate (STAMNES_TABLE *table,
		       double meas_zenith, double meas_ratio,
		       double *column)
{
  int counter=0, status=0;
  double *a0=NULL, *a1=NULL, *a2=NULL, *a3=NULL;
  double *b0=NULL, *b1=NULL, *b2=NULL, *b3=NULL;
  double lower=0, upper=0;

  /* reset ozone column */
  *column = -1;

  /* look for neighbour points in  array of ratios */
  while (counter < table->n_ratio && table->ratio[counter] < meas_ratio)
    counter++;

  /* check if meas_ratio within Stamnes table limits */ 
  if (counter<=0 || counter >= table->n_ratio)
    return RATIO_OUT_OF_RANGE;


  /* calculate spline coefficients for ozone as function of zenith */
  status = spline_coeffc (table->zenith, table->ozone[counter-1], table->n_zenith, 
			  &a0, &a1, &a2, &a3);

  if (status!=0)  {
    free(a0);
    free(a1);
    free(a2);
    free(a3);
    return ERROR_SPLINING_ZENITH;
  }
    

  status = spline_coeffc (table->zenith, table->ozone[counter], table->n_zenith, 
			  &b0, &b1, &b2, &b3);
  
  if (status!=0)  {
    free(a0);
    free(a1);
    free(a2);
    free(a3);
    free(b0);
    free(b1);
    free(b2);
    free(b3);
    return ERROR_SPLINING_ZENITH;
  }

  /* calculate interpolated data */
  status =  calc_splined_value (meas_zenith, &lower, 
				table->zenith, table->n_zenith, 
				a0, a1, a2, a3);

  free(a0);
  free(a1);
  free(a2);
  free(a3);
  

  if (status!=0)  {
    free(b0);
    free(b1);
    free(b2);
    free(b3);
    return ERROR_SPLINING_ZENITH;
  }

  status =  calc_splined_value (meas_zenith, &upper, 
				table->zenith, table->n_zenith, 
				b0, b1, b2, b3);

  free(b0);
  free(b1);
  free(b2);
  free(b3);


  if (status!=0) 
    return ERROR_SPLINING_ZENITH;


  /* now interpolate linearly for ratio */
  *column = lower + (upper - lower) / (table->ratio[counter] - table->ratio[counter-1]) *
    (meas_ratio - table->ratio[counter-1]);

  return 0;  /* if o.k. */
}

double cloud_clear_ratio (double **table,
			  double meas_zenith, double meas_ratio, int rows, double *cratio)
{
  int status = 0, i;
  double *a0=NULL, *a1=NULL, *a2=NULL, *a3=NULL;
  double *x=NULL, *y=NULL;

  if (meas_zenith < table[0][0] || meas_zenith > table[rows-1][0]) {
    fprintf (stderr, "Solar zenith angle %f out of range [%f - %f]\n", 
	     meas_zenith, table[0][0], table[rows-1][0]);
    return -1;
  }

  x = calloc (rows, sizeof(double));
  y = calloc (rows, sizeof(double));
  for (i=0;i<rows;i++) {
    x[i] = table[i][0];
    y[i] = table[i][1];
  }


  /* calculate spline coefficients for clear sky irradiance as function of zenith */
  status = spline_coeffc (x, y, rows, 
			  &a0, &a1, &a2, &a3);

  if (status!=0)  {
    free(a0);
    free(a1);
    free(a2);
    free(a3);
    return -1;
  }

  /* calculate interpolated data */
  status =  calc_splined_value (meas_zenith, cratio, 
				x, rows, 
				a0, a1, a2, a3);
  if (status!=0)  {
    free(a0);
    free(a1);
    free(a2);
    free(a3);
    return -1;
  }
  
  free(a0);
  free(a1);
  free(a2);
  free(a3);
  free(x);
  free(y);

  return 0;  /* if o.k. */
}

void panic()   /* ;-) */
{
  fprintf (stderr, "Help!\n");
  exit(-1);
}

int main(int argc, char **argv) 
{
  int status=0;
  int c=0;
  int print =1, rat=0, zen=0, clear=0;
  int day_of_year = 0;
  char filename[FILENAME_MAX+200] = "";
  char filename_clear[FILENAME_MAX+200] = "";
  char *dummy=NULL;

  double meas_ratio=0.0;  /* measured ratio 340 / 305 */ 
  double meas_zenith=0.0; /* solar zenith angle       */ 
  double column=0.0;
  double eccentricity_factor=1.0, eccentricity_day_one=1.0;
  int  crows=0, cmin_columns=0, cmax_columns=0;
  double cratio=0;
  double** clear_table=NULL;
	
  STAMNES_TABLE *table = NULL;

    /* check number of remaining command line arguments */

  while ((c=getopt (argc, argv, "cd:pr:s:?h")) != EOF)  {
    switch(c)  {
    case 'c':  /* cloud/clear sky factor */  
      clear = 1;
      break;   
    case 'd':  /* day of year */  
      day_of_year  = atoi(optarg);
      rat = 1;
      break;   
    case 'p':  /* no print messages */  
      print = 0;
      break;   
    case 'r':  /* no print messages */  
      meas_ratio  = strtod (optarg, &dummy);
      rat = 1;
      break;   
    case 's':  /* no print messages */  
      meas_zenith = strtod (optarg, &dummy);
      zen = 1;
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
  if (argc - optind != 1)  {
    fprintf(stderr,"No Stamnes table file specified\n"); 
    print_usage (argv[0]);
    return -1;
  }
  else if (!zen)  {
    fprintf(stderr,"No solar zenith angle given\n"); 
    print_usage (argv[0]);
    return -1;
  }
  else if (!rat)  {
    fprintf(stderr,"No measured irradiance given\n"); 
    print_usage (argv[0]);
    return -1;
  }
  
  /* Correct for Earth-Sun distance, not used for ozone, but
     important for water optical thickness */
 
  if (day_of_year > 0) {
    eccentricity_day_one = eccentricity (1);   /* The wc optical depth table is generated for day 1. */
    eccentricity_factor = eccentricity (day_of_year); 
    /*    fprintf(stderr,"ecc: %f %f %f %f %f\n", 
	   eccentricity_day_one, eccentricity_factor, eccentricity_day_one/eccentricity_factor,
	   meas_ratio, meas_ratio*eccentricity_day_one/eccentricity_factor); */
    meas_ratio = meas_ratio*eccentricity_day_one/eccentricity_factor;
  }

  strcpy (filename, argv[optind]);
  status = read_stamnes_table (filename, &table);
  /*  status = read_stamnes_table ("stamnes.dat", &table);*/


  if (status!=0) {
    switch (status) {
    case ASCIIFILE_NOT_FOUND:
      fprintf(stderr,"Stamnes table file not found %s\n",filename);
      break;
    default:
      fprintf(stderr,"Problem reading Stamnes table file %s\n",filename);
      break;
    }
    exit(-1);
  }

  status =  stamnes_calculate (table, meas_zenith, meas_ratio,
			       &column);

  if (status!=0)
    panic();  

  if (clear) {
    strcpy (filename_clear, filename);
    strcat(filename_clear, ".CLEAR");
    /* read clear sky file */
    status = ASCII_file2double (filename_clear,
				&crows, &cmax_columns, &cmin_columns, 
				&clear_table);
    status = cloud_clear_ratio (clear_table, meas_zenith, meas_ratio, crows, &cratio);
    if (status!=0)
      panic();  
  }


  if (print)
    fprintf (stdout, "The value interpolated from the table is %g\n", column);

  if (clear)
    fprintf (stdout, "%g %g\n", column, meas_ratio/cratio);
  else
    fprintf (stdout, "%g\n", column);

  free_stamnes_table (table);

  return 0;
}
