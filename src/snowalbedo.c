/*--------------------------------------------------------------------
 * $Id: snowalbedo.c 2792 2012-08-31 08:41:47Z svn-kylling $
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
  @17c@                
 */

/*
<lpdoc>
\subsection{Calculate albedo of snow - \codeidx{Gen\_snow\_tab}, \codeidx{snowalbedo}}

The \code{Gen\_snow\_tab.pl} script and the \code{snowalbedo} program 
may be used to calculate the diffuse and direct albedo of snow as 
formulated by \citet{Warren1980}. First a table of various snow optical 
properties must be generated. This is done by the Perl\code{Gen\_snow\_tab.pl} 
script. The resulting tables will be read by the \code{snowalbedo} 
program which will calculate the wanted surface albedo quantities.

Generating the tables by the \code{Gen\_snow\_tab.pl} script is straightforward 
as the script only takes one argument, namely the name of the file body (It will
also print a small help message if \code{--help} is given to it). The script
will generate three files with extensions .gg, .qext and .ssa.
  \begin{Verbatim}[fontsize=\footnotesize]
    perl Gen_snow_tab.pl --file <name>
  \end{Verbatim}

The generated tables is read by the \code{snowalbedo} program which requires 
the following options:
\begin{description}
   \item[-l]            Equivalent depth of liquid water in snowpack (g cm -2)
   \item[-r]            mean grain radius ($\mu$m)
   \item[-u]            cosine of solar zenith angle
\end{description}
The options below are optional

\begin{description}
   \item[-a]            albedo of underlying surface, default 0.03
   \item[-p]            turn of printing of messages
   \item[-h] Print help message.
\end{description}

A typical usage of \code{snowalbedo} is 
(\code{Gen\_snow\_tab.pl --file ../examples/MIE\_ICE\_TAB} has been executed first)
  \begin{Verbatim}[fontsize=\footnotesize]
      snowalbedo ../examples/MIE_ICE_TAB -l 0.05 -r 50 -u 0.5 -p
  \end{Verbatim}
This will produce the following output (only two first output lines shown)
  \begin{Verbatim}[fontsize=\footnotesize]
      290.0 2.00893  0.9999776000 0.88037   0.9728   0.9689
      291.0 2.01212  0.9999782400 0.88064   0.9731   0.9693
  \end{Verbatim}
Here, the various columns have the following content
\begin{enumerate}
   \item  wavelength (nm)
   \item  Q\_ext
   \item  Single scattering albedo
   \item  Asymmetry parameter
   \item  Direct albedo
   \item  Diffuse albedo
\end{enumerate}


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

#include "snow.h"
#include "ascii.h"
#include "numeric.h"
#include "table.h"
#include "f77-uscore.h"

#define PROGRAM "snowalbedo"
#define VERSION "1.0"

static void print_usage (char *name)
{
  fprintf (stderr, "\nusage:   %s MIE_ICE.tab [options]\n\n", name);
  fprintf (stderr, "%s  calculates the diffuse and direct albedo of snow as\n", name);
  fprintf (stderr, "formulated by Wiscombe and Warren, Journal of the\n");
  fprintf (stderr, "Atmospheric Sciences, vol, 37, 2712-2733, 1980.\n\n");
  fprintf (stderr, "Required arguments:\n");
  fprintf (stderr, "  -l            Equivalent depth of liquid water in snowpack (g cm -2)\n");
  fprintf (stderr, "  -r            mean grain radius (um)\n");
  fprintf (stderr, "  -u            cosine of solar zenith angle\n");
  fprintf (stderr, "Optional arguments:\n");
  fprintf (stderr, "  -a            albedo of underlying surface, default 0.03\n");
  fprintf (stderr, "  -h            prints this message\n");
  fprintf (stderr, "  -?            prints this message\n");
  fprintf (stderr, "  -p            turn of printing of messages\n");
  fprintf (stderr, "\n");
  fprintf (stderr, "The MIE_ICE.tab.?? files must first be generated using\n");
  fprintf (stderr, "the Gen_snow_tab.pl perl script. The use of Gen_snow_tab.pl\n");
  fprintf (stderr, "is shown when executing perl Gen_snow_tab.pl\n");
  fprintf (stderr, "Note that Gen_snow_tab.pl generates three files\n");
  fprintf (stderr, "with extensions .gg, .qext and .ssa.\n");
  fprintf (stderr, "snowalbedo should just be fed the name of the common \n");
  fprintf (stderr, "body of the three files, see example below \n");
  fprintf (stderr, "\n");
  fprintf (stderr, "For example  %s -p ../examples/MIE_ICE.tab -l 0.05 -r 50 -u 0.5 -p \n", PROGRAM);
  fprintf (stderr, "will create: (only two first output lines shown)\n");
  fprintf (stderr, "290.0 2.00893  0.9999776000 0.88037   0.9728   0.9689\n");
  fprintf (stderr, "291.0 2.01212  0.9999782400 0.88064   0.9731   0.9693\n");
  fprintf (stderr, "\n");
  fprintf (stderr, "Here,\n");
  fprintf (stderr, "1. Column: wavelength (nm)\n");
  fprintf (stderr, "2. Column: Q_ext\n");
  fprintf (stderr, "3. Column: Single scattering albedo\n");
  fprintf (stderr, "4. Column: Asymmetry parameter\n");
  fprintf (stderr, "5. Column: Direct albedo\n");
  fprintf (stderr, "6. Column: Diffuse albedo\n");
}




int main (int argc, char **argv)
{
  char filename[FILENAME_MAX+200] = "";
  char *file_qext, *file_ssa, *file_gg;
  int print =1;
  int status=0;
  int iv=0;
  int c;
  char *dummy=NULL;
  double omega, tau, gg, lambda; 
  double surface_albedo, umu0; 
  double albedo_diffuse;
  double albedo_direct;
  double LWD, Qext, r_mean, rho_ice;

  TABLE *table_qext = NULL, *table_ssa = NULL, *table_gg = NULL;

  rho_ice = 0.917; /* The density of ice, in g cm-3. */
  LWD     = 10.0;  /* Equivalent depth of liquid water in snowpack, g cm-2. */
  Qext    = 2.01;
  r_mean  = 0.050;  /* Mean snow grain radius. */
  surface_albedo = 0.03;
  umu0 = 0.5;

  while ((c=getopt (argc, argv, "a:c?hl:pr:u:")) != EOF)  {
    switch(c)  {
    case 'a':   
      surface_albedo = strtod (optarg, &dummy);
      break;
    case 'l':   
      LWD = strtod (optarg, &dummy);
      break;
    case 'r':   
      r_mean = strtod (optarg, &dummy);
      break;
    case 'u':   
      umu0 = strtod (optarg, &dummy);
      break;
    case 'p':  /* no print messages */  
      print = 0;
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
      
  /* check number of remaining command line arguments */
  if (argc - optind != 1)  {
    fprintf(stderr,"No table file specified\n"); 
    print_usage (argv[0]);
    return -1;
  }
  
  strcpy (filename, argv[optind]);
  file_qext = (char *) calloc (strlen(filename)+6, sizeof (char));
  strcpy (file_qext, filename);
  file_qext = strncat (file_qext, ".qext", 5);

  file_ssa  = (char *) calloc (strlen(filename)+5, sizeof (char));
  strcpy (file_ssa, filename);
  file_ssa  = strncat (file_ssa, ".ssa", 4);

  file_gg   = (char *) calloc (strlen(filename)+4, sizeof (char));
  strcpy (file_gg, filename);
  file_gg   = strncat (file_gg, ".gg", 3);

  status = read_table (file_qext, &table_qext);
  if (status!=0) {
    switch (status) {
    case ASCIIFILE_NOT_FOUND:
      fprintf(stderr,"Table file not found %s\n",file_qext);
      break;
    default:
      fprintf(stderr,"Problem reading table file %s\n",file_qext);
      break;
    }
    exit(-1);
  }
  status = read_table (file_ssa, &table_ssa);
  if (status!=0) {
    switch (status) {
    case ASCIIFILE_NOT_FOUND:
      fprintf(stderr,"Table file not found %s\n",file_ssa);
      break;
    default:
      fprintf(stderr,"Problem reading table file %s\n",file_ssa);
      break;
    }
    exit(-1);
  }
  status = read_table (file_gg, &table_gg);
  if (status!=0) {
    switch (status) {
    case ASCIIFILE_NOT_FOUND:
      fprintf(stderr,"Table file not found %s\n",file_gg);
      break;
    default:
      fprintf(stderr,"Problem reading table file %s\n",file_gg);
      break;
    }
    exit(-1);
  }

  /* Output albedo for all wavelengths */
  
  for (iv=0;iv<(table_qext)->n_yy;iv++) {
    lambda = (table_qext)->yy[iv];
    status =  table_calculate (table_qext, r_mean, lambda, &Qext);
    status =  table_calculate (table_ssa,  r_mean, lambda, &omega);
    status =  table_calculate (table_gg,   r_mean, lambda, &gg);

    /* Factor 1.e-06 converts from um to cm for r_mean */
    tau   = (3.*LWD*Qext)/(4.*r_mean*1.e-06*rho_ice);
    
    snowalbedo (omega, tau, gg, surface_albedo, umu0, &albedo_diffuse, &albedo_direct);
    
    if (print)
      fprintf (stderr, "Direct albedo:  %8.4f\nDiffuse albedo: %8.4f\n",
	       albedo_direct, albedo_diffuse);
    else
      printf ("%7.1f %7.5f %13.10f %7.5f %8.4f %8.4f\n",
	      lambda, Qext, omega, gg, albedo_direct, albedo_diffuse);

  }
  return status;  
}
