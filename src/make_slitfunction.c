/*--------------------------------------------------------------------
 * $Id: make_slitfunction.c 2623 2011-12-23 10:52:38Z robert.buras $
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
#include <math.h>
#include "getopt.h"

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif


static void usage();

int main(int argc, char **argv) 
{

  int c=0,nfwhm=4;
  double resolution=0.0, fwhm=0.0;

  double lambda, lambda_range, func, func_step, ln2, fact1, fact2;
  int i, n_points=0, type=1;

  while ((c=getopt (argc, argv, "hr:f:n:t:")) != EOF)  {
    switch(c)  {
    case 'r': 
      resolution = atof(optarg);
      break;
    case 't': 
      type = atoi(optarg);
      break;
    case 'f': 
      fwhm = atof(optarg);
      break;
    case 'n': 
      nfwhm = atof(optarg);
      break;
    case 'h':  /* Guess what */
      usage();
      return (-1);
      break;
    default:
      usage();
      return (-1);
    }
  }
   
  switch (type) {
  case 1: 
    lambda_range = 2*fwhm;
    n_points     = (int) (lambda_range/resolution) + 1;
    func_step    = 1./(n_points/2);
    for (i=0;i<n_points/2;i++) {
      lambda = -lambda_range/2. + i*resolution;
      func   = (double) i * func_step;
      printf("%8.3f %10.6f\n", lambda, func);
    }
    for (i=0;i<=n_points/2;i++) {
      lambda = i*resolution;
      func   = (double) (n_points/2 - i) * func_step;
      printf("%8.3f %10.6f\n", lambda, func);
    }
    break;
  case 2: 
    lambda_range = 2*fwhm;
    n_points     = (int) (lambda_range/resolution) + 1;
    func_step    = 1./(n_points/2);
    for (i=0;i<n_points/2;i++) {
      lambda = -lambda_range/2. + i*resolution;
      func   = (double) 1.0;
      printf("%8.3f %10.6f\n", lambda, func);
    }
    for (i=0;i<=n_points/2;i++) {
      lambda = i*resolution;
      func   = (double) 1.0;
      printf("%8.3f %10.6f\n", lambda, func);
    }
    break;
  case 3: 
    lambda_range = nfwhm*fwhm;
    n_points     = (int) (lambda_range/resolution) + 1;
    func_step    = 1./(n_points/2);
    ln2 = log(2.0);
    fact1 = 1.0; //2.0*fwhm*sqrt(ln2/PI);
    fact2 = fwhm*fwhm / (4.0*ln2);
    for (i=0;i<n_points/2;i++) {
      lambda = -lambda_range/2. + i*resolution;
      func   = fact1 * exp(-(lambda*lambda/fact2));
      printf("%8.3f %10.6f\n", lambda, func);
    }
    for (i=0;i<=n_points/2;i++) {
      lambda = i*resolution;
      func   = fact1 * exp(-(lambda*lambda/fact2));
      printf("%8.3f %10.6f\n", lambda, func);
    }
    break;
  }
  return 0;   /* if o.k. */
}

/*************************************************************/
/* print usage information                                   */
/*************************************************************/
/*
  Documentation in latex for the User's Guide. Extracted by sdoc2.awk 
  and input in the doc/tools.tex.
*/
/*
  <lpdoc>
  \subsection{Slit function generator - \codeidx{make\_slitfunction}}

  To generate standard slit functions to be used by \code{uvspec} the 
  \code{make\_slitfunction} tool may be used. For a given set of input
  it outputs to stdout in two column format the wavelength and corresponding value for the
  wanted slit function.

   The \code{make\_slitfunction} tool is invoked on the command line as
\begin{Verbatim}[fontsize=\footnotesize]
    make_angresfunc [-hrtf]
\end{Verbatim}
where the various options are
\begin{description}
   \item[-t] type of slitfunction
      \begin{enumerate}
         \item triangular (default)
         \item rectangular
         \item Gaussian
      \end{enumerate}
   \item[-f] full width at half maximum, in nm
   \item[-r] resolution, in nm
   \item[-n] number of fwhm (in nm) spanned by the slit function.
             Only applicable with Gaussian (type 3) slit function.
             Default value is 4.
   \item[-h] Print help message.
\end{description}

The following invocation of \code{make\_slitfunction} calculates the
a triangular slit function with FWHM of 0.75 nm and a resolution of
0.01 nm. The output is found in the 
\file{examples/TRI_SLIT.DAT}.

\begin{Verbatim}[fontsize=\footnotesize]
    make_slitfunction -f 0.75 -r 0.01 -t 1
\end{Verbatim}


  </lpdoc>

*/

static void usage()
{
  fprintf (stderr, "\nOutputs to stdout the wavelength and \n");
  fprintf (stderr, "corresponding value for the wanted slit function.\n");
  fprintf (stderr, "All options must be set.\n\n");
  fprintf (stderr, "No error checking yet, so be careful.\n\n");

  fprintf (stderr, "Usage: make_slitfunction [-hrft]\n");
  fprintf (stderr, " -t type of slitfunction\n");
  fprintf (stderr, "    1: triangular (default)\n");
  fprintf (stderr, "    2: rectangular\n");
  fprintf (stderr, "    3: gaussian\n");
  fprintf (stderr, " -f full width at half maximum, in nm\n");
  fprintf (stderr, " -r resolution, in nm\n");
  fprintf (stderr, " -n number of fwhm (in nm) spanned by the slit function.\n");
  fprintf (stderr, "    Only applicable with Gaussian (type 3) slit function.\n");
  fprintf (stderr, "    Default value is 4.\n");
  fprintf (stderr, " -h    Print this message.\n");  
}
