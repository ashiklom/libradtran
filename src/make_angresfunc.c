/*--------------------------------------------------------------------
 * $Id: make_angresfunc.c 2623 2011-12-23 10:52:38Z robert.buras $
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

  int c=0;
  double resolution=0.0;

  double ang, ang_start, ang_end, func=0, x;
  int i, n_points=0, type=1, angunit=1;

  while ((c=getopt (argc, argv, "hr:a:t:")) != EOF)  {
    switch(c)  {
    case 'r': 
      resolution = atof(optarg);
      break;
    case 't': 
      type = atoi(optarg);
      break;
    case 'a': 
      angunit = atoi(optarg);
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
   
  ang_start = 180;
  ang_end   = 0;
  n_points     = (int) (abs(ang_end-ang_start)/resolution) + 1;
  for (i=0;i<n_points;i++) {
    ang    = ang_start - i*resolution;
    switch (type) {
    case 1: 
      func   = (double) cos(PI*ang/180.);
      break;
    case 2: 
    case 3: 
      func   = (double) 1;
      break;
    }
    if (ang < 90 && type != 3) func = 0;
    x = ang;
    if (angunit ==2) 
      x=cos(PI*ang/180.);
    printf("%10.6f %10.6f\n", x, fabs(func));
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
  \subsection{Angular response function - \codeidx{make\_angresfunc}}

   The \code{make\_angresfunc} tool calculates various angular response
   functions to be used by for example the \code{angres} tool. All
   output is to stdout in two column format. The first column is the angle and
   the second column contains the corresponding value for a given angular response.
   The output angles follow \code{disort} conventions.

   The \code{make\_angresfunc} tool is invoked on the command line as

\begin{Verbatim}[fontsize=\footnotesize]
    make_angresfunc [-hart]
\end{Verbatim}
where the various options are
\begin{description}
   \item[-t] type of angular response
      \begin{enumerate}
         \item cosine (default)
         \item 2pi actinic flux
         \item 4pi actinic flux
      \end{enumerate}
   \item[-a] angular output format
      \begin{enumerate}
         \item angles (default)
         \item cosine of angle
      \end{enumerate}
   \item[-r] resolution, in degrees
   \item[-h] Print help message.
\end{description}

The following invocation of \code{make\_angresfunc} calculates the
angular response for a perfect cosine detector. The output is found in the 
\file{examples/ANGRES_1_ANG.DAT}.

\begin{Verbatim}[fontsize=\footnotesize]
    make_angresfunc -t 1 -r 1
\end{Verbatim}


  </lpdoc>

*/
static void usage()
{
  fprintf (stderr, "\nOutputs to stdout the angle and \n");
  fprintf (stderr, "corresponding value for a given angular response.\n");
  fprintf (stderr, "The output angles follow disort conventions.\n");
  fprintf (stderr, "\n");
  fprintf (stderr, "All options must be set.\n\n");
  fprintf (stderr, "No error checking yet, so be careful.\n\n");

  fprintf (stderr, "Usage: make_angresfunc [-hart]\n");
  fprintf (stderr, " -t type of angular response\n");
  fprintf (stderr, "    1: cosine (default)\n");
  fprintf (stderr, "    2: 2pi actinic flux\n");
  fprintf (stderr, "    3: 4pi actinic flux\n");
  fprintf (stderr, " -a angular output format\n");
  fprintf (stderr, "    1: angles (default)\n");
  fprintf (stderr, "    2: cosine of angle\n");
  fprintf (stderr, " -r resolution, in degrees\n");
  fprintf (stderr, " -h    Print this message.\n");  
}
