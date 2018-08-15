/*--------------------------------------------------------------------
 * $Id: snow.c 2792 2012-08-31 08:41:47Z svn-kylling $
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
#include <float.h>

#include <getopt.h>
#include "ascii.h"
#include "numeric.h"
#include "table.h"
#include "f77-uscore.h"
#include "solver.h"


/***********************************************************************************/
/* Function: snowalbedo                                                   @69_30i@ */
/* Description:                                                                    */
/*  Calculate the diffuse and direct albedo as formulated by Wiscombe and Warren,  */
/*  Journal of the Atmospheric Sciences, vol, 37, 2712-2733, 1980. Equation        */
/*  numbers below refer to equations in their paper.                               */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double omega:      snow single scattering albedo                               */
/*  double tau:        snow optical depth                                          */
/*  double gg:         snow assymmetry factor                                      */
/*  double surface_albedo: albedo of underlying surface                            */
/*  double umu0:       cosine of solar zenith angle                                */
/*  double albedo_diffuse: diffuse snow albedo, Eq. 6                              */
/*  double albedo_direct: direct snow albedo, Eq. 3                                */
/*                                                                                 */
/* Return value:                                                                   */
/*  None                                                                           */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling, arve.kylling@@nilu.no                                     */
/*                                                                        @i69_30@ */
/***********************************************************************************/

void snowalbedo (double omega, double tau, double gg, double surface_albedo, double umu0,
		 double* albedo_diffuse, double* albedo_direct) 
{

  double F77_FUNC (dei, DEI) (double* x);
  double omegas, taus, ggs;
  double as, bs, chi, PP, gamma, QQP, QQM, QQ;
  double arg, dei1, dei2, dei3;

  /* Delta-Eddington scaling */
  ggs       = gg/(1.+gg);                           /* Eq. 2a */
  omegas    = ((1.-gg*gg)*omega)/(1.-gg*gg*omega);  /* Eq. 2b */
  taus      = (1.-gg*gg*omega)*tau;                 /* Eq. 2c */

  /* Definitions following Eq. 3. */
  as        = 1-omegas*ggs;
  bs        = ggs/as;
  chi       = sqrt(3*as*(1.-omegas));
  PP        = 2.*chi/(3*as);
  gamma     = (1.-surface_albedo)/(1.+surface_albedo);
  QQM       = (gamma-PP)*exp(-chi*taus);
  QQP       = (gamma+PP)*exp(+chi*taus);
  QQ        = (1.+PP)*QQP - (1.-PP)*QQM;
    

  /* For some combinations of input variables (large -l and -r) QQP and QQ overflow. */
  /* Set to maximum double value. Arve Kylling 11062010.                             */
  if (QQP >= DBL_MAX) QQP = DBL_MAX*1e-100;
  if (QQ  >= DBL_MAX) QQ  = DBL_MAX*1e-100;
  //aky  if (QQP >= HUGE_VAL) QQP = DBL_MAX*1e-100;
  //aky  if (QQ  >= HUGE_VAL) QQ  = DBL_MAX*1e-100;

  /* Direct snow albedo, Eq. 3 and 4*/
  if ( chi*taus > 7.0 ) 
    *albedo_direct = (omegas/(1+PP)) * ((1.-bs*chi*umu0)/(1.+chi*umu0));  
  else 
    *albedo_direct = (2*(PP*(1.-gamma+omegas*bs) + 
			 omegas*(1.+bs)*(gamma*chi*umu0-PP)/(1.-chi*chi*umu0*umu0))*exp(-taus/umu0) 
		      - omegas*bs*(QQP - QQM)
		      + omegas*(1.+bs)*(QQP/(1.+chi*umu0) - QQM/(1.-chi*umu0) )
		      )/QQ;


  /* Diffuse snow albedo, Eq. 6 */
  arg = -taus;    

  /* Added tests on all args function dei to avoid numerical overflow problems. Arve Kylling 11062010. */
  if (arg < -700 )                 dei1 = 0.0;
  else
    dei1 = F77_FUNC (dei, DEI) (&arg);

  arg = -(1.+chi)*taus;
  if (arg < -700 )                 dei2 = 0.0;
  else
    dei2 = F77_FUNC (dei, DEI) (&arg);

  arg = -(1.-chi)*taus;
  if (arg < -700 )                 dei3 = 0.0;
  else
    dei3 = F77_FUNC (dei, DEI) (&arg);

  /* Sometimes dei3 overflow.                                                        */
  /* Set to maximum double value. Arve Kylling 11062010.                             */
  if (dei3 >= HUGE_VAL) dei3 = DBL_MAX;

  if ( chi*taus > 7.0 ) 
    /* Diffuse snow albedo, semi-infinite approximation. Eq. 7 */
    *albedo_diffuse = (2.*omegas/(1+PP))*(((1.+bs)/(chi*chi))*(chi-log(1+chi)) -bs/2.);
  else 
    /* Diffuse snow albedo, Eq. 6 */
    *albedo_diffuse = (2*PP*( (1.-gamma+omegas*bs)*(1.-taus)-(gamma*omegas)*(1.+bs)/(1.-omegas))*exp(-taus)
		       - 2*PP*( omegas*(1.+bs)*(2./(chi*chi)+(gamma*taus)/(1-omegas))+
				(1.-gamma+omegas*bs)*taus*taus)*dei1
		       + 2.*omegas*(1.+bs)/(chi*chi)*(QQP*(dei2 + chi - log(1.+chi)) -
						      QQM*(dei3 - chi - log(fabs(1.-chi))))
		       - omegas*bs*(QQP-QQM)
		       )/QQ;



}


