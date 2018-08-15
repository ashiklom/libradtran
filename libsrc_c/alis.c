/************************************************************************
 * $Id: alis.c 3072 2014-10-06 10:58:52Z Claudia.Emde $
 *
 * MYSTIC - Monte Carlo code for the physically correct tracing of
 *          photons in cloudy atmospheres.
 *
 * Copyright (c) 2000-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * Correspondence: bernhard.mayer@lmu.de
 *
 ************************************************************************/

#include <math.h>

#include "alis.h"
#include <stdlib.h>

/***********************************************************************************/
/* Function: spectral_is_weight                                           @62_30i@ */
/* Description: Calculate total absorption weight for spectral importance          */
/*              sampling (corresponing to Eq. 14 of ALIS paper (Emde et al., 2011).*/
/*                                                                                 */
/*                                                                                 */
/* Parameters:  Output:                                                            */
/*              totweight_spectral     spectral weight                             */
/*                                                                                 */ 
/*              Input:                                                             */
/*              p                      photon structure of local estimate photon   */
/*              photon                 photon structure of "real" photon           */
/*              atmos                  atmosphere structure                        */
/* Known bugs:                                                                     */
/* Author: Claudia Emde                                                            */
/* Date: 2011-07-20                                                                */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int spectral_is_weight ( double           **totweight_spectral,
			 photon_struct     *p,
			 photon_struct     *photon,
			 atmosphere_struct *atmos )
{
  double static *tau_spectral=NULL;
  int iv=0, kc=0; 
  
  if (tau_spectral==NULL)
    tau_spectral = calloc((size_t) atmos->nlambda_abs, sizeof(double));

  for (iv=0; iv<atmos->nlambda_abs; iv++){
    tau_spectral[iv]=0.0;
    
    for (kc=0; kc<atmos->Nz; kc++)
      /* Spectral absorption optical thickness */
      tau_spectral[iv]+= p->pathlength_per_layer[kc] *
        ( atmos->kabs_spectral[MCCAOTH_TOT][iv][kc] - 
          (atmos->kabs->prof [MCCAOTH_TOT])[kc]
	  + atmos->ksca_spectral[MCCAOTH_TOT][iv][kc] - 
          (atmos->ksca[MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->prof [MCCAOTH_TOT])[kc]);
  }
  
  /* XXX RPB: this has nothing to do with NEW_REFLECT, but the correct thing
     here is p instead of photon. However, I did not dare to commit without
     more testing,*/ 
  /* ???? CE: I think only for albedo weight p is needed */
  for (iv=0; iv<atmos->nlambda_abs; iv++)
    (*totweight_spectral)[iv] = exp(-tau_spectral[iv])*photon->q_spectral[iv]*
      photon->q2_spectral[iv] 
      * p->q_albedo_spectral[iv];
  
  return 0; 
}

/***********************************************************************************/
/* Function: concentration_is_weight                                      @62_30i@ */
/* Description: Calculate total absorption weight for spectral importance          */
/*              sampling (corresponing to Eq. 14 of ALIS paper (Emde et al., 2011).*/
/*                                                                                 */
/*                                                                                 */
/* Parameters:  Output:                                                            */
/*              totweight_concentration  concentration weight                      */
/*                                                                                 */ 
/*              Input:                                                             */
/*              p                      photon structure of local estimate photon   */
/*              photon                 photon structure of "real" photon           */
/*              atmos                  atmosphere structure                        */
/* Known bugs:                                                                     */
/* Author: Claudia Emde and Marius Schmidl                                         */
/* Date: 2012-03-06                                                                */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int concentration_is_weight ( double           **totweight_concentration,
			      photon_struct     *p,
			      photon_struct     *photon, 
			      atmosphere_struct *atmos )
{
  double tauabs_concentration=0.0, tausca_concentration=0.0; 
  int ic=0, kc=0; 

  for (ic=0; ic<atmos->Nc; ic++){
    tauabs_concentration=0.0;
    tausca_concentration=0.0;
    
    for (kc=0; kc<atmos->Nz; kc++){
      
      /* Concentration absorption optical thickness */
      tauabs_concentration+= p->pathlength_per_layer[kc] *
        ( - (atmos->kabs->prof [MCCAOTH_AER])[kc] + atmos->kabs_scaled[ic][kc]);
      
      /* Concentration optical thickness, only applied for aerosol scattering */
      tausca_concentration+= p->pathlength_per_layer[kc] *
        (  -(atmos->ksca[MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL] ->prof [MCCAOTH_AER])[kc]
           + atmos->ksca_scaled[ic][kc]);
    }
    
    
    /* Total weight for concentration importance sampling: */
    
   /*  fprintf(stderr, "ic %d tausca_concentration %g tauabs_concentration %g *photon->q_concentration[ic] %g *photon->q2_concentration[ic] %g\n", ic, tausca_concentration, tauabs_concentration, photon->q_concentration[ic], p->q2_concentration[ic]); */
    
    (*totweight_concentration)[ic] = 
      exp(-tauabs_concentration-tausca_concentration)*photon->q_concentration[ic]*p->q2_concentration[ic];
    
  }
  return 0; 
}
