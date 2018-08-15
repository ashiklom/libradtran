/************************************************************************
 * $Id: alis.h 2772 2012-08-20 15:36:02Z robert.buras $
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

#ifndef _ALIS_H
#define _ALIS_H 1

#include "mystic.h"

int spectral_is_weight ( double           **totweight_spectral,
			 photon_struct     *p, 
			 photon_struct     *photon,
			 atmosphere_struct *atmos );

int concentration_is_weight ( double           **totweight_concentration,
			      photon_struct     *p,
			      photon_struct     *photon, 
			      atmosphere_struct *atmos );



#endif /* _ALIS_H */

