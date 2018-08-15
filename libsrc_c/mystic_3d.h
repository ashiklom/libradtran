/************************************************************************
 * $Id: mystic_3d.h 2731 2012-07-10 11:29:51Z robert.buras $
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

#ifndef _MYSTIC_3D_H
#define _MYSTIC_3D_H 1

#include "mystic.h"

int step3D ( photon_struct     *p,
	     atmosphere_struct *atmos,
	     double             step,
	     int                bcond,
	     int                photonpath,
	     int                spherical3D,
	     int                visualize );

int intersection3D ( photon_struct     *p,
		     atmosphere_struct *atmos,
		     double             tau,
		     double             tausca,
		     double            *step);



#endif /* _MYSTIC_3D_H */

