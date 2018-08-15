/************************************************************************
 * $Id: sofi.h 3313 2017-12-11 14:27:51Z Claudia.Emde $
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

#ifndef __sofi_h
#define __sofi_h
#include "mystic.h"
int sample_photons_sofi (float wvnmlo, float wvnmhi, 
			  double *pd);

void generate_photon_sofi ( atmosphere_struct* atmos, photon_struct* p, 
                            double sza, double phi, double* pd );

#endif
