/************************************************************************
 * $Id: ambralsfor.h 2623 2011-12-23 10:52:38Z robert.buras $
 ************************************************************************/

#ifndef __ambralsfor_h__
#define __ambralsfor_h__

#include <math.h>
#include "f77-uscore.h"

#define D2R(p) M_PI*(p)/180.    /* degrees->radians */
#define BR 1.0                  /* LiSparse b/r */
#define HB 2.0                  /* LiSparse h/b */

double ambrals_brdf (double iso, double vol, double geo, 
		     double mu1, double mu2, double phi);

void F77_FCN (ambralsfort) (float *iso, float *vol, float *geo,
			    float *mu1, float *mu2, float *phi, float *bdref);


#endif
