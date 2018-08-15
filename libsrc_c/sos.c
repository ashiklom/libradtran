/*--------------------------------------------------------------------
 * $Id: sos.c 2626 2011-12-28 15:22:24Z svn-kylling $
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
#include <time.h>
#include "ascii.h"
#include "f77-uscore.h"
#include "fortran_and_c.h"

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif

void F77_FUNC (qgausn, QGAUSN) (int* n, float* cmu, float*cwt);
void F77_FUNC (lepolys, LEPOLYS) (int* nn, int* mazim, int* mxcmu, 
		       int* nstr, float* cmu, float* ylmc);

float chapmanfac(int lc, int *nfac, float **fac, float *dtauc ) {

  /*     Calculates the Chapman-factor, Eq. B2 (DS).                     */

  /* I n p u t       v a r i a b l e s:                                  */
   
  /*      lc        : Computational layer                                */
  /*      dtauc     : Optical thickness of layer lc (un-delta-m-scaled)  */
  /*      fac       : Geomtrical correction factor, calculated by routine*/
  /*                  geocorfac.                                         */
  /*      nfac      : Number of terms in Eq. B2 (DS)                     */
   
  /* O u t p u t      v a r i a b l e s:                                 */
   
  /*      chapman   : Chapman-factor. In a pseudo-spherical atmosphere   */
  /*                  replace EXP( -tau/umu0 ) by EXP( -chapman ) in the */
  /*                  beam source.                                       */
   
  /* I n t e r n a l     v a r i a b l e s:                              */
   
  /*      fact      : =1 for first sum in Eq. B2 (DS)                    */
  /*                  =2 for second sum in Eq. B2 (DS)                   */
   
  float fact, sum=0;
  int j;  

  sum = 0.0;
  if ( nfac[lc] < 0 ) { sum = 1.e+20; return sum; }
  else {
    for (j=0;j<=nfac[lc];j++) {
      /* Include factor of 2 for zenang .gt. 90, second sum in Eq. B2 (DS) */
      fact = 1.0;
      if( j > lc ) {fact = 2.0;}
      sum += dtauc[j]*fact*fac[lc][j];
    }
    if ( nfac[lc] > lc ) {
      sum += dtauc[lc] * fac[lc][j+1];
    }
  }
  return sum;
}

int geocorfac(float **fac, int *nfac, float z_lay, int nlyr,
	      float *zd, float zenang, float r )  
{
  /* Calculates the geometric correction factor needed for spherical      */
  /* geometry. See Appendix B, (DS).                                      */

  /* I n p u t       v a r i a b l e s:                                   */

  /*      nlyr      : Number of layers in atmospheric model               */
  /*      zd(lc)    : lc = 0, nlyr. zd(lc) is distance from bottom        */
  /*                  surface to top of layer lc. zd(nlyr) = 0.0 km       */
  /*      zenang    : Solar zenith angle as seen from bottom surface      */ 
  /*      r         : Radius of earth. NOTE: Use the same dimension as zd,*/
  /*                  for instance both in km.                            */
  /*      z_lay     : Where in the layer the Chapman function is to be    */
  /*                  computed. E.g. 0.0, bottom of layer, 0.5, middle    */ 
  /*                  of layer and 1.0 top of layer.                      */

  /* O u t p u t      v a r i a b l e s:                                  */
     
  /*      fac        : geometrical correction factor                      */
  /*      nfac       : Number of terms in Eq. B2 (DS)                     */

  /* I n t e r n a l     v a r i a b l e s:                               */

  /*      dhj       : delta-h-sub-j in Eq. B2 (DS)                        */
  /*      dsj       : delta-s-sub-j in Eq. B2 (DS)                        */
  /*      fact      : =1 for first sum in Eq. B2 (DS)                     */
  /*                  =2 for second sum in Eq. B2 (DS)                    */
  /*      rj        : r-sub-j in Eq. B1 (DS)                              */
  /*      rjp1      : r-sub-j+1 in Eq. B1 (DS)                            */
  /*      xpsinz    : The length of the line OG in Fig. 1, (DS)           */

  float dsj, dhj, fact, rj, rjp1, xp, xpsinz, zenrad;
  int id, j, lc, status =0;

  zenrad = zenang * PI / 180.0;
  for (lc=0;lc<nlyr;lc++) {
    nfac[lc] = 0;
  }

  for (lc=0;lc<nlyr;lc++) {
    xp     = r +  zd[lc+1] + (zd[lc] - zd[lc+1] ) * z_lay;
    xpsinz = xp * sin( zenrad );
    if( (zenang > 90.0) && (xpsinz < r) ) {
      nfac[lc] = -1;
    }
    else {
      /*  Find index of layer in which the screening height lies */
      id = lc;
      if ( zenang > 90.0 ) {
	for (j=lc;j<nlyr;j++) {
	  if( (xpsinz < ( zd[j] + r ) ) &&
	      (xpsinz >= ( zd[j+1] + r )) ) {id = j;}
	}
      }
      for (j=0;j<=id;j++) {
	fact = 1.0;
	/* Include factor of 2 for zenang .gt. 90, second sum in Eq. B2 (DS)*/
	if( j > lc ) { fact = 2.0;}
	if(j==id && id==lc && zenang>90.0) { fact = -1.0;}
	rj = r + zd[j];
	rjp1 = r + zd[j+1];
	if(j==lc && id==lc) {rjp1 = xp;}
	dhj = zd[j] -zd[j+1];
	if (id>lc && j==id) { dsj = sqrt(rj*rj - xpsinz*xpsinz );}
	else { dsj = sqrt( rj*rj - xpsinz*xpsinz ) -
		 fact * sqrt( rjp1*rjp1 - xpsinz*xpsinz );}
	fac[lc][j] = dsj / dhj;
      }
      nfac[lc] = id;
     
      /* Third term in Eq. B2 (DS) */
      if( id > lc ) {
	dhj = zd[lc] -zd[lc+1];
	dsj = sqrt( xp*xp - xpsinz*xpsinz ) -
	  sqrt( (zd[lc+1]+r)*(zd[lc+1]+r) - xpsinz*xpsinz );
	fac[lc][j+1] = dsj / dhj;
      }
    }
  }
  return status;
}

void trans( int nlyr, float *chtau, float *trs) {
  int lc;
  for (lc=0;lc<=nlyr;lc++) {
    trs[lc] = exp(-chtau[lc]);
  }
}

void trans_double( int nlyr, double *chtau, double *trs) {
  int lc;
  for (lc=0;lc<=nlyr;lc++) {
    trs[lc] = exp(-chtau[lc]);
  }
}

void sinsca(float albedo, float* ssalb, float* dtauc, float umu0, float* trs,
	    float **pmom, float* cmu, int nlyr, int nstr, int numu, float* umu,
	    float **ylmc, float **ylm0, float **ylmu, float ***u0uc, float ***u0uu) {
  int iq, is, iu=0, lc, l, lu;
  float phase=0;

  is = 0;  /* First order scattering, stored at index 0 */

  /* First downwelling radiation */

  for (iq=0;iq<nstr/2;iq++) {
    lu = 0;
    u0uc[is][lu][iq] = 0;
    for (lc=0;lc<nlyr;lc++) {
      lu    = lc+1;
      phase = 0;
      for (l=0;l<nstr;l++) {
	phase += (2*l+1) * ylmc[iq][l] * ylm0[0][l] * pmom[lc][l];
      }
      /*      sum = (1/(4*PI))*(ssalb[lc]*umu0*phase/(-cmu[iq]-umu0))*trs[lc]*  */
      /*	(exp(dtauc[lc]/cmu[iq])-exp(-dtauc[lc]/umu0));                  */
      u0uc[is][lu][iq] = 
	u0uc[is][lu-1][iq]*exp(dtauc[lc]/cmu[iq])
	+(1/(4*PI))*(ssalb[lc]*umu0*phase/(-cmu[iq]-umu0))*trs[lc]*
	(exp(dtauc[lc]/cmu[iq])-exp(-dtauc[lc]/umu0));
      /*      printf("trs %2d %2d %7.5f %7.5f %7.5f %7.5f %8.5f %8.5f %7.5f \n", */
      /*	     iq, lc, sum, u0uc[is][lu][iq], trs[lc], dtauc[lc],          */
      /*	     cmu[iq], umu0, phase);                                      */
    }
  }

  /* Second upwelling radiation */

  for (iq=nstr/2;iq<nstr;iq++) {
    lu = nlyr;
    u0uc[is][lu][iq] = albedo*umu0*trs[nlyr]/PI;
    for (lc=nlyr-1;lc>=0;lc--) {
      lu    = lc;
      phase = 0;
      for (l=0;l<nstr;l++) {
	phase += (2*l+1) * ylmc[iq][l] * ylm0[0][l] * pmom[lc][l];
      }
      /*      sum = (1/(4*PI))*(ssalb[lc]*umu0*phase/(cmu[iq]+umu0))*trs[lc]* */
      /*	(1-exp(-dtauc[lc]*(1/umu0+1/cmu[iq])));                       */

      u0uc[is][lu][iq] = 
	u0uc[is][lu+1][iq]*exp(-dtauc[lc]/cmu[iq])
	+(1/(4*PI))*(ssalb[lc]*umu0*phase/(cmu[iq]+umu0))*trs[lc]*
	(1-exp(-dtauc[lc]*(1/umu0+1/cmu[iq])));

      /*      printf("trs %2d %2d %7.5f %7.5f %7.5f %7.5f %8.5f %8.5f %7.5f \n", */
      /*	     iq, lc, sum, u0uc[is][lu][iq], trs[lc], dtauc[lc],          */
      /*	     cmu[iq], umu0, phase);                                      */
    }
  }

  if ( numu > 0 ) {
    for (iu=0;iu<numu;iu++) {
      if ( umu[iu] < 0 ) { /* Downwelling user defined angles */
	lu = 0;
	u0uu[is][lu][iu] = 0;
	for (lc=0;lc<nlyr;lc++) {
	  lu    = lc+1;
	  phase = 0;
	  for (l=0;l<nstr;l++) {
	    phase += (2*l+1) * ylmu[iu][l] * ylm0[0][l] * pmom[lc][l];
	  }
	  u0uu[is][lu][iu] = 
	    u0uu[is][lu-1][iu]*exp(dtauc[lc]/umu[iu])
	    +(1/(4*PI))*(ssalb[lc]*umu0*phase/(-umu[iu]-umu0))*trs[lc]*
	    (exp(dtauc[lc]/umu[iu])-exp(-dtauc[lc]/umu0));
	}
      }
      else {              /* Upwelling user defined angles */
	lu = nlyr;
	u0uu[is][lu][iu] = albedo*umu0*trs[nlyr]/PI;
	for (lc=nlyr-1;lc>=0;lc--) {
	  lu    = lc;
	  phase = 0;
	  for (l=0;l<nstr;l++) {
	    phase += (2*l+1) * ylmu[iu][l] * ylm0[0][l] * pmom[lc][l];
	  }

	  u0uu[is][lu][iu] = 
	    u0uu[is][lu+1][iu]*exp(-dtauc[lc]/umu[iu])
	    +(1/(4*PI))*(ssalb[lc]*umu0*phase/(umu[iu]+umu0))*trs[lc]*
	    (1-exp(-dtauc[lc]*(1/umu0+1/umu[iu])));
	}
      }
    }    
  }


}

void morsca(float albedo, float* ssalb, float* dtauc, float *cmu, float *cwt, 
	     int nlyr, int nstr, int numu, float *umu, float **sf, float **sfu,
	    int is, float ***u0uc, float ***u0uu) {
  int iq, iu=0, lc, lu;
  float E_surface=0;

  /* First downwelling radiation */

  for (iq=0;iq<nstr/2;iq++) {
    lu = 0;
    u0uc[is][lu][iq] = 0;
    for (lc=0;lc<nlyr;lc++) {
      lu    = lc+1;
      u0uc[is][lu][iq] = u0uc[is][lu-1][iq]*exp(dtauc[lc]/cmu[iq])
	+0.5*ssalb[lc]*sf[lc][iq]*(1-exp(dtauc[lc]/cmu[iq]));

      /*      printf("trs %2d %2d %7.5f %7.5f %7.5f %7.5f %8.5f %8.5f %7.5f \n", */
      /*	     iq, lc, sum, u0uc[is][lu][iq], trs[lc], dtauc[lc],          */
      /*	     cmu[iq], umu0, phase);                                      */
    }
  }

  /* Second upwelling radiation */

  /* Need downward flux at surface from previous scattering event */
  E_surface = 0;
  for (iq=0;iq<nstr/2;iq++) {
    E_surface += cwt[iq]*(-cmu[iq])*u0uc[is-1][nlyr][iq];
  }
  E_surface *= 2*PI;

  for (iq=nstr/2;iq<nstr;iq++) {
    lu = nlyr;
    u0uc[is][lu][iq] = albedo*E_surface/PI;
    for (lc=nlyr-1;lc>=0;lc--) {
      lu    = lc;

      u0uc[is][lu][iq] = u0uc[is][lu+1][iq]*exp(-dtauc[lc]/cmu[iq])
	+0.5*ssalb[lc]*sf[lc][iq]*(1-exp(-dtauc[lc]/cmu[iq]));

      /*      printf("trs %2d %2d %7.5f %7.5f %7.5f %7.5f %8.5f %8.5f %7.5f \n", */
      /*	     iq, lc, sum, u0uc[is][lu][iq], trs[lc], dtauc[lc],          */
      /*	     cmu[iq], umu0, phase);                                      */
    }
  }

  if ( numu > 0 ) {
    for (iu=0;iu<numu;iu++) {
      if ( umu[iu] < 0 ) { /* Downwelling user defined angles */
	lu = 0;
	u0uu[is][lu][iu] = 0;
	for (lc=0;lc<nlyr;lc++) {
	  lu    = lc+1;
	  u0uu[is][lu][iu] = u0uu[is][lu-1][iu]*exp(dtauc[lc]/umu[iu])
	    +0.5*ssalb[lc]*sfu[lc][iu]*(1-exp(dtauc[lc]/umu[iu]));
	}
      }
      else {              /* Upwelling user defined angles */
	lu = nlyr;
	u0uu[is][lu][iu] = albedo*E_surface/PI;
	for (lc=nlyr-1;lc>=0;lc--) {
	  lu    = lc;
	  u0uu[is][lu][iu] = u0uu[is][lu+1][iu]*exp(-dtauc[lc]/umu[iu])
	    +0.5*ssalb[lc]*sfu[lc][iu]*(1-exp(-dtauc[lc]/umu[iu]));	  
	}
      }
    }
  }
}

int source_function(float ** pmom, float* cwt, int is, int nlyr, 
		     int nstr, int numu, float **ylmc, float **ylmu,
		    float ***u0uc, float **sf, float **sfu) {
  int iq, iu=0, jq, lc, l, status=0;
  float sum=0, phase=0;
  float **u0um=NULL;

  /* Get average radiation field for the layer */
  if ((status = ASCII_calloc_float (&u0um, nlyr, nstr)) != 0)
    return status;
  for (lc=0;lc<nlyr;lc++) {
    for (iq=0;iq<nstr;iq++) {
      u0um[lc][iq] = (u0uc[is][lc][iq]+u0uc[is][lc+1][iq])/2;
    }
  }

  /* Calculate the source function */
  for (lc=0;lc<nlyr;lc++) {
    for (iq=0;iq<nstr;iq++) {
      sum   = 0;
      for (jq=0;jq<nstr;jq++) {
	phase = 0;
	for (l=0;l<nstr;l++) {
	  phase += (2*l+1) * ylmc[iq][l] * ylmc[jq][l] * pmom[lc][l];
	}      
	sum   += phase*cwt[jq]*u0um[lc][jq];
      }
      sf[lc][iq] = sum;
    }
  }

  /* Calculate the source function at user angles*/
  if ( numu > 0 ) {
    for (lc=0;lc<nlyr;lc++) {
      for (iu=0;iu<numu;iu++) {
	sum   = 0;
	for (jq=0;jq<nstr;jq++) {
	  phase = 0;
	  for (l=0;l<nstr;l++) {
	    phase += (2*l+1) * ylmu[iu][l] * ylmc[jq][l] * pmom[lc][l];
	  }      
	  sum   += phase*cwt[jq]*u0um[lc][jq];
	}
	sfu[lc][iu] = sum;
      }
    }
  }

  if ( (status=ASCII_free_float (u0um, nlyr)) !=0 )
    return status;
  return status;
}



int sos (int nlyr, int newgeo, int nstr, int nscat, float albedo,
	 float radius, float *zd, float *ssalb, float **pmom,
	 float *dtauc, float zenang, int ntau, 
	 int numu, float *umu, float *utau,
	 float *rfldir, float *rfldn,  float *flup,
	 float *uavgso, float *uavgdn, float *uavgup, float **u0u)
{
  float cumtau=0;
  float umu0, **fac=NULL, *chtau=NULL, *tauc=NULL;
  float *cmu=NULL, *cwt=NULL, *tmpylmc, *tmpylmu,
    **ylmc=NULL, **ylmu=NULL, **ylm0=NULL, ***u0uc, ***u0uu, sgn, 
    **sf=NULL, **sfu=NULL;
  float *tmp;
  float *trs=NULL;
  int *layru=NULL, *nfac=NULL;
  int status=0, mazim, twonm1;
  int iq, is, iu=0, l, lc, lu, nn;

  /* Do some initialization stuff */
  umu0 = cos(zenang/180.0*PI);
  if ( newgeo ) {
    /* Calculate geometric correction factor needed for chapman function */
    /* Only need to do this if the geometry changes, that is, it is not  */
    /* wavelength dependent. */
    if ( (layru = (int *) calloc ((size_t) (ntau), sizeof (int))) == NULL )  
      return ASCII_NO_MEMORY;
    if ( (nfac = (int *) calloc ((size_t) (nlyr), sizeof (int))) == NULL )  
      return ASCII_NO_MEMORY;
    if ( (tauc = (float *) calloc ((size_t) (nlyr+1), sizeof (float))) == NULL )
      return ASCII_NO_MEMORY;
    if ((status = ASCII_calloc_float (&fac, nlyr, 2*nlyr)) != 0) return status;
    status = geocorfac( fac, nfac, 0.0, nlyr, zd, zenang, radius ) ;
    if (status!=0) {
      fprintf (stderr, "error %d during geocorfac\n", status);
      return status; 
    }
    /* Allocate for nscat so that index is=0 is first order of scattering etc. */
    if ((status = ASCII_calloc_float_3D (&u0uc, nscat, nlyr+1, nstr)) != 0)
      return status;
    if ((status = ASCII_calloc_float_3D (&u0uu, nscat, nlyr+1, numu)) != 0)
      return status;
    chtau = calloc ((size_t) (nlyr+1), sizeof(float));
    if ( (cmu = (float *) calloc ((size_t) (nstr), sizeof (float))) == NULL )
      return ASCII_NO_MEMORY;
    if ( (cwt = (float *) calloc ((size_t) (nstr), sizeof (float))) == NULL )
      return ASCII_NO_MEMORY;
    nn = nstr/2;
    F77_FUNC (qgausn, QGAUSN) (&nn,cmu,cwt);
    /* Rearrange to get cmu in ascending order */
    for (iq=0;iq<nn;iq++) { cmu[nn+iq] = cmu[iq]; cwt[nn+iq] = cwt[iq]; }
    for (iq=0;iq<nn;iq++) { cmu[iq] = -cmu[nstr-1-iq]; cwt[iq] = cwt[nstr-1-iq];}
    
    if ( (tmpylmc = (float *) calloc ((size_t) ((nstr+1)*nstr), sizeof (float))) == NULL ) 
      return ASCII_NO_MEMORY;
    mazim=0;
    twonm1=nstr-1;
    nn=1;
    if ( (tmp = (float *) calloc ((size_t) (1), sizeof (float))) == NULL )
      return ASCII_NO_MEMORY;
    tmp[0] = -umu0;
    F77_FUNC (lepolys, LEPOLYS) ( &nn, &mazim, &nstr, &twonm1, tmp, tmpylmc );
    free(tmp);
    ylm0 = fortran2c_2D_float_ary(nstr, nstr+1, tmpylmc);
    nn=nstr/2;
    F77_FUNC (lepolys, LEPOLYS) ( &nn, &mazim, &nstr, &twonm1, cmu, tmpylmc );
    ylmc = fortran2c_2D_float_ary(nstr, nstr+1, tmpylmc);
    free(tmpylmc);

    /* Evaluate Legendre polynomials with negative -cmu- from those with*/
    /* positive -cmu-;  Dave Armstrong Eq. (15) */
    sgn  = -1.0;
    for(l=0;l<nstr;l++) {
      sgn = -sgn;
      for (iq=nn;iq<nstr;iq++) {
	ylmc[iq][l] = sgn*ylmc[nstr-1-iq][l];
      }
    }    

    if ( numu > 0 ) {
      if ( (tmpylmu = (float *) calloc ((size_t) ((nstr+1)*numu), sizeof (float))) == NULL ) 
	return ASCII_NO_MEMORY;
      F77_FUNC (lepolys, LEPOLYS) ( &numu, &mazim, &numu, &twonm1, umu, tmpylmu );
      ylmu = fortran2c_2D_float_ary(numu, numu+1, tmpylmu);
      free(tmpylmu);
    }

  }

  /* Chapman factor */
  cumtau = 0;
  tauc[0]=0;
  chtau[0]=0;
  for (lc=0;lc<nlyr;lc++) {
    chtau[lc+1] = chapmanfac(lc, nfac, fac, dtauc);
    cumtau += dtauc[lc];
    tauc[lc+1] = cumtau;
  }

  /* Set arrays defining location of user output levels */
  for (lu=0;lu<ntau;lu++) {
    for (lc=0;lc<nlyr;lc++) {
      if ( utau[lu]>=tauc[lc] && utau[lu]<tauc[lc+1] )
	goto done; /* All good programs need at least one goto statement:))*/
    }
    lc = nlyr;
  done:
    layru[lu] = lc;
  }

  /* The real stuff for user defined order of scattering */
  
  /* Transmittance of the atmosphere  */
  if ( (trs = (float *) calloc ((size_t) (nlyr+1), sizeof (float))) == NULL )
    return ASCII_NO_MEMORY;
  trans( nlyr, chtau, trs );

  /* Single scattering */
  if (nscat >=1) {
    is = 0;
    sinsca(albedo, ssalb, dtauc, umu0, trs, pmom, cmu, 
	   nlyr, nstr, numu, umu, ylmc, ylm0, ylmu, u0uc, u0uu);
  }

  /* Higher order scattering */
  if (nscat > 1) {
    for (is=1;is<nscat;is++) {
      if (is == 1) {
	if ((status = ASCII_calloc_float (&sf, nlyr, nstr)) != 0)
	  return status;
	if ((status = ASCII_calloc_float (&sfu, nlyr, numu)) != 0)
	  return status;
      } 

      /* Calculate the source function, done for one order of scattering less */
      /* than the one we are interested in.                                   */
      if ((status = source_function( pmom, cwt, is-1, nlyr, nstr, numu,
				     ylmc, ylmu, u0uc, sf, sfu)) !=0)
	return status;

      /* Calculate higher orders of scattering */
      morsca(albedo, ssalb, dtauc, cmu, cwt, nlyr, nstr, numu, umu, 
	     sf, sfu, is, u0uc, u0uu);
    }
  }
  if ( nscat > 1) {
    free(sf);
  }


  /* Output at user levels */
  for (lu=0;lu<ntau;lu++) {
    rfldir[lu] = 0;
    rfldn[lu]  = 0;
    flup[lu]   = 0;
    uavgso[lu] = 0;
    uavgdn[lu] = 0;
    uavgup[lu] = 0;
    for (iu=0;iu<numu;iu++) {
      u0u[lu][iu] = 0;
    }
    for (is=0;is<nscat;is++) {
      for (iq=0;iq<nstr/2;iq++) {
	rfldn[lu]  += cwt[iq]*(-cmu[iq])*u0uc[is][layru[lu]][iq];
	uavgdn[lu] += cwt[iq]*u0uc[is][layru[lu]][iq];
	/*	printf("rfldn %d %d %d %d %f %f %e\n",                                 */
	/*	       is,lu,iq,layru[lu],uavgdn[lu],cwt[iq],u0uc[is][layru[lu]][iq]); */
      }
      for (iq=nstr/2;iq<nstr;iq++) {
	flup[lu]   += cwt[iq]*cmu[iq]*u0uc[is][layru[lu]][iq];
	uavgup[lu] += cwt[iq]*u0uc[is][layru[lu]][iq];
	/*		printf("flup %d %d %d %d %f %f %f %f\n",    */
	/*		       is,lu,iq,layru[lu],flup[lu],cmu[iq], */
	/*                     cwt[iq],u0uc[is][layru[lu]][iq]);    */
      }
      for (iu=0;iu<numu;iu++) {
	u0u[lu][iu]  += u0uu[is][layru[lu]][iu];
      }
    }
    rfldn[lu] *= 2*PI;
    flup[lu]  *= 2*PI;
    rfldir[lu] = umu0 * trs[layru[lu]];
    uavgso[lu] = trs[layru[lu]]/(4*PI);
    uavgdn[lu] *= 0.5;
    uavgup[lu] *= 0.5;

    /*    printf("rfldn %d %f %f\n",is,rfldn[lu], flup[lu]); */
    /*    printf("rfldir %f\n",rfldir[lu]);                  */
  }
  free(trs);
  return status;
}
