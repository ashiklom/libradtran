/************************************************************************/
/* solver.h                                                             */
/*                                                                      */
/* C header file for the Fortran RTE solvers (disort, sdisort, twostr)  */
/*                                                                      */
/* Author: Bernhard Mayer,                                              */
/*         NCAR, bmayer@ucar.edu                                        */
/*                                                                      */
/*----------------------------------------------------------------------*/
/* Copyright (C) 1999 Bernhard Mayer                                    */
/*                                                                      */
/* This program is free software; you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation; either version 1, or (at your option)  */
/* any later version.                                                   */
/*                                                                      */
/* This program is distributed in the hope that it will be useful,      */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the         */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* To obtain a copy of the GNU General Public License write to the      */
/* Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139,   */
/* USA.                                                                 */
/*----------------------------------------------------------------------*/
/************************************************************************/

#ifndef __solver_h
#define __solver_h

#if defined (__cplusplus)
extern "C" {
#endif

#include "f77-uscore.h"

typedef struct {
  double re;
  double im;
} pol_complex;


/* prototypes */
void F77_FUNC  (disort, DISORT)(int* nlyr, float* dtauc, float* ssalb, float* pmom, 
		       float* temper, float* wvnmlo, 
		       float* wvnmhi, int* usrtau, int* ntau, float* utau,
		       int* nstr, int* usrang, int* numu, float* umu, 
		       int* nphi, float* phi, int* ibcnd, float* fbeam, 
		       float* umu0, float* phi0, float* fisot, int* lamber, 
		       float* albedo, float* hl, float* btemp, float* ttemp, 
		       float* temis, int* deltam, int* planck,
		       int* onlyfl, float* accur, 
		       int* prndis, const char* header, int* maxcly, 
		       int* nzout, int* maxumu, int* maxcmu,
		       int* maxphi, float*  rfldir, float* rfldn, 
		       float* flup, float* dfdt, float*  uavg, 
		       float* uu, float* u0u, float* albmed, float* trnmed, 
		       float* uavgdn, float* uavgso, float* uavgup, int* quiet );

void F77_FUNC  (disort2, DISORT2)(int* nlyr, float* dtauc, float* ssalb, int* nmom, float* pmom,
				  int* nphas, float* phaso, double* mup,
				  float* temper, float* wvnmlo, float* wvnmhi, int* usrtau, 
				  int* ntau, float* utau, int* nstr, int* usrang, int* numu, 
				  float* umu, int* nphi, float* phi, int* ibcnd, float* fbeam, 
				  float* umu0, float* phi0, float* fisot, int* lamber, 
				  float* albedo, float* btemp, float* ttemp, 
				  float* temis, int* planck,
				  int* onlyfl, float* accur, 
				  int* prndis, const char* header, int* maxcly, 
				  int* nzout, int* maxumu, int* maxphi, int* maxmom,
				  float*  rfldir, float* rfldn, 
				  float* flup, float* dfdt, float*  uavg, 
				  float* uu, float* u0u, float* albmed, float* trnmed, 
				  float* uavgdn, float* uavgso, float* uavgup,
				  int* brdftype, float* rho0, float* k, float* theta, 
				  float *bsigma, float *bt1, float *bt2, float* scale, 
				  float *biso, float *bvol, float *bgeo, 
				  float* u10, float* pcl, float* xsal,
				  int* intcor, int* oldintcor,
				  int* quiet );
  
void F77_FUNC  (polradtran, POLRADTRAN)(int* nstokes, int* nummu, int* aziorder, 
			   double* max_delta_tau, int* src_code,
			   const char* quad_type, const char* deltam,
			   double* direct_flux, double* direct_mu,
			   double* ground_temp, const char* ground_type,
			   double* ground_albedo, pol_complex *ground_index,
			   double* sky_temp, double* wavelength,
			   int* num_layers, double* height,
			   double* temperatures, double* gas_extinct,
			   const char* scat_files, int* noutlevels,
			   int* outlevels, double* mu_values,
			   double* up_flux, double* down_flux,
			   double* up_rad, double* down_rad 
			   );
  
void F77_FUNC  (radtran, RADTRAN)(int* nstokes, int* nummu, int* aziorder, 
			double* max_delta_tau, int* src_code,
			const char* quad_type, const char* deltam,
			double* direct_flux, double* direct_mu,
			double* ground_temp, const char* ground_type,
			double* ground_albedo, pol_complex *ground_index,
			double* sky_temp, double* wavelength,
			int* num_layers, double* height, 
			double* temperatures, double* gas_extinct,
			const char* scat_files, int* noutlevels,
			int * outlevels, double* mu_values, double* up_flux, 
			double* down_flux, double* up_rad, double* down_rad );

void F77_FUNC  (sdisort, SDISORT)(int* nlyr, float* dtauc, float* ssalb, float* pmom, 
			float* temper, float* wvnmlo, 
			float* wvnmhi, int* usrtau, int* ntau, float* utau,
			int* nstr, int* usrang, int* numu, float* umu, 
			int* nphi, float* phi, float* fbeam, float* beta, int* nil,
			float* umu0, float* phi0, int* newgeo, float* zd,
			int* spher, float* radius, float* fisot,
			float* albedo, float* btemp, float* ttemp, 
			float* temis, int* deltam, int* planck,
			int* onlyfl, float* accur, int* quiet, int* ierror_d, 
			int* prndis, const char* header, int* maxcly, 
			int* nzout, int* maxumu, int* maxcmu,
			int* maxphi, float*  rfldir, float* rfldn, 
			float* flup, float* dfdt, float*  uavg, 
			float* uu, float* u0u, int* nscat,
			float* uavgdn, float* uavgso, float* uavgup,
			int* nrefrac, int* ichap, float* refind,
			int* ndenssza, float* denssza, float* denstab_sig,
			float* denstab, float* dtauc_mb);

void F77_FUNC  (spsdisort, SPSDISORT)(int* nlyr, float* dtauc, float* ssalb, float* pmom, 
			float* temper, float* wvnmlo, 
			float* wvnmhi, int* usrtau, int* ntau, float* utau,
			int* nstr, int* usrang, int* numu, float* umu, 
			int* nphi, float* phi, float* fbeam, float* beta, int* nil,
			float* umu0, float* phi0, int* newgeo, float* zd,
			int* spher, float* radius, float* fisot,
			float* albedo, float* btemp, float* ttemp, 
			float* temis, int* deltam, int* planck,
			int* onlyfl, float* accur, int* quiet, int* ierror_d, 
			int* prndis, const char* header, int* maxcly, 
			int* nzout, int* maxumu, int* maxcmu,
			int* maxphi, float*  rfldir, float* rfldn, 
			float* flup, float* dfdt, float*  uavg, 
			float* uu, float* u0u, 
			float* uavgdn, float* uavgso, float* uavgup );

void F77_FUNC  (twostr, TWOSTR)(float*  albedo, float* btemp, int* deltam, float* dtauc,
		       float* fbeam, float* fisot, float* gg, const char*header,
		       int* ierror_t, int* maxcly, int* maxulv, int* newgeo, 
		       int* nlyr, int* planck, int* ntau, int* prnt, int* quiet,
		       float* radius, int* spher, float* ssalb, float* temis,  
		       float* temper, float* ttemp, float* umu0,  int* usrtau, 
		       float* utau, float* wvnmlo, float* wvnmhi, float* zd, 
		       float* dfdt, float* flup, float* rfldir, float* rfldn, 
		       float* uavg );




void F77_FUNC  (tzs, TZS)(int *nlyr, float *dtauc, float *ssalb, float *temper, float *wvnmlo, float *wvnmhi, 
		      int *usrtau, int *ntau, float *utau, int *usrang, 
		      int *numu, float *umu, int *nphi, float *phi, 
		      float *albedo, float *btemp, float *ttemp, float *temis, int *plank, 
		      int *prnt, const char *header, int *maxcly, int *maxulv, int *maxumu, int *maxphi, 
		      float *rfldir, float *rfldn, float *flup, float *dfdt, float *uavg, float *uu, 
		      float *albmed, float *trnmed, float *uavgdn, float *uavgso, float *uavgup);

void F77_FUNC  (sss, SSS)(int* nlyr, float* dtauc, float* ssalb, int* nmom, float* pmom, 
		    float* temper, float* wvnmlo, float* wvnmhi, int* usrtau, 
		    int* ntau, float* utau, int* nstr, int* usrang, int* numu, 
		    float* umu, int* nphi, float* phi, int* ibcnd, float* fbeam, 
		    float* umu0, float* phi0, float* fisot, int* lamber, 
		    float* albedo, float* btemp, float* ttemp, 
		    float* temis, int* planck,
		    int* onlyfl, float* accur, 
		    int* prndis, const char* header, int* maxcly, 
		    int* nzout, int* maxumu, int* maxphi, int* maxmom,
		    float*  rfldir, float* rfldn, 
		    float* flup, float* dfdt, float*  uavg, 
		    float* uu, float* albmed, float* trnmed, 
		    float* uavgdn, float* uavgso, float* uavgup,
		    int* brdftype, float* rho0, float* k, float* theta, 
		    float* u10, float* pcl, float* xsal );
  
void F77_FUNC  (sssi, SSSI)(int* nlyr, float* dtauc, float* ssalb, int* nmom, float* pmom, 
		     float* temper, float* wvnmlo, float* wvnmhi, int* usrtau, 
		     int* ntau, float* utau, int* nstr, int* usrang, int* numu, 
		     float* umu, int* nphi, float* phi, int* ibcnd, float* fbeam, 
		     float* umu0, float* phi0, float* fisot, int* lamber, 
		     float* albedo, float* btemp, float* ttemp, 
		     float* temis, int* planck,
		     int* onlyfl, float* accur, 
		     int* prndis, const char* header, int* maxcly, 
		     int* nzout, int* maxumu, int* maxphi, int* maxmom,
		     float*  rfldir, float* rfldn, 
		     float* flup, float* dfdt, float*  uavg, 
		     float* uu, float* albmed, float* trnmed, 
		     float* uavgdn, float* uavgso, float* uavgup,
		     int* brdftype, float* rho0, float* k, float* theta, 
		     float* u10, float* pcl, float* xsal, float *ref, int *lctop);
  

void F77_FUNC  (cplkavg, CPLKAVG)(float* wvnmlo, float *wvnmhi, float *t, float *plkavg);


int   F77_FUNC  (dcheck, DCHECK) (int *nlyr, int *ntau, int *nstr, int *numu, int *nphi, int *optimize, int *delta);
int   F77_FUNC  (tcheck, TCHECK) (int *nlyr, int *ntau, int *nstr, int *numu, int *nphi, int *optimize, int *delta);

void F77_FUNC  (setout, SETOUT) (float *dtauc, int *nlyr, int *ntau, float *utau, 
			float *z, float *zout);
void F77_FUNC  (setoutd, SETOUTD) (float *dtauc, int *nlyr, int *ntau, float *utau, 
			float *z, float *zout);

void F77_FUNC  (qgausn, QGAUSN) (int* nn, float *cmu, float *cwt);

double F77_FUNC (dei, DEI) (double* x);

#if defined (__cplusplus)
}
#endif

#endif






