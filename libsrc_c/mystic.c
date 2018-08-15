/************************************************************************
 * $Id: mystic.c 3328 2017-12-19 10:54:25Z Claudia.Emde $
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
/*#include <values.h>*/

#include "uvspecrandom.h"
#include "wcloud3d.h"
#include "mystic.h"
#if HAVE_MYSTIC3D
  #include "mystic_3d.h"
#endif
#if HAVE_ALIS
  #include "alis.h"
#endif
#if HAVE_OPENGL
  #include "GLmystic.h"
#endif 
#include "yang56.h"
#include "sofi.h"
#include "errors.h"

#include "ascii.h"
#include "miecalc.h"
#include "phasetable.h"
#include "numeric.h"
#include "ambralsfor.h"
#include "ocean.h"
#include "f77-uscore.h"
#include "solver.h"
#include "rayleigh.h"
#include "locate.h"
#include "cdisort.h"

#ifdef HAVE_LIBTENSTREAM
#include "tenstream.h"
#endif

#if HAVE_LIBGSL 
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h> 
#endif

#if HAVE_LIBNETCDF
#include <netcdf.h>
#endif

#if HAVE_VROOM
/* include vroom, and evtl. also lidar stuff */
#include "vroom.h"
#if HAVE_LIDAR
#include "lidar.h"
#endif
#endif

/* will be deleted soon */
#define CLDPRP

/* spherical 3D, not for distribution */
#undef HAVE_SPHER
/* #define HAVE_SPHER */

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif

/* This is a fix needed at MIM, need to test whether this also works on DLR */
#if COMPILING_CONDOR_AT_MIM
#define fscanf __fscanf
#endif


/* CE: This has been included as option since to simulate limb sounding a */
/* variable earth radius is required CHECK!!! */ 
/* #define R_EARTH 6371211.0 */

/* for impatient users. If set, touching the file "fastexit"   */
/* in the working directly will let the code end prematurely   */
/* results are summarized as usual, only the number of simul   */
/* ated photons is reduced                                     */
#define ALLOW_FILE_CHECKS

/* this should probably be on for lidar, FODDIS is bad */
#undef NOFODDIS

/* experimental stuff for Lidar, not yet debugged */
#undef NEWRISQIDD
#undef NEWVISQIDD

#undef SPIKEOUT
#undef ESCAPES
#undef TIMINGS

/* absorption optical depth along path is calculated separately */
/* for different absorbers; if ONLY_TOT is true, only the total */
/* optical depth is calculated                                  */
#define ONLY_TOT      1

/* if PRINT_PHOTON_PATH is defined, the location of the photon at   */
/* the start and end of its journey as well as at each scattering   */
/* event will be written to stdout. */
/*WARNING: Will produce stack overflows for very long photon paths!*/
/* #define PRINT_PHOTON_PATH 1 */

/* experimental option to calculate full-sky radiances with "cone"  */
/* sampling, BM 27.7.2016                                           */
#define NTHETA 180
#define NPHI 180

#include "raytracing.h"

/* ??? Attention, this is only to reduce noise ??? CHECK!!! */
/* ??? but it actually throws away part of the ??? */
/* ??? signal. To get correct results, set     ??? */
/* ??? ACTINIC_CUTOFF to 0                     ??? */
#define ACTINIC_CUTOFF 0.01 

/* if MYSTIC_DEBUG is switched on, some extra tests are performed; */
/* might be useful for testing MYSTIC on new machines              */

/* #define MYSTIC_DEBUG  1 */

/* angle from which on the phase function will be tabulated in high resolution */
#define MUHIGH 0.99

/* global variable: random number generator */
#if HAVE_LIBGSL 
  gsl_rng *uvspecrng; 
#endif

/* prototypes of internal functions */
static int mystic_check_consistency ( sample_struct    *sample,
				      elevation_struct *elev,
				      surftemp_struct  *surftemp,
				      double            delta_scaling_mucut,
				      double            truncate,
				      int               source,
				      char             *umufilename,
				      int               ipa );

static int mystic_set (/* input */
		       int            nlyr,
		       int           *nmomaer,
		       int          **nthetaaer,
		       int            source,
		       int            absorption,
		       int            visualize,
		       int           *loaddata,
		       int            nzoutorg,
		       float         *zoutorg,
		       float          rpv_rho0,
		       float          rpv_k,
		       float          rpv_theta, 
		       float          hapke_h,
		       float          hapke_b0,
		       float          hapke_w,
		       float          rossli_iso,
		       float          rossli_vol,
		       float          rossli_geo,
		       float          u10,
		       float          pcl,
		       float          xsal,
		       float          bpdf_u10, 
		       char          *albedo_filename,
		       char          *alb_type_filename,
		       char          *rpvfilename,
		       char          *ambralsfilename,
		       char          *rosslifilename,
		       int            reflectalways,
		       /* throughput */
		       int           *first,
		       sample_struct *sample,
		       albedo_struct **albedo, 
		       float         *delta_scaling_mucut,
		       float         *truncate,
		       /* output */
		       int           *aerosol_scatter_type,
		       long int      *nphotons,
		       int           *nzout,
		       float        **zout,
		       float         *rfldir,
		       float         *rfldn,
		       float         *flup,
		       float         *uavgso,
		       float         *uavgdn,
		       float         *uavgup,
		       int            quiet );

static int mystic_verbose ( int            quiet,
			    sample_struct *sample,
			    albedo_struct *albedo,
			    int            aerosol_scatter_type,
			    long int       nphotons,
			    int            ipa,
			    double         delta_scaling_mucut,
			    double         truncate );

#ifdef ALLOW_FILE_CHECKS
static int check_premature_exit ( char          *filename,
				  sample_struct *sample );

static void check_progress ( long int       photoncounter,
			     long int       nphotons,
			     sample_struct *sample,
			     float          wavelength );
#endif

#if HAVE_LIBGSL 
static inline int read_and_write_random_status ( sample_struct *sample,
						 char          *randomstatusfilename,
						 int            photoncounter,
						 int           *readrandomstatus );
#endif

#if HAVE_OPENGL
static void initialize_plot_photon_path ( photon_struct *p,
					  sample_struct *sample,
					  result_struct *result,
					  int            photoncounter,
					  int           *plotpathcounter,
					  int           *plotpath,
					  int           *plotresultcounter,
					  int           *plotresult );
#endif

static int cleanup_after_photon_journey ( photon_struct *p,
					  sample_struct *sample,
					  result_struct *result,
					  int            n_caoth,
					  int            chouf );

static void free_mystic ( result_struct     *result,
			  atmosphere_struct *atmos,
			  sample_struct     *sample,
			  float             *zout,
			  int                absorption,
			  int                thermal_heating_method); /* **CK 2013.09.23 added thermal_heating_method */

/* static int istrue(double p); */
inline double random_tau();
static int reflection (double albedo);
static int emission  (double albedo);
static void random_direction (direction *dir);
/* now elsewhere: static void new_direction (double mu, double phi, direction *dir); */
static ddprofile3D *calloc_ddprofile3D (int Nz, int Nx, int Ny, int *threed);
static radiation_field *calloc_radiation_field (int Nx, int Ny, int Nd, int Nr, int Nt, int Np, int std, int cldprp);
static radiation_field_t *calloc_radiation_field_t(int Nc, int Nv, int Nx, int Ny, int Np);
static profile *calloc_profile (int ntypes, int n);

static profile3D *calloc_profile3D ( int   n_caoth,
				     int  *tocalloc,
				     int   Nz,
				     int   Nx,
				     int   Ny,
				     int **threed );

static void free_profile   (profile *);
static void free_profile3D (profile3D *);
static void free_ddprofile3D (ddprofile3D *);
static void free_radiation_field (radiation_field *res, int Nx, int Nr, int Nd, int Np, int std);

static int cp_1D_to_3D (profile3D *p3D, profile *p);

static void cp_optical_depth (optical_depth *dest, optical_depth *source);

static inline void mc_add_optical_depth (atmosphere_struct *atmos,
					 photon_struct     *p,
                                         double             length,
					 int                calc_tauext );

static void mc_add_absorption_3D (profile3D *kabs3D, profile *kabs,
                                  int kc, int ic, int jc, 
                                  double length, double totabs, int std,
                                  ddprofile3D *absorption3D, ddprofile3D *absorption3D2,
                                  struct tnode **wtree);

static double elevation (surface surf, double x, double y);

static int cross_bilinear (surface surf, double *a, double *b, double x0, double y0, double *solution);

static int setup_caoth3D ( caoth3d_out_struct  *caoth3d,
			   int                  nxcld,
			   int                  nycld,
			   int                  nzcld,
			   double               dxcld,
			   double               dycld,
			   float                wavelength,
			   atmosphere_struct   *atmos,
			   char                *datapath,
			   int                  delta_scaling,
			   int                  spherical3D,
			   int                  spherical3D_scene,
			   double               spherical3D_scene_lon_min,
			   double               spherical3D_scene_lon_max,
			   double               spherical3D_scene_lat_min,
			   double               spherical3D_scene_lat_max,
			   float                r_earth,
			   float                sza,
			   int                  setup_tipa,
			   int                  visualize,
			   int                  sample_cldprp,
			   int                  quiet );

static int setup_thermal (float         ***temper,
			  int              temper3d,
			  int              cldNx,
			  int              cldNy,  
			  int              nlyr,
			  float            wvnmlo,
			  float            wvnmhi,
                          profile         *kabs,
			  profile3D       *kabs3D,
			  float           *Z, 
                          int              source,
			  int              absorption,
			  int              thermal_heating_method,
                          profile3D      **Bplanck,
			  double          *maxemis,
			  double          *Watm,
			  surftemp_struct *surftemp,
			  int              backward,
			  float            zstart,
			  int             *backemis_kc,
                          int              nlambda_abs,
                          float           *lambda,
                          double        ***Bplanck_spectral,
                          int              quiet,
			  long int         nphotons,
			  sample_struct    *sample);  /* **CK 2013.08.27 Added nphotons and sample_struct, added cldNx and CldNy */

static int setup_albedo2D (char *albedo_filename, char *albedo_spectralfilename, char *rpv_filename, char *rossli_filename, int isAmbralsFile,
                           atmosphere_struct *atmos,
                           albedo_struct *albedo,
                           double alb, float *alb_type,
                           float *rpv_rho0, float *rpv_k, float *rpv_theta, 
                           float *rpv_scale, float *rpv_sigma,
			   float *rpv_t1, float *rpv_t2, 
                           char **rpv_labels, int rpv_nlabels,
			   float *hapke_h,float *hapke_b0, float *hapke_w,
			   float *rossli_iso, float *rossli_vol, float *rossli_geo, 
			   int rossli_hotspot,
                           float u10, float pcl, float xsal, float uphi, int solar_wind,
			   float bpdf_u10, 
                           int spectral, float wavelength, int polarisation,
			   int spherical3D, int spherical3D_scene,
			   double spherical3D_scene_lon_min,
			   double spherical3D_scene_lon_max,
			   double spherical3D_scene_lat_min,
			   double spherical3D_scene_lat_max,
			   int nlambda_abs,
			   double *alis_albedo,
			   double **alis_alb_type,
			   int ixmin,
			   int ixmax,
			   int iymin,
			   int iymax,
                           int quiet);

static int setup_umu2D (char *umu_filename, sample_struct *sample, 
			atmosphere_struct *atmos, int quiet);

static int setup_sample2D (atmosphere_struct *atmos,
                           float *zout, int nzout,
                           float *sza, float *phi0,
                           float *wavelength,
                           sample_struct *sample,
                           int quiet);

static int setup_mc_result (atmosphere_struct *atmos, sample_struct *sample,
                            result_struct *result, int absorption, int thermal_heating_method,
			    int quiet);

static int setup_Legendre_table_aerosol (double **moment, int nmom, int nphamat, double truncate, 
                                         pft *phase, int quiet);

static inline double get_Lambertian_phase(double mu);

static double sc_interp_mu (double reff, 
                            phase_function_table *phase, int scaled);

static int incone (double *dx, radang rad);

static double HG (double g, double mu);
static double HG2 (double g1, double g2, double ff, double mu);

static double HGint (double g, double mu);
static double HG2int (double g1, double g2, double ff, double mu);

static void count_photon (radiation_field *res, 
                          photon_struct *p, 
                          sample_struct *sample, 
                          double *totweight, 
                          int surfaceparallel, int incoming);

static void count_thermal_backward_photon (result_struct *result, 
                                           photon_struct *p, 
                                           atmosphere_struct *atmos,
                                           surftemp_struct *surftemp,
                                           sample_struct *sample, 
                                           double btemp_plkavg,
                                           double albedo, 
                                           int quiet, int boundary);

static int forward2backward (atmosphere_struct *atmos,
                             float **zout, int *nzout,
                             float *sza, float *phi0,
                             sample_struct *sample,
                             result_struct * result,
                             int quiet);

static inline int perform_periodic_bc (photon_struct *p, atmosphere_struct *atmos);

/* static photon_struct *calloc_photon (sample_struct *sample); */

static photon_struct *generate_photon ( int                source,
					int                photoncounter,
					atmosphere_struct *atmos,
					albedo_struct     *albedo,
					elevation_struct  *elev,
					sample_struct     *sample,
					float              wvnmlo,
					float              wvnmhi,
					float              wavelength,
					double             sza,
					double             phi,
					int                ipa,
					long int           nphotons);   /* **CK 2013.08.28 Include nphotons */

void generate_photon_backward_photon_direction_panorama ( sample_struct* sample, 
                                                          photon_struct* p,
                                                          double phi, 
                                                          double sza, 
                                                          double cossza,
                                                          double sinsza, 
                                                          double cosphi, 
                                                          double sinphi );

int generate_photon_backward_vertical_position_heat (sample_struct* sample, 
                                                     atmosphere_struct* atmos, 
                                                     photon_struct* p, 
                                                     int* hit, 
                                                     long int nphotons, 
                                                     int* counter1 );

static photon_struct* generate_photon_backward_photon_direction (sample_struct* sample, 
                                                                 atmosphere_struct* atmos, 
                                                                 photon_struct* p,
                                                                 elevation_struct* elev,
                                                                 int* hit, 
                                                                 long int nphotons, 
                                                                 int* counter1,
                                                                 double sza,
                                                                 double phi,
                                                                 double phase,
                                                                 int source,
                                                                 double* norm);

static int generate_photon_backward_vertical_position (sample_struct* sample, 
                                                       atmosphere_struct* atmos, 
                                                       photon_struct* p, 
                                                       int* hit, 
                                                       long int nphotons, 
                                                       int* counter1);

static photon_struct* generate_photon_solar_thermal_backward (sample_struct* sample, 
                                                              atmosphere_struct* atmos, 
                                                              photon_struct* p, 
                                                              elevation_struct* elev,
                                                              int* hit, 
                                                              long int nphotons, 
                                                              int* counter1,
                                                              double sza,
                                                              double phi,
                                                              double phase,
                                                              int source,
                                                              int sofi,
                                                              double* norm,
                                                              double* pd);

static void hunt_modified (float *xx, int n, float x, int *jlo, int *hit);

static double slt2hrz (elevation_struct *elev, photon_struct *p, int surfaceparallel, int incoming);

static void area_average (radiation_field *res, double ***back,
                          int Nx, int Ny,
                          int islower, int isupper, int jslower, int jsupper, 
                          int   *ndir, int   *ndn, int   *nup,
                          float *edir, float *edn, float *eup,
                          float *fdir, float *fdn, float *fup,
                          int backward, double **surface_area, int **surface_area_counter,
			  int elev2D, int surfaceparallel);

static int summarize_result (result_struct *result, sample_struct *sample,
                             atmosphere_struct *atmos, albedo_struct *albedo,
                             int source, long int mcphotons, int mcsimulations, double sza,
                             int absorption, int write_files, int write_output_as_netcdf,
                             int tenstream,
                             char *basename,
                             float *rfldir, float *rfldn,  float *flup,
                             float *uavgso, float *uavgdn, float *uavgup,
                             float ***rfldir3d, float ***rfldn3d, float ***flup3d, 
                             float ***uavgso3d, float ***uavgdn3d, float ***uavgup3d, 
                             float ****radiance3d, float ***absback3d, 
                             float ***abs3d, float ******radiance3s_is,
                             float ***rfldir3d_var, float ***rfldn3d_var, float ***flup3d_var, 
                             float ***uavgso3d_var, float ***uavgdn3d_var, float ***uavgup3d_var, 
                             float ****radiance3d_var, float ***absback3d_var, 
                             float ***abs3d_var, 
			     int escape, int elev2D, double wavelength, int quiet);

static int reflect (photon_struct *p, albedo_struct *albedo, elevation_struct *elev,
                    atmosphere_struct *atmos, sample_struct *sample, 
                    result_struct *result,
                    float wvnmlo, float wvnmhi, float *refind);

static double reflection_probability_tot (albedo_struct *albedo, photon_struct *p,
                                          double murad, double mu1, double deltaphirad, 
                                          int ia, int ja, int il, int DoLE,
                                          float wvnmlo, float wvnmhi, double paw, 
					  int spherical3D, 
					  int spectral_is, int nlambda,
					  int *status);

static int scattering ( atmosphere_struct *atmos,
			photon_struct     *p,
			sample_struct     *sample,
			result_struct     *result,
			elevation_struct  *elev,
			albedo_struct     *albedo,
			float              wvnmlo,
			float              wvnmhi, 
			float             *refind );

static inline int random_scatter_type (atmosphere_struct *atmos, photon_struct *p);

static int setup_mystic (caoth3d_out_struct *caoth3d,
                         atmosphere_struct  *atmos,
			 result_struct      *result,
                         albedo_struct      *albedo,
                         surftemp_struct    *surftemp,
                         sample_struct      *sample,
			 float               delta_scaling_mucut,
                         int                 absorption, 
			 int                 thermal_heating_method,
                         int                 nlyr,
			 int                 n_caoth,
                         float             **dt_s,
			 float             **om_s,
			 float             **g1_s,
			 float             **g2_s,
			 float             **f_s,
                         float             **ds_s,
                         alis_struct        *alis,
                         float               rayleigh_depol,
                         float             **re_s,
                         float           ***temper,
			 int                temper3d, 
                         float              *zprof,
                         float              *sza,
			 float              *phi0,
			 float               sza_spher,
			 float               phi0_spher,
                         int                 aerosol_scatter_type,
			 float            ***momaer,
			 int                *nmomaer, 
			 int               **nthetaaer,
			 float            ***thetaaer,
			 double           ***muaer,
			 float            ***phaseaer,
                         int                 nphamataer,
                         float               alb,
			 float               *alb_type,
			 float               *rpv_rho0,
			 float               *rpv_k,
			 float               *rpv_theta, 
			 float               *rpv_scale,
			 float               *rpv_sigma,
			 float               *rpv_t1,
			 float               *rpv_t2, 
			 float               *hapke_h,
			 float               *hapke_b0,
			 float               *hapke_w,
			 float               *rossli_iso,
			 float               *rossli_vol,
			 float               *rossli_geo,
			 int                  rossli_hotspot,
                         float               *u10,
			 float               *pcl,
			 float               *xsal,
			 float               *uphi,
			 int                 *solar_wind,
                         float               *bpdf_u10,
                         int                  source,
			 float               *wvnmlo,
			 float               *wvnmhi,
			 float               *wavelength,
                         float               *zout,
			 int                  nzout,
                         int                  nxcld,
			 int                  nycld,
			 int                  nzcld,
			 double               dxcld,
			 double               dycld,
                         char                *umufilename,
                         char                *sunshape_filename,
                         char                *albedo_filename,
                         char                *alb_type_filename,
                         char                *rpvfilename,
			 char               **rpv_labels,
			 int                  rpv_nlabels,
                         char                *ambralsfilename,
                         char                *rosslifilename,
                         double               truncate,
                         int                  loaddata, 
                         char                *datapath,
			 int                  spectral,
			 double               r_earth,
                         int                  visualize,
			 int                  quiet, 
			 long int             nphotons); /* **CK 2013.08.27 Added nphotons */

static int setup_sunshape ( sample_struct *sample,
			    char          *filename,
			    float          wvnmlo,
			    float          wvnmhi );

static inline int calc_surface_emission ( albedo_struct *albedo,
					  int            source,
					  float          btemp,
					  float         *wvnmlo,
					  float         *wvnmhi,
					  float         *bplkavg,
					  double        *Wsurf );

static inline int setup_aerosol1D ( pft     **phase_aer, 
				    int       nlyr,
				    float    *dtaer,
				    float    *omaer,
				    float  ***momaer,
				    int      *nmomaer, 
				    int     **nthetaaer,
				    double ***muaer,
				    float  ***phaseaer,
				    int       nphamataer,
				    double    truncate,
				    int       quiet );

static void free_mystic_for_load (albedo_struct *albedo, 
				  sample_struct *sample);

inline double sc_Isotropic_phi(void);

static void free_albedo    (albedo_struct *albedo);
static void free_rpv       (albedo_struct *albedo);
static void free_rossli    (albedo_struct *albedo);

static void free_atmosphere (atmosphere_struct *atmos);

static void free_result (result_struct *result, int thermal_heating_method,  /* **CK 2013.09.23 added thermal_heating_method */
                         atmosphere_struct *atmos, sample_struct *sample, 
                         int absorption);

static int register_at_start (photon_struct *p, sample_struct *sample,
                              atmosphere_struct *atmos, result_struct *result,
                              int source, int absorption);

static photon_struct *create_escape_photon ( photon_struct     *photon, 
					     atmosphere_struct *atmos,
					     sample_struct     *sample,
					     int                id );

static int escape_probability (photon_struct *photon,
                               atmosphere_struct *atmos,
                               sample_struct *sample, result_struct *result,
                               elevation_struct *elev, int id, int il, int it,
                               int surface, float *refind);

static void count_escape (photon_struct *p, sample_struct *sample,
			  atmosphere_struct *atmos,
			  radiation_field *result_alt, int id, int il,
			  int backward, double ***result_back, double ***result_back2,
                          double *****result_back_t, 
			  double **circcontr,
                          jacobian_result_field *jacobian,
                          radiation_field_t *result_rad_t,
			  double tauabs, double slant2horz, 
                          double *totweight_spectral,
                          double *totweight_concentration,
                          double *pathlength_per_layer_tot,
                          mishchenko_cb_struct* mish);

static int equal (double x, double y);

static inline void sample_coord  (photon_struct *p, sample_struct *sample, int *is, int *js);
static int  sample_radius (photon_struct *p, sample_struct *sample, int *ir);
static int  sample_time   (photon_struct *p, sample_struct *sample, int *it);

static inline void calc_normal_vector ( photon_struct     *p,
					elevation_struct  *elev,
					sample_struct     *sample,
					atmosphere_struct *atmos,
					int                use_elev,
					int                use_sensordirection,
					double            *norm );

static inline void elev_coord_normal_vector (photon_struct *p, elevation_struct *elev, double *dx);

static inline void albedo_coord (photon_struct *p, albedo_struct *albedo, int *ia, int *ja);

static inline void surftemp_coord (photon_struct *p, surftemp_struct *surftemp, int *it, int *jt);

static int calculate_emission (atmosphere_struct *atmos, result_struct *result);
static int calculate_backward_emission (atmosphere_struct *atmos, sample_struct *sample, int kc, double *emis);


#ifdef PRINT_PHOTON_PATH
  static void add_to_photonpath (struct photon_path **p, double *x);
  static void print_and_clear_photonpath (struct photon_path *p);
#endif

static int setup_profiles1D ( int                n_caoth,
			      float            **dt_s,
			      float            **om_s,
			      float            **g1_s,
			      float            **g2_s,
			      float            **f_s,
			      float            **ds_s,
                              float            **re_s,
			      float             *zprof,
			      int                nlyr,
			      sample_struct     *sample, 
			      atmosphere_struct *atmos,
                              alis_struct    *alis
                              );

static int calloc_1D_atmos ( atmosphere_struct *atmos,
			     int                nlyr,
			     int                n_caoth,
			     int                nlambda_abs,
			     int                nlambda_ref );

static inline double get_ksca_spectral (atmosphere_struct *atmos, photon_struct *p, int isp, int iv);
static inline double get_kabs_tot (atmosphere_struct *atmos, photon_struct *p);
static inline double get_g1       (atmosphere_struct *atmos, photon_struct *p, int isp);
static inline double get_g2       (atmosphere_struct *atmos, photon_struct *p, int isp);
static inline double get_ff       (atmosphere_struct *atmos, photon_struct *p, int isp);
static inline double get_reff     (atmosphere_struct *atmos, photon_struct *p, int isp);
static inline double get_dxlwc    (atmosphere_struct *atmos, photon_struct *p, int isp);
static inline double get_dylwc    (atmosphere_struct *atmos, photon_struct *p, int isp);
static inline double get_dzlwc    (atmosphere_struct *atmos, photon_struct *p, int isp);
static inline double get_kext     (atmosphere_struct *atmos, photon_struct *p);

static void stokes_vector_sca(double *stokes_vector,
                              double *P,
                              double *pdir_sca,
                              double *pdir_inc,
                              double  phi_sun,
                              double  phi_det);

static int get_phase_matrix_caoth ( atmosphere_struct *atmos,
				    photon_struct     *p,
				    int                isp,
				    double             mu,
				    int                np,
				    /* Output */
				    double            *phase_matrix );

static inline int calc_stokes_vector_for_escape ( sample_struct     *sample,
						  atmosphere_struct *atmos,
						  double            *phase_matrix,
						  double             mu2,
						  direction          dir_inc,
						  double             phi_esc,
						  direction          dir_esc,
						  int                polmat,
						  /* In-/Output */
						  photon_struct     *p );

int get_reflection_probability_matrix ( albedo_struct *albedo,
					photon_struct *p,
					int            nstokes,
					int            backward,
					double         mu_inc,
					double         phi_inc,
					double         mu_esc,
					double         phi_esc, 
					int            ia,
					int            ja,
					int            il,
					float          wvnmlo,
					float          wvnmhi,
					int            iv,
					/* Output */
					double      ***refl_mat );

static double reflection_polarized ( albedo_struct *albedo,
				     photon_struct *p,
				     int            nstokes,
				     int            backward,
				     double        *dir_inc,
				     double        *dx_out,
				     double        *n_hor,
				     double         mu_inc,
				     double         phi_inc,
				     double         mu_esc,
				     double         phi_esc, 
				     double         phi_target,
				     int            ia,
				     int            ja,
				     int            il,
				     float          wvnmlo,
				     float          wvnmhi,
				     int            escape,
				     int            spherical3D, 
				     int            spectral_is,
				     int            nlambda,
				     float         *lambda,
				     int           *status );

void phase_matr_rot ( double **Z_in, double cos2alpha, double sin2alpha, double **Z_out );
void rot_phase_matr ( double **Z_in, double cos2alpha, double sin2alpha, double **Z_out );

/*static inline void dtauabs_spectral_calc(atmosphere_struct *atmos, photon_struct *p*/
/*  double step); */

static int correct_escape_dir_refract(double *dir_sensor, 
                                      photon_struct *photon, 
                                      atmosphere_struct *atmos, 
                                      sample_struct *sample, 
                                      float *refind);

static inline double std_noisy (double var, double avg, double dmcphotons);

int set_panorama_dir ( sample_struct *sample, float *sza, float *phi0 );

double Fresnel ( double nr,
		 double ni,
		 double coschi );
void index_water ( double wl,
		   double xsal,
		   double *nr,
		   double *ni );

static double find_max (double x, double y, double z);  /* **CK 2013.08.27  */

static double find_min (double x, double y, double z); /* **CK 2013.08.27 */

static double pst (double mu, int i, double dtau); /* **CK 2013.08.27 */

static int emabsopt_location (photon_struct *p, atmosphere_struct *atmos, sample_struct *sample, long int photons); /* **CK 2013.08.27 */

static int denet_weight (photon_struct *p, atmosphere_struct *atmos, sample_struct *sample, long int nphotons, int *counter1); /* **CK 2013.08.27 */

static int denet_set_direction (photon_struct *p, atmosphere_struct *atmos, sample_struct *sample, long int nphotons, int *counter1, double *norm); /* **CK 2013.08.27 */

static int dist_photon_direction_by_phasemax(photon_struct *p, sample_struct *sample, 
    double *sza, double *cossza, double *sinsza); /* Bernhard Reinhardt 2013.12.23  */

#ifdef CLDPRP
static int summarize_result_cldprp (sample_struct* sample, atmosphere_struct* atmos,
                                     result_struct* result, char* basename, int elev2D,
                                     int islower, int isupper, int jslower, int jsupper);
#endif
inline static void write_backward_irradiance(sample_struct* sample, double result,
                                             FILE* backfile, int is, int js, int nn);

static int summarize_result_forward_surface(sample_struct* sample, 
                                            atmosphere_struct* atmos,
                                            result_struct* result,
                                            char* basename,
                                            double dmcphotons,
                                            double factor,
                                            double incident,
                                            int write_files,
                                            int elev2D,
                                            char* flxext, char* flxext2,
                                            char* radext, char* radext2,
                                            char* rplext, char* rplext2);

static int summarize_result_backward_irradiance(sample_struct* sample, 
                                                atmosphere_struct* atmos,
                                                result_struct* result,
                                                char* basename,
                                                char* specradfilename,
                                                char* bacext,
                                                char* bacext2,
                                                double dmcphotons,
                                                double factor,
                                                double wavelength,
                                                int mcsimulations,
                                                int write_files,
                                                int elev2D,
                                                int islower,
                                                int isupper,
                                                int jslower,
                                                int jsupper);

static int summarize_result_backward_radiance(sample_struct* sample, 
                                              atmosphere_struct* atmos,
                                              result_struct* result,
                                              char* basename,
                                              char* specradfilename,
                                              char* bacext,
                                              char* bacext2,
                                              char* picext,
                                              char* bacjacext,
                                              char* bacjacext2,
                                              double dmcphotons,
                                              double factor,
                                              double sza,
                                              int source,
                                              int mcsimulations,
                                              int write_files,
                                              int elev2D,
                                              int islower,
                                              int isupper,
                                              int jslower,
                                              int jsupper);

static int summarize_result_forward_altitude_radiance(sample_struct* sample, 
                                                      atmosphere_struct* atmos,
                                                      result_struct* result,
                                                      char* basename,
                                                      char* specradfilename,
                                                      char* radext,
                                                      char* radext2,
                                                      double dmcphotons,
                                                      double factor,
                                                      int write_files,
                                                      int kc);

static int summarize_result_passback3D(sample_struct* sample, 
                                       atmosphere_struct* atmos,
                                       result_struct* result,
                                       float ***rfldir3d, float ***rfldn3d, float ***flup3d, 
                                       float ***uavgso3d, float ***uavgdn3d, float ***uavgup3d, 
                                       float ****radiance3d, float ***absback3d,
                                       float ***rfldir3d_var, float ***rfldn3d_var, float ***flup3d_var, 
                                       float ***uavgso3d_var, float ***uavgdn3d_var, float ***uavgup3d_var, 
                                       float ****radiance3d_var, float ***absback3d_var,
                                       float ***abs3d, float ******radiance3d_is,
                                       int escape,
                                       int islower,
                                       int isupper,
                                       int jslower,
                                       int jsupper,
                                       int* lu3D);

static int summarize_result_altitude_rpl(sample_struct* sample, 
                                         atmosphere_struct* atmos,
                                         result_struct* result,
                                         char* basename,
                                         char* rplext,
                                         char* rplext2,
                                         double dmcphotons,
                                         double incident,
                                         int write_files,
                                         int kc);

static int summarize_result_absorption(sample_struct* sample, 
                                       atmosphere_struct* atmos,
                                       result_struct* result,
                                       float ***abs3d, 
                                       float ***abs3d_var, 
                                       char* basename,
                                       char* absext,
                                       char* absext2,
                                       double dmcphotons,
                                       double incident,
                                       int write_files,
                                       int absorption);

static int summarize_result_forward_altitude_irradiance (sample_struct* sample, 
                                                         atmosphere_struct* atmos,
                                                         result_struct* result,
                                                         char* basename,
                                                         double dmcphotons,
                                                         double factor,
                                                         int write_files,
                                                         int kc,
                                                         char* flxext, char* flxext2);

static int summarize_result_boxairmass (sample_struct* sample, 
                                        atmosphere_struct* atmos,
                                        result_struct* result,
                                        char* basename,
					int islower,
                                        int isupper,
                                        int jslower,
                                        int jsupper);

int mystic (int                  *nlyr,
	    int                   n_caoth,
            float               **dt_s,
	    float               **om_s,
	    float               **g1_s,
	    float               **g2_s,
	    float               **f_s,
            float               **ds_s,
            alis_struct          *alis,
            float                *refind,
            float                 r_earth,
            float                *rayleigh_depol,
            caoth3d_out_struct   *caoth3d,
            float               **re_s,
            float              ***temper,
	    int                   temper3d, 
            float                *zprof,
            float              ***momaer,
	    int                  *nmomaer,
	    int                 **nthetaaer,
	    float              ***thetaaer,
	    double             ***muaer,
	    float              ***phaseaer,
	    int                   nphamataer,
	    int                  *spectral,
            float                *alb,
            float                *alb_type,
	    float                 sza,
	    float                 phi0,
	    float                 sza_spher,
	    float                 phi0_spher,
	    long int              nphotons,
            int                  *source,
	    float                *wvnmlo,
	    float                *wvnmhi,
	    float                *wavelength,
            float                *zoutorg,
	    int                  *nzoutorg,
            float                *rpv_rho0,
	    float                *rpv_k,
	    float                *rpv_theta, 
            float                *rpv_scale,
	    float                *rpv_sigma,
	    float                *rpv_t1,
	    float                *rpv_t2,
	    float                *hapke_h,
	    float                *hapke_b0,
	    float                *hapke_w,
	    float                *rossli_iso,
	    float                *rossli_vol,
	    float                *rossli_geo,
	    int                   rossli_hotspot,
            float                *u10,
	    float                *pcl,
	    float                *xsal,
	    float                *uphi,
	    int                  *solar_wind,
            float                *bpdf_u10, 
            int                  *tenstream,
	    char                 *tenstream_options,
            int                  *ipa,
	    int                  *absorption,
	    int                  *thermal_heating_method,
	    int                  *loaddata,
            sample_struct        *sample,
            elevation_struct     *elev,
            surftemp_struct      *surftemp,
            char                 *basename, 
            char                 *umufilename, 
            char                 *sunshape_filename, 
            char                 *albedo_filename,
            char                 *alb_type_filename,
            char                 *rpvfilename,
	    char                **rpv_labels,
	    int                  *rpv_nlabels,
            char                 *ambralsfilename,
            char                 *rosslifilename,
            float                *delta_scaling_mucut,
            float                *truncate,
            int                  *reflectalways,
            int                   quiet,
            float                *rfldir,
	    float                *rfldn,
	    float                *flup,
            float                *uavgso,
	    float                *uavgdn,
	    float                *uavgup,
            float              ***rfldir3d,
	    float              ***rfldn3d,
	    float              ***flup3d, 
            float              ***uavgso3d,
	    float              ***uavgdn3d,
	    float              ***uavgup3d, 
            float             ****radiance3d,
	    float              ***absback3d,
            float              ***abs3d,
            float           ******radiance3d_is,
            float              ***rfldir3d_var,
	    float              ***rfldn3d_var,
	    float              ***flup3d_var, 
            float              ***uavgso3d_var,
	    float              ***uavgdn3d_var,
	    float              ***uavgup3d_var, 
            float             ****radiance3d_var,
	    float              ***absback3d_var,
            float              ***abs3d_var, 
	    int                   nxcld,
	    int                   nycld,
	    int                   nzcld,
	    double                dxcld,
	    double                dycld,
            char                 *datapath,
            char                 *randomstatusfilename,
	    int                   readrandomstatus,
	    int                   write_output_as_netcdf,
            int                   write_files,
	    int                   visualize )
{
  /* local variables */
  static int first=1;

  int status=0;

  long int photoncounter=0;
  long int totalphotoncounter=0;
  int totalsimulations=0;
  int outofdomaincounter=0;

  static albedo_struct *albedo = NULL;

  result_struct     *result = calloc (1, sizeof(result_struct));
  atmosphere_struct *atmos  = calloc (1, sizeof(atmosphere_struct));
  
  photon_struct *p = NULL;

  int starttime=0, endtime=0;

  int aerosol_scatter_type=0;
  
  float *zout=NULL; 
  int dummy=0;
  int *nzout=&dummy;

  int plotpath=0;

#ifdef ALLOW_FILE_CHECKS
  long int fileteststarttime=0;
  long int newfileteststarttime=0;
#endif

#if HAVE_OPENGL
  int plotpathcounter=0, plotresult=0, plotresultcounter=0;
#endif

#ifdef TIMINGS
  int clockstarttime=0;
#endif

  int chouf=0; /* needed by MUCHOUT etc */

  /****************************************************/
  /* 1. return if mystic does not need to do anything */
  /****************************************************/

  /* do nothing if SZA >= 90 degrees, unless in spherical geometry */
  if (!sample->spherical && !sample->spherical3D && sza>=90.0)
    return 0;

  /* do nothing if number of photons smaller or equal zero */
  if (*absorption!=MCFORWARD_ABS_EMISSION && nphotons<=0)
    return 0;

  /****************************************************/
  /* 2. check on errors                               */
  /****************************************************/
  status = mystic_check_consistency ( sample, elev, surftemp,
				      *delta_scaling_mucut, *truncate,
				      *source, umufilename, *ipa );
  if (status)
    return fct_err_out (status, "mystic_check_consistency", ERROR_POSITION);

  /****************************************************/
  /* 3. set things                                    */
  /****************************************************/
  status = mystic_set ( *nlyr, nmomaer, nthetaaer, *source,
			*absorption, visualize, loaddata,
			*nzoutorg, zoutorg,
			*rpv_rho0, *rpv_k, *rpv_theta, 
			*hapke_h, *hapke_b0, *hapke_w,
			*rossli_iso, *rossli_vol, *rossli_geo,
			*u10, *pcl, *xsal, *bpdf_u10, 
			albedo_filename, alb_type_filename, rpvfilename,
			ambralsfilename, rosslifilename,
			*reflectalways, &first, sample, &albedo,
			delta_scaling_mucut, truncate,
			&aerosol_scatter_type, &nphotons, nzout, &zout,
			rfldir, rfldn, flup, uavgso, uavgdn, uavgup,
			quiet );
  if (status)
    return err_out ("Error, mystic_set returned status %d\n", status);

  /****************************************************/
  /* 4. talk to user                                  */
  /****************************************************/
  status = mystic_verbose ( quiet, sample, albedo,
			    aerosol_scatter_type, nphotons, *ipa,
			    *delta_scaling_mucut, *truncate);
  if (status)
    return err_out ("Error, mystic_verbose returned status %d\n", status);
  
  /* initial setup: read 3D caoths, 2D albedo, and elevation fields; */
  /* allocate memory and initialize arrays.                          */
  status = setup_mystic (caoth3d, atmos, result, albedo, surftemp, sample,
			 *delta_scaling_mucut,
                         *absorption, *thermal_heating_method,
                         *nlyr,
			 n_caoth,
                         dt_s, om_s, g1_s, g2_s, f_s, ds_s,
                         alis,
                         *rayleigh_depol,
                         re_s,
                         temper,
			 temper3d,
                         zprof,
                         &sza, &phi0,
                         sza_spher, phi0_spher,
                         aerosol_scatter_type, momaer, nmomaer,
			 nthetaaer, thetaaer, muaer, phaseaer,
			 nphamataer,
                         *alb, alb_type, rpv_rho0, rpv_k, rpv_theta,
			 rpv_scale, rpv_sigma, rpv_t1, rpv_t2,
			 hapke_h, hapke_b0, hapke_w,
			 rossli_iso, rossli_vol, rossli_geo, rossli_hotspot,
                         u10, pcl, xsal, uphi, solar_wind, bpdf_u10,
                         *source, wvnmlo, wvnmhi, wavelength,
                         zout, *nzout,
                         nxcld, nycld, nzcld, dxcld, dycld,
                         umufilename,
                         sunshape_filename,
                         albedo_filename, alb_type_filename, 
                         rpvfilename, rpv_labels, *rpv_nlabels, 
                         ambralsfilename, rosslifilename,
                         *truncate,
                         *loaddata, 
                         datapath, *spectral,
			 r_earth,
                         visualize, quiet, nphotons);  /* **CK  2013.08.27 Added nphotons */

  if (status)
    return err_out ("Error, setup_mystic returned status %d\n", status);

  /* in case of atmospheric thermal radiation:                          */
  /* if maximum emission equals zero, we are done (no thermal emission) */
  if (*source == MCSRC_THERMAL_ATMOSPHERE)
    if (atmos->maxemis<=0)
      return 0;

  /* if MCFORWARD_ABS_EMISSION was specified, we only need to calculate */
  /* the emission of each grid cell.                                    */
  if (*absorption == MCFORWARD_ABS_EMISSION) {
    status = calculate_emission (atmos, result);
    if (status)
      return err_out ("Error %d returned by calculate_emission()\n", status);
  }

#if HAVE_OPENGL
  if (visualize) {
    GLmystic_setdomain (GLMYSTIC_PHOTONS,
                        atmos->X[0], atmos->Y[0], atmos->Z[0], 
                        atmos->X[atmos->Nx], 
                        atmos->Y[atmos->Ny], 
                        atmos->Z[atmos->Nz]);
  }
#endif

  /***********************/
  /* end of preparations */ 
  /***********************/
  if (*tenstream && nphotons>0 ) {
#ifdef HAVE_LIBTENSTREAM
    status = tenstream_like_mystic( 
        *source,
        dxcld,dycld,atmos->Z,
        albedo->albedo, phi0, sza,
        atmos,sample,surftemp,
        temper,*wvnmlo, *wvnmhi,
        rfldir,rfldn,flup,
        rfldir3d,rfldn3d,flup3d,
        *absorption, abs3d );

    if (status!=0) {
      fprintf (stderr, "Error %d returned by tenstream_like_mystic()\n", status);
      return status;
    }

    nphotons = 0;
#else

    fprintf (stderr, "Error you tried running the TENSTREAM solver with input option mc_tenstream -- your installation however is not configured to run it. Check your installation");
    return -1;
#endif
  }

  if (!quiet)
    fprintf (stderr, " ... starting first photon!\n");
  
  starttime = (int) time(NULL);
  
  /*****************************************/
  /* loops over lasers and lidar detectors */
  /*****************************************/
  for (sample->ili = 0; sample->ili < sample->Nli; sample->ili++) {

    if (!quiet && sample->LidarLocEst && (sample->Nli > 1))
      fprintf(stderr,"starting laser-detector-system %d ...\n",sample->ili);

    /*****************************************************************************/
    /* loop over backward sample pixels; in forward mode, upper and lower pixels */
    /* are identical and the photon loop is executed only once, as it should be  */
    /*****************************************************************************/
    for (sample->backward_is = sample->backward_islower;
	 sample->backward_is <= sample->backward_isupper; sample->backward_is++) {
    for (sample->backward_js = sample->backward_jslower;
	 sample->backward_js <= sample->backward_jsupper; sample->backward_js++) {

      if (!quiet && sample->backward)
        fprintf (stderr, " ... sample pixel %d %d, %ld photons\n",
		 sample->backward_is, sample->backward_js, nphotons);
          
      if (!sample->backward && (sample->backward_is>0 || sample->backward_js>0)) {
	fprintf (stderr, "\n");
        fprintf (stderr, "*** Warning: loop over backward pixels initiated although forward\n");
        fprintf (stderr, "*** calculation. This is not supposed to happen!\n");
	fprintf (stderr, "\n");
      }

#ifdef ALLOW_FILE_CHECKS
      newfileteststarttime=(long int) clock()/10000;
      /* wait at least one second before next file test */
      if (newfileteststarttime - fileteststarttime > 1000)
	if ( check_premature_exit ( "slowexit", sample ) )
	  break;
#endif

#ifdef TIMINGS
      clockstarttime=clock();
#endif

      if (sample->panorama) {
	status = set_panorama_dir ( sample, &sza, &phi0 );
	if (status<0)
	  return err_out ("Error %d returned by set_panorama_dir()\n", status);
	if (status==1)
	  continue;
      }
      
      /* calculate emission for backward heating rate calculation */


      if (sample->backward == MCBACKWARD_EMIS)  {
	status = calculate_backward_emission ( atmos, sample, sample->backemis_kc,
					       &(result->backemis[sample->backward_is][sample->backward_js]) );
	if (status)
	  return err_out ("Error %d returned by calculate_backward_emission()\n", status);
      }
      



      /*********************/
      /* loop over photons */ 
      /*********************/

      photoncounter=0;
      
      while(photoncounter++ < nphotons) {

#ifdef ALLOW_FILE_CHECKS
	newfileteststarttime=(long int) clock()/10000;
	/* wait at least one second before next file test */
	if (newfileteststarttime - fileteststarttime > 1000) {
	  check_progress ( photoncounter, nphotons, sample, *wavelength );
	  if ( check_premature_exit ( "fastexit", sample ) )
	    break;
	  fileteststarttime=newfileteststarttime;
	}
#endif

#if HAVE_LIBGSL
	status = read_and_write_random_status ( sample,
						randomstatusfilename,
						photoncounter,
						&readrandomstatus );
	if (status)
	  return err_out ("Error %d returned by read_and_write_random_status()\n", status);
#endif

	/* for VROOM paper, Fig. 4 */
	/*	fprintf(stderr,"%d %e %e\n",photoncounter,
		result->back[sample->backward_is][sample->backward_js][0]
		/((double) photoncounter),
		result->back2[sample->backward_is][sample->backward_js][0]
		/((double) photoncounter));*/

        /* generate a new photon */

        if ((p = generate_photon (*source, photoncounter,
                                  atmos, albedo, elev, sample,
                                  *wvnmlo, *wvnmhi, *wavelength,
                                  sza, phi0, *ipa, nphotons)) == NULL) {  /* **CK  2013.08.27 Add nphotons */
          fprintf (stderr, "Error, failed to generate a photon\n");
          return -1;
        }

        /* don't trace photons with zero weight */
        if (p->weight == 0.0) {
          destroy_photon (p, atmos->n_caoth);
          continue;
        }


#if HAVE_OPENGL
        if (visualize)
	  initialize_plot_photon_path ( p, sample, result, photoncounter,
					&plotpathcounter, &plotpath,
					&plotresultcounter, &plotresult );
#endif


/* ============================Start Ris Factor================================ */
/* RIS-Factor works only with lidar and radar.                                  */
/* There are 3 variables:                                                       */
/* sample->ris_factor: Value from config, does not change, not used after this. */
/* sample->ris_optical_depth: Like sample->ris_factor                           */
/* atmos->ris_factor: Is set below and used everywhere else, must not be 0,     */
/*                    defaults to 1!                                            */
/*                                   Cases:                                     */
/* A. If sample->ris_optical_depth != 0 then find_ris_factor() is called and    */
/*    atmos->ris_factor is determined and set in that function                  */
/* B. If sample->ris_factor != 0 then atmos->ris_factor = sample->ris_factor.   */
/*                                                                              */
/* The following code is executed once per lidar/radar shot and ris factor is   */
/* set anew for each shot. Otherwise it's just executed once before the first   */
/* photon starts its journey.                                                   */
/*                                                                              */

        if (photoncounter == 1) {
          if (!sample->LidarLocEst && (sample->ris_optical_depth != 0. || sample->ris_factor != 0.)) {
            fprintf(stdout, "Error, use ris factor only with lidar or radar!\n");
            return -1;
          }
          if (sample->ris_optical_depth != 0. && sample->ris_factor != 0.) {
            fprintf(stdout, "Error, use only one ris factor method at a time!\n");
            return -1;
          }
      /* Case A */
          else if (sample->ris_optical_depth != 0.) {
            find_ris_factor(p,atmos,sample,elev,refind);
            if (!quiet)
              fprintf(stderr," ... RIS-FACTOR: %e\n",atmos->ris_factor);
          }
      /* Case B */
          else if (sample->ris_factor != 0.) {
            atmos->ris_factor = sample->ris_factor;
            if (!quiet)
              fprintf(stderr," ... RIS-FACTOR: %e\n",atmos->ris_factor);
          }
        }

/* ===========================End Ris Factor=================================== */

	/* calculate direct contribution to radiation */
	status = direct_radiation (sample, atmos, p, result, elev,
				   source, refind);

        if (status)
          return err_out ("Error %d returned by direct_radiation()\n", status);

        /* register photon at its start position */
        status = register_at_start (p, sample, atmos, result,
                                    *source, *absorption);

        if (status)
          return err_out ("Error %d returned by register_at_start()\n", status);

#if HAVE_LIDAR
	if (sample->LLE_taumax != 0.0 && p->photoncounter == 1) {
	  status = mc_lidar_fix_taumax ( sample,  atmos, result,  elev, p, refind, quiet );
	  if (status>0)
	    break;
	  if (status<0)
	    return err_out ("Error %d returned by mc_lidar_fix_taumax()\n", status);
	}
#endif

#ifdef MUCHOUT
	/*      chouf=1;*/
	if (chouf) {
	  fprintf(stderr,"choufing\n");
	  p->muchoutcounter=MUCHOUTPHOTON;
	}
#endif

        /*********************/
        /* start the journey */
        /*********************/

	status = photon_journey ( p,
				  atmos,
				  sample,
				  result,
				  elev,
				  albedo,
				  surftemp,
				  absorption,
				  source, 
				  plotpath,
				  visualize,
				  wvnmlo,
				  wvnmhi,
				  refind, 
				  quiet );
	if (status<0) {
	  fprintf(stderr,
		  "Error %d occuring in photon_journey() for photon %d, exiting...\n",
		  status, p->photoncounter);
	  return -1;
	}

        outofdomaincounter += (p->photon_status == MCSTATUS_OUTOFDOMAIN);

	status = cleanup_after_photon_journey ( p, sample, result, atmos->n_caoth, chouf );
	if (status)
	  return err_out ("Error %d returned by cleanup_after_photon_journey()\n", status);

	if (sample->pan_quicklook) {
	  photoncounter++;
	  break;
	}

      } /* end of photon loop */

      totalphotoncounter += photoncounter-1;
      totalsimulations++;

#ifdef TIMINGS
      fprintf(stderr,"cputime %d %d %d\n",sample->backward_is,sample->backward_js,
	      ((int) clock()-clockstarttime)/10000);
#endif

    }  /* end of sample->backward_js loop */
    }  /* end of sample->backward_is loop */

/* Chris - special output for CB stuff, must be within ili-loop */
    if (sample->coherent_backscatter) {
#if HAVE_LIDAR
      if ( sample->LidarLocEst ) {
        summarize_result_lidar_cb (result->mish, result->lidcb, sample, 
                                   atmos->ris_factor, nphotons, basename );
      }
      else
        summarize_result_mishchenko_cb (result->mish, basename );
#else
      summarize_result_mishchenko_cb (result->mish, basename );
#endif
    }

  }  /* end of sample->ili loop */
  
  endtime = (int) time(NULL);

  if (!quiet)
    fprintf (stderr, " ... intended to trace %ld photons, traced %ld in %d simulations\n", 
        nphotons, totalphotoncounter, totalsimulations);

  if (!quiet)
    fprintf (stderr, " ... summarizing\n");


#ifdef CLDPRP
  /* summarized effective radius has to be stored for each wavelength seperatly */
  result->wavelength = 0.5*1e4*(1.0/(*wvnmlo) + 1.0/(*wvnmhi));
#endif

  /**************************************************************/
  /* calculate result arrays, write data to files and to stderr */
  /**************************************************************/

  status = summarize_result (result, sample,
      atmos, albedo,
      *source, totalphotoncounter, totalsimulations, sza,
      *absorption, write_files, write_output_as_netcdf,
      *tenstream,
      basename,
      rfldir, rfldn,  flup,
      uavgso, uavgdn, uavgup,
      rfldir3d, rfldn3d, flup3d, 
      uavgso3d, uavgdn3d, uavgup3d, 
      radiance3d, absback3d, abs3d, radiance3d_is,
      rfldir3d_var, rfldn3d_var, flup3d_var, 
      uavgso3d_var, uavgdn3d_var, uavgup3d_var, 
      radiance3d_var, absback3d_var, abs3d_var, 
      sample->escape, elev->elev2D,
      2.0/((*wvnmlo+*wvnmhi))*1e7, quiet);

  if (status)
    return err_out ("Error %d returned by summarize_result()\n", status);

  if (!quiet) 
    fprintf (stderr, " ... back from summarize_result\n");
  fflush (stderr);

  free_mystic ( result, atmos, sample, zout, *absorption, *thermal_heating_method ); /* **CK 2013.09.23 added thermal_heating_method */

  /* other stuff, unsorted */

  if (!quiet) {
    if (endtime - starttime > 0) 
      fprintf (stderr, " ... used %d seconds for tracing; %d photons per second\n", 
             endtime - starttime, 
               (int) ((double) totalphotoncounter / (double) (endtime-starttime)));
    else 
      fprintf (stderr, " ... used less than 1 second for tracing\n");
  }

  if (sample->bcond==MCBCOND_ABSORB && !sample->LidarLocEst && outofdomaincounter>0) {
    fprintf (stderr, "ATTENTION! %d/%ld = %.1f%% photons tried to leave the domain sideways\n",
             outofdomaincounter, totalphotoncounter,
	     (double) outofdomaincounter / (double) totalphotoncounter * 100.0);
    fprintf (stderr, "           or were absorbed due to numerical problems \n");
    fprintf (stderr, "           and were eliminated - the result is affected by that!\n");
  }
  
#if HAVE_OPENGL
  if (visualize) {
    fprintf (stderr, " ... waiting for 3D graphics to terminate\n");
    wait(&status);
    wait(&status);

    /* free shared memory */
    GLmystic_free_shared (GLMYSTIC_PHOTONS);
    GLmystic_free_shared (GLMYSTIC_ELEVATION);
  }
#endif

  /***** I3RC CASE 7 VERSION ******/
  /* fprintf (stderr, "ATTENTION! THIS IS A SPECIAL I3RC CASE 7 VERSION!\n"); */
  /* fprintf (stderr, "DO NOT USE FOR ANY OTHER PURPOSE\n");                  */
  /***** I3RC CASE 7 VERSION ******/

  /* restore original directions - they might be needed again */ 
  /* if MYSTIC called more than once                          */
  if (sample->Nd > 0) {
    sample->rad[0].theta = sample->forward_vza;
    sample->rad[0].phi   = sample->forward_phi;
  }

  return 0;
}


/***********************************************************************************/
/* Function: mystic_check_consistency                                     @62_30i@ */
/* Description:                                                                    */
/*  check consistency of settings for mystic, and exit if inconsistent             */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: RPB                                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int mystic_check_consistency ( sample_struct    *sample,
				      elevation_struct *elev,
				      surftemp_struct  *surftemp,
				      double            delta_scaling_mucut,
				      double            truncate,
				      int               source,
				      char             *umufilename,
				      int               ipa )
{

  if (delta_scaling_mucut>-1.0) {
    if (truncate>-1.0) {
      fprintf (stderr, "Error! Combination of delta-scaling AND additional truncation");
      fprintf (stderr, " using 'mc_truncate' not allowed. \n");
      return -1;
    } 
  }

  if (ipa!=0 && sample->reference_to_NN) {
      fprintf (stderr, "Error! Combination of reference_to_nn AND mc_ipa not allowed. \n");
      return -1;
  }

  /* adjust for periodic boundary conditions */
  switch (sample->bcond) {
  case MCBCOND_PERIODIC:
  case MCBCOND_ABSORB:
    break;

  case MCBCOND_MIRROR: 
    fprintf (stderr, "Error, mirroring boundary conditions not implemented! The reason is\n");
    fprintf (stderr, "that is does not bring any advantage to handle mirroring\n");
    fprintf (stderr, "boundary conditions internally. If you want those, please\n");
    fprintf (stderr, "mirror your input data manually!\n");
    return -1;

  default:
    fprintf (stderr, "Error, unknown boundary conditions %d\n", sample->bcond);
    return -1;
  }

  if (elev->elev2D && source==MCSRC_THERMAL_SURFACE) {
    fprintf (stderr, "Error, emission from inclined surfaces not yet included!\n");
    return -1;
  }

  if (sample->spherical3D && !(sample->backward && sample->panorama && source==MCSRC_SOLAR)) {
    fprintf (stderr, "DANGER! DANGER! DANGER! DANGER! DANGER! DANGER! DANGER! DANGER!\n");
    fprintf (stderr, " 3D spherical geometry only works in combination with backward,\n");
    fprintf (stderr, " panorama, and solar! Your photons will be entering a dimension\n");
    fprintf (stderr, " where no photon has gone before!\n");
  }

  if ( sample->spherical3D && strlen(umufilename) != 0 ) {
    fprintf (stderr, "Error, 3D spherical geometry and umu2d can not be used together!\n");
    return -1;
  }

  if (elev->elev2D && sample->spherical3D) {
    fprintf (stderr, "Error, elevation in spherical geometry not yet implemented!\n");
    return -1;
  }

  if (elev->elev2D && sample->spherical) {
    fprintf (stderr, "Error, elevation in spherical geometry not yet implemented!\n");
    return -1;
  }

  /* ??? currently 2D surface temperature is only available in backward mode; CHECK!!! */
  /* ??? it's easy to implement in forward mode as well.                      */
  if (surftemp->surf2D && source!=MCSRC_THERMAL_BACKWARD) {
    fprintf (stderr, "Error, 2D surface temperature currently only implemented for\n");
    fprintf (stderr, "backward thermal calcultions. Sorry!\n");
    return -1;
  }

//  if (sample->vroom && (sample->spectral_is || sample->concentration_is) ) {
//    fprintf (stderr, "Error, vroom probably does not work together with ALIS and AerIS (spectral_is, aerosol_is)!\n");
//    return -1;
//  }


  return 0;
}


/***********************************************************************************/
/* Function: mystic_set                                                   @62_30i@ */
/* Description:                                                                    */
/*  basic settings for mystic                                                      */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: RPB                                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int mystic_set (/* input */
		       int            nlyr,
		       int           *nmomaer,
		       int          **nthetaaer,
		       int            source,
		       int            absorption,
		       int            visualize,
		       int           *loaddata,
		       int            nzoutorg,
		       float         *zoutorg,
		       float          rpv_rho0,
		       float          rpv_k,
		       float          rpv_theta, 
		       float          hapke_h,
		       float          hapke_b0,
		       float          hapke_w,
		       float          rossli_iso,
		       float          rossli_vol,
		       float          rossli_geo,
		       float          u10,
		       float          pcl,
		       float          xsal,
		       float          bpdf_u10, 
		       char          *albedo_filename,
		       char          *alb_type_filename,
		       char          *rpvfilename,
		       char          *ambralsfilename,
		       char          *rosslifilename,
		       int            reflectalways,
		       /* throughput */
		       int           *first,
		       sample_struct *sample,
		       albedo_struct **albedo, 
		       float         *delta_scaling_mucut,
		       float         *truncate,
		       /* output */
		       int           *aerosol_scatter_type,
		       long int      *nphotons,
		       int           *nzout,
		       float        **zout,
		       float         *rfldir,
		       float         *rfldn,
		       float         *flup,
		       float         *uavgso,
		       float         *uavgdn,
		       float         *uavgup,
		       int            quiet )
{
  int kc=0, lu=0;
  int sumaer=0;

  /* HG aerosol or user-defined aerosol? */
  for (kc=0; kc<nlyr; kc++)
    sumaer += nmomaer[kc] + nthetaaer[kc][0];

  if (sumaer>0)
    *aerosol_scatter_type=MCSCAT_AER;
  else
    *aerosol_scatter_type=MCSCAT_HG1;

  if (absorption==MCFORWARD_ABS_EMISSION || sample->backward == MCBACKWARD_EMIS) {
    if (!quiet)
      fprintf (stderr, " ... calculating atmospheric emission; no need to trace any photons!\n");
    *nphotons=0;
  }


  /* **CK , 2013.08.27 - Adujsting number of photons to a multiple of 12 for backward thermal heating rates */
  if (sample->backward == MCBACKWARD_HEAT) {
    if (*nphotons % 4 != 0) 
      *nphotons = 4*(*nphotons/4+1);
    if (!quiet)
      fprintf (stderr, " ... calculating backward thermal heating rates; adjusting photon number to %ld \n", *nphotons);   
  } /* **CK  2013.08.27 End MCBACKWARD_HEAT */
  
 

  /* move this to setup_sample_grid in caoth3d.c ? */
  if (sample->escape && source==MCSRC_THERMAL_BACKWARD) {
    if (!quiet) {
      fprintf (stderr, "\n");
      fprintf (stderr, "*** Warning, backward thermal and mc_escape doesn't work!\n");
      fprintf (stderr, "*** Switching off mc_escape!\n");
      fprintf (stderr, "\n");
    }
    sample->escape=0;
    sample->vroom=0;
    if (*delta_scaling_mucut>-1.0 || *truncate>-1.0) { 
      fprintf (stderr, "%s %s\n %s %s\n",
	       "... ATTENTION, Combination of delta-scaling",
	       "or other truncation of phase function",
	       "... doesn't make much sense with backward thermal mode.",
	       "Switching off truncation!");
      *truncate=-1.0;
      *delta_scaling_mucut=-1.0;
      sample->delta_scaling = -1;
    }
  }

  /* memory allocation for static structures */
  if (*first) {
    
#if HAVE_OPENGL
    if (visualize) {
      GLmystic_init();
      GLmystic_calloc_shared (GLMYSTIC_PHOTONS, 10000);
      GLmystic_calloc_shared (GLMYSTIC_RESULT, sample->Nx*sample->Ny);
    }
#endif

    /* first call to mystic(); need to load the data in any case */
    *loaddata=1;
    *first=0;

    /* allocate memory for structs */
    *albedo = calloc (1, sizeof(albedo_struct));
  }
  else {

    /* mystic() has been called before; need to free memory if new */
    /* data are to be loaded                                       */

    if (*loaddata!=0) {
      free_mystic_for_load (*albedo, sample);

      /* allocate memory for structs */
      *albedo = calloc (1, sizeof(albedo_struct));
    }
    else
      if (!quiet)
	fprintf (stderr, " ... no need to load data again\n");
  }

  /* make a copy of zout because zout[] may be changed by mystic */
  /* and we want the original values for the next call.          */
  if (*zout!=NULL) {
    free (*zout);
    *zout=NULL;
  }
  *nzout = 0;
  if (nzoutorg>0) { /* XXX this is probably bug!!! * is missing? */
    *nzout = nzoutorg;
    *zout = calloc (*nzout, sizeof(float));
    for (lu=0;lu<*nzout;lu++)
      (*zout)[lu] = zoutorg[lu];
  }

  /* reset output */
  for (lu=0; lu<*nzout; lu++) {
    rfldir [lu] = 0.0;
    rfldn  [lu] = 0.0;
    flup   [lu] = 0.0;
    uavgso [lu] = 0.0;
    uavgdn [lu] = 0.0;
    uavgup [lu] = 0.0;
    
  }    

  /**********************************************/
  /* default settings for elevation and albedo; */
  /* can be changed by the user but may be      */
  /* overwritten by the mystic() function       */
  /* parameters                                 */
  /**********************************************/

  /* albedo->method  = MCALB_LAM2D; */      /* 2D Lambertian surface albedo  */
  (*albedo)->method  = MCALB_LAM;              /* homogeneous Lambertian albedo */

  /* ????? need more checking to avoid contraditionary specifications ????? */ 
  /* override albedo->method if RPV specified CHECK!!! */ 
  /* BCA: this is nasty, rpv_k is normally a vector, improve if statement */
  if (rpv_k!=0 || rpv_rho0!=0 || rpv_theta!=0)
    (*albedo)->method  = MCALB_RPV;
    
  /* override albedo->method if HAPKE specified */
  if (hapke_h!=0 || hapke_b0!=0 || hapke_w!=0)
    (*albedo)->method  = MCALB_HAPKE;

  /* override albedo->method if ROSSLI specified */
  if (rossli_iso!=0 || rossli_vol!=0 || rossli_geo!=0)
    (*albedo)->method  = MCALB_ROSSLI;

  /* override albedo->method if Cox and Munk specified */
  if (u10>=0 || pcl>=0 || xsal>=0)
    (*albedo)->method  = MCALB_COXANDMUNK;

  /* override albedo->method if Tsang BPDF is specified */
  if (bpdf_u10>0)
    (*albedo)->method  = MCALB_TSANG;

  /* override albedo->method if 2D albedo file specified */
  if (*albedo_filename)
    (*albedo)->method  = MCALB_LAM2D;

  /* override albedo->method if 2D spectral albedo file specified */
  if (*alb_type_filename)
    (*albedo)->method  = MCALB_LAM2D_SPECTRAL;

  /* override albedo->method if 2D RPV file specified;      */
  /* it is important that this comes AFTER the Cox_and_Munk */
  /* check because we allow MCALB_RPV2D_SPECTRAL with ocean pixels.  */
  if (*rpvfilename)
    (*albedo)->method  = MCALB_RPV2D_SPECTRAL;

  /* override albedo->method if 2D ROSSLI file specified */
  if (*ambralsfilename)
    (*albedo)->method  = MCALB_ROSSLI2D;

  if (*rosslifilename)
    (*albedo)->method  = MCALB_ROSSLI2D;

  /* reflect each photon and weight with albedo or reflect only a fraction of "albedo" */
  (*albedo)->reflectalways = reflectalways;

  /* CP - Hotfix to avoid rare error from mc_vroom_DDIS_weight_and_prep_stuff()
     over reflect() despite vroom being off and albedo being zero */
#if HAVE_LIDAR
  if (sample->LidarLocEst == MCRADAR && !(sample->vroom) && (*albedo)->albedo == 0.)
    (*albedo)->reflectalways = 0;
#endif

  switch ((*albedo)->method) {
  case MCALB_LAM:
  case MCALB_LAM2D:
  case MCALB_LAM2D_SPECTRAL:
    break;

  case MCALB_RPV:
  case MCALB_RPV2D_SPECTRAL:
  case MCALB_COXANDMUNK: 
  case MCALB_HAPKE:
  case MCALB_ROSSLI:
  case MCALB_ROSSLI2D:
  case MCALB_TSANG:
    if (!(*albedo)->reflectalways) {
      (*albedo)->reflectalways=1;
      if (!quiet) 
        fprintf (stderr, " ... switching on reflectalways\n");
    }
    break;
  default:
    fprintf (stderr, "Fatal error, unknown albedo type %d\n", (*albedo)->method);
    return -1;
  }

  return 0;
}


/***********************************************************************************/
/* Function: mystic_verbose                                               @62_30i@ */
/* Description:                                                                    */
/*  talk to user: tell user what is being done                                     */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: RPB                                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int mystic_verbose ( int            quiet,
			    sample_struct *sample,
			    albedo_struct *albedo,
			    int            aerosol_scatter_type,
			    long int       nphotons,
			    int            ipa,
			    double         delta_scaling_mucut,
			    double         truncate )
{
  static int firstwarning=1;

  if (!quiet) {
    if (ipa)
      fprintf (stderr, "\n ... start independent pixel approximation!\n");
    else
      fprintf (stderr, "\n ... start three-dimensional calculation\n");

    if (aerosol_scatter_type==MCSCAT_AER)
      fprintf (stderr, " ... user-defined aerosol properties\n");


    /***** I3RC CASE 7 VERSION *****/
    /* fprintf (stderr, "ATTENTION! THIS IS A SPECIAL I3RC CASE 7 VERSION!\n"); */
    /* fprintf (stderr, "DO NOT USE FOR ANY OTHER PURPOSE\n");                  */
    /***** I3RC CASE 7 VERSION ******/

    if (delta_scaling_mucut>-1.0) {
      fprintf (stderr, " ... delta-scaling ... truncating phase function at mu = %f.\n",
	       delta_scaling_mucut);
      fprintf (stderr, "     USE ONLY IF YOU REALLY, REALLY KNOW WHAT YOU ARE DOING!\n");
    }

    switch (sample->bcond) {
    case MCBCOND_PERIODIC:
      fprintf (stderr, " ... periodic boundary conditions\n");
      break;

    case MCBCOND_MIRROR: 
      fprintf (stderr, " ... mirroring boundary conditions\n");
      return -1;

    case MCBCOND_ABSORB:
      fprintf (stderr, " ... absorbing boundary conditions\n");
      break;

    default:
      fprintf (stderr, "Error, unknown boundary conditions %d\n", sample->bcond);
      return -1;
    }

#if HAVE_VROOM
    if (sample->escape && sample->escape_eps_ddis_upf)
      fprintf(stderr,"DDISUPF for escape (%e)",sample->escape_eps_ddis_upf);
#endif

    fprintf (stderr, " ... running %ld photons\n", nphotons);

  } /* end if quiet */

  /* if truncate is true, the phase function will be truncated     */
  /* but only for the calculation of the escape probability after  */
  /* the first scattering event; this helps reducing hot spots     */
  if (truncate>-1.0) {
    fprintf (stderr, " ... truncating phase function at mu = %f.\n", truncate);
    fprintf (stderr, " ... USE ONLY IF YOU KNOW WHAT YOU ARE DOING!\n");
  }    

  if (firstwarning==1) {

    if ((sample->spherical3D || sample->spherical3D_scene) && !quiet) {
      fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      fprintf (stderr, "!!! Warning! Spherical3D options are still very       !!!\n");
      fprintf (stderr, "!!! experimental! It is not sure whether what you are !!!\n");
      fprintf (stderr, "!!! trying to calculate will work! Please be extremely!!!\n");
      fprintf (stderr, "!!! careful with interpreting your results, and talk  !!!\n");
      fprintf (stderr, "!!! to Claudia Emde about whether what you are        !!!\n");
      fprintf (stderr, "!!! calculating can work at all!                      !!!\n");
      fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    }

    if (sample->LidarLocEst && sample->polarisation) {
      fprintf (stderr, "Warning, polarisation not yet validated for Lidar! Be careful with results!\n");
    }

    if (sample->LidarLocEst) {
      switch (albedo->method) {
      case MCALB_LAM:
      case MCALB_LAM2D:
      case MCALB_LAM2D_SPECTRAL:
	break;

      case MCALB_RPV:
      case MCALB_RPV2D_SPECTRAL:
      case MCALB_COXANDMUNK: 
      case MCALB_HAPKE:
      case MCALB_ROSSLI:
      case MCALB_ROSSLI2D:
      case MCALB_TSANG:
	fprintf (stderr, "Warning, BRDF not yet validated for Lidar! Be careful with results!\n");
	break;
      default:
	fprintf (stderr, "Fatal error, unknown albedo type %d\n", albedo->method);
	break;
      }
    }

    if (sample->vroom && sample->spherical3D && !quiet) {
      fprintf (stderr, "Warning, vroom not yet validated for spherical geometry! Be careful with results!\n");
    }

    if (sample->vroom && sample->spherical && !quiet) {
      fprintf (stderr, "Warning, vroom not yet validated for spherical geometry! Be careful with results!\n");
    }

    if (sample->refraction) {
      fprintf (stderr, "Warning, refraction not yet working for reflection! Be careful with results!\n");
    }

    if (!sample->pan_distr_photons_over_pixel && sample->pan_with_direct_rad && !quiet){
      fprintf (stderr, "\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
      fprintf (stderr,   " !!!! ATTENTION, it may be a good idea to use   !!!!!!\n");
      fprintf (stderr,   " !!!! mc_panorama_distr_photons_over_pixel when !!!!!!\n");
      fprintf (stderr,   " !!!! using mc_panorama_with_direct_rad.        !!!!!!\n");
      fprintf (stderr,   " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
    }  

    if (sample->planar_switch) {
      fprintf (stderr, "Warning, you are using experimental options! They can be wrong! Be careful!\n");
    }

    if (sample->spherical && !quiet) {
      fprintf (stderr, "ATTENTION: For spherical MonteCarlo, the photons are started\n");
      fprintf (stderr, "           at the central point of the domain (not in the central\n");
      fprintf (stderr, "           sample pixel, as usual).\n");
    }

    firstwarning=0;
  }

  return 0;
}


/***********************************************************************************/
/* Function: check_premature_exit                                         @62_30i@ */
/* Description:                                                                    */
/*  check whether user has decided to exit prematurely with output of stuff so far */
/*  calculated.                                                                    */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: RPB                                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/

#ifdef ALLOW_FILE_CHECKS
static int check_premature_exit ( char          *filename,
				  sample_struct *sample )
{
  FILE *file;

  file = fopen(filename,"r");
  if (file!=NULL) {
    fclose(file);

    sample->ili = sample->Nli;
    sample->backward_is = sample->backward_isupper;
    sample->backward_js = sample->backward_jsupper;
    fprintf(stderr,"Exiting prematurely...\n");
    return 1;
  }
  return 0;
}
#endif


/***********************************************************************************/
/* Function: check_premature_exit                                         @62_30i@ */
/* Description:                                                                    */
/*  check whether user has decided to exit prematurely with output of stuff so far */
/*  calculated.                                                                    */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: RPB                                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/

#ifdef ALLOW_FILE_CHECKS
static void check_progress ( long int       photoncounter,
			     long int       nphotons,
			     sample_struct *sample,
			     float          wavelength )
{
  FILE *file;

  file = fopen("progress","r");
  if (file!=NULL) {
    fclose(file);
    fprintf(stderr,
	    "%s %e, is: %d of [%d-%d], js: %d of [%d-%d], ili %d of %d, nphots %ld of %ld\n",
	    "Progress reached: wvl: ",
	    wavelength,
	    sample->backward_is, sample->backward_islower, sample->backward_isupper,
	    sample->backward_js, sample->backward_jslower, sample->backward_jsupper,
	    sample->ili, sample->Nli,
	    photoncounter, nphotons);
  }
  return;
}
#endif


/***********************************************************************************/
/* Function: read_and_write_random_status                                 @62_30i@ */
/* Description:                                                                    */
/*  check whether to read/write random status.                                     */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: RPB                                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/

#if HAVE_LIBGSL 
static inline int read_and_write_random_status ( sample_struct *sample,
						 char          *randomstatusfilename,
						 int            photoncounter,
						 int           *readrandomstatus )
{
  FILE *randomstatusfile=NULL;
  int status=0;

  /************************************+**************************/
  /* starting with the randomseed is nice, but starting with the */
  /* actual photon that produced the problem is even nicer...    */
  /* here we read and write the status of the random number      */
  /* generator                                                   */
  /* this only works with the gsl !!!                            */
  /************************************+**************************/

  /* read randomstatusfile if asked for */
  if (*readrandomstatus && photoncounter==1) {
    if ( ( randomstatusfile = fopen(randomstatusfilename, "r")) != NULL) {
      status = gsl_rng_fread (randomstatusfile, uvspecrng);
      if (status)
	return err_out ("Error %d returned by gsl_rng_fread()\n", status);
      if ( ! fscanf(randomstatusfile,"%d %d %d",&sample->ili,&sample->backward_is,&sample->backward_js) ) {
	fprintf (stderr, "Error reading file %s\n", randomstatusfilename);
	return -1;
      }
      (void) fclose (randomstatusfile);
      fprintf(stderr,"reading random status for photon...\n");
    }
    else
      return err_out ("Error %d, could not open randomstatusfile\n", -1);
  }

#ifdef WRITERANDOMSTATUS
  /* write status every 10 photons, starting with photon 1 */
  if (photoncounter % 10 == 1 && !(*readrandomstatus)) {
    if ( ( randomstatusfile = fopen(randomstatusfilename, "w")) != NULL) {
      status = gsl_rng_fwrite (randomstatusfile, uvspecrng);
      if (status)
	return err_out ("Error %d returned by gsl_rng_fwrite()\n", status);
      fprintf(randomstatusfile,"%d %d %d\n",sample->ili,sample->backward_is,sample->backward_js);
      (void) fclose (randomstatusfile);
    }
    else
      return err_out ("Error %d, could not open randomstatusfile\n", -1);
  }
#endif
  *readrandomstatus=0;

  return 0;
}
#endif


/***********************************************************************************/
/* Function: initialize_plot_photon_path                                  @62_30i@ */
/* Description:                                                                    */
/*  check whether to plot photon path of current photon.                           */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: RPB                                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/

#if HAVE_OPENGL
static void initialize_plot_photon_path ( photon_struct *p,
					  sample_struct *sample,
					  result_struct *result,
					  int            photoncounter,
					  int           *plotpathcounter,
					  int           *plotpath,
					  int           *plotresultcounter,
					  int           *plotresult )
{

  /* start new photon path if first photon or a multiple of 10000 */
  *plotpath=0;
  (*plotpathcounter)++;
  if ( *plotpathcounter==10000 || photoncounter==1 ) {
    *plotpathcounter=0;
    *plotpath=1;
    
    GLmystic_start_photon_path (p);
  }

  /* copy global surface irradiance to result field for plotting, */
  /* either for the first photon or for every 100000 photons      */
  *plotresult=0;
  (*plotresultcounter)++;
  if ( *plotresultcounter==100000|| photoncounter==1 ) {
    *plotresultcounter=0;
    *plotresult=1;

    GLmystic_write_global (p, sample, result);
  }
}
#endif


/***********************************************************************************/
/* Function: cleanup_after_photon_journey                                 @62_30i@ */
/* Description:                                                                    */
/*  save and free stuff after photon journey.                                      */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: RPB                                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int cleanup_after_photon_journey ( photon_struct *p,
					  sample_struct *sample,
					  result_struct *result,
					  int            n_caoth,
					  int            chouf )
{
#if defined MUCHOUT || defined SPIKEOUT || defined ESCAPES
#if 0
  static double oldres=0.0;

  /*	if (chouf) return -1; */
  double churn=result->back[sample->backward_is][sample->backward_js][0]-oldres;
  /*	if (sample->backward_is==76)*/
  fprintf(stderr,"delta %30.20f %d %d\n",churn,sample->backward_is,p->photoncounter);
#if defined SPIKEOUT
  if (result->back[sample->backward_is][sample->backward_js][0]-oldres
      > 145220. && chouf ) {
    fprintf(stderr,"found spike...\n");
    return -1;
  }
  if (churn==7.11312911017375881784)
    chouf=1;
#endif
  oldres=result->back[sample->backward_is][sample->backward_js][0];
#endif
#endif

  treeprint_stddev (p->wtree);  /* calculate standard deviation */
  free(p->wtree);
  p->wtree=NULL;
#if HAVE_LIDAR
  /* special treatment for standard deviation of jacobians */
  if (sample->LLE_jacobian_std)
    lidar_collect_jacobian_stddev (result->lidar, sample->LLE_Nt, sample->ili, n_caoth);
  if (sample->abs_jacobian_std)
    abs_collect_jacobian_stddev (result->jacobian, sample->backward_is, sample->backward_js);
#endif

#ifdef PRINT_PHOTON_PATH
  /* print photon path to stdout and free memory */
  print_and_clear_photonpath (p->path);

  /* photon path separator */
  fprintf (stdout, "! x x\n");
#endif

  destroy_photon (p, n_caoth);

  return 0;
}


/***********************************************************************************/
/* Function: free_mystic                                                  @62_30i@ */
/* Description:                                                                    */
/*  Deallocate memory after mystic has ended.                                      */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: RPB                                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static void free_mystic ( result_struct     *result,
			  atmosphere_struct *atmos,
			  sample_struct     *sample,
			  float             *zout,
			  int                absorption,
			  int                thermal_heating_method) /* **CK 2013.09.23 added thermal_heating_method */
{
  int kc=0;

  /* Free the result structure; this structure is allocated and  */
  /* zeroed during each call of mystic, other than the other     */
  /* structures which are allocated only if loaddata!=0. Reason: */
  /* freeing and reallocating is a convenient way to reset the   */
  /* structures to 0 without forgetting anything.                */

  free_result (result, thermal_heating_method, atmos, sample, absorption); /* **CK 2013.09.23 added thermal_heating_method */

  free(zout);

  if ( sample->escape_eps_ddis_upf || sample->LLE_D_DIS || sample->n_phase_max!=0) {
    free_iphase(sample->phase_max[0]);
    free(sample->phase_max[0]);
    free(sample->phase_max);
    free_hybrid3D_field (&(atmos->spiky_box), atmos);
  }

  if (atmos->phase_aer!=NULL) {
    for (kc=0; kc<atmos->Nz; kc++)
      free_pft (&((atmos->phase_aer)[kc]));
    free (atmos->phase_aer);
  }

  /* free stuff in functions */
  get_phase_matrix_total ( atmos, NULL, 0.0, 0, 0, 0, 0, 0, 1, NULL, NULL, NULL, NULL );
  mc_vroom_scattering_calc_phases_and_jacobians ( NULL, NULL, atmos, 0.0, 1, NULL, NULL, NULL );

  free_atmosphere (atmos);
  free (atmos);
}


/***********************************************************************************/
/* Function: direct_radiation                                             @62_30i@ */
/* Description:                                                                    */
/*  Calculate direct contribution to radiation.                                    */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: RPB                                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int direct_radiation (sample_struct     *sample,
		      atmosphere_struct *atmos,
		      photon_struct     *p,
		      result_struct     *result,
		      elevation_struct  *elev,
		      int               *source, 
		      float             *refind)
{
  int id=0, status=0, j=0;
  double norm[3] = {0,0,0};
  double cotheta=0, slope=0.0, alpha2=0.0;
  photon_struct *p_escape=NULL;

  /* radiance contribution of the new photon at its start position */
  if (sample->escape && p->scattercounter >= sample->minscatters) {

    switch (sample->backward) {
    case MCBACKWARD_EGLOB:
    case MCBACKWARD_EDIR:
    case MCBACKWARD_FDIR:

      /* in case of direct irradiance we only count the escape     */
      /* radiance of the starting photon; here we need to consider */
      /* the projection of the surface area perpendicular to the   */
      /* photon direction                                          */

      /* determine upward normal to the surface      */
      calc_normal_vector ( p, elev, sample, atmos,
			   sample->surfaceparallel && sample->zstart==-999.0,
			   sample->sensordirection && sample->zstart==-999.0,
			   norm );
         
      v_mult ( norm, sample->rad[0].dir.dx, &cotheta );
          
      /* no contribution to direct irradiance if sensor looking downward */
      /* CE: not for spherical atmosphere! CHECK!!!*/ 
      /* RPB: think! If norm is calculated correctly, ok for spherical atmosphere! */
      if (cotheta > 0 || sample->spherical || sample->spherical3D ){

	p_escape = create_escape_photon ( p, atmos, sample, 0 );

	p_escape->pdir[0] = 4.0 * fabs(cotheta);
              
	status = escape_probability (p_escape, atmos,
				     sample, result,
				     elev, 0, -1, -1,
				     1, refind);

	if (status!=0)
	  return err_out ("Error, escape_probability() returned status %d\n", status);

	/* free memory */
	destroy_photon (p_escape, atmos->n_caoth);
      }
            
      break;
          
    default:
      break;
    }

    switch (*source) {
    case MCSRC_THERMAL_ATMOSPHERE:
      for (id=0; id<sample->Nd; id++) {
	p_escape = create_escape_photon ( p, atmos, sample, id );

	p_escape->pdir[id] = 1.0;

	status = escape_probability (p_escape, atmos,
				     sample, result,
				     elev, id, -1, -1,
				     0, refind);

	if (status!=0)
	  return err_out ("Error, escape_probability() returned status %d\n", status);

	/* free memory */
	destroy_photon (p_escape, atmos->n_caoth);
      }
      break;
      
    case MCSRC_THERMAL_SURFACE:
      for (id = 0; id < sample->Nd; id++) {
	p_escape = create_escape_photon ( p, atmos, sample, id );
        
	/* ??? the following statement is not correct     ??? CHECK!!! */ 
	/* ??? for thermal emission from inclined surface ??? */
	if (sample->rad[id].dir.dx[2] > 0) /* upward   */
	  p_escape->pdir[id] = 4.0 * sample->rad[id].dir.dx[2]; 
	else                               /* downward */      
	  p_escape->pdir[id] = 0; 
          
	status = escape_probability (p_escape, atmos,
				     sample, result,
				     elev, id, -1, -1,
				     1, refind);

	if (status!=0)
	  return err_out ("Error, escape_probability() returned status %d\n", status);

	/* free memory */
	destroy_photon (p_escape, atmos->n_caoth);
      }
      break;

    case MCSRC_NONE:
      break;
    case MCSRC_SOLAR:
      /* Do a 0-order local estimate (no scattering is taking place) */
      /* to account for DIRECT RADIATION IN PANORAMA */
      if (sample->pan_with_direct_rad) {
	/* At this point we should still be at the location where the */
	/* photon is generated. */

	if ( sample->sample_backward_sunshape ) {

	  if (sample->pan_alignment!=MCPAN_ALIGNMENT_SUN) {
	    fprintf (stderr, "Error, pan_with_direct_rad with sunshape does not work without pan_alignment sun! Exiting...\n");
	    return -1;
	  }

	  v_mult_mu ( p->dir.dx, sample->rad_dir0.dx, &cotheta );

	  /* Check if photon direction is pointing somewhere inside sun */
	  /* disc and if above horizon                                  */
	  if ( ( cotheta >= cosd(sample->sun_radius) ) && p->dir.dx[2] >= 0 ){
	    p_escape = create_escape_photon ( p, atmos, sample, 0 );

	    /* alpha2: relative position on sun disc. scalar product of */
	    /* direction to sun center and photon direction gives the   */
	    /* cosine of the angle between them. division by sunradius  */
	    /* gives the relative position on sun disc                  */
	    alpha2 = acosd(cotheta) / sample->sun_radius;

	    /* find nearest alpha */
	    j = locate ( sample->sample_backward_sunshape_alpha,
			 sample->sample_backward_sunshape_n_lines,
			 alpha2 );

	    /* slope is the slope of the linear interpolation */
	    slope = sample->sample_backward_sunshape_slope [j];

	    /* Calc. probability of photon originating from this spot on the sun */  
	    /* Interpolate sunshape to actual alpha */

	    p_escape->pdir[0] = sample->sample_backward_sunshape_p[j] +
	      slope * ( alpha2 - sample->sample_backward_sunshape_alpha[j] );

	    /* We want radiances, so account for the solid angle of the sun */

	    p_escape->pdir[0] *= 4.0 * PI / ( 2. * PI  * ( 1.0 - cosd ( sample->sun_radius ) ) );

	    /* In count_escape the weight is divided by              */
	    /* sample->rad[id].dir.cotheta but instead it should be  */
	    /* divided by the cotheta of the direction to sun center */
	    /* (sample->rad_dir0.cotheta) we account for that by	   */
	    /* modifying p->pdir[0] accordingly                      */

	    id=0;
	    p_escape->pdir[0] *= sample->rad[id].dir.cotheta / sample->rad_dir0.cotheta;

	    status = escape_probability (p_escape, atmos, sample, result, elev,
					 0, -1, 0, 0, refind);
	    if (status!=0)
	      return err_out ("Error, escape_probability() returned status %d\n", status);

	    /* free memory */
	    destroy_photon (p_escape, atmos->n_caoth);
	  }
	}
	else {
	  p_escape = create_escape_photon ( p, atmos, sample, 0 );

	  p_escape->pdir[0]=1.0;

	  if (sample->pan_weight_with_cos) {
	    v_mult_mu( p->dir00.dx, sample->rad_dir0.dx , &cotheta );
	    p_escape->pdir[0] *= 4.0*PI*cotheta;
	  }

	  status = escape_probability (p_escape, atmos, sample, result, elev,
				       0, -1, 0, 0, refind);
	  if (status!=0)
	    return err_out ("Error, escape_probability() returned status %d\n", status);

	  /* free memory */
	  destroy_photon (p_escape, atmos->n_caoth);
	}
      }
      break;
    case MCSRC_BLITZ:
      if (sample->pan_with_direct_rad) {
	p_escape = create_escape_photon ( p, atmos, sample, 0 );
	p_escape->pdir[0] = 1.0;
	status = escape_probability (p_escape, atmos, sample, result, elev,
				     0, -1, 0, 0, refind);
	if (status!=0)
	  return err_out ("Error, escape_probability() returned status %d\n", status);
	/* free memory */
	destroy_photon (p_escape, atmos->n_caoth);
      }
      break;
    default:
      fprintf (stderr, "Error, unknown source type\n");
      return -1;
    }
  }

  if (sample->pan_weight_with_cos)
    p->weight*=-2.0*PI*p->weight1;

  return 0;
}


/***********************************************************************************/
/* Function: photon_journey                                               @62_30i@ */
/* Description:                                                                    */
/*  Handle the journey of one photon.                                              */
/*  The function loops infinitely over the switch, which does something according  */
/*  to the current status of the photon. Each switch must change the status of     */
/*  the photon to a different status, else the loop does not work.                 */
/*  Two cases lead to an exit from the loop, OUTOFDOMAIN and PURGE. PURGE should   */
/*  be the standard way to exit the loop, OUTOFDOMAIN is counted outside photon_   */
/*  journey and leads to a warning.                                                */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: RPB                                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int photon_journey ( photon_struct     *p,
		     atmosphere_struct *atmos,
		     sample_struct     *sample,
		     result_struct     *result,
		     elevation_struct  *elev,
		     albedo_struct     *albedo,
		     surftemp_struct   *surftemp,
		     int               *absorption,
		     int               *source, 
		     int                plotpath,
		     int                visualize,
		     float             *wvnmlo,
		     float             *wvnmhi,
		     float             *refind,
		     int                quiet )
{
  int status=0, ie=0, je=0, ip=0, old_photon_status=MCSTATUS_DEFAULT;
  double tau=0.0, zsurf=0.0, slant2horz=0.0;
  double totweight[4];

  /* nothing to do if direct irradiance is to be calculated backwards */
  /* in this case only the escape radiance from the start counts      */
  if ( sample->maxscatters == 0 ||
       sample->backward == MCBACKWARD_EDIR ||
       sample->backward == MCBACKWARD_FDIR ||
       sample->backward == MCBACKWARD_EMIS)
    return 0;

  while (1==1)   {

#ifdef MUCHOUT
    if (p->muchoutcounter==MUCHOUTPHOTON) {
      fprintf(stderr,"counters %d %d %d %d --- photonstatus %d isclone %d \n",
	      p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
	      p->photon_status,p->isclone);
    }
#endif

    old_photon_status = p->photon_status;

    /**************************************************/
    /* process photon according to its current status */
    /**************************************************/

    switch(p->photon_status) {
    case MCSTATUS_TRAVEL:
      
      /* random step, in optical depth units        */
      if (sample->pan_quicklook)
	tau = 0.1;
      else
	tau = random_tau();

      /********************************/
      /* now travel optical depth tau */
      /********************************/

      status = travel_tau (p, atmos,
			   sample, result,
			   elev, tau, *absorption, *source, plotpath,
			   visualize, refind, quiet);

      if (status<0)
	return err_out ("Error %d returned by travel_tau()\n", status);
      
      /***** I3RC CASE 7 VERSION ******/
      /* if (p->pathlength > 15360) */
      /*   break;                   */
      /***** I3RC CASE 7 VERSION ******/

#ifdef PRINT_PHOTON_PATH
      /* add photon to photon path */
      add_to_photonpath (&(p->path), p->x);
#endif

#if HAVE_OPENGL
      if (visualize)
	if (plotpath)
	  GLmystic_add_node_to_photon_path (p->x[0], p->x[1], p->x[2]);
#endif

      p->photon_status = status;

      break;

    case MCSTATUS_SCATTER:
      /**************************************/
      /* ... normal scattering, no boundary */
      /**************************************/

#if HAVE_LIDAR
      /* save tauext for CB */
      if (sample->coherent_backscatter && sample->LidarLocEst)
        add_tauext_to_cohebasca ( p, sample, result->lidcb );
#endif

      status = scattering (atmos, p, sample, result, elev,
			   albedo, *wvnmlo, *wvnmhi, refind);

      if (status<0)
	return err_out ("Error, scattering() returned status %d\n", status);

      if (status==MCSTATUS_VIRTUALSCATTER)
	p->photon_status=MCSTATUS_TRAVEL;
      else
	p->photon_status = status;
      break;
      
    case MCSTATUS_ABSORB:
      /***********************************************/
      /* ... real absorption event in the atmosphere */
      /***********************************************/
      /*TZ bt...*/
      
      /* contribution is Bplanck(T) in this case */ 
      count_thermal_backward_photon (result, p, atmos, surftemp, sample, 
				     0.0, albedo->albedo,
				     quiet, p->photon_status);

      /* now kill photon */
      p->photon_status=MCSTATUS_PURGE;

      break;

    case MCSTATUS_SURFACE:
      /*********************************************************************/
      /* ... 2D surface hit preparation:                                   */
      /* adjust x[2] if photon hit the 2D surface to avoid roundoff errors */
      /*********************************************************************/
      
      /* elevation coordinates */
      elev_coord (p, elev, &ie, &je);

      /* final consistency check, if the photon really hit the surface */
      zsurf = elevation (elev->surf[ie][je],
			 p->x[0] - (double) ie * elev->delX, 
			 p->x[1] - (double) je * elev->delY);
      
      /* here we test if the altitude deviates more than      */
      /* MC_EPSILON * maximum altitude from the surface;      */
      /* previously, MC_EPSILON was used as threshold, but    */
      /* this occasionally produced FATAL errors; are we      */
      /* still conservative enough?                           */

      if (fabs(p->x[2] - zsurf) > MC_EPSILON * elev->surfmax) {
	fprintf (stderr, "FATAL error:\n");
	fprintf (stderr, "cross_surface() found that photon hit the surface at\n");
	fprintf (stderr, "x = %g, y = %g, z = %.10f\n", p->x[0], p->x[1], p->x[2]);
	fprintf (stderr, "while the surface elevation at this xy-location\n");
	fprintf (stderr, "is actually %.10f.\n", zsurf);
        
	return -1;
      }

      /* and adjust x[2] to the "true" surface altitude */
      p->x[2] = zsurf;

      /* no break */

    case MCSTATUS_BOUNDARY_LOWER:
      /************************************************************/
      /* ... if photon arrived at the surface or lower boundary   */
      /************************************************************/

      /* register photon at the lower boundary */
      if (*source!=MCSRC_THERMAL_BACKWARD ) {
          
	/* conversion factor to correct for slant surface */
	slant2horz = slt2hrz (elev, p, sample->surfaceparallel, 1);
	if (slant2horz<0.0)
	  return fct_err_out (-1, "slt2hrz", ERROR_POSITION);

        for (ip=0; ip<sample->nstokes; ip++)
          totweight[ip] = p->weight * p->stokes[ip] *
            exp(-p->tauabs.tot) * slant2horz; 
        
        /* register incoming photon */
	count_photon (result->surf, p, sample, 
		      totweight, sample->surfaceparallel, 1);
      }

      /* reflect photon at the lower boundary */
        
      /**********************/
      /* surface reflection */
      /**********************/
      status = reflect (p, albedo, elev, atmos, sample, result,
			*wvnmlo, *wvnmhi, refind);              

      p->reflectcounter=1;

      if (status<0)
	return err_out ("FATAL error! reflect() returned status %d\n", status);

      /* don't do anything more if photon was destroyed */
      if (status==MCSTATUS_PURGE) {
	p->photon_status = status;
	break;
      }

      /* if photon was absorbed by surface: */
      if (status==MCSTATUS_ABSORB) {
	/* count contribution of surface emission, if no reflection */
        if (*source ==MCSRC_THERMAL_BACKWARD) {    /*TZ bt ...*/
          count_thermal_backward_photon (result, p, atmos, surftemp, sample,
					 albedo->bplkavg, albedo->albedo,
					 quiet, p->photon_status);
	  /* contribution is here Bplanck(T) */
	}

	/* forward+solar backward: kill photon if no reflection */
	p->photon_status=MCSTATUS_PURGE;
	break;

      }/*... TZ bt*/ 

      if (*source!=MCSRC_THERMAL_BACKWARD) {/*TZ bt*/ 
	/* conversion factor to correct for slant surface */
	slant2horz = slt2hrz (elev, p, sample->surfaceparallel, 0);
	if (slant2horz<0.0)
	  return fct_err_out (-1, "slt2hrz", ERROR_POSITION);
              
	/* register reflected photon */
        for (ip=0; ip<sample->nstokes; ip++)
          totweight[ip] = p->weight * p->stokes[ip] * exp(-p->tauabs.tot)
            * slant2horz;
        
	count_photon (result->surf, p, sample, 
		      totweight, sample->surfaceparallel, 0);
              
	/* if the reflection was off the 1D lower boundary */
	/* register also into altitude profile             */
	if (p->photon_status == MCSTATUS_BOUNDARY_LOWER && sample->sample[0]!=0)
	  count_photon (result->alt[0], p, sample, 
			totweight, 0, 0);
      }/*TZ bt*/

      p->photon_status=MCSTATUS_INCREMENT;

      break;

    case MCSTATUS_BOUNDARY_UPPER:
      /**************************************************/
      /* ... if photon arrived at the upper boundary    */
      /**************************************************/
  
      /* leave TOA */
      if (*source==MCSRC_THERMAL_BACKWARD) {/*TZ bt ...*/
	count_thermal_backward_photon (result, p, atmos, surftemp, sample,
				       0.0, albedo->albedo, 
				       quiet, p->photon_status);
	/* contribution is here Bplanck*emissivity. T(TOA)=0 => Bplanck=0*/
      }/*... TZ bt*/

      /* now kill photon */
      p->photon_status=MCSTATUS_PURGE;

      break;

    case MCSTATUS_INCREMENT:
      /*********************************************************************************/
      /* photon has finished scattering/reflecting, now increment scattercounters etc. */
      /*********************************************************************************/
      /* is not done if virtual scatter (see VIS-FOD) */
      /* surface reflection is also counted           */

      /* if escapescattercounter was smaller or equal than scattercounter then escape was done */
      if (p->escapescattercounter <= p->scattercounter)
	p->escapescattercounter++;

      if (p->escapescattercounter <= p->scattercounter) {
	fprintf(stderr,"Error! escape seems to have skipped a von-Neumann-order! %d %d\n",p->escapescattercounter, p->scattercounter);
	return -1;
      }

      p->scattercounter++;

      if (p->isclone)
	p->clonescattercounter++;

      p->direct=0;
      if (sample->tipadir==3) /* should be TIPA_DIR3D */
	p->ipa=1;

      if ( sample->delta_scaling>0 && (p->scattercounter >= sample->delta_scaling) )
	p->SC_mode = MCSC_MODE_DELTA_SCALE;
      if ( sample->delta_scaling>0 && (p->scattercounter >= sample->delta_scaling) )
	p->DDIS_SC_mode = MCSC_MODE_DELTA_SCALE;

      if (sample->maxscatters > 0)
	if (p->scattercounter >= sample->maxscatters) {
	  p->photon_status=MCSTATUS_PURGE;
          break;
        }

      if (sample->pan_quicklook && p->scattercounter>0) {
	p->photon_status=MCSTATUS_PURGE;
	break;
      }

      p->photon_status=MCSTATUS_TRAVEL;

      /* go to cloning if VROOM is on */
      if (sample->ntupelLE)
	p->photon_status=MCSTATUS_CLONE;
      else
	p->photon_status=MCSTATUS_TRAVEL;

      break;

#if HAVE_VROOM
    case MCSTATUS_CLONE:

      status = mc_vroom_cloning (p, atmos, sample, result, elev, albedo,
				 surftemp, absorption, source, plotpath, visualize,
				 wvnmlo, wvnmhi, refind, quiet);
      if (status<0)
	return err_out ("FATAL error! mc_vroom_cloning() returned status %d\n", status);

      p->photon_status = MCSTATUS_SPLIT;
      break;

    case MCSTATUS_SPLIT:


      status = mc_vroom_splitting_and_rr (p, atmos, sample, result, elev, albedo,
					  surftemp, absorption, source, plotpath, visualize,
					  wvnmlo, wvnmhi, refind, quiet);
      if (status<0)
	return err_out ("FATAL error! mc_vroom_splitting_and_rr() returned status %d\n", status);

      p->photon_status = status;

      /* after split comes travel */
      if (p->photon_status == MCSTATUS_DEFAULT)
	p->photon_status = MCSTATUS_TRAVEL;

      break;
#endif

    case MCSTATUS_DEFAULT:     /* should not happen here */
    case MCSTATUS_OUTOFDOMAIN: /* should not happen here */
    case MCSTATUS_PURGE:       /* should not happen here */
    case MCSTATUS_VIRTUALSCATTER: /* uld not happen here */
    default:
      fprintf(stderr,"Error in photon_journey() loop, no valid photon_status %d\n",p->photon_status);
      return -1;
    }

    /* if the photon status did not change in the switch, something is wrong */
    /* the photon will be in the same status in the next switch, and might   */
    /* stay so infinitely. Hence, exit is better.                            */
    if (old_photon_status == p->photon_status) {
      fprintf(stderr,"Error, photon_status %d did not change, this potentially means infinite loop!\n",p->photon_status);
      return -1;
    }

    /* ============ Don't trace photons with zero weight! ============== */
    /* Testing exp(-p->tauabs.tot) == 0 should not be done because the   */
    /* exp function is too expensive to calculate, therefore we check    */
    /* p->tauabs.tot >= 709. since exp(-709) = 1.217e-308, which is near */
    /* the lowest possible double precision number, so it's save to kill */
    /* photons if they ever reach that point                             */
    if (p->weight == 0. || p->stokes[0] == 0. || p->tauabs.tot >= 709. ) 
      p->photon_status = MCSTATUS_PURGE;

    /**************************************************************/
    /* stop if photon absorbed, in space, or due to other reasons */
    /**************************************************************/
    if (p->photon_status == MCSTATUS_PURGE)
      break;

    /* exit loop if photon out of domain */
    if (p->photon_status == MCSTATUS_OUTOFDOMAIN)
      break;

  } /* end while */

  return 0;
}


/***********************************************************************************/
/* Function: free_albedo                                                  @62_30i@ */
/* Description:                                                                    */
/*  Free memory of struct albedo_struct.                                           */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static void free_albedo (albedo_struct *albedo) 
{
  int ia=0;

  fprintf (stderr, " ... freeing albedo structure\n");

  for (ia=0; ia<albedo->Nx; ia++) {
    fprintf (stderr, " ... ia=%d\n", ia);
    fflush (stderr);
    free(albedo->albedo2D[ia]);
  }
    
  free(albedo->albedo2D);

  free(albedo);
}


/***********************************************************************************/
/* Function: free_rpv                                                     @62_30i@ */
/* Description:                                                                    */
/*  Free memory of struct albedo_struct.                                           */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static void free_rpv (albedo_struct *albedo) 
{
  int ia=0;

  for (ia=0; ia<albedo->Nx; ia++)
    free (albedo->rpv_index[ia]);
    
  free (albedo->rpv_index);
  free (albedo->rpv_count);
  free (albedo->rpv);
  free (albedo);
}


/***********************************************************************************/
/* Function: free_rossli                                                  @62_30i@ */
/* Description:                                                                    */
/*  Free memory of struct albedo_struct.                                           */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static void free_rossli (albedo_struct *albedo) 
{
  int ix=0;

  for (ix=0; ix<albedo->Nx; ix++)
    free (albedo->rossli2D[ix]);
  free (albedo->rossli2D);
  free (albedo);
}


/***********************************************************************************/
/* Function: free_atmosphere                                              @62_30i@ */
/* Description:                                                                    */
/*  Free memory of struct atmosphere_struct.                                       */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static void free_atmosphere (atmosphere_struct *atmos)
{
  int ids=0, iris=0, ivs=0, iv=0, ic=0, isp=0;

  free (atmos->X);
  free (atmos->Y);
  free (atmos->Z);

  for (isp=0; isp<=atmos->n_caoth; isp++)
    free(atmos->threed[isp]);
  free (atmos->threed);
  
  if (atmos->ddis_eps!=NULL)
    free (atmos->ddis_eps);

  free_profile (atmos->kabs);

  for (ids=0; ids<atmos->nscaDS; ids++)
    for (iris=0; iris<atmos->nscaRIS; iris++)
      free_profile (atmos->ksca[ids][iris]);
  for (ids=0; ids<atmos->nscaDS; ids++)
    free (atmos->ksca[ids]);
  free (atmos->ksca);

  for (ids=0; ids<atmos->nscaDS; ids++)
    for (iris=0; iris<atmos->nscaRIS; iris++)
      for (ivs=0; ivs<atmos->nscaVIS; ivs++)
	free_profile (atmos->kext[ids][iris][ivs]);
  for (ids=0; ids<atmos->nscaDS; ids++)
    for (iris=0; iris<atmos->nscaRIS; iris++)
      free (atmos->kext[ids][iris]);
  for (ids=0; ids<atmos->nscaDS; ids++)
    free (atmos->kext[ids]);
  free (atmos->kext);

  free_profile (atmos->g1);
  free_profile (atmos->g2);
  free_profile (atmos->ff);

  free_profile (atmos->reff);

  free_profile3D (atmos->Bplanck);

  for (ids=0; ids<atmos->nscaDS; ids++)
    for (iris=0; iris<atmos->nscaRIS; iris++)
      free_profile3D (atmos->ksca3D[ids][iris]);
  for (ids=0; ids<atmos->nscaDS; ids++)
    free (atmos->ksca3D[ids]);
  free (atmos->ksca3D);

  for (ids=0; ids<atmos->nscaDS; ids++)
    for (iris=0; iris<atmos->nscaRIS; iris++)
      for (ivs=0; ivs<atmos->nscaVIS; ivs++)
	free_profile3D (atmos->kext3D[ids][iris][ivs]);
  for (ids=0; ids<atmos->nscaDS; ids++)
    for (iris=0; iris<atmos->nscaRIS; iris++)
      free (atmos->kext3D[ids][iris]);
  for (ids=0; ids<atmos->nscaDS; ids++)
    free (atmos->kext3D[ids]);
  free (atmos->kext3D);

  free_profile3D (atmos->kabs3D);

  free_profile3D (atmos->g1_3D);
  free_profile3D (atmos->g2_3D);
  free_profile3D (atmos->ff_3D);

  free_profile3D (atmos->reff_3D);

#ifdef CLDPRP

  if (atmos->sample_cldprp){
    free_profile3D (atmos->dxlwc_3D);
    free_profile3D (atmos->dylwc_3D);
    free_profile3D (atmos->dzlwc_3D);
  }
  
#endif

  /* FIXCE better to use sample->spectral_is as criterion CHECK!!! */ 
  if(atmos->nlambda_abs>1){
    free (atmos->lambda);
    for (isp=0; isp<=atmos->n_caoth; isp++) {
      for (iv=0; iv<atmos->nlambda_abs; iv++){
        free( atmos->kabs_spectral[isp][iv] ); 
        free( atmos->ksca_spectral[isp][iv] ); 
      }
      free( atmos->kabs_spectral[isp] ); 
      free( atmos->ksca_spectral[isp] ); 
    }
    free( atmos->kabs_spectral ); 
    free( atmos->ksca_spectral ); 
    
    for (iv=0; iv<atmos->nlambda_abs; iv++) 
      free((atmos->Bplanck_spectral)[iv]);
    free(atmos->Bplanck_spectral);
  }
    
  if (atmos->Nc > 1){
    for (ic=0; ic<atmos->Nc; ic++){
      free( atmos->kabs_scaled[ic] ); 
      free( atmos->ksca_scaled[ic] );  
    }
    free( atmos->kabs_scaled ); 
    free( atmos->ksca_scaled );
  } 
  
}


/***********************************************************************************/
/* Function: free_result                                                  @62_30i@ */
/* Description:                                                                    */
/*  Free memory of struct result_struct.                                           */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static void free_result (result_struct *result, int thermal_heating_method,
                         atmosphere_struct *atmos, sample_struct *sample, int absorption) /* **CK 2013.09.23 added thermal_heating_method */
{
  int is=0, js=0, kc=0, ic=0, il=0, ip=0;
  
  /* free memory of result arrays */
  free_radiation_field(result->surf, sample->Nx, sample->Nr, sample->Nd, sample->Ny, sample->std);

#if HAVE_LIDAR
  if (sample->LidarLocEst)
    free_lidar_result_field(result->lidar, sample->Nli, sample->LLE_Nt, sample->std);
  if (sample->abs_jacobian)
    free_jacobian_result_field (result->jacobian);
#endif
  
  /* CE: This should be moved to subroutine CHECK!!! */ 
  if (sample->spectral_is || sample->concentration_is){
    
    for (ic=0; ic<atmos->Nc; ic++){
      for (is=0; is<sample->Nx; is++) {
	for (js=0; js<sample->Ny; js++){
	  for (ip=0; ip<sample->nstokes; ip++){
	    free(result->surf_t->rad_t[ic][is][js][ip]);
	  }
	  free(result->surf_t->rad_t[ic][is][js]);
	}
	free(result->surf_t->rad_t[ic][is]);
      }
      free(result->surf_t->rad_t[ic]);
    }
    free(result->surf_t->rad_t);
    
    for (kc=0; kc<=atmos->Nz; kc++){
      if (sample->sample[kc]!=0){
        for (ic=0; ic<atmos->Nc; ic++){
	  for (is=0; is<sample->Nx; is++) {
	    for (js=0; js<sample->Ny; js++){
	      for (ip=0; ip<sample->nstokes; ip++){
                free(result->rad_t[kc]->rad_t[ic][is][js][ip]);
              }
              free(result->rad_t[kc]->rad_t[ic][is][js]);
            }
            free(result->rad_t[kc]->rad_t[ic][is]);
          }
          free(result->rad_t[kc]->rad_t[ic]);
        }
        free(result->rad_t[kc]->rad_t);
      }
    }
  }
  
  
  for (kc=0; kc<=atmos->Nz; kc++)
    if (sample->sample[kc]!=0)
      free_radiation_field (result->alt[kc], sample->Nx, sample->Nr, sample->Nd, 
                            sample->Ny, sample->std);
  
  free(result->alt);

  if (sample->backward) {
    for (is=0; is< sample->Nx; is++) {
      for (js=0; js< sample->Ny; js++) {
        free (result->back[is][js]);
        free (result->back2[is][js]);
	if (sample->backward == MCBACKWARD_HEAT) /* **CK 2013.09.23 */
	  if (thermal_heating_method == MCBACKWARD_HEAT_DENET || thermal_heating_method == MCBACKWARD_HEAT_HYBRID) {
	    free (result->back_dEnet[is][js]);
	    free (result->back_dEnet2[is][js]);
	  }
      }
      free (result->back[is]);
      free (result->back2[is]);
      free (result->backemis[is]);
      if (sample->backward == MCBACKWARD_HEAT) /* **CK 2013.09.23 */
	if (thermal_heating_method == MCBACKWARD_HEAT_DENET || thermal_heating_method == MCBACKWARD_HEAT_HYBRID) {
	  free (result->back_dEnet[is]);
	  free (result->back_dEnet2[is]);
	}    
    }
    free (result->back);
    free (result->back2);
    free (result->backemis);
    if (sample->backward == MCBACKWARD_HEAT) /* **CK 2013.09.23 */
      if (thermal_heating_method == MCBACKWARD_HEAT_DENET || thermal_heating_method == MCBACKWARD_HEAT_HYBRID) {
	free (result->back_dEnet);
 	free (result->back_dEnet2);    
      }    

    if (sample->spectral_is || sample->concentration_is){
      for (ic=0; ic < atmos->Nc; ic++){
	for (is=0; is< sample->Nx; is++) {
	  for (js=0; js< sample->Ny; js++) {
            for (ip=0; ip < sample->nstokes; ip++){
              free (result->back_t[ic][is][js][ip]);
            }
            free (result->back_t[ic][is][js]);
          }
          free (result->back_t[ic][is]);
        }
        free (result->back_t[ic]);
      }
      free (result->back_t); 
    }  
  }
  
  if (absorption) {
    free_ddprofile3D (result->absorption3D);
    if (sample->std)
      free_ddprofile3D (result->absorption3D2);
  }

  if (sample->boxairmass)
    free(result->pathlength_per_layer_tot);
  
  if (sample->ncirc) {
    /* QUICK FIX*/
    for (il=0;il<10;il++)
      free(result->circcontr[il]);
    free(result->circcontr);
  }
 
  free (result->mish);
#if HAVE_LIDAR
  free (result->lidcb);
#endif

  free (result);
}


/***********************************************************************************/
/* Function: setup_mystic                                                 @62_30i@ */
/* Description:                                                                    */
/*  Allocate memory, create and initialize 1D and 3D arrays, etc.                  */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int setup_mystic (caoth3d_out_struct *caoth3d,
                         atmosphere_struct  *atmos,
			 result_struct      *result,
                         albedo_struct      *albedo,
                         surftemp_struct    *surftemp,
                         sample_struct      *sample,
			 float               delta_scaling_mucut,
                         int                 absorption, 
			 int                 thermal_heating_method,
                         int                 nlyr,
			 int                 n_caoth,
                         float             **dt_s,
			 float             **om_s,
			 float             **g1_s,
			 float             **g2_s,
			 float             **f_s,
                         float             **ds_s,
                         alis_struct        *alis, 
                         float               rayleigh_depol,
                         float             **re_s,
                         float            ***temper,
			 int                 temper3d, 
                         float              *zprof,
                         float              *sza,
			 float              *phi0,
			 float               sza_spher,
			 float               phi0_spher,
                         int                 aerosol_scatter_type,
			 float            ***momaer,
			 int                *nmomaer, 
			 int               **nthetaaer,
			 float            ***thetaaer,
			 double           ***muaer,
			 float            ***phaseaer,
                         int                 nphamataer,
                         float               alb,
			 float               *alb_type,
			 float               *rpv_rho0,
			 float               *rpv_k,
			 float               *rpv_theta, 
			 float               *rpv_scale,
			 float               *rpv_sigma,
			 float               *rpv_t1,
			 float               *rpv_t2, 
			 float               *hapke_h,
			 float               *hapke_b0,
			 float               *hapke_w,
			 float               *rossli_iso,
			 float               *rossli_vol,
			 float               *rossli_geo,
			 int                  rossli_hotspot,
                         float               *u10,
			 float               *pcl,
			 float               *xsal,
			 float               *uphi,
			 int                 *solar_wind,
                         float               *bpdf_u10,
                         int                  source,
			 float               *wvnmlo,
			 float               *wvnmhi,
			 float               *wavelength,
                         float               *zout,
			 int                  nzout,
                         int                  nxcld,
			 int                  nycld,
			 int                  nzcld,
			 double               dxcld,
			 double               dycld,
                         char                *umufilename,
                         char                *sunshape_filename,
                         char                *albedo_filename,
                         char                *alb_type_filename,
                         char                *rpvfilename,
			 char               **rpv_labels,
			 int                  rpv_nlabels,
                         char                *ambralsfilename,
                         char                *rosslifilename,
                         double               truncate,
                         int                  loaddata, 
                         char                *datapath,
			 int                  spectral,
			 double               r_earth,
                         int                  visualize,
			 int                  quiet,
			 long int             nphotons) /* **CK  2013.08.27 Added nphotons */
{
  int status=0;
  int doddis=0;
  int isp=0;
  int ispo=0;
  int isAmbralsFile=0;

  if ( sample->spherical3D ) {
    *sza  =  90.0 - sza_spher;
    *phi0 = -90.0 - phi0_spher;
  }

  /* define dimensions of scatter profiles concerning different scatter modes */
  atmos->nscaDS = MCSC_MODE_NORMAL+1;
  if (sample->delta_scaling!=-1)
    atmos->nscaDS = MCSC_MODE_DELTA_SCALE+1;

  /* define dimensions of scatter profiles concerning different RIS modes */
  atmos->nscaRIS = MCRIS_MODE_NORMAL+1;
#if HAVE_LIDAR
  if (sample->LLE_RIS_MAS || sample->RIS_MS)
    atmos->nscaRIS = MCRIS_MODE_MAS+1;
#ifdef NEWRISQIDD
    atmos->nscaRIS++;
#endif
#endif

  /* define dimensions of scatter profiles concerning different VIS modes */
  atmos->nscaVIS = MCVIS_MODE_NORMAL+1;
#if HAVE_LIDAR
  if (sample->LLE_VIS_QIDD)
    atmos->nscaVIS = MCVIS_MODE_QIDD+1;
#endif
  
  /* Rayleigh depolarisation */
  atmos->rayleigh_depol = rayleigh_depol;
  
  status = setup_profiles1D (n_caoth,
			     dt_s, om_s, g1_s, g2_s, f_s, ds_s,
                             re_s,
                             zprof, nlyr,
			     sample,
                             atmos, 
                             alis);
  
  if (status!=0)
    return err_out ("Error %d returned by setup_profiles1D()\n", status);

  if (!quiet)
    fprintf (stderr, " ... setting up 3D grid\n");

  status = setup_caoth3D (caoth3d,
			  nxcld,
			  nycld,
			  nzcld,
			  dxcld,
			  dycld,
                          *wavelength,
                          atmos,
			  datapath,
                          sample->delta_scaling,
			  sample->spherical3D,
			  sample->spherical3D_scene,
			  sample->spherical3D_scene_lon_min,
			  sample->spherical3D_scene_lon_max,
			  sample->spherical3D_scene_lat_min,
			  sample->spherical3D_scene_lat_max,
			  r_earth,
			  *sza,
			  sample->tipadir,
			  visualize,
			  sample->cldprp,
			  quiet);
  
  if (status!=0)
    return err_out ("Error %d returned by setup_caoth3D()\n", status);

  /* read sampling information */
  status = setup_sample2D (atmos,
                           zout, nzout,
                           sza, phi0,
                           wavelength,
                           sample, 
                           quiet);
  
  if (status!=0)
      return err_out ("Error %d returned by setup_sample2D()\n", status);
  
  if (loaddata) {

    if (strlen(umufilename) != 0 ) {
      /* setup umu data */
      status = setup_umu2D (umufilename, sample, atmos, quiet);
      
      if (status!=0)
	return err_out ("Error %d returned by read_2D_umu()\n", status);
    }

    /*  Check if either a sunshape file or a sun radius is given, Bernhard Reinhardt */
    if ( (strlen(sunshape_filename) != 0) || (sample->sun_radius > 0.0) ){	
      status = setup_sunshape ( sample, sunshape_filename, *wvnmlo, *wvnmhi);

      if (status!=0)
	return err_out ("Error %d returned by setup_sunshape()\n", status);
    }

    /* setup albedo data */

    if (*albedo_filename)
      strcpy (albedo->filename, albedo_filename);
    else 
      strcpy (albedo->filename, "./albedo2D.dat");
  
    if (*alb_type_filename)
      strcpy (albedo->spectralfilename, alb_type_filename);
    else 
      strcpy (albedo->spectralfilename, "./albedospectral2D.dat");
  
    if (*rpvfilename)
      strcpy (albedo->rpv_filename, rpvfilename);
    else
      strcpy (albedo->rpv_filename, "./rpv2D.dat");
  
    if (*rosslifilename)
      strcpy (albedo->rossli_filename, rosslifilename);
    else
      strcpy (albedo->rossli_filename, "./rossli2D.dat");
  
    if (*ambralsfilename) {
      isAmbralsFile=1;
      strcpy (albedo->rossli_filename, ambralsfilename);
    }else
      strcpy (albedo->rossli_filename, "./ambrals2D.dat");
  

    if (!quiet) {
      switch (albedo->method) {
      case MCALB_LAM:
	fprintf (stderr, " ... homogeneous Lambertian surface albedo\n");
	break;
      case MCALB_LAM2D:
	fprintf (stderr, " ... reading 2D albedo data from %s\n", albedo->filename);
	break;
      case MCALB_LAM2D_SPECTRAL:
	fprintf (stderr, " ... reading 2D albedo data from %s\n", albedo->spectralfilename);
	break;
      case MCALB_RPV:
	fprintf (stderr, " ... homogeneous RPV BRDF\n");
	break;
      case MCALB_RPV2D_SPECTRAL:
	fprintf (stderr, " ... reading RPV BRDF from %s\n", albedo->rpv_filename);
	break;
      case MCALB_COXANDMUNK:
	fprintf (stderr, " ... homogeneous Cox and Munk ocean BRDF\n");
	break;
      case MCALB_HAPKE:
	fprintf (stderr, " ... homogeneous Hapke BRDF\n");
	break;
      case MCALB_ROSSLI:
	fprintf (stderr, " ... homogeneous Ross-Li BRDF\n");
	break;
      case MCALB_ROSSLI2D:
	fprintf (stderr, " ... reading Ross-Li BRDF from %s\n", albedo->rossli_filename);
	if (isAmbralsFile){
	  fprintf(stderr, " ... using AMBRALS (Lucht2000) implementation of Ross-Li\n");
	}else{
	  fprintf(stderr, " ... using Lin(2015) implementation of Ross-Li\n");
	}
	break;
      case MCALB_TSANG:
	fprintf (stderr, " ... homogeneous Tsang ocean BPDF\n");
	break;
      default:
	fprintf (stderr, "Fatal error, unknown albedo type %d\n", albedo->method);
	break;
      }
    }

    status = setup_albedo2D (albedo->filename, albedo->spectralfilename, albedo->rpv_filename, 
                             albedo->rossli_filename, isAmbralsFile, atmos, albedo,
                             alb, alb_type, rpv_rho0, rpv_k, rpv_theta, rpv_scale,
			     rpv_sigma, rpv_t1, rpv_t2, rpv_labels,
                             rpv_nlabels,
			     hapke_h, hapke_b0, hapke_w,
			     rossli_iso, rossli_vol, rossli_geo, rossli_hotspot,
			     *u10, *pcl, *xsal, *uphi, 
                             *solar_wind, *bpdf_u10, spectral, *wavelength,
                             sample->polarisation,
			     sample->spherical3D, sample->spherical3D_scene,
			     sample->spherical3D_scene_lon_min,
			     sample->spherical3D_scene_lon_max,
			     sample->spherical3D_scene_lat_min,
			     sample->spherical3D_scene_lat_max,
			     alis->nlambda_abs, alis->albedo, alis->alb_type,
			     sample->ixmin,
			     sample->ixmax,
			     sample->iymin,
			     sample->iymax,
			     quiet);
   
    if (status)
      return err_out ("Error %d returned by setup_albedo2D()\n", status);
  } /* endif loaddata */

  /* formerly, the aerosol setup was performed only if loaddata */
  /* now, phase_aer etc. is deallocated at the end of MYSTIC, so */
  /* at the next call, the aerosol must always be setup! */
  /* As soon as MCSCAT_AER has been merged with CAOTH, all this will be */
  /* obsolete */
  if ( atmos->tscatot[MCCAOTH_AER] > 0.0 )
    atmos->scatter_type[MCCAOTH_AER]=aerosol_scatter_type;
  else
    atmos->scatter_type[MCCAOTH_AER]=MCSCAT_SKIP;

  if (atmos->scatter_type[MCCAOTH_AER]==MCSCAT_AER) {
    status = setup_aerosol1D ( &(atmos->phase_aer), nlyr,
			       dt_s[CAOTH_AER], om_s[CAOTH_AER], momaer, nmomaer,
			       nthetaaer, muaer, phaseaer, nphamataer,
			       truncate, quiet );
    if (status)
      return err_out ("Error %d returned by setup_aerosol1D()\n", status);
  }
#if HAVE_VROOM
  /* following should be standing somewhere else !!! */
  for (isp=MCCAOTH_AER;isp<=atmos->n_caoth;isp++) /* molecular should not be tested here */
    /* doddis is only zero if pure rayleigh and no scattering in any other caoth */
    if ( atmos->tscatot[isp] > 0.0 )
      doddis++;
    else {
      ispo=isp-MCCAOTH_FIR;
      if (ispo>=0)
	if (caoth3d[ispo].nthreed>=1)
	  doddis++;
    }

  /* prepare vroom, in particular define  representative phase function */
  if (sample->vroom!=0 || sample->LidarLocEst) {
    status = mc_vroom_prepare (sample, atmos, delta_scaling_mucut, nlyr, atmos->tscatot[MCCAOTH_AER]>0.0, quiet);

    if (status!=0)
      return err_out ("Error %d returned by mc_vroom_prepare()\n", status);
  }

  if (!sample->LidarLocEst) {
    status = set_vroom_settings (sample->vroom, sample, quiet);
    if (status!=0)
      return err_out ("Error %d returned by set_vroom_settings()\n", status);
  }

  /* if vroom then enforce reflectalways */
  if (sample->vroomreflectalways)
    albedo->reflectalways = 1;

  status = mc_vroom_check_and_verbose (sample, quiet, doddis );
  if (status!=0)
    return err_out ("Error, mc_vroom_check_and_verbose returned status %d\n", status);
#endif

  status = setup_mc_result (atmos, sample, result, absorption, thermal_heating_method, quiet);

  if (status!=0)
    return err_out ("Error %d returned by setup_mc_result()\n", status);

  /* calculate energy emitted from the surface */
  status = calc_surface_emission ( albedo, source, surftemp->btemp, wvnmlo, wvnmhi,
				   &(albedo->bplkavg), &(albedo->Wsurf) );
  if (status)
    return err_out ("Error %d returned by calc_surface_emission()\n", status);



  status = setup_thermal (temper, temper3d, atmos->Nx, atmos->Ny, atmos->Nz, *wvnmlo, *wvnmhi, 
                          atmos->kabs, atmos->kabs3D, atmos->Z, source, absorption, thermal_heating_method,
                          &(atmos->Bplanck), &(atmos->maxemis), &(atmos->Watm), surftemp,
			  sample->backward, sample->zstart, &(sample->backemis_kc), 
                          alis->nlambda_abs, alis->lambda, &(atmos->Bplanck_spectral), quiet, nphotons, sample); /* **CK  2013.08.27 Added nphotons and sample; added atmos->Nx, atmos->Ny */

  if (status!=0)
    return err_out ("Error %d returned by setup_thermal()\n", status);

#if HAVE_LIDAR
  check_lidar_in_atmos (sample, atmos);

  if (status!=0)
    return err_out ("Error %d returned by check_lidar_in_atmos()\n", status);
#endif

  return 0;
}


/***********************************************************************************/
/* Function: setup_sunshape                                               @62_30i@ */
/* Description:                                                                    */
/*  Setup sunshape.                                                                */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int setup_sunshape ( sample_struct *sample,
			    char          *filename,
			    float          wvnmlo,
			    float          wvnmhi )
{
  double sunshape_norm_factor = 0.0;
  double **mu_new=NULL;
  float *temp=NULL;
  int i=0, status=0;

  /*  Use default spectral sunshape if sunshape filename is "default" or not specified */   
  if (strcmp("default", filename) == 0 || strcmp("", filename) == 0 ){
    sample->sample_backward_sunshape_n_lines = 200; 
    sample->sample_backward_sunshape_p =
      (float*) calloc(sample->sample_backward_sunshape_n_lines, sizeof(float));
    sample->sample_backward_sunshape_alpha =
      (double*) calloc(sample->sample_backward_sunshape_n_lines, sizeof(double));
    sample->sample_backward_sunshape_F =
      (double*) calloc(sample->sample_backward_sunshape_n_lines, sizeof(double));
    /* generate spectral sunshape according to Koepke2001 */
    gen_default_sunshape(wvnmlo, wvnmhi, sample->sample_backward_sunshape_n_lines,
			 sample->sample_backward_sunshape_p,
			 sample->sample_backward_sunshape_alpha);
  }
  else  {
    /*  Read sun shape file, Bernhard Reinhardt */

    /*  Do a first read of the sun shape file. This is to allocate sunshape_F */
    /*  and to read alpha as doubles. */
    status = read_2c_file( filename,
			   &sample->sample_backward_sunshape_F,
			   &sample->sample_backward_sunshape_alpha,
			   &sample->sample_backward_sunshape_n_lines );
    if (status!=0)
      return err_out ("Error %d returned by read_2c_file() when reading sun shape input file\n",
		      status);
        
    /* Now do a second read to get p as float. we read the second */
    /* column as float into temp to throw it away                 */

    status = read_2c_file_float( filename,
				 &sample->sample_backward_sunshape_p,
				 &temp,
				 &sample->sample_backward_sunshape_n_lines );
    if (status!=0)
      return err_out ("Error %d returned by read_2c_file_float() when reading sun shape input file\n", status);
    free(temp);
  }

  /* Now we have either read in or generated a sun shape. Let??s do */
  /* some postprocessing                                           */

  sample->sample_backward_sunshape = 1;
  sample->sample_backward_sunshape_p2 =
    (float*) calloc(sample->sample_backward_sunshape_n_lines, sizeof(float));
  status = ASCII_calloc_double(&mu_new,1,sample->sample_backward_sunshape_n_lines);

  /* We need the temp variable mu_new in matrix notation to call   */
  /* normalize_phase. Weight the PDF (sunshape_p2) with sine to    */
  /* make sure that we do not send to many photons to the center   */
  /* of the sun. sunshape_F is the cumulative distribution	       */
  /* function of sunshape_p2. sunshape_F is used in the routines   */
  /* for the diffuse radiation, while direct radiation is computed */
  /* using sunshape_p (not weighted with sine)                     */

  for (i=0; i < sample->sample_backward_sunshape_n_lines; i++){
    mu_new[0][i] = sample->sample_backward_sunshape_alpha[i] ;
    sample->sample_backward_sunshape_p2[i]
      = sample->sample_backward_sunshape_p[i];
    sample->sample_backward_sunshape_p2[i]
      *= sind ( sample->sample_backward_sunshape_alpha[i] * sample->sun_radius );
  }

  normalize_phase( mu_new,
		   &sample->sample_backward_sunshape_p2,
		   sample->sample_backward_sunshape_F,
		   &sample->sample_backward_sunshape_n_lines, 1, 0 );

  /*normalize to 1 instead of 2 */
  for (i=0; i < sample->sample_backward_sunshape_n_lines; i++){
    sample->sample_backward_sunshape_F[i] /= 2.0;
  }
      
  /* sunshape_norm_factor is the normalizing factor for the       */
  /* sunshape pdf. It is calculated as			      */
  /* integral(sin(theta))_0..sunrad /			      */
  /* integral(p(theta)*sin(theta))_0..sunrad. This works probably */
  /* only for geometric distribution of the photons               */

  /* First calc the denominator of the term: */
  for (i=1; i < sample->sample_backward_sunshape_n_lines; i++){
    sunshape_norm_factor
      += ( cosd( sample->sample_backward_sunshape_alpha[i-1] * sample->sun_radius ) -
	   cosd( sample->sample_backward_sunshape_alpha[i]   * sample->sun_radius ) 
	   ) *
      ( sample->sample_backward_sunshape_p[i] +
	sample->sample_backward_sunshape_p[i-1]
	) / 2.;
  }

  /* Now the numerator. Here the integral is easy to solve: */
  sunshape_norm_factor = ( 1. - cosd(sample->sun_radius) ) / sunshape_norm_factor;
  for (i=0; i < sample->sample_backward_sunshape_n_lines; i++){
    sample->sample_backward_sunshape_p[i] *= sunshape_norm_factor;
  }

  ASCII_free_double(mu_new,1);

  /* calculate sunshape slope */
  sample->sample_backward_sunshape_slope
    = (double*) calloc(sample->sample_backward_sunshape_n_lines-1, sizeof(double));
  for (i=0; i<sample->sample_backward_sunshape_n_lines-1; i++)
    sample->sample_backward_sunshape_slope [i] =
      ( sample->sample_backward_sunshape_p [i+1]
	- sample->sample_backward_sunshape_p [i] ) /
      ( sample->sample_backward_sunshape_alpha [i+1]
	- sample->sample_backward_sunshape_alpha [i] );

  return 0;
}


/***********************************************************************************/
/* Function: calc_surface_emission                                        @62_30i@ */
/* Description:                                                                    */
/*  Calculate albedo->bplkavg and albedo->Wsurf                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static inline int calc_surface_emission ( albedo_struct *albedo,
					  int           source,
					  float          btemp,
					  float         *wvnmlo,
					  float         *wvnmhi,
					  float         *bplkavg,
					  double        *Wsurf)
{
  double avg_albedo=0.0;
  int ia=0, ja=0;

  /* calculate Planck emission for surface temperature */
  switch (source) {
  case MCSRC_THERMAL_SURFACE:
    switch (albedo->method) {
    case MCALB_LAM2D:
    case MCALB_LAM2D_SPECTRAL:
      /* ??? need to check if this averaging procedure gives correct results ??? CHECK!!! */
      if (albedo->method == MCALB_LAM2D) {
	avg_albedo = 0;
	for (ia=0; ia<albedo->Nx; ia++) 
	  for (ja=0; ja<albedo->Ny; ja++) {
	    if (albedo->method == MCALB_LAM2D)
	      avg_albedo += albedo->albedo2D[ia][ja];
	    else
	      avg_albedo += albedo->alb_type[(int) albedo->rpv_index[ia][ja]];
	  }
	avg_albedo /= (double) (albedo->Nx * albedo->Ny);
      }
      break;

    case MCALB_LAM:
      avg_albedo = albedo->albedo;
      break;

    case MCALB_RPV:
    case MCALB_RPV2D_SPECTRAL:
    case MCALB_COXANDMUNK: 
    case MCALB_HAPKE:
    case MCALB_ROSSLI:
    case MCALB_ROSSLI2D:
      fprintf (stderr, "**** hmm, bi-directional reflectance specified together\n");
      fprintf (stderr, "**** with thermal source; ignoring BRDF and using albedo %g\n", 
	       albedo->albedo);
      fprintf (stderr, "**** instead!\n"); 
      return -1;

    default:
      fprintf (stderr, "Error, albedo->method %d not yet implemented!\n", albedo->method);
      return -1;
    }

    /* no break */

  case MCSRC_THERMAL_ATMOSPHERE:
  case MCSRC_THERMAL_BACKWARD:
    if (btemp>0)
      /*  F77_FUNC (cplkavg, CPLKAVG) (wvnmlo, wvnmhi, &btemp, bplkavg); */
      *bplkavg = c_planck_func1(*wvnmlo, *wvnmhi, btemp);
    else 
      *bplkavg=0;

    if ( source == MCSRC_THERMAL_SURFACE )
      *Wsurf = PI * (1.0 - avg_albedo) * *bplkavg;
    break;

  case MCSRC_NONE:
  case MCSRC_SOLAR:
  case MCSRC_LIDAR:
  case MCSRC_BLITZ:
    break;

  default: 
    fprintf (stderr, "Error, unknown source %d\n", source);
  }

  return 0;
}


/***********************************************************************************/
/* Function: setup_aerosol1D                                              @62_30i@ */
/* Description:                                                                    */
/*  Setup the structures for 1D aerosols.                                          */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: RPB                                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static inline int setup_aerosol1D ( pft     **phase_aer, 
				    int       nlyr,
				    float    *dtaer,
				    float    *omaer,
				    float  ***momaer,
				    int      *nmomaer, 
				    int     **nthetaaer,
				    double ***muaer,
				    float  ***phaseaer,
				    int       nphamataer,
				    double    truncate,
				    int       quiet )
{
  int sumphases=0, kc=0, i=0, ip=0, status=0;

  double **moment;

  /* aerosol phase function table CHECK!!! */ 
  /* ??? this is not a good solution; the phase function table is */
  /* ??? calculated separately for each layer, although it is     */
  /* ??? very likely that several layers use the same table;      */
  /* ??? how to improve:                                          */
  /* ???  read each phase function file only once in aerosol.c,   */
  /* ???  read_optprop_files(); create table that stores          */
  /* ???  the different moment arrays and an index table that     */
  /* ???  just stores which array is to be used for each layer    */
  /* ???  (similar to the handling of 2D BRDFs                    */

  (*phase_aer) = calloc (nlyr, sizeof(pft)); 

  for (kc=0; kc<nlyr; kc++)
    sumphases += nthetaaer[kc][0];

  if (!sumphases) {
    for (kc=0; kc<nlyr; kc++) {
      if (dtaer[kc] > 0.0){
	moment = calloc (nphamataer, sizeof(double *));
          
	for (ip=0; ip<nphamataer;ip++){
	  moment[ip] = calloc ((size_t) nmomaer[kc], sizeof(double));
            
	  for (i=0; i<nmomaer[kc]; i++)
	    moment[ip][i] = momaer[kc][ip][i];
	}
	if (!quiet)
	  fprintf (stderr, " ... creating aerosol scattering table for layer %d, using %d moments, phase matrix elements %d\n", 
		   kc, nmomaer[kc], nphamataer);
          
	status = setup_Legendre_table_aerosol
	  (moment, nmomaer[kc], nphamataer,  truncate,
	   &((*phase_aer)[kc]), quiet);
          
	if (status!=0)
	  return err_out ("Error %d returned by setup_Legendre_table_aerosol()\n", status); 
	for (ip=0; ip<nphamataer;ip++)
	  free(moment[ip]);
	free(moment);
      }
    }
  }
  else {
    for (kc=0; kc<nlyr; kc++)
      if (dtaer[kc] > 0.0){
	/* calculate cumulative table */
	status = calc_cumulative_table (muaer[kc], phaseaer[kc], nthetaaer[kc], nphamataer,
					truncate, &((*phase_aer)[kc]), -999.99, quiet);
  
	if (status!=0)
	  return err_out ("Error %d returned by calc_cumulative_table()\n", status);
      }
  }

  return 0;
}


/***********************************************************************************/
/* Function: free_mystic_for_load                                         @62_30i@ */
/* Description:                                                                    */
/*  Free the MYSTIC memory structures.                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static void free_mystic_for_load (albedo_struct *albedo, 
				  sample_struct *sample)
{

  switch (albedo->method) {
  case MCALB_LAM:
  case MCALB_RPV:
  case MCALB_HAPKE:
  case MCALB_ROSSLI:
  case MCALB_TSANG:
  case MCALB_COXANDMUNK: 
    /* don't need to do anything because 1D albedo */
    break;
    
  case MCALB_LAM2D:
    free_albedo(albedo);
    break;
      
  case MCALB_RPV2D_SPECTRAL:
  case MCALB_LAM2D_SPECTRAL:
    free_rpv(albedo);
    break;
      
  case MCALB_ROSSLI2D:
    free_rossli(albedo);
    break;
      
  default:
    fprintf (stderr, "Fatal error, unknown albedo type %d\n", albedo->method);
    break;
  }
  
  /* ??????????  CHECK!!!
  if (nonHGaer)
    for (kc=0; kc<nlyr; kc++)
      free_pft (phase_aer[kc]);
      ?????????? */

}


/***********************************************************************************/
/* Function: sc_Isotropic_phi                                             @62_30i@ */
/* Description:                                                                    */
/*   Calculate a random azimuth angle for an azimutally isotropic phase function.  */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

inline double sc_Isotropic_phi()
{
  return 360.0 * uvspec_random();
}


/***********************************************************************************/
/* Function: sc_Rayleigh_mu                                               @62_30i@ */
/* Description:                                                                    */
/*   Calculate a random scattering angle for the Rayleigh phase function.          */
/*   Effect of depolarization is not included, see sc_Rayleigh_mu_depol()          */
/*   Returns mu=cos(theta), for more speed!                                        */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: ...., last big modification by Robert Buras                             */
/*                                                                        @i62_30@ */
/***********************************************************************************/

double sc_Rayleigh_mu()
{
  double P = uvspec_random();
  double q = 8.0 * P - 4.0;
  /* ??? cbrt() is a GNU extension ??? CHECK!!! */ 
  double u = cbrt (-q/2.0 + sqrt(1.0 + q*q/4.0));
  double v = -1.0 / u;
  
  double mu = u + v;
  
  /* ATTENTION: No error check, but very rude correction */

  if (mu > 1.0) {
    if (mu > 1.0 + MC_EPSILON) {
      fprintf(stderr,"Error! mu = %e in sc_Rayleigh_mu() is unphysical!\n",mu);
      return NOT_A_NUMBER;
    }
    mu=1.0;
  }
  
  if (mu < -1.0) {
    if (mu < -1.0 - MC_EPSILON) {
      fprintf(stderr,"Error! mu = %e in sc_Rayleigh_mu() is unphysical!\n",mu);
      return NOT_A_NUMBER;
    }
    mu=-1.0;
  }

  return mu;
}


/***********************************************************************************/
/* Function: sc_Rayleigh_mu_depol                                         @62_30i@ */
/* Description:                                                                    */
/*   Calculate a random scattering angle for the Rayleigh phase function           */
/*   including depolarisation.                                                     */
/*   Returns mu=cos(theta), for more speed!                                        */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: ...., last big modification by Robert Buras                             */
/*                                                                        @i62_30@ */
/***********************************************************************************/

double sc_Rayleigh_mu_depol (double depol)
{
  double gamma=depol/(2.0-depol);
  double P = uvspec_random();

  double p = (3.0*gamma+1.0)/(1.0-gamma);
  double q = (1.0+2.0*gamma)/(1.0-gamma)*(2.0 - 4.0 * P);

  double D = sqrt(q*q + p*p*p);

  /* ??? cbrt() is a GNU extension ??? CHECK!!! */ 
  double u = cbrt(-q+D);
  double v = -p / u;
  double mu = u + v;

  /* ATTENTION: No error check, but very rude correction */

  if (mu > 1.0) {
    if (mu > 1.0 + MC_EPSILON) {
      fprintf(stderr,"Error! mu = %e in sc_Rayleigh_mu_depol() is unphysical!\n",mu);
      return NOT_A_NUMBER;
    }
    mu=1.0;
  }
  
  if (mu < -1.0) {
    if (mu < -1.0 - MC_EPSILON) {
      fprintf(stderr,"Error! mu = %e in sc_Rayleigh_mu_depol() is unphysical!\n",mu);
      return NOT_A_NUMBER;
    }
    mu=-1.0;
  }

  return mu;
}


/***********************************************************************************/
/* Function: sc_Isotropic_mu                                              @62_30i@ */
/* Description:                                                                    */
/*   Calculate a random scattering angle for an isotropic phase function.          */
/*   Returns mu=cos(theta), for more speed!                                        */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: ...., last big modification by Robert Buras                             */
/*                                                                        @i62_30@ */
/***********************************************************************************/

inline double sc_Isotropic_mu()
{
  double P = uvspec_random();

  return (1.0 - 2.0*P);
}


/***********************************************************************************/
/* Function: sc_Isotropic_upward_mu                                       @62_30i@ */
/* Description:                                                                    */
/*   Calculate an upward random scattering angle for an isotropic phase function.  */
/*   Returns mu=cos(theta), for more speed!                                        */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: ...., last big modification by Robert Buras                             */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static inline double sc_Isotropic_upward_mu()
{
  return uvspec_random();
}

/***********************************************************************************/
/* Function: sc_Isotropic_downward_mu                                     @62_30i@ */
/* Description:                                                                    */
/*   Calculate a downward random scattering angle for an isotropic phase function. */
/*   Returns mu=cos(theta), for more speed!                                        */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: ...., last big modification by Robert Buras                             */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static inline double sc_Isotropic_downward_mu()
{
  return -uvspec_random();
}



/***********************************************************************************/
/* Function: sc_Lambertian_mu                                             @62_30i@ */
/* Description:                                                                    */
/*   Calculate a random polar angle for a Lambertian source.                       */
/*   Returns mu=cos(theta), for more speed!                                        */
/*   Is not exactly equivalent to old version, should be sqrt(1-P) for exact equiv */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: ...., last big modification by Robert Buras                             */
/*                                                                        @i62_30@ */
/***********************************************************************************/

inline double sc_Lambertian_mu()
{
  double P = uvspec_random();
  return sqrt(P); /* would be equivalent to old version if sqrt(1-P) */
}


/***********************************************************************************/
/* Function: sc_HG_mu                                                     @62_30i@ */
/* Description:                                                                    */
/*  Calculate a random polar angle for the Henyey-Greenstein phase function.       */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: ...., last big modification by Robert Buras                             */
/*                                                                        @i62_30@ */
/***********************************************************************************/

double sc_HG_mu(double g)
{
  double P = uvspec_random();

  double mu=0.0;

  /* ATTENTION: No error check, but very rude correction */

  if (g<-1.0) {
    fprintf (stderr, " ... g = %g out of range, setting to -1\n", g);
    g=-1;
  }

  if (g>1.0) {
    fprintf (stderr, " ... g = %g out of range, setting to 1\n", g);
    g=1;
  }
    

  if (g==0)
    mu = 1.0 - 2.0*P;
  else 
    mu = 1.0/2.0/g*(1.0+g*g-(g*g-1.0)*(g*g-1.0) / 
                         ((2.0*g*P-g-1.0)*(2.0*g*P-g-1.0)));
  
  /* ATTENTION: No error check, but very rude correction */

  if (mu > 1.0) {
    if (mu > 1.0 + MC_EPSILON) {
      fprintf(stderr,"Error! mu = %e in sc_HG_mu() is unphysical!\n",mu);
      return NOT_A_NUMBER;
    }
    mu=1.0;
  }
  
  if (mu < -1.0) {
    if (mu < -1.0 - MC_EPSILON) {
      fprintf(stderr,"Error! mu = %e in sc_HG_mu() is unphysical!\n",mu);
      return NOT_A_NUMBER;
    }
    mu=-1.0;
  }

  return mu;
}


/***********************************************************************************/
/* Function: HG                                                           @62_30i@ */
/* Description:                                                                    */
/*  Evaluate the Henyey-Greenstein phase function at cosine of polar angle         */
/*  mu = cos(theta) for asymmetry parameter g. The phase function                  */
/*  p(mu) = p(cos(theta)) is normalized to 2.                                      */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static double HG (double g, double mu) 
{
  double temp = 1.0 + g*g - 2.0*g*mu;

  return (1.0-g*g) / (temp * sqrt(temp));
}


/***********************************************************************************/
/* Function: HG2                                                          @62_30i@ */
/* Description:                                                                    */
/*  Evaluate the double-Henyey-Greenstein phase function at cosine of polar angle  */
/*  mu = cos(theta) for asymmetry parameter g. The phase function                  */
/*  p(mu) = p(cos(theta)) is normalized to 2.                                      */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static double HG2 (double g1, double g2, double ff, double mu) 
{
  return ff * HG(g1,mu) + (1.0-ff) * HG(g2,mu);
}


/***********************************************************************************/
/* Function: setup_Legendre_table_aerosol                                 @62_30i@ */
/* Description:                                                                    */
/*  Calculate lookup-tables of the cumulative probability distribution and the     */
/*  phase function from the moments of the phase function                          */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int setup_Legendre_table_aerosol (double **moment, int nmom, 
                                         int nphamat, 
                                         double truncate, 
                                         pft *phase, int quiet)
{
  int status=0, ip=0;
  double **mu;
  float **p;
  int *n;
  
  /* calculate phase function p(mu) */
  status = calc_Legendre_phase (moment, nmom, nphamat, &mu, &p, &n);
  
  if (status!=0)
    return err_out ("Error %d returned by calc_Legendre_phase()\n", status);

  /* calculate cumulative table */
  status = calc_cumulative_table (mu, p, n, nphamat, truncate, phase, -999.99, quiet); /*TZ ds: no delta-scaling for aerosol*/
  
  if (status!=0)
    return err_out ("Error %d returned by calc_cumulative_table()\n", status);
    
  for(ip=0; ip<nphamat; ip++){
    free(mu[ip]);
    free(p[ip]);
  }
      
  free(mu); free(p); free(n);

  return 0;
}


/***********************************************************************************/
/* Function: get_F                                                        @62_30i@ */
/* Description:                                                                    */
/*  Extract the integrated phase function for a given polar angle mu from a        */
/*  phase function table.                                                          */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

double get_F (pft *iphase, double mu, int scaled)
{
  int it=0, ip=0;
  double F=0; 

  /* if inside the truncation range */
  if (iphase->truncate && mu > iphase->mutrnc) {
    fprintf(stderr, "ERROR! truncate not implemented for get_F!\n");
    return NOT_A_NUMBER;
  }

  /* make sure that mu is physical, this should be checked directly at mu calculation site! */
  /*  if (mu > 1.0) mu=1.0;
      if (mu < -1.0) mu=-1.0;
  */
  
  ip = 0;
  if (iphase->is_linear) {
    /* this is a new version of getting it, it is much quicker, but less precise, CHECK */
    it = (int) ((acos(mu) - iphase->theta_0)*iphase->dthetainv);
    /* careful with rounding errors at the array boundaries */
    if (it>iphase->n[ip]-1 || it<0) {
      fprintf(stderr, "ERROR! mu is outside delta-scaled phasetable! Probably only infinitesimally, fix this in get_F\n");
      fprintf(stderr, "it %d nit %d mu %e rit %e \n",it,iphase->n[ip]-1,mu,((acos(mu) - iphase->theta_0)*iphase->dthetainv));
      return NOT_A_NUMBER;
    }
  }
  else
    it = locate (iphase->mu[ip], iphase->n[ip], mu);

  if (it==iphase->n[ip]-1)
    F = iphase->F[scaled][iphase->n[ip]-1];
  else
    /* bug fix, old version was unstable for spiky phase functions */
    F = iphase->F[scaled][it] + iphase->p[scaled][0][it] * ( mu - iphase->mu[scaled][it] )
      + iphase->A[scaled][it] * ( mu - iphase->mu[scaled][it] ) * ( mu - iphase->mu[scaled][it] );
    /*    F = iphase->C[scaled][it] + mu * iphase->B[scaled][it] + mu * mu * 
	  iphase->A[scaled][it]; */

  if (F<0.0) F=0.0;

  return F;
}


/***********************************************************************************/
/* Function: get_Lambertian_phase                                         @62_30i@ */
/* Description:                                                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: ...., last big modification by Robert Buras                             */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static inline double get_Lambertian_phase(double mu)
{
  return 4.*mu;
}


/***********************************************************************************/
/* Function: sc_mu                                                        @62_30i@ */
/* Description:                                                                    */
/*  Calculate a random scattering angle mu = cos(theta) using a pre-calculated     */
/*  table, see description of calc_cumulative_table(). If newrand is TRUE, a new   */
/*  random number is selected, otherwise the old one is recycled. CAREFUL: THIS    */
/*  FUNCTION IS NOT THREAD-SAFE AS IT STORES THE RANDOM NUMBER LOCALLY!            */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

double sc_mu (pft *iphase, int newrand, int scaled, double F_max, double F_min)
{
  int it=0;
  double mu=0.0, mu1=0.0, mu2=0.0, d=0.0, norm=2.0;
  static double U=0;

  /* limit F to within F_min and F_max */
  if (F_max != 0.0)
    norm = F_max - F_min;

  if (newrand!=0)
    U = norm * uvspec_random() + F_min;

  it = locate (iphase->F[scaled], iphase->n[0], U);

  /* special case: forward peak */
  if (it==iphase->n[0]-1)
    return 1.0;

  /* special case: phase function constant within interval */
  if (iphase->A[scaled][it]==0)
    return (U-iphase->C[scaled][it]) / iphase->B[scaled][it];
  
  /* CE 20100705: This calculation for d is unstable for plates CHECK!!! */ 
  /*  d = sqrt(iphase->B[scaled][it]*iphase->B[scaled][it] */
  /* 	   -4.0*iphase->A[scaled][it]*(iphase->C[scaled][it]-U)); */
  /* CE: The following seems to be more stable */
  d = sqrt( iphase->p[scaled][0][it]*iphase->p[scaled][0][it] - 
            4.*iphase->A[scaled][it] * (iphase->F[scaled][it]-U) );
  
  if (iphase->B[scaled][it]<0.0) {
    mu1 = (2.0 * (iphase->C[scaled][it] - U)) / (-iphase->B[scaled][it] + d);
    mu2 = (-iphase->B[scaled][it] + d) / (2.0 * iphase->A[scaled][it]);
  }
  else {
    mu1 = (-iphase->B[scaled][it] - d) / (2.0 * iphase->A[scaled][it]);
    mu2 = (2.0 * (iphase->C[scaled][it] - U)) / (-iphase->B[scaled][it] - d);
  }    

  if (iphase->p[scaled][0][it] < iphase->p[scaled][0][it+1])
    mu = (mu1>mu2?mu1:mu2);
  else
    mu = (mu1<mu2?mu1:mu2);

  if (mu > 1.0) {
    if (mu > 1.0 + MC_EPSILON) {
      fprintf(stderr,"Error! mu= %e in sc_mu() is unphysical!\n",mu);
      return NOT_A_NUMBER;
    }
    /* fprintf (stderr, " ... sc_mu(): cutting mu=%f to 1\n", mu); */
    mu = 1.0;
  }

  if (mu < -1.0) {
    if (mu < -1.0 - MC_EPSILON) {
      fprintf(stderr,"Error! mu= %e in sc_mu() is unphysical!\n",mu);
      return NOT_A_NUMBER;
    }
    /* fprintf (stderr, " ... sc_mu(): cutting mu=%f to -1\n", mu); */
    mu = -1.0;
  }

  return mu;
}


/***********************************************************************************/
/* Function: sc_interp_mu                                                 @62_30i@ */
/* Description: Calculate a random scattering angle mu = cos(theta) using a        */
/* pre-calculated table, see description of calc_cumulative_table().               */
/* Interpolation between to tables is needed here since reff lies between the      */
/* reff's of two table entries.                                                    */
/*                                                                                 */
/* CAUTION!!!!                                                                     */
/* Note that the "interpolation" in reff has been done on purpose!                 */
/* The former version was inconsistent with the way the phase function is          */
/* calculated for a given mu. However, methods like VROOM rely on the fact that    */
/* randomly chosing alpha from the phase function (sc_interp_mu) and deriving the  */
/* phase function from mu (get_phase_matrix_pft_interpol_reff) are done with       */
/* exactly the same phase function.                                                */
/*                                                                                 */
/* So: If you change anything in this subroutine, you must do so consistently also */
/*     in get_phase_matrix_pft_interpol_reff!!!                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static double sc_interp_mu (double reff, phase_function_table *phase, int scaled)
{
  int i1=0,i2=0;
  double temp=0, r1=0, r2=0;
  double mu=0;
  double weight=0.0;

  /* determine lower and upper reff's */
  temp = (reff - phase->r0) / phase->dr;
  /* this is a new version of getting i1 and i2, it is much quicker, but less precise */
  i1 = (int) (temp);
  i2 = i1;
  if ( (float) i1 != (float) temp ) i2++;
  
  /* careful with rounding errors at the array boundaries */
  i1 = (i1>phase->n-1?phase->n-1:i1);
  i2 = (i2>phase->n-1?phase->n-1:i2);
  
  i1 = (i1<0?0:i1);
  i2 = (i2<0?0:i2);

  /* adjacent radii */
  r1 = phase->r0 + (double) i1 * phase->dr;
  r2 = phase->r0 + (double) i2 * phase->dr;

  if (i2>i1)
    weight = ( r2 - reff ) * phase->sca[i1] /
      ( ( r2 - reff ) * phase->sca[i1] + ( reff - r1 ) * phase->sca[i2] );
  else
    weight = 1.00;

  if (uvspec_random() < weight)
    mu = sc_mu (phase->iphase[i1], 1, scaled, 0., 0.);
  else
    mu = sc_mu (phase->iphase[i2], 1, scaled, 0., 0.);

  return mu;

}


/***********************************************************************************/
/* Function: reflection                                                   @62_30i@ */
/* Description:                                                                    */
/*  Return 1 with a probability of albedo and 0 with a probability of (1-albedo).  */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static inline int reflection(double albedo) 
{
  if (uvspec_random() < albedo)
    return 1;

  return 0;
}
  

/***********************************************************************************/
/* Function: emission                                                     @62_30i@ */
/* Description:                                                                    */
/*  Return 1 with a probability of (1-albedo) and 0 with a probability of albedo.  */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static inline int emission(double albedo) 
{
  if (uvspec_random() < albedo)
    return 0;

  return 1;
}


/***********************************************************************************/
/* Function: random_tau                                                   @62_30i@ */
/* Description:                                                                    */
/*  Calculate a random optical thickness with a pdf defined                        */
/*  by Lambert-Beer's law.                                                         */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

inline double random_tau ()
{
  double temp = uvspec_random();

  if (temp>0)    /* careful with 0 */
    return (-log(temp));
  else 
    return DBL_MAX;
}


/***********************************************************************************/
/* Function: random_direction                                             @62_30i@ */
/* Description:                                                                    */
/*  Calculate a random direction.                                                  */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void random_direction (direction *dir)
{
  double phi     = sc_Isotropic_phi();
  double cotheta = sc_Isotropic_mu();
  
  double sitheta = sqrt( 1.0 - cotheta*cotheta );
  
  dir->dx[0] = cosd(phi)*sitheta;
  dir->dx[1] = sind(phi)*sitheta;
  dir->dx[2] = cotheta;

  dir->cotheta = fabs(dir->dx[2]);
  
  hitflag(dir);
}


/***********************************************************************************/
/* Function: random_Isotropic_normal                                      @62_30i@ */
/* Description:                                                                    */
/*  Calculate a random upward direction relative to normal vector.                 */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void random_Isotropic_normal (direction *dir, double *norm)
{     
  int i=0;
  double phi=0, mu=0;

  /* surface/elevation upward normal vector */
  for (i=0; i<3; i++)
    dir->dx[i] = norm[i];

  /* Lambertian reflection */
  phi = sc_Isotropic_phi();
  mu  = sc_Isotropic_upward_mu();

  /* calculate new direction */
  new_direction (mu, phi, dir, 0.);
}


/***********************************************************************************/
/* Function: random_variation_above                                       @62_30i@ */
/* Description:                                                                    */
/*  Calculate a random value within a certain interval above the original value    */
/*  A more sophisticated version would include a distribution function             */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Bernhard Reinhardt                                                      */
/*                                                                        @i62_30@ */
/***********************************************************************************/

double random_variation_above (double orig_value, double interval_width)
{
  return orig_value + interval_width * uvspec_random();
}


/***********************************************************************************/
/* Function: random_theta_above                                           @62_30i@ */
/* Description:                                                                    */
/*  Calculate a random theta within a certain interval above the original value    */
/*  Distribute according to p(theta)=sin(theta) 			           */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Bernhard Reinhardt                                                      */
/*                                                                        @i62_30@ */
/***********************************************************************************/

double random_theta_above (double orig_value, double interval_width)
{
  double norm_factor;	
  norm_factor = cosd(orig_value) - cosd(orig_value + interval_width);
  return acosd(cosd(orig_value) - norm_factor*uvspec_random());
}


/***********************************************************************************/
/* Function: random_cone_mu                                               @62_30i@ */
/* Description:                                                                    */
/*  Calculate a random mu within a certain cone around mu = 1                      */
/*  A more sophisticated version would include a detector sensitivity function     */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

double random_cone_mu (double *cotheta)
{
  return 1. - ( 1. - *cotheta ) * uvspec_random();
}


/***********************************************************************************/
/* Function: random_cone_mu_gauss                                         @62_30i@ */
/* Description:                                                                    */
/*  Calculate a random mu within a certain cone around mu = 1                      */
/*  A more sophisticated version would include a detector sensitivity function     */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

double random_cone_mu_gauss (double *cotheta)
{
  double alpha = uvspec_random();

  if (alpha > 1.0-1e-7)
    alpha=1.0-1e-7;
  
  double result = 1.0 / sqrt(1.0 - ( 1. / (*cotheta * *cotheta) - 1.0) * log(1.0-alpha));
  
  if (result > 1.0)
    result = 1.0;
  
  return result;
}


/***********************************************************************************/
/* Function: random_Lambertian_normal                                     @62_30i@ */
/* Description:                                                                    */
/*  Random (Lambertian) reflection at the surface.                                 */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void random_Lambertian_normal (direction *dir, double *norm)
{     
  int i=0;
  double phi=0, mu=0;

  /* surface/elevation upward normal vector */
  for (i=0; i<3; i++)
    dir->dx[i] = norm[i];

  /* Lambertian reflection */
  phi = sc_Isotropic_phi();
  mu  = sc_Lambertian_mu();

  /* calculate new direction */
  new_direction (mu, phi, dir, 0.);
}


/***********************************************************************************/
/* Function: new_direction                                                @62_30i@ */
/* Description:                                                                    */
/*  Calculate new direction at angles (mu, phi) to the original direction.         */
/*  This is an optimized version, thus being much faster than the old version      */
/*  It produces the same results within 1e-12                                      */
/*  In the case of mu=-1, it even produces the correct result! (mu=-1)             */
/*  The old version of new_direction produced mu=+1 in that case                   */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void new_direction (double mu, double phi,
                    direction *dir, double phi_inc)
{
  
  int i=0;

  double u[3], v[3], e[3];
  double siphi=0, cophi=0, sitheta=0, cotheta=0;
  double u2inv=0;

  /********************************/
  /* u,v are vectors of length 1; */
  /* u,v,dx form a left-handed    */
  /* orthogonal system            */
  /********************************/

  u[2] = + sqrt( dir->dx[0] * dir->dx[0] + dir->dx[1] * dir->dx[1] );

  if (u[2] == 0.) {
    v[0] = -cosd(phi_inc)*dir->dx[2];
    v[1] = sind(phi_inc)*dir->dx[2];
    v[2] = 0.;
  } else {
    u2inv = 1./u[2];

    v[0] = - dir->dx[1] * u2inv;
    v[1] = + dir->dx[0] * u2inv;
    v[2] =   0.;
  }

  u[0] = - dir->dx[2] * v[1];
  u[1] = + dir->dx[2] * v[0];

  /****************************/
  /* create a new vector from */
  /* u, v, and dx             */
  /****************************/

  cotheta = mu;
  sitheta = sqrt( 1.0 - mu*mu );

  cophi = cosd(phi);
  siphi = sind(phi);

  e[0] = sitheta * siphi;
  e[1] = sitheta * cophi;
  e[2] = cotheta;

  for (i=0; i<3; i++)
    dir->dx[i] =  e[0] * u[i] + e[1] * v[i] + e[2] * dir->dx[i];

  /* cosine of solar zenith angle for radiance calculation */
  dir->cotheta = fabs (dir->dx[2]);

  hitflag (dir);
}




/***********************************************************************************/
/* Function: ddprofile3D                                                  @62_30i@ */
/* Description:                                                                    */
/*  Allocate memory for struct ddprofile3D and initialize structure.               */ 
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static ddprofile3D *calloc_ddprofile3D (int Nz, int Nx, int Ny, int *threed)
{
  int kc=0, ic=0;
  ddprofile3D *p = calloc(1, sizeof(ddprofile3D));
  if (p==NULL)
    return NULL;

  
  p->Nz = Nz;
  p->Nx = Nx;
  p->Ny = Ny;

  p->threed = calloc((size_t) Nz, sizeof(int));
  if (p->threed==NULL)
    return NULL;
  
  for (kc=0; kc<Nz; kc++) {
    p->threed[kc] = threed[kc];
    if (p->threed[kc]>=1)
      p->nthreed++;
  }

  /*
  p->aer = calloc((size_t) Nz, sizeof(double **));
  if (p->aer==NULL)
    return NULL;
  */

  /*
  p->ozo = calloc((size_t) Nz, sizeof(double **));
  if (p->ozo==NULL)
    return NULL;
  */

  /*
  p->ice = calloc((size_t) Nz, sizeof(double **));
  if (p->ice==NULL)
    return NULL;
  */

  p->cld = calloc((size_t) Nz, sizeof(double **));
  if (p->cld==NULL)
    return NULL;

  p->tot = calloc((size_t) Nz, sizeof(double **));
  if (p->tot==NULL)
    return NULL;
  
  for (kc=0; kc<Nz; kc++)
    if (p->threed[kc]>=1) {

      /*
      p->aer[kc] = calloc(Nx, sizeof(double *));
      for (ic=0; ic<Nx; ic++)
        p->aer[kc][ic] = calloc((size_t) Ny, sizeof(double));
      */

      /*
      p->ozo[kc] = calloc((size_t) Nx, sizeof(double *));
      for (ic=0; ic<Nx; ic++)
        p->ozo[kc][ic] = calloc((size_t) Ny, sizeof(double));
      */

      p->cld[kc] = calloc((size_t) Nx, sizeof(double *));
      for (ic=0; ic<Nx; ic++)
        p->cld[kc][ic] = calloc((size_t) Ny, sizeof(double));
      
      /*
      p->ice[kc] = calloc((size_t) Nx, sizeof(double *));
      for (ic=0; ic<Nx; ic++)
        p->ice[kc][ic] = calloc((size_t) Ny, sizeof(double));
      */

      p->tot[kc] = calloc((size_t) Nx, sizeof(double *));
      for (ic=0; ic<Nx; ic++)
        p->tot[kc][ic] = calloc((size_t) Ny, sizeof(double));
    }

  return p; 
}


/***********************************************************************************/
/* Function: cp_direction                                                 @62_30i@ */
/* Description:                                                                    */
/*  Copy struct direction.                                                         */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void cp_direction (direction *dest, direction *source)
{
  int i=0;

  for (i=0; i<3; i++)
    dest->dx[i] = source->dx[i];

  for (i=0; i<3; i++)
    dest->hit[i] = source->hit[i];

  dest->cotheta = source->cotheta; 
}


/***********************************************************************************/
/* Function: cp_optical_depth                                             @62_30i@ */
/* Description:                                                                    */
/*  Copy struct optical_depth                                                      */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static void cp_optical_depth (optical_depth *dest, optical_depth *source)
{
  dest->mol = source->mol;
  dest->aer = source->aer;
  dest->cld = source->cld;
  dest->ice = source->ice;
  dest->tot = source->tot;
}


/***********************************************************************************/
/* Function: cp_photon_struct                                             @62_30i@ */
/* Description:                                                                    */
/*  Copy struct photon_struct                                                      */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void cp_photon_struct ( photon_struct *dest,
			photon_struct *source, 
			sample_struct *sample,
			int            n_caoth )
{
  int i=0, j=0, k=0, l=0, id=0, ip=0, iv=0, ic=0;
  
  for (i=0; i<3; i++) 
    dest->x[i] = source->x[i];

  cp_direction (&(dest->dir),    &(source->dir));
  cp_direction (&(dest->dir00),   &(source->dir00));

  dest->ic = source->ic;
  dest->jc = source->jc;
  dest->kc = source->kc;
  dest->jcg = source->jcg;
  dest->iv_alis = source->iv_alis;

  dest->photon_status = source->photon_status;
  for (i=0; i<3; i++) 
    dest->reallyhit[i] = source->reallyhit[i];

  dest->SC_mode  = source->SC_mode;
  dest->DDIS_SC_mode  = source->DDIS_SC_mode;
  dest->RIS_mode = source->RIS_mode;
  dest->VIS_mode = source->VIS_mode;

  dest->direct         = source->direct;

  dest->scattercounter  = source->scattercounter;
  dest->reflectcounter  = source->reflectcounter;
  dest->spikewarningcounter  = source->spikewarningcounter;
  dest->photoncounter   = source->photoncounter;
  dest->backward_is     = source->backward_is;
  dest->backward_js     = source->backward_js;
  dest->rayleighcounter = source->rayleighcounter;
  dest->muchoutcounter  = source->muchoutcounter;

  dest->ncircos=source->ncircos;
  for (ic=0;ic<dest->ncircos;ic++)
    dest->tocirco[ic] = source->tocirco[ic];

  dest->update_atmos_coord = source->update_atmos_coord;

  dest->pathlength = source->pathlength;
  cp_optical_depth (&(dest->tauabs), &(source->tauabs));
  dest->tauris = source->tauris;
  dest->p_norm = source->p_norm;
  
  if (sample->spectral_is || sample->concentration_is){
    dest->nlambda = source->nlambda;
    dest->Nc = source->Nc; 
    for (iv=0; iv<dest->nlambda; iv++){
      dest->q_spectral[iv] = source->q_spectral[iv];
      dest->q2_spectral[iv] = source->q2_spectral[iv];
      dest->q_albedo_spectral[iv]=source->q_albedo_spectral[iv];
      /* dest->dtauabs_spectral[iv] =  source->dtauabs_spectral[iv];*/
    }
    for (ic=0; ic<dest->Nc; ic++){
      dest->q_concentration[ic] = source->q_concentration[ic];
      dest->q2_concentration[ic] = source->q2_concentration[ic];
    }
    
    dest->Nz_alis = source->Nz_alis;
    for (k=0; k<dest->Nz_alis; k++) 
      dest->pathlength_per_layer[k] =  source->pathlength_per_layer[k]; 
  }    

  if(sample->boxairmass){
    dest->Nz_alis = source->Nz_alis;
    for (k=0; k<dest->Nz_alis; k++)
      dest->pathlength_per_layer[k] =  source->pathlength_per_layer[k];
  }

  dest->weight     = source->weight;
  dest->phi0       = source->phi0;
  dest->fw_phi0    = source->fw_phi0;
  dest->fw_phi     = source->fw_phi;
  for(ip=0; ip<sample->nstokes; ip++)
    dest->stokes0[ip]     = source->stokes0[ip];
  for(ip=0; ip<sample->nstokes; ip++)
    dest->stokes[ip]     = source->stokes[ip];

  for (i=0; i<sample->nstokes; i++)
    for (j=0; j<sample->nstokes; j++)
      dest->phamat[i][j] = source->phamat[i][j];

#if HAVE_VROOM
  dest->q_isoene   = source->q_isoene;
  dest->special_weight = source->special_weight;

  if (dest->q_jacobian != NULL) { /* FIJ */
    for (j=0; j<2; j++)
      for (k=0; k<n_caoth; k++)
	for (l=0; l<dest->Nz_jac; l++)
	  dest->q_jacobian[j][k][l] = source->q_jacobian[j][k][l];
    for (l=0; l<dest->Nz_jac; l++)
      dest->r_jacobian[l] = source->r_jacobian[l];
    dest->Nz_jac = source->Nz_jac;
  }

  for (id=0; id<sample->Nd; id++)
    dest->pdir[id] = source->pdir[id];
  
  cp_locest (&(dest->lest), &(source->lest), sample, n_caoth);

#if HAVE_LIDAR
  /* CB - copy cb and pss */
  if (sample->coherent_backscatter && sample->LidarLocEst) {
    copy_cohebasca (&(dest->cb), &(source->cb));
    copy_pss(&(dest->pss),&(source->pss));
  }
  dest->tauext_tot = source->tauext_tot;
#endif

  dest->isclone = source->isclone;
  dest->clonescattercounter = source->clonescattercounter;
  dest->escapescattercounter = source->escapescattercounter;
#endif

  /* these two should not be in photon structure !!! CHECK!!! */ 
  dest->ipa = source->ipa;
  dest->maxpathlength = source->maxpathlength;

  /* do not copy photon weight tree !!!??? CHECK!!! */ 
  /* do not copy photon_path !!!??? */

  dest->vis_beta = source->vis_beta;
  dest->risqidd_beta = source->risqidd_beta;

  /* it is quite useful to know the parent of the photon! Used e.g. in escape_probability */
  dest->parent_photon = source;
  
#ifdef CLDPRP
  if (sample->cldprp){
    dest->cldprp.reff_wc  = source->cldprp.reff_wc;
    dest->cldprp.reff_ic  = source->cldprp.reff_ic;
    dest->cldprp.rhit_wc  = source->cldprp.rhit_wc;
    dest->cldprp.rhit_ic  = source->cldprp.rhit_ic;
    dest->cldprp.tau_wc   = source->cldprp.tau_wc;
    dest->cldprp.tau_ic   = source->cldprp.tau_ic;
    dest->cldprp.dxlwc    = source->cldprp.dxlwc;
    dest->cldprp.dylwc    = source->cldprp.dylwc;
    dest->cldprp.dzlwc    = source->cldprp.dzlwc;
    dest->cldprp.dxiwc    = source->cldprp.dxiwc;
    dest->cldprp.dyiwc    = source->cldprp.dyiwc;
    dest->cldprp.dziwc    = source->cldprp.dziwc;
  }
#endif
  
}


/***********************************************************************************/
/* Function: calloc_radiation_field                                       @62_30i@ */
/* Description:                                                                    */
/*  Allocate memory for struct radiation_field and initialize structure.           */ 
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static radiation_field *calloc_radiation_field (int Nx, int Ny, int Nd, int Nr, int Nt, int Np, int std, int cldprp)
{
  int id=0, ir=0, is=0, iy=0;
  
  radiation_field *temp=calloc(1, sizeof(radiation_field));
  

  temp->ndir = calloc((size_t) Nx, sizeof(int *));
  temp->ndn  = calloc((size_t) Nx, sizeof(int *));
  temp->nup  = calloc((size_t) Nx, sizeof(int *));

  temp->edir = calloc((size_t) Nx, sizeof(double *));
  temp->edn  = calloc((size_t) Nx, sizeof(double *));
  temp->eup  = calloc((size_t) Nx, sizeof(double *));

  temp->fdir = calloc((size_t) Nx, sizeof(double *));
  temp->fdn  = calloc((size_t) Nx, sizeof(double *));
  temp->fup  = calloc((size_t) Nx, sizeof(double *));

  if (std) {  /* standard deviations */
    temp->edir2 = calloc((size_t) Nx, sizeof(double *));
    temp->edn2  = calloc((size_t) Nx, sizeof(double *));
    temp->eup2  = calloc((size_t) Nx, sizeof(double *));
    
    temp->fdir2 = calloc((size_t) Nx, sizeof(double *));
    temp->fdn2  = calloc((size_t) Nx, sizeof(double *));
    temp->fup2  = calloc((size_t) Nx, sizeof(double *));
  }

  if (Nd>0) {
    temp->raddir  = calloc((size_t) Nd, sizeof(double ***));
    temp->raddif  = calloc((size_t) Nd, sizeof(double ***));
    temp->radesc  = calloc((size_t) Nd, sizeof(double ***));
    temp->radpat  = calloc((size_t) Nd, sizeof(double ***));
    temp->radpes  = calloc((size_t) Nd, sizeof(double ***));

    if (std) {
      temp->raddir2  = calloc((size_t) Nd, sizeof(double ***));
      temp->raddif2  = calloc((size_t) Nd, sizeof(double ***));
      temp->radesc2  = calloc((size_t) Nd, sizeof(double ***));
      temp->radpat2  = calloc((size_t) Nd, sizeof(double ***));
      temp->radpes2  = calloc((size_t) Nd, sizeof(double ***));

      temp->nraddir  = calloc((size_t) Nd, sizeof(int ***));
      temp->nraddif  = calloc((size_t) Nd, sizeof(int ***));
      temp->nradesc  = calloc((size_t) Nd, sizeof(int ***));
      temp->nradpat  = calloc((size_t) Nd, sizeof(int ***));
      temp->nradpes  = calloc((size_t) Nd, sizeof(int ***));
    }
  }

  for (is=0; is<Nx; is++) {
    temp->ndir[is] = calloc((size_t) Ny, sizeof(int));
    temp->ndn [is] = calloc((size_t) Ny, sizeof(int));
    temp->nup [is] = calloc((size_t) Ny, sizeof(int));

    temp->edir[is] = calloc((size_t) Ny, sizeof(double));
    temp->edn [is] = calloc((size_t) Ny, sizeof(double));
    temp->eup [is] = calloc((size_t) Ny, sizeof(double));
    temp->fdir[is] = calloc((size_t) Ny, sizeof(double));
    temp->fdn [is] = calloc((size_t) Ny, sizeof(double));
    temp->fup [is] = calloc((size_t) Ny, sizeof(double));

    if (std) {
      temp->edir2[is] = calloc((size_t) Ny, sizeof(double));
      temp->edn2 [is] = calloc((size_t) Ny, sizeof(double));
      temp->eup2 [is] = calloc((size_t) Ny, sizeof(double));
      
      temp->fdir2[is] = calloc((size_t) Ny, sizeof(double));
      temp->fdn2 [is] = calloc((size_t) Ny, sizeof(double));
      temp->fup2 [is] = calloc((size_t) Ny, sizeof(double));
    }
  }
  
  for (id=0; id<Nd; id++) {
    temp->raddir[id] = calloc((size_t) Nx, sizeof(double **));
    temp->raddif[id] = calloc((size_t) Nx, sizeof(double **));
    temp->radesc[id] = calloc((size_t) Nx, sizeof(double **));

    temp->radpat[id] = calloc((size_t) Nr, sizeof(double **));
    temp->radpes[id] = calloc((size_t) Nr, sizeof(double **));

    if (std) {
      temp->raddir2[id] = calloc((size_t) Nx, sizeof(double **));
      temp->raddif2[id] = calloc((size_t) Nx, sizeof(double **));
      temp->radesc2[id] = calloc((size_t) Nx, sizeof(double **));
      temp->radpat2[id] = calloc((size_t) Nr, sizeof(double **));
      temp->radpes2[id] = calloc((size_t) Nr, sizeof(double **));
      
      temp->nraddir[id] = calloc((size_t) Nx, sizeof(int **));
      temp->nraddif[id] = calloc((size_t) Nx, sizeof(int **));
      temp->nradesc[id] = calloc((size_t) Nx, sizeof(int **));
      temp->nradpat[id] = calloc((size_t) Nr, sizeof(int **));
      temp->nradpes[id] = calloc((size_t) Nr, sizeof(int **));
    }

    for (is=0; is<Nx; is++) {
      temp->raddir[id][is] = calloc((size_t) Ny, sizeof(double *));
      temp->raddif[id][is] = calloc((size_t) Ny, sizeof(double *));
      temp->radesc[id][is] = calloc((size_t) Ny, sizeof(double *));
      
      if (std) {
        temp->raddir2[id][is] = calloc((size_t) Ny, sizeof(double *));
        temp->raddif2[id][is] = calloc((size_t) Ny, sizeof(double *));
        temp->radesc2[id][is] = calloc((size_t) Ny, sizeof(double *));

        temp->nraddir[id][is] = calloc((size_t) Ny, sizeof(int *));
        temp->nraddif[id][is] = calloc((size_t) Ny, sizeof(int *));
        temp->nradesc[id][is] = calloc((size_t) Ny, sizeof(int *));
      }
      
      for (iy=0; iy<Ny; iy++){
        temp->raddir[id][is][iy] = calloc((size_t) Np, sizeof(double));
        temp->raddif[id][is][iy] = calloc((size_t) Np, sizeof(double));
        temp->radesc[id][is][iy] = calloc((size_t) Np, sizeof(double));
        
        if (std) {
          temp->raddir2[id][is][iy] = calloc((size_t) Np, sizeof(double));
          temp->raddif2[id][is][iy] = calloc((size_t) Np, sizeof(double));
          temp->radesc2[id][is][iy] = calloc((size_t) Np, sizeof(double));
          
          temp->nraddir[id][is][iy] = calloc((size_t) Np, sizeof(int));
          temp->nraddif[id][is][iy] = calloc((size_t) Np, sizeof(int));
          temp->nradesc[id][is][iy] = calloc((size_t) Np, sizeof(int));
        }
      }  
    }
    for (ir=0; ir<Nr; ir++) {
      temp->radpat[id][ir] = calloc((size_t) Nt, sizeof(double *));
      temp->radpes[id][ir] = calloc((size_t) Nt, sizeof(double *));
      
      if (std) {
        temp->radpat2[id][ir] = calloc((size_t) Nt, sizeof(double *));
        temp->radpes2[id][ir] = calloc((size_t) Nt, sizeof(double *));
        
        temp->nradpat[id][ir] = calloc((size_t) Nt, sizeof(int *));
        temp->nradpes[id][ir] = calloc((size_t) Nt, sizeof(int *));
      }
      for (iy=0; iy<Np; iy++){
        temp->radpat[id][ir][iy] = calloc((size_t) Np, sizeof(double));
        temp->radpes[id][ir][iy] = calloc((size_t) Np, sizeof(double));
       
        if (std) {
          temp->radpat2[id][ir][iy] = calloc((size_t) Np, sizeof(double));
          temp->radpes2[id][ir][iy] = calloc((size_t) Np, sizeof(double));
          
          temp->nradpat[id][ir][iy] = calloc((size_t) Np, sizeof(int));
          temp->nradpes[id][ir][iy] = calloc((size_t) Np, sizeof(int));
        } 
        
      }
    }
  }

#ifdef CLDPRP
  
  temp->cldprp.sample_cldprp = cldprp;
  
  if (temp->cldprp.sample_cldprp){
    temp->cldprp.reff_wc = (double **)calloc((size_t) Nx,sizeof(double*));
    temp->cldprp.reff_ic = (double **)calloc((size_t) Nx,sizeof(double*));
    temp->cldprp.rhit_wc = (double **)calloc((size_t) Nx,sizeof(double*));
    temp->cldprp.rhit_ic = (double **)calloc((size_t) Nx,sizeof(double*));
    temp->cldprp.tau_wc  = (double **)calloc((size_t) Nx,sizeof(double*));
    temp->cldprp.tau_ic  = (double **)calloc((size_t) Nx,sizeof(double*));
    temp->cldprp.dxlwc      = (double **)calloc((size_t) Nx,sizeof(double*));
    temp->cldprp.dylwc      = (double **)calloc((size_t) Nx,sizeof(double*));
    temp->cldprp.dzlwc      = (double **)calloc((size_t) Nx,sizeof(double*));
    temp->cldprp.dxiwc      = (double **)calloc((size_t) Nx,sizeof(double*));
    temp->cldprp.dyiwc      = (double **)calloc((size_t) Nx,sizeof(double*));
    temp->cldprp.dziwc      = (double **)calloc((size_t) Nx,sizeof(double*));
    temp->cldprp.totweights = (double **)calloc((size_t) Nx,sizeof(double*));
    temp->cldprp.wc_weights = (double **)calloc((size_t) Nx,sizeof(double*));
    temp->cldprp.ic_weights = (double **)calloc((size_t) Nx,sizeof(double*));
    
    for(is = 0; is<(size_t) Nx; is++){
      temp->cldprp.reff_wc[is] = (double*)calloc((size_t) Ny,sizeof(double));
      temp->cldprp.reff_ic[is] = (double*)calloc((size_t) Ny,sizeof(double));
      temp->cldprp.rhit_wc[is] = (double*)calloc((size_t) Ny,sizeof(double));
      temp->cldprp.rhit_ic[is] = (double*)calloc((size_t) Ny,sizeof(double));
      temp->cldprp.tau_wc[is]  = (double*)calloc((size_t) Ny,sizeof(double));
      temp->cldprp.tau_ic[is]  = (double*)calloc((size_t) Ny,sizeof(double));
      temp->cldprp.dxlwc[is]      = (double*)calloc((size_t) Ny,sizeof(double));
      temp->cldprp.dylwc[is]      = (double*)calloc((size_t) Ny,sizeof(double));
      temp->cldprp.dzlwc[is]      = (double*)calloc((size_t) Ny,sizeof(double));
      temp->cldprp.dxiwc[is]      = (double*)calloc((size_t) Ny,sizeof(double));
      temp->cldprp.dyiwc[is]      = (double*)calloc((size_t) Ny,sizeof(double));
      temp->cldprp.dziwc[is]      = (double*)calloc((size_t) Ny,sizeof(double));
      temp->cldprp.totweights[is] = (double*)calloc((size_t) Ny,sizeof(double));
      temp->cldprp.wc_weights[is] = (double*)calloc((size_t) Ny,sizeof(double));
      temp->cldprp.ic_weights[is] = (double*)calloc((size_t) Ny,sizeof(double));

    }
  }
  
#endif  
  
  return temp;
}

/***********************************************************************************/
/* Function: calloc_radiation_field_t                                              */
/* Description:                                                                    */
/*  Allocate spectral radiation field (for calculation of high resolution          */ 
/*  spectra using importance sampling).                                            */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Claudia Emde                                                            */
/*                                                                                 */
/***********************************************************************************/
static radiation_field_t *calloc_radiation_field_t(int Nc, int Nx, int Ny, int Np, int Nv)
{
  int ic=0, ip=0, is=0, js=0;

  radiation_field_t *temp=calloc(1, sizeof(radiation_field_t));

  temp->Nc=Nc;
  temp->Nv=Nv;
  temp->Np=Np; 

  temp->rad_t = calloc((size_t) Nc, sizeof(double ****));
  for (ic=0; ic<Nc; ic++){
    temp->rad_t[ic]=calloc((size_t) Nx, sizeof(double ***));
    for (is=0; is<Nx; is++){
      temp->rad_t[ic][is]=calloc((size_t) Ny, sizeof(double **));
      for (js=0; js<Ny; js++){
	temp->rad_t[ic][is][js]=calloc((size_t) Np, sizeof(double *));
	for (ip=0; ip<Np; ip++){
	  temp->rad_t[ic][is][js][ip] = calloc((size_t) Nv, sizeof(double));
        }
      }
    }
  }
    
  return temp;

}

/***********************************************************************************/
/* Function: free_profile                                                 @62_30i@ */
/* Description:                                                                    */
/*  Free memory of struct profile.                                                 */ 
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static void free_profile (profile *p)
{
  int isp=0;

  for (isp=0; isp<=p->n_caoth; isp++)
    free (p->prof[isp]);

  free (p->prof);
  free (p);
}


/***********************************************************************************/
/* Function: free_profile3D                                               @62_30i@ */
/* Description:                                                                    */
/*  Free memory of struct profile3D.                                               */ 
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static void free_profile3D (profile3D *p)
{
  int isp=0, kc=0, ic=0;

  for (isp=0; isp<=p->n_caoth; isp++) {
    if (p->tocalloc[isp]) {
      for (kc=0; kc<p->Nz; kc++)
	if (p->threed[isp][kc]>=1) {
	  for (ic=0; ic<p->Nx; ic++)
	    free(p->prof [isp][kc][ic]);
	  free(p->prof [isp][kc]);
	}
      free(p->prof [isp]);
    }
    free (p->threed[isp]);
  }

  free (p->tocalloc);
  free (p->prof);
  free (p->threed);
  free (p->nthreed);
  free (p);
}


/***********************************************************************************/
/* Function: free_ddprofile3D                                             @62_30i@ */
/* Description:                                                                    */
/*  Free memory of struct ddprofile3D.                                             */ 
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static void free_ddprofile3D (ddprofile3D *p)
{
  int kc=0, ic=0;

  for (kc=0; kc<p->Nz; kc++)
    if (p->threed[kc]>=1) {

      /*
      for (ic=0; ic<p->Nx; ic++)
        free(p->aer[kc][ic]);
      free(p->aer[kc]);
      */

      /*
      for (ic=0; ic<p->Nx; ic++)
        free(p->ozo[kc][ic]);
      free(p->ozo[kc]);
      */

      /*
      for (ic=0; ic<p->Nx; ic++)
        free(p->ice[kc][ic]);
      free(p->ice[kc]);
      */

      for (ic=0; ic<p->Nx; ic++)
        free(p->cld[kc][ic]);
      free(p->cld[kc]);

      for (ic=0; ic<p->Nx; ic++)
        free(p->tot[kc][ic]);
      free(p->tot[kc]);
    }
  
  
  /* free (p->aer); */
  /* free (p->ozo); */
  /* free (p->ice); */
  free (p->cld);
  free (p->tot);

  free (p->threed);

  free (p);
}


/***********************************************************************************/
/* Function: free_radiation_field                                         @62_30i@ */
/* Description:                                                                    */
/*  Free memory of struct radiation_field.                                         */ 
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static void free_radiation_field (radiation_field *res, int Nx, int Nr, int Nd, int Ny, int std)
{
  int id=0, ir=0, is=0, iy=0;

  for (is=0; is<Nx; is++) {
    free(res->ndir[is]);
    free(res->ndn [is]);
    free(res->nup [is]);
    free(res->edir[is]);
    free(res->edn [is]);
    free(res->eup [is]);
    free(res->fdir[is]);
    free(res->fdn [is]);
    free(res->fup [is]);

    if (std) {
      free(res->edir2[is]);
      free(res->edn2 [is]);
      free(res->eup2 [is]);
      free(res->fdir2[is]);
      free(res->fdn2 [is]);
      free(res->fup2 [is]);
    }      
  }

  for (id=0; id<Nd; id++) {
    for (is=0; is<Nx; is++) {
      for (iy=0; iy<Ny; iy++){
        free (res->raddir[id][is][iy]);
        free (res->raddif[id][is][iy]);
        free (res->radesc[id][is][iy]);
        
        if (std) {
          free (res->raddir2[id][is][iy]);
          free (res->raddif2[id][is][iy]);
          free (res->radesc2[id][is][iy]);
          
          free (res->nraddir[id][is][iy]);
          free (res->nraddif[id][is][iy]);
          free (res->nradesc[id][is][iy]);
        }
      }
      
      free (res->raddir[id][is]);
      free (res->raddif[id][is]);
      free (res->radesc[id][is]);
      
      if (std) {
        free (res->raddir2[id][is]);
        free (res->raddif2[id][is]);
        free (res->radesc2[id][is]);
        
        free (res->nraddir[id][is]);
        free (res->nraddif[id][is]);
        free (res->nradesc[id][is]);
      }
    }
    
    for (ir=0; ir<Nr; ir++) {
      for (iy=0; iy<Ny; iy++){
        free (res->radpat[id][ir][iy]);
        free (res->radpes[id][ir][iy]);
        
        if (std) {
          free (res->radpat2[id][ir][iy]);
          free (res->radpes2[id][ir][iy]);
          
          free (res->nradpat[id][ir][iy]);
          free (res->nradpes[id][ir][iy]);
        }
      }
      
      free (res->radpat[id][ir]);
      free (res->radpes[id][ir]);
      
      if (std) {
        free (res->radpat2[id][ir]);
        free (res->radpes2[id][ir]);
        
        free (res->nradpat[id][ir]);
        free (res->nradpes[id][ir]);
      }
    }
    
    free (res->raddir[id]);
    free (res->raddif[id]);
    free (res->radesc[id]);
    free (res->radpat[id]);
    free (res->radpes[id]);

    if (std) {
      free (res->raddir2[id]);
      free (res->raddif2[id]);
      free (res->radesc2[id]);
      free (res->radpat2[id]);
      free (res->radpes2[id]);

      free (res->nraddir[id]);
      free (res->nraddif[id]);
      free (res->nradesc[id]);
      free (res->nradpat[id]);
      free (res->nradpes[id]);
    }
  }

  free (res->ndir);
  free (res->ndn);
  free (res->nup);
  free (res->edir);
  free (res->edn);
  free (res->eup);
  free (res->fdir);
  free (res->fdn);
  free (res->fup);

  free (res->raddir);
  free (res->raddif);
  free (res->radesc);
  free (res->radpat);
  free (res->radpes);

  if (std) {
    free (res->edir2);
    free (res->edn2);
    free (res->eup2);
    free (res->fdir2);
    free (res->fdn2);
    free (res->fup2);

    free (res->raddir2);
    free (res->raddif2);
    free (res->radesc2);
    free (res->radpat2);
    free (res->radpes2);

    free (res->nraddir);
    free (res->nraddif);
    free (res->nradesc);
    free (res->nradpat);
    free (res->nradpes);
  }

#ifdef CLDPRP
  
  if (res->cldprp.sample_cldprp){
    for (is=0; is< Nx; is++) {
      free (res->cldprp.totweights[is]);
      free (res->cldprp.wc_weights[is]);
      free (res->cldprp.ic_weights[is]);
      free (res->cldprp.reff_wc[is]);
      free (res->cldprp.reff_ic[is]);
      free (res->cldprp.rhit_wc[is]);
      free (res->cldprp.rhit_ic[is]);
      free (res->cldprp.tau_wc[is]);
      free (res->cldprp.tau_ic[is]);
      free (res->cldprp.dxlwc[is]);
      free (res->cldprp.dylwc[is]);
      free (res->cldprp.dzlwc[is]);
      free (res->cldprp.dxiwc[is]);
      free (res->cldprp.dyiwc[is]);
      free (res->cldprp.dziwc[is]);
    }
    free (res->cldprp.totweights);
    free (res->cldprp.wc_weights);
    free (res->cldprp.ic_weights);
    free (res->cldprp.reff_wc);
    free (res->cldprp.reff_ic);
    free (res->cldprp.rhit_wc);
    free (res->cldprp.rhit_ic);
    free (res->cldprp.tau_wc);
    free (res->cldprp.tau_ic);
    free (res->cldprp.dxlwc);
    free (res->cldprp.dylwc);
    free (res->cldprp.dzlwc);
    free (res->cldprp.dxiwc);
    free (res->cldprp.dyiwc);
    free (res->cldprp.dziwc);
  }
  
#endif
  
  free (res);
}


/***********************************************************************************/
/* Function: compare_levels                                               @62_30i@ */
/* Description:                                                                    */
/*  Compare two altitude grids.                                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int compare_levels (float *z1, int n1,
                           float *z2, int n2)
{
  int kc=0;

  if (n1!=n2) {
    fprintf (stderr, "Error, number of layers differing, %d vs %d\n", n1, n2);
    return -1;
  }

  for (kc=0; kc<=n1; kc++)
    if (fabs(z1[kc] - z2[kc]) > MC_EPSILON * z1[kc]) {
      fprintf (stderr, "Error, difference at level %d, altitude %g vs. %g\n",   
               kc, z1[kc], z2[kc]);
      return -1;
    }

  return 0;  /* if o.k. */
}
      

/***********************************************************************************/
/* Function: mc_add_optical_depth                                         @62_30i@ */
/* Description:                                                                    */
/*  Add an optical depth to tauabs, using the extinction coefficient kabs          */
/*  of layer kc and the pathlength.                                                */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: ...., last big modification by Robert Buras                             */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static inline void mc_add_optical_depth (atmosphere_struct *atmos,
					 photon_struct     *p,
                                         double             length,
					 int                calc_tauext )
{
  /* commented four fields out to save computational time */

  /* 1D */
  /* p->tauabs.aer += length * (atmos->kabs->aer)[p->kc]; */
  /* p->tauabs.mol += length * (atmos->kabs->ozo)[p->kc]; */
  /* p->tauabs.ice += length * (atmos->kabs->ice)[p->kc]; */

  /* old version */
  /* if (is_threed==1) {/\* 3D layer *\/ */
  /*   /\* p->tauabs.cld += length * (atmos->kabs3D->cld)[p->kc][p->ic][p->jc]; *\/ */
  /*   p->tauabs.tot += length * atmos->kabs3D->prof [MCCAOTH_TOT][p->kc][p->ic][p->jc]; */
  /* } */
  /* else { */
  /*   /\* p->tauabs.cld += length * (atmos->kabs->cld)[p->kc]; *\/ */
  /*   p->tauabs.tot += length * atmos->kabs->prof [MCCAOTH_TOT][p->kc]; */
  /* } */

  p->tauabs.tot += length * get_kabs(atmos, p, MCCAOTH_TOT);
  p->cb.tausca  += length * get_kscaIS(atmos, p, MCCAOTH_TOT);
  if (calc_tauext)
    p->tauext_tot += length * (get_ksca(atmos, p, MCCAOTH_TOT) + get_kabs(atmos,p, MCCAOTH_TOT));

}


/***********************************************************************************/
/* Function: mc_add_absorption_3D                                         @62_30i@ */
/* Description:                                                                    */
/*  Calculate the fraction of radiation absorbed in pixel (kc, ic,jc) and add      */
/*  it to the total amount absorbed.                                               */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static void mc_add_absorption_3D (profile3D *kabs3D, profile *kabs, 
                                  int kc, int ic, int jc, 
                                  double length, double totweight, int std,
                                  ddprofile3D *absorption3D, ddprofile3D *absorption3D2, 
                                  struct tnode **wtree)
{
  double weight=0;

  /* commented four fields out to save computational time */
  
  
  /* 1D */
  /* absorption3D->aer[kc][ic][jc] += i0 * (1.0 - exp(- length * kabs->prof [MCPROF_AER][kc])); */
  /* absorption3D->ozo[kc][ic][jc] += i0 * (1.0 - exp(- length * kabs->prof [MCPROF_MOL][kc])); */

  /* 3D */
  /* absorption3D->cld[kc][ic][jc] += i0 * (1.0 - exp(- length * kabs3D->prof [MCPROF_WC ][kc][ic][jc])); */
  /* absorption3D->ice[kc][ic][jc] += i0 * (1.0 - exp(- length * kabs3D->prof [MCPROF_IC ][kc][ic][jc])); */

  weight = totweight * (1.0 - exp ( - length * kabs3D->prof [MCCAOTH_TOT][kc][ic][jc] ) );
  absorption3D->tot[kc][ic][jc] += weight;

  if (std)
    *wtree = addtree_stddev (*wtree, &(absorption3D2->tot[kc][ic][jc]), weight);
}


/***********************************************************************************/
/* Function: elevation                                                    @62_30i@ */
/* Description:                                                                    */
/*  Calculate surface elevation using bi-linear interpolation.                     */
/*  ATTENTION: x and y are trucated coordinates, that is, relative to the left     */
/*  and lower edges of the pixel.                                                  */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static double elevation (surface surf, double x, double y)
{
  return surf.a*x + surf.b*y + surf.c*x*y + surf.d;
}
  

/***********************************************************************************/
/* Function: cross_bilinear                                               @62_30i@ */
/* Description:                                                                    */
/*  Calculate crossing points between a straight line and the bi-linear surface    */
/*  for a given pixel.                                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int cross_bilinear (surface surf, double *a, double *b, 
                           double x0, double y0, 
                           double *solution)
{
  double alpha=0, beta=0, gamma=0;
  double disc=0, sqdisc=0;
  double r1=0, r2=0;

  
  /* reset solution */
  solution[0] = -1;
  solution[1] = -1;


  gamma =   surf.a * a[0]
          - surf.a * x0
          + surf.b * a[1]
          - surf.b * y0
          + surf.c * a[0] * a[1]
          + surf.c * x0   * y0
          - surf.c * a[0] * y0
          - surf.c * a[1] * x0
          + surf.d 
          - a[2];

  beta =    surf.a * b[0]
          + surf.b * b[1]
          + surf.c * (
                         a[0] * b[1]
                       + a[1] * b[0]
                       - b[0] * y0
                       - b[1] * x0)
          - b[2];

  alpha =   surf.c * b[0] * b[1];

  if (alpha==0) {
    if (beta==0)
      return 0;
    else {
      solution[0] = -gamma/beta;
      return 1;
    }
  }
  else {
    /* solve quadratic equation */
    disc = beta*beta - 4.0*alpha*gamma;
    
    if (disc<0)   /* no solution */
      return 0;
    else {
      if (disc==0) {
        solution[0] = -beta/2.0/alpha;
        return 1;
      }
      else {
        sqdisc = sqrt(disc);
      
        if (beta<0)
          r1 = 2.0*gamma / (-beta + sqdisc);
        else
          r1 = (-beta - sqdisc)/2.0/alpha;
          

        if (beta<0)
          r2 = (-beta + sqdisc)/2.0/alpha;
        else 
          r2 = 2.0*gamma / (-beta - sqdisc);


        /* sort in ascending order */
        if (r1<r2) {
          solution[0] = r1;
          solution[1] = r2;
        }
        else {
          solution[0] = r2;
          solution[1] = r1;
        }

        return 2;
      }
    }
  }

 
  fprintf (stderr, "\nFATAL error, arrived at the end of cross_bilinear()\n");
  return -1;   /* this should not happen */
}


/***********************************************************************************/
/* Function: cross_surface                                                @62_30i@ */
/* Description:                                                                    */
/*  Determine crossing point of a photon with the surface.                         */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int cross_surface (photon_struct *p, double step,
		   elevation_struct *elev,
		   int bcond,
		   double *gamma)
{
  int i=0, status=0;

  int ie=0, je=0;
  double X_elev=0, Y_elev=0;

  double alpha[2]  = {0,0};
  int reallyhit[2] = {0,0};
  double xp[2] = {0,0};
  double stepact = 0, hop=0;

  double xact[3]  = {p->x[0], p->x[1], p->x[2]};
  double n[3]     = {0,0,0};

  double X=0, Y=0;
  
  int used2nd=0, heureka=0;

  double solution[2] = {0,0};

#ifdef MUCHOUT
  if(p->muchoutcounter==MUCHOUTPHOTON)
    fprintf(stderr,"counters %d %d %d %d --- cross_surface %e %e %e\n",
	    p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,p->x[0],p->x[1],p->x[2]);
#endif

  /* photon cannot cross the surface if the whole path is above */
  /* the maximum surface elevation                              */
  
  if (p->x[2] > elev->surfmax + MC_EPSILON &&  
      p->x[2] + step * p->dir.dx[2] > elev->surfmax + MC_EPSILON)
    return 0;

#ifdef MUCHOUT
  if(p->muchoutcounter==MUCHOUTPHOTON)
    fprintf(stderr,"counters %d %d %d %d --- could cross_surface\n",
	    p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode);
#endif

  /* determine elevation pixel */
  elev_coord (p, elev, &ie, &je);
  
  /* hmm, in the version before March 2, 2008 the right and upper edges were treated differently: */
  /* when the photon was at the right/upper boundary, X_elev/Y_elev was set to the right/upper    */
  /* boundary and afterwards ie/je was reduced by one. Now we reduce before (in elev_coord) and   */
  /* X_elev/Y_elev is therefore to the left/lower of the pixel which seems correct                */

  X_elev = (double) ie * elev->delX;
  Y_elev = (double) je * elev->delY;

  while (1) { 
    status = cross_bilinear (elev->surf[ie][je], xact, p->dir.dx, 
                             X_elev,
                             Y_elev,
                             solution);
    
    if (status<0) {
      fprintf (stderr, "FATAL error:\n");
      fprintf (stderr, "cross_bilinear returned -1!\n");
      return -1;
    }
    
    used2nd=0;

    if (status>0) {    /* if there is a solution */
      *gamma = solution[0];  
      
      /* use 2nd solution if 1st is smaller than 0 */
      if (status==2 && *gamma<0) {
        *gamma = solution[1];
        used2nd = 1;
      }

      /* proceed only if a positive solution exists */
      if (*gamma >= 0) {
        
        if (-MC_EPSILON < *gamma && stepact + *gamma < step+MC_EPSILON) {
          xp[0] = xact[0] + *gamma * p->dir.dx[0];
          xp[1] = xact[1] + *gamma * p->dir.dx[1];
          
          if (X_elev <= xp[0] && xp[0] <= X_elev + elev->delX &&
              Y_elev <= xp[1] && xp[1] <= Y_elev + elev->delY) {
            
            /* check if the photon is going towards the surface or */
            /* away from the surface.                              */
            
            /* upward normal */
            n[0] = -elev->surf[ie][je].a - elev->surf[ie][je].c * (xp[1] - Y_elev);
            n[1] = -elev->surf[ie][je].b - elev->surf[ie][je].c * (xp[0] - X_elev);
            n[2] = 1.0;
            
            /* only if the scalar product between surface normal and  */
            /* photon direction is larger than 0, the photon hits the */
            /* surface from outside.                                  */
            
            if ( n[0]*p->dir.dx[0] + n[1]*p->dir.dx[1] + n[2]*p->dir.dx[2] < 0 ) {
              heureka=1;
              *gamma += stepact;
              
              break;   /* done */
            }
            else {

              if (used2nd==0) {
                *gamma = solution[1];

                if (-MC_EPSILON < *gamma && stepact + *gamma < step+MC_EPSILON) {
                  xp[0] = xact[0] + *gamma * p->dir.dx[0];
                  xp[1] = xact[1] + *gamma * p->dir.dx[1];
                  
                  if (X_elev <= xp[0] && xp[0] <= X_elev + elev->delX &&
                      Y_elev <= xp[1] && xp[1] <= Y_elev + elev->delY) {
            
                    /* check if the photon is going towards the surface or */
                    /* away from the surface.                              */
                    
                    /* upward normal */
                    n[0] = -elev->surf[ie][je].a - elev->surf[ie][je].c * (xp[1] - Y_elev);
                    n[1] = -elev->surf[ie][je].b - elev->surf[ie][je].c * (xp[0] - X_elev);
                    n[2] = 1.0;
            
                    /* only if the scalar product between surface normal and  */
                    /* photon direction is larger than 0, the photon hits the */
                    /* surface from outside.                                  */
            
                    if ( n[0]*p->dir.dx[0] + n[1]*p->dir.dx[1] + n[2]*p->dir.dx[2] < 0 ) {
                      heureka=1;
                      *gamma += stepact;
              
                      break;   /* done */
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    /* stop if solution found */
    if (heureka==1)
      break;


    /* no need to check any further if vertical movement */
    if (p->dir.hit[0]==-1 && p->dir.hit[1]==-1) 
      break;


    /* advance to next 2D elevation grid point */

    /* reset reallyhit */
    reallyhit[0] = 0;
    reallyhit[1] = 0;

    if (p->dir.hit[0] >= 0) {   /* x-axis */
      X = (double) (ie + p->dir.hit[0]) * elev->delX;
      alpha[0] = (X - xact[0]) / p->dir.dx[0];
      
      xp[1] = xact[1] + alpha[0] * p->dir.dx[1];
      
      if (Y_elev <= xp[1] && xp[1] <= Y_elev + elev->delY)
        reallyhit[0] = 1;
      else 
        reallyhit[0] = 0;
    }
    
    if (p->dir.hit[1] >= 0) {   /* y-axis */
      Y = (double) (je + p->dir.hit[1]) * elev->delY;
      alpha[1] = (Y - xact[1]) / p->dir.dx[1];
    
      xp[0] = xact[0] + alpha[1] * p->dir.dx[0];
      
      if (X_elev <= xp[0] && xp[0] <= X_elev + elev->delX)
        reallyhit[1] = 1;
      else 
        reallyhit[1] = 0;
    }
              
              
    /* adjust indices */
    if (reallyhit[0]) {
      hop = alpha[0];
      
      if (p->dir.hit[0]==0)
        ie--;
      else
        ie++;

      X_elev = (double) ie * elev->delX;
    }
  
    if (reallyhit[1]) {
      hop = alpha[1];

      if (p->dir.hit[1]==0)
        je--;
      else
        je++;

      Y_elev = (double) je * elev->delY;
    }

    for (i=0; i<3; i++)
      xact[i] += hop * p->dir.dx[i];
    

    /* finally, adjust data to grid to avoid roundoff errors */
    if (reallyhit[0]) {
      if (fabs(xact[0] - X) > MC_EPSILON) {
        fprintf (stderr, "PANIC, expected xact[0] = %g to be X = %g in cross_surface().\n",
                 xact[0], X);
        return -1;
      }
      else
        xact[0] = X;
    }

    if (reallyhit[1]) {
      if (fabs(xact[1] - Y) > MC_EPSILON) {
        fprintf (stderr, "PANIC, expected xact[1] = %g to equal Y = %g in cross_surface().\n",
                 xact[1], Y);
        return -1;
      }
      else
        xact[1] = Y;
    }


    /* stop if step is exceeded */
    stepact += hop;
    if (stepact>step)
      break;
    

    /* adjust for periodic boundary conditions */
    switch (bcond) {
    case MCBCOND_PERIODIC:
      
      if (ie==-1) {
        ie = elev->Nx-1;
        X_elev = (double) ie * elev->delX;
        xact[0] += (double) elev->Nx * elev->delX;
      }
      else {
        if (ie==elev->Nx) {
          ie = 0;
          X_elev = (double) ie * elev->delX;
          xact[0] -= (double) elev->Nx * elev->delX;
        }
      }
            
      if (je==-1) {
        je = elev->Ny-1;
        Y_elev = (double) je * elev->delY;
        xact[1] += (double) elev->Ny * elev->delY;
      }
      else {
        if (je==elev->Ny) {
          je = 0;
          Y_elev = (double) je * elev->delY;
          xact[1] -= (double) elev->Ny * elev->delY;
        }
      }
      break;

    case MCBCOND_MIRROR: 
      fprintf (stderr, "Error, mirroring boundary conditions not implemented!\n");
      return -1;

    case MCBCOND_ABSORB: 
      if (ie==-1 || ie==elev->Nx || je==-1 || je==elev->Ny) {
	/*fprintf(stderr,"absorbing surface photon \n");*/
        p->photon_status = MCSTATUS_OUTOFDOMAIN;
	p->weight = 0.0;
	return 0;
      }

      break;
    default:
      fprintf (stderr, "Error, unknown boundary conditions %d\n", bcond);
      return -1;
    }
  }

#ifdef MUCHOUT
  if(p->muchoutcounter==MUCHOUTPHOTON)
    fprintf(stderr,"counters %d %d %d %d --- result of cross_surface: heureka %d gamma %e\n",
	    p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
	    heureka,*gamma);
#endif

  return heureka;
}


/***********************************************************************************/
/* Function: cp_1D_to_3D                                                  @62_30i@ */
/* Description:                                                                    */
/*  Copy 1D profile to 3D profile, assuming constant optical properties within     */
/*  each layer.                                                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int cp_1D_to_3D (profile3D *p3D, profile *p)
{
  int isp=0, kc=0, ic=0, jc=0;

  if (p3D->n_caoth > p->n_caoth) {
    fprintf (stderr,"ERROR, cannot copy profile, dimensions different! %d %d\n",
	     p3D->n_caoth, p->n_caoth);
    return -1;
  }

  for (isp=0; isp<=p3D->n_caoth; isp++)
    if (p3D->tocalloc[isp])
      for (kc=0; kc<p3D->Nz; kc++)
	if (p3D->threed[isp][kc]>=1)
	  for (ic=0; ic<p3D->Nx; ic++)
	    for (jc=0; jc<p3D->Ny; jc++)
	      p3D->prof [isp][kc][ic][jc] = p->prof [isp][kc];

  return 0;
}


/***********************************************************************************/
/* Function: setup_thermal                                                @62_30i@ */
/* Description:                                                                    */
/*  Initialize thermal emission; calculate emission of each box                    */
/*  (emissivity * planck function), determine maximum emission and the total       */
/*  energy emitted by the atmosphere.                                              */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int setup_thermal (float         ***temper,
			  int              temper3d,
			  int              cldNx,
			  int              cldNy,
			  int              nlyr,
			  float            wvnmlo,
			  float            wvnmhi,
                          profile         *kabs,
			  profile3D       *kabs3D,
			  float           *Z, 
                          int              source,
			  int              absorption,
			  int              thermal_heating_method,
                          profile3D      **Bplanck,
                          double          *maxemis,
			  double          *Watm,
			  surftemp_struct *surftemp,
			  int              backward,
			  float            zstart,
			  int             *backemis_kc,
                          int              nlambda_abs,
                          float           *lambda,
                          double        ***Bplanck_spectral,
                          int              quiet,
			  long int         nphotons,
			  sample_struct    *sample)   /* **CK  2013.08.27 Added nphotons and sample_struct, added cldNx and CldNy*/
{
  int ic=0, jc=0, lc=0;
  int it=0, jt=0;
 /* **CK 2013.08.27 Add new variables for thermal backward heating rate calculation */
  int i=0;  
  int phot_fourth = 0, n_miss=0;
  int use_dEnet = 0;  /* checking weather dEnet is used */
  static int isNx = 0;  /*store number of x gridboxes of actual layer */
  double a_xy= 0.0, a_xz= 0.0, a_yz= 0.0, a = 0.0; /* area fraction of grid box face and sum of the three fractions */
  double v = 0.0;  /* grid box volume */
  double rho = 0.0; /* random number */
  double taux = 0.0, tauy = 0.0, tauz = 0.0, taumin = 0.0;  /* optical thickness in x,y and z direction and minimum optical thickness of those three */
  double tau_crossover = 0.0; /* crossover optical thickness to switch methods */
 /* **CK END */

  double emis=0;
  double **Ws=0, add=0;
  

  float temp=0; 
  /*  float plkavg=0; */

  int iv=0;
  
  int *tocalloc=NULL;
  int **threed=NULL;
  int isp=0; 

  /* reset output data */
  *maxemis=0;
  *Watm=0;

  /* 3D temperature is defined for full atmosphere, therefore all levels will be allocated for profile3D BPlanck*/ 
  tocalloc=calloc(2, sizeof(int));
  ASCII_calloc_int(&threed, 2, nlyr+1); 
  for (isp=0; isp<=1; isp++){
    tocalloc[isp]=1;
    for (lc=0; lc<nlyr+1; lc++){
      threed[isp][lc]=1;
    }
  }
  
  /* need nlyr+1 points because emission is defined */ 
  /* at levels, not at layers         */
  *Bplanck = calloc_profile3D (1, tocalloc, nlyr+1, kabs3D->Nx, 
			       kabs3D->Ny, threed);
  if (*Bplanck==NULL) {
    fprintf (stderr, "Error allocating memory for Bplanck\n");
      return -1;
  }
  
  /* 3DAbs Planck needs to be extended to 3D */
  if (nlambda_abs > 0){
    *Bplanck_spectral=calloc(nlambda_abs, sizeof(double*));
    for (iv=0; iv<nlambda_abs; iv++)
      (*Bplanck_spectral)[iv] = calloc(nlyr+1, sizeof(double)); 
  }
    
  if (source == MCSRC_THERMAL_ATMOSPHERE || source == MCSRC_THERMAL_BACKWARD) {/*TZ bt*/
    
    /* thermal emission:                                         */
    /* calculate wavelength integrated Planck function at levels */
    for (lc=0; lc<nlyr+1; lc++) {
      for (ic=0; ic<kabs3D->Nx; ic++){
	for (jc=0; jc<kabs3D->Ny; jc++){ 
	  /*  F77_FUNC (cplkavg, CPLKAVG) (&wvnmlo, &wvnmhi, &temper[lc], &plkavg); */
	  /*       (*Bplanck)->prof [MCCAOTH_TOT][lc] = plkavg; */
	  if(temper3d)
	    (*Bplanck)->prof [MCCAOTH_TOT][lc][ic][jc] = c_planck_func1(wvnmlo, wvnmhi, temper[ic][jc][lc]);
	  else
	    (*Bplanck)->prof [MCCAOTH_TOT][lc][ic][jc] = c_planck_func1(wvnmlo, wvnmhi, temper[0][0][lc]);
	}
      }
      if (nlambda_abs > 1){
	for (iv=0; iv<nlambda_abs; iv++){
          (*Bplanck_spectral)[iv][lc]=c_planck_func1(1e7/lambda[iv]-0.5, 1e7/lambda[iv]+0.5, temper[ic][jc][lc]);
        }
      }
    } 
    
    /* determine average total column emissivity: required for */
    /* absolute calibration of irradiance and radiance         */
    
    Ws = calloc((size_t) kabs3D->Nx, sizeof(double *));
    for (ic=0; ic<kabs3D->Nx; ic++)
      Ws[ic] = calloc((size_t) kabs3D->Ny, sizeof(double));
    
    /* allocate memory for column integrated emitted flux */
    for (lc=0; lc<nlyr; lc++) {
      if (kabs3D->threed[MCCAOTH_TOT][lc]>=1) {
        for (ic=0; ic<kabs3D->Nx; ic++) {
          for (jc=0; jc<kabs3D->Ny; jc++) {
            Ws[ic][jc] += 0.5 * ( (*Bplanck)->prof [MCCAOTH_TOT][lc][ic][jc]
				  + (*Bplanck)->prof [MCCAOTH_TOT][lc+1][ic][jc] ) *
              kabs3D->prof [MCCAOTH_TOT][lc][ic][jc] * (Z[lc+1] - Z[lc]);
          }
        }
      }
      else {
        add = 0.5 * ( (*Bplanck)->prof [MCCAOTH_TOT][lc][0][0]
		      + (*Bplanck)->prof [MCCAOTH_TOT][lc+1][0][0] ) *
          kabs->prof [MCCAOTH_TOT][lc] * (Z[lc+1] - Z[lc]);
        
        for (ic=0; ic<kabs3D->Nx; ic++)
          for (jc=0; jc<kabs3D->Ny; jc++)
            Ws[ic][jc] += add;
      }
    }
    
    /* average Ws over all columns                                         */
    /* ??? need to check if this gives the correct number for 3D calcs ??? CHECK!!! */ 
    for (ic=0; ic<kabs3D->Nx; ic++)
      for (jc=0; jc<kabs3D->Ny; jc++) 
        *Watm += Ws[ic][jc];
    
    *Watm /= (double) (kabs3D->Nx * kabs3D->Ny);
    *Watm *= (4.0*PI);  /* integrate radiance over 4 pi */
    
    /* free memory */
    for (ic=0; ic<kabs3D->Nx; ic++)
      free(Ws[ic]);
    
    free(Ws);
    
    /* determine maximum emission (= emissivity * Bplanck) */
    for (lc=0; lc<nlyr; lc++) {
      
      if (kabs3D->threed[MCCAOTH_TOT][lc]>=1) {
        for (ic=0; ic<kabs3D->Nx; ic++) {
          for (jc=0; jc<kabs3D->Ny; jc++) {
            emis = ( (*Bplanck)->prof [MCCAOTH_TOT][lc+1][ic][jc] > (*Bplanck)->prof [MCCAOTH_TOT][lc][ic][jc] ? 
                     (*Bplanck)->prof [MCCAOTH_TOT][lc+1][ic][jc] : (*Bplanck)->prof [MCCAOTH_TOT][lc][ic][jc] )
	           * kabs3D->prof [MCCAOTH_TOT][lc][ic][jc];
            if (emis > *maxemis)
              *maxemis = emis;
          }
        }
      }
      else {
        emis = ( (*Bplanck)->prof [MCCAOTH_TOT][lc+1][0][0] > (*Bplanck)->prof [MCCAOTH_TOT][lc][0][0] ? 
                 (*Bplanck)->prof [MCCAOTH_TOT][lc+1][0][0] : (*Bplanck)->prof [MCCAOTH_TOT][lc][0][0] )
	  * kabs->prof [MCCAOTH_TOT][lc];
        if (emis > *maxemis)
          *maxemis = emis;
      }
    }  
    
   if (source!=MCSRC_THERMAL_BACKWARD) {/*TZ bt*/
     /* normalize Bplanck to maximum probability;    */
     /* attention: after this normalization, Bplanck */
     /* is of course no longer the Planck function   */
     
     if (*maxemis>0) 
       for (lc=0; lc<=nlyr; lc++){
	 for (ic=0; ic<kabs3D->Nx; ic++) {
          for (jc=0; jc<kabs3D->Ny; jc++) {
	    (*Bplanck)->prof [MCCAOTH_TOT][lc][ic][jc] /= *maxemis;
	  }
	 }
       }
   }/*TZ bt*/
  }
  
  
  if (absorption == MCFORWARD_ABS_EMISSION) {
    
    /* thermal emission:                                         */
    /* calculate wavelength integrated Planck function at levels */
    
    for (lc=0; lc<nlyr+1; lc++) {
      for (ic=0; ic<kabs3D->Nx; ic++) {
	for (jc=0; jc<kabs3D->Ny; jc++) {
	  /* F77_FUNC (cplkavg, CPLKAVG) (&wvnmlo, &wvnmhi, &temper[lc], &plkavg); */
	  /*       (*Bplanck)->prof [MCCAOTH_TOT][lc] = plkavg; */
	  (*Bplanck)->prof [MCCAOTH_TOT][lc][ic][jc] = c_planck_func1(wvnmlo, wvnmhi, temper[ic][jc][lc]);
	}
      }
    }
  }

  /* 2D surface temperature distribution */
  if (surftemp->surf2D)
    for (it=0; it<surftemp->Nx; it++)
      for (jt=0; jt<surftemp->Ny; jt++) {
        temp=surftemp->temp2D[it][jt];
        /*  F77_FUNC (cplkavg, CPLKAVG) (&wvnmlo, &wvnmhi, &temp, &plkavg); */
        /*         surftemp->plkavg2D[it][jt] = plkavg; */
        surftemp->plkavg2D[it][jt] = c_planck_func1(wvnmlo, wvnmhi, temp);
      }

  if (backward == MCBACKWARD_EMIS || backward == MCBACKWARD_HEAT)  {

    /* determine vertical box index kc */
    /* ??? this could be made faster! CHECK!!! */ 
    if (zstart > 0)
      temp = zstart; /* zstart contains zout[0] (swap in setup_sample2D)*/
    else 
      temp = Z[0];

    *backemis_kc=-1;
    /* Determine the internal z-index of the zstart given in input file */
    while (Z[++(*backemis_kc)] < temp);

    if (!quiet)
      fprintf (stderr, 
	       " ... calculating emission for backward heating rate calculations in layer %d\n",
	       *backemis_kc);
  }


  
  if (backward == MCBACKWARD_HEAT)  {
    /* ********************************************************************************* */
    /* **CK 2013.08.27 Setting flag for thermal backward heating rate calculation method */
    /* ********************************************************************************* */
 
    /* ** free heat_flag if already allocated and allocate new! */
      
    if (sample->heat_flag != NULL){
      for (ic=0; ic<isNx; ic++){
	free (sample->heat_flag[ic]);
      }
      free (sample->heat_flag);
    }
      
    sample->heat_flag = calloc(sample->Nx, sizeof(int *));
    for (ic=0; ic<sample->Nx; ic++)
      sample->heat_flag[ic] = calloc(sample->Ny, sizeof(int));


    /* ******************************************************************************************************** */
    /* **CK 2013.08.28 Check if sample grid and cloud grid have the same size. Optimized methods only work then */
    /* ******************************************************************************************************** */

    if (thermal_heating_method != MCBACKWARD_HEAT_EMABS) {  

      if ((sample->Nx!=kabs3D->Nx)||(sample->Ny!=kabs3D->Ny)){ 
	fprintf (stderr, " ***********************************************************************\n");
	fprintf (stderr, "  Atmospheric grid and the sample grid do not have the same size! \n");
	fprintf (stderr, "  Therefore, standard EMABS method is used, instead of optimized. \n");
	fprintf (stderr, " ***********************************************************************\n");
	for (ic=0; ic<sample->Nx; ic++) 
	  for (jc=0; jc<sample->Ny; jc++) 	  
	    sample->heat_flag[ic][jc] = MCBACKWARD_HEAT_EMABS;
	return 0;
      }
    }

    /* ************************************************* */
    /* **CK 2013.08.27 EMABS is used. Set flag to EMABS  */
    /* ************************************************* */

    switch (thermal_heating_method) {
    case MCBACKWARD_HEAT_EMABS:
      if (!quiet)
	fprintf (stderr, " EMABS-method (unoptimized) is choosen. \n");

      for (ic=0; ic<sample->Nx; ic++) 
	for (jc=0; jc<sample->Ny; jc++) 	  
	  sample->heat_flag[ic][jc] = MCBACKWARD_HEAT_EMABS;
      
      break;

      /* **************************************************************************** */
      /* **CK 2013.08.27 Check which optimized method is choosen and setup variables. */
      /*            If   EMABS_OPT is used. Set flag to EMABSOPT.                     */
      /*                 DENET is used. Set flag to DENET.                            */
      /*                 HYBRID is used. Set flag to HYBRID. Hybrid is also default.  */
      /* **************************************************************************** */

    case MCBACKWARD_HEAT_EMABSOPT:
      if (!quiet)
	fprintf (stderr, " EMABS_OPT method is choosen. \n");

	
      for (ic=0; ic<kabs3D->Nx; ic++) {
	for (jc=0; jc<kabs3D->Ny; jc++) { 
	      
	  if (kabs3D->threed[MCCAOTH_TOT][*backemis_kc]>=1) { /* if 3d */
	  
	    taux = kabs3D->prof [MCCAOTH_TOT][*backemis_kc][ic][jc] * sample->delX;
	    tauy = kabs3D->prof [MCCAOTH_TOT][*backemis_kc][ic][jc] * sample->delY;
	    tauz = kabs3D->prof [MCCAOTH_TOT][*backemis_kc][ic][jc] * (Z[*backemis_kc+1] - Z[*backemis_kc]);
	    taumin = find_min(taux,tauy,tauz);

	    if (taux*tauy*tauz > 100000000) {
	      fprintf (stderr, " ***********************************************************************\n");
	      fprintf (stderr, "  The optical thickness is very large. Large optical thicknesses        \n");
	      fprintf (stderr, "  need huge amounts of memory if EMABS_OPT is chosen. Please            \n");
	      fprintf (stderr, "  use DENET or HYBRID method and restart!                               \n");
	      fprintf (stderr, " ***********************************************************************\n");
	      return -1;
	    }
	  }
	  if(taumin<=2) /*make sure that grid box can be split in sub cubes, if not, use EMABS*/
	    sample->heat_flag[ic][jc] = MCBACKWARD_HEAT_EMABS;
	  else
	    sample->heat_flag[ic][jc] = MCBACKWARD_HEAT_EMABSOPT;
	}	
      }

      break;

    case MCBACKWARD_HEAT_DENET:
      if (!quiet)
	fprintf (stderr, " DENET method is choosen. \n");

      for (ic=0; ic<sample->Nx; ic++) 
	for (jc=0; jc<sample->Ny; jc++) 	  
	  sample->heat_flag[ic][jc] = MCBACKWARD_HEAT_DENET;
      use_dEnet =1;

      break;

    case MCBACKWARD_HEAT_HYBRID:
      if (!quiet)
	fprintf (stderr, " HYBRID method is choosen. \n");

      /* *** set flag with crossover optical thickness - using EMABS_OPT or DENET Method *** */

      for (ic=0; ic<kabs3D->Nx; ic++) {
	for (jc=0; jc<kabs3D->Ny; jc++) {

	  if (kabs3D->threed[MCCAOTH_TOT][*backemis_kc]>=1) { /* if 3d */
	      
	    taux = kabs3D->prof [MCCAOTH_TOT][*backemis_kc][ic][jc] * sample->delX;
	    tauy = kabs3D->prof [MCCAOTH_TOT][*backemis_kc][ic][jc] * sample->delY;
	    tauz = kabs3D->prof [MCCAOTH_TOT][*backemis_kc][ic][jc] * (Z[*backemis_kc+1] - Z[*backemis_kc]);
	      
	    taumin = find_min(taux,tauy,tauz);
	  }
	  else
	    taumin = kabs->prof [MCCAOTH_TOT][*backemis_kc];
	    
	  tau_crossover = 5.0;	  
	    
	  if (taumin < tau_crossover) {
	    // change 14.9.2017, CK and BM: MCBACKWARD_HEAT_EMABSOPT causes problems with "flat" grid boxes (dx >> dz)
	    // where a lot of memory is allocated in generate_photon_backward_vertical_heat(); therefore
	    // using EMABS instead of EMABS_OPT in accordance with the conclusions in Klinger and Mayer (2014)
	    
	    // sample->heat_flag[ic][jc] = MCBACKWARD_HEAT_EMABSOPT;
	    // if(taumin<=2) /*make sure that grid box can be split in sub cubes, if not, use EMABS*/
	    //   sample->heat_flag[ic][jc] = MCBACKWARD_HEAT_EMABS;
	    
	    sample->heat_flag[ic][jc] = MCBACKWARD_HEAT_EMABS;
	  }
	  else {
	    sample->heat_flag[ic][jc] = MCBACKWARD_HEAT_DENET;
	    use_dEnet =1;
	  }
	}
      }
      break; /* end HYBRID */

    default:
      fprintf (stderr, " ***********************************************************************\n");
      fprintf (stderr, "  Error! No calculation method for thermal backward heating rates is set! \n");
      fprintf (stderr, " ***********************************************************************\n");
      return -1;
    }
    

    /* *********************************************************************************** */
    /* if DENET method is used, or DENET is used within HYBRID, make following calculation */
    /* Set area fraction for photon distribution. Photons are distirbuted according to the */
    /* face area of a grid box which may vary if the grid box is non-cubic.                */
    /* *********************************************************************************** */

    if (kabs3D->threed[MCCAOTH_TOT][*backemis_kc]>=1) {  /* if 3d */ 
      if (use_dEnet) {

	phot_fourth = nphotons/4;
   
	a_xy = sample->delX * sample->delY;
	a_xz = sample->delX * (Z[*backemis_kc+1] - Z[*backemis_kc]);
	a_yz = sample->delY * (Z[*backemis_kc+1] - Z[*backemis_kc]);
    
	v= sample->delX * sample->delY * (Z[*backemis_kc+1] - Z[*backemis_kc]);
    
	a = a_xy + a_xz + a_yz;
    
	sample->weight_heat = (Z[*backemis_kc+1] - Z[*backemis_kc]) / v * a;
    
	a_xy=a_xy/a;
	a_xz=a_xz/a;
	a_yz=a_yz/a;
    
	sample->n_xy = (int) (a_xy*phot_fourth);
	sample->n_xz = (int) (a_xz*phot_fourth);
	sample->n_yz = (int) (a_yz*phot_fourth);
    
	n_miss = phot_fourth-sample->n_xy-sample->n_xz-sample->n_yz;

	for (i=0; i<n_miss; i++) {
      
	  rho=uvspec_random();
      
	  if (rho < ((a_xy*phot_fourth - sample->n_xy)/(phot_fourth-sample->n_xy-sample->n_xz-sample->n_yz)))
	    sample->n_xy++;
	  else if (rho < 1-((a_yz*phot_fourth - sample->n_yz)/(phot_fourth-sample->n_xy-sample->n_xz-sample->n_yz)))
	    sample->n_xz++;
	  else
	    sample->n_yz++;
      
	}

	sample->n_xy=sample->n_xy*4;
	sample->n_xz=sample->n_xz*4;
	sample->n_yz=sample->n_yz*4;
    
	if (sample->n_xy+sample->n_xz+sample->n_yz!=nphotons){
	  fprintf(stderr,"Error, number of photons not consistent!! %d %ld \n", sample->n_xy+sample->n_xz+sample->n_yz ,nphotons);
	  return 0;
	}
      } /* close use_dEnet */

      isNx=kabs3D->Nx;  /* set isNx to number of x-gridboxes; this simulation is performed in 3d mode -- necessary to free sample->heat_flag for next layer */
    } /* end 3d */
    

    else {      /* if 1d */

      if (use_dEnet) {

	phot_fourth = nphotons/4;
    
	sample->n_xy = (int) phot_fourth;

    
	sample->n_xy=sample->n_xy*4;
    
	if (sample->n_xy!=nphotons){
	  fprintf(stderr,"Error, number of photons not consistent!! %d %ld \n", sample->n_xy ,nphotons);
	  return 0;
	}
      } /* close use_dEnet */


      isNx=1; /* set isNx to number of x-gridboxes; this simulation is performed in 3d mode -- necessary to free sample->heat_flag for next layer */
    } /* end 1d */
  } /* close backward  = MCBACKWARD_HEAT*/

  return 0;
}
  

/***********************************************************************************/
/* Function: setup_caoth3D                                                @62_30i@ */
/* Description:                                                                    */
/*  Initialize the 3D part of the atmosphere_struct with the caoth data read from  */
/*  filename.                                                                      */ 
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int setup_caoth3D ( caoth3d_out_struct  *caoth3d,
			   int                  nxcld,
			   int                  nycld,
			   int                  nzcld,
			   double               dxcld,
			   double               dycld,
			   float                wavelength,
			   atmosphere_struct   *atmos,
			   char                *datapath,
			   int                  delta_scaling,
			   int                  spherical3D,
			   int                  spherical3D_scene,
			   double               spherical3D_scene_lon_min,
			   double               spherical3D_scene_lon_max,
			   double               spherical3D_scene_lat_min,
			   double               spherical3D_scene_lat_max,
			   float                r_earth,
			   float                sza,
			   int                  setup_tipa,
			   int                  visualize,
			   int                  sample_cldprp,
			   int                  quiet )
{
  int status=0, i=0, ic=0, jc=0, kc=0, ids=0, iris=0, ivs=0, isp=0, ispo=0, n_caoth=0;
  int found_raytracing_prop=0;
  int *tocalloc;

  /* copy 3D caoth data to final destination */
  atmos->Nx = nxcld;
  atmos->Ny = nycld;

  atmos->delX = dxcld;
  atmos->delY = dycld;

  atmos->scatter_type = calloc ((size_t) atmos->n_caoth+1, sizeof(int));
  
  /* allocate indices array to locate water and ice profiles for subsequent lookups */
  atmos->i_wc = calloc ((size_t) atmos->n_caoth+1, sizeof(int));
  atmos->i_ic = calloc ((size_t) atmos->n_caoth+1, sizeof(int));

  atmos->r_earth = (double) r_earth;

  /* molecular scatter type */
  atmos->scatter_type[MCCAOTH_MOL] = MCSCAT_MOL;

  /* aerosol scatter type is set in setup_aerosol1D */

  /* scatter types of other caoth */
  for (isp=MCCAOTH_FIR; isp<=atmos->n_caoth; isp++) {
    ispo=isp-MCCAOTH_FIR;
    if (caoth3d[ispo].nonHG)
      atmos->scatter_type[isp] = MCSCAT_PFT;
    else
      atmos->scatter_type[isp] = MCSCAT_HG2;
    
    /* transfer indices of water and ice from coath3d to atmos structure 
     for subsequent lookups since caoth arrangement can be arbitrary */
    if ( strncasecmp(caoth3d[ispo].name, "wc", 2 ) == 0 )
      atmos->i_wc[isp] = 1;
    if ( strncasecmp(caoth3d[ispo].name, "ic", 2 ) == 0 )
      atmos->i_ic[isp] = 1;
  }

  atmos->phase = calloc ((size_t) atmos->n_caoth+1, sizeof(phase_function_table *));
  for (isp=MCCAOTH_FIR; isp<=atmos->n_caoth; isp++) {
    ispo=isp-MCCAOTH_FIR;
    atmos->phase[isp] = &(caoth3d[ispo].phase);
  }

  atmos->caoth_name = calloc ((size_t) atmos->n_caoth+1, sizeof(char *));
  atmos->caoth_name[MCCAOTH_TOT] = (char *) calloc (strlen("total")+1, sizeof (char));
  strcpy (atmos->caoth_name[MCCAOTH_TOT], "total");
  atmos->caoth_name[MCCAOTH_MOL] = (char *) calloc (strlen("molecular")+1, sizeof (char));
  strcpy (atmos->caoth_name[MCCAOTH_MOL], "molecular");
  atmos->caoth_name[MCCAOTH_AER] = (char *) calloc (strlen("aerosol")+1, sizeof (char));
  strcpy (atmos->caoth_name[MCCAOTH_AER], "aerosol");
  for (isp=MCCAOTH_FIR; isp<=atmos->n_caoth; isp++) {
    ispo=isp-MCCAOTH_FIR;
    atmos->caoth_name[isp] = (char *) calloc (strlen(caoth3d[ispo].fullname)+1,
					      sizeof (char));
    strcpy (atmos->caoth_name[isp], caoth3d[ispo].fullname);
  }

  /* Set scatter type to Rayleigh for molecular_3d. */
  for (isp=MCCAOTH_FIR; isp<=atmos->n_caoth; isp++)
    if (strcmp(atmos->caoth_name[isp], "profile molecular_3d") == 0 )
      atmos->scatter_type[isp] = MCSCAT_MOL;
  
  if (!quiet) {
    for (isp=MCCAOTH_FIR; isp<=atmos->n_caoth; isp++) {
      if (atmos->scatter_type[isp] == MCSCAT_MOL)
	fprintf(stderr, " ... %s considered with Rayleigh phase function \n",
		atmos->caoth_name[isp]);
      else if (atmos->scatter_type[isp]==MCSCAT_PFT)
	fprintf (stderr, " ... %s considered with explicit phase function\n",
		 atmos->caoth_name[isp]);
      else 
	fprintf (stderr, " ... %s considered with Henyey-Greenstein phase function\n",
		 atmos->caoth_name[isp]);
    }
  }

  /* XXX part of the stuff above here should be somewhere else */

  for (isp=MCCAOTH_FIR; isp<=atmos->n_caoth; isp++) {
    ispo=isp-MCCAOTH_FIR;
    /* check if water caoth file levels agree with the 1D atmospheric profile */
    if (caoth3d[ispo].nthreed>=1) {
      if (!quiet)
	fprintf (stderr, " ... checking 3D %s levels\n",
		 caoth3d[ispo].fullname );
    
      /* make sure that 1D and 3D properties are defined on the same altitude grid */
      status = compare_levels (caoth3d[ispo].zd, caoth3d[ispo].nlyr,
			       atmos->Z, atmos->Nz);
      if (status!=0) {
	fprintf (stderr, "\nFATAL error: z-grid of 1D atmosphere and 3D caoths are different\n");
	return status;
      }
    }
  }

  /********************************************************/
  /* determine threed profile: threed is true for a layer */
  /* if it contains caoth                                 */
  /********************************************************/

  atmos->nthreed = calloc ((size_t) atmos->n_caoth+1, sizeof(int));
  atmos->threed = calloc ((size_t) atmos->n_caoth+1, sizeof(int *));
  for (isp=0; isp<=atmos->n_caoth; isp++)
    atmos->threed[isp] = calloc ((size_t) atmos->Nz, sizeof(int));

  for (kc=0; kc<atmos->Nz; kc++) {
    atmos->threed[MCCAOTH_TOT][kc] = 0;

    for (isp=MCCAOTH_FIR; isp<=atmos->n_caoth; isp++) {
      ispo=isp-MCCAOTH_FIR;
      if ( caoth3d[ispo].nthreed >= 1 ) {
	atmos->threed[isp][kc]         = caoth3d[ispo].threed[kc];
	atmos->threed[MCCAOTH_TOT][kc] += caoth3d[ispo].threed[kc];
      }
    }

    for (isp=0; isp<=atmos->n_caoth; isp++)
      if (atmos->threed[isp][kc]>=1)
	atmos->nthreed[isp]++;
  }

  if (!quiet)
    fprintf (stderr, " ... found %d 3D layers\n", atmos->nthreed[MCCAOTH_TOT]);

  /************************************************************/
  /* Now check if we have 1D and 3D caoths in the same layer. */
  /* We don't want that because it may cause confusion!       */
  /* Exception is in case of spherical3D-scenes, then 1d is   */
  /*   outside the scene                                      */
  /************************************************************/

  status=0;
  for (isp=MCCAOTH_FIR; isp<=atmos->n_caoth; isp++) {
    ispo=isp-MCCAOTH_FIR;
    if (caoth3d[ispo].nthreed>=1)
      for (kc=0; kc<atmos->Nz; kc++) 
	if (caoth3d[ispo].threed[kc]>=1) 
	  if ( atmos->ksca [MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->prof [isp][kc]
	       + atmos->kabs->prof [isp][kc] > 0.0 ) {
	    fprintf (stderr, "Error, found 1D AND 3D %s in layer %d.\n",
		     caoth3d[ispo].fullname, kc);
	    fprintf (stderr, "Don't know how to handle that - must be programming error!\n");
	    status++;
	  }
  }

  if (status!=0)
    return status;

  /* allocate memory for 3D profiles */
  tocalloc = calloc( atmos->n_caoth+1, sizeof(int) );

  /* absorption and extinction coefficients are needed for total only */
  for (isp=0; isp<=atmos->n_caoth; isp++)
    tocalloc[isp]=0;
  if (atmos->nthreed[MCCAOTH_TOT]>=1)
    tocalloc[MCCAOTH_TOT]=1;
  
  n_caoth = 0;
  
#ifdef CLDPRP
  
  if (sample_cldprp) { 
    for (isp=MCCAOTH_FIR; isp<=atmos->n_caoth; isp++){
      if (atmos->i_wc[isp] || atmos->i_ic[isp])
        tocalloc[isp] = 1;
    }
    n_caoth = atmos->n_caoth;
  }
  
#endif

  atmos->kabs3D = calloc_profile3D ( n_caoth,
				     tocalloc,
				     atmos->Nz,
				     atmos->Nx,
				     atmos->Ny,
				     atmos->threed );
  if (atmos->kabs3D==NULL)
    return mem_err_out ( "atmos->kabs3D", ERROR_POSITION );

  atmos->kext3D   = calloc( atmos->nscaDS, sizeof(profile3D **) );
  if (atmos->kext3D==NULL)
    return mem_err_out ( "atmos->kext3D", ERROR_POSITION );
  for (ids=0; ids<atmos->nscaDS; ids++) {
    atmos->kext3D[ids] = calloc( atmos->nscaRIS, sizeof(profile3D *) );
    if (atmos->kext3D[ids]==NULL)
      return mem_err_out ( "atmos->kext3D[]", ERROR_POSITION );
  }
  for (ids=0; ids<atmos->nscaDS; ids++)
    for (iris=0; iris<atmos->nscaRIS; iris++) {
      atmos->kext3D[ids][iris] = calloc( atmos->nscaVIS, sizeof(profile3D) );
      if (atmos->kext3D[ids][iris]==NULL)
	return mem_err_out ( "atmos->kext3D[][]", ERROR_POSITION );
    }
  for (ids=0; ids<atmos->nscaDS; ids++)
    for (iris=0; iris<atmos->nscaRIS; iris++)
      for (ivs=0; ivs<atmos->nscaVIS; ivs++) {
	atmos->kext3D[ids][iris][ivs] = calloc_profile3D ( 0,
							  tocalloc,
							  atmos->Nz,
							  atmos->Nx,
							  atmos->Ny,
							  atmos->threed );
	if (atmos->kext3D[ids][iris][ivs]==NULL)
	  return mem_err_out ( "atmos->kext3D[][][]", ERROR_POSITION );
      }

  /* scattering coefficients are needed for all cases in which threed layers exist */
  for (isp=0; isp<=atmos->n_caoth; isp++)
    tocalloc[isp]=0;
  for (isp=0; isp<=atmos->n_caoth; isp++)
    if (atmos->nthreed[isp]>=1)
      tocalloc[isp]=1;

  atmos->ksca3D = calloc( atmos->nscaDS, sizeof(profile3D *) );
  if (atmos->ksca3D==NULL)
    return mem_err_out ( "atmos->ksca3D", ERROR_POSITION );
  for (ids=0; ids<atmos->nscaDS; ids++) {
    atmos->ksca3D[ids] = calloc( atmos->nscaRIS, sizeof(profile3D) );
    if (atmos->ksca3D [ids]==NULL)
      return mem_err_out ( "atmos->ksca3D[]", ERROR_POSITION );
  }
  for (ids=0; ids<atmos->nscaDS; ids++)
    for (iris=0; iris<atmos->nscaRIS; iris++) {
      atmos->ksca3D [ids][iris] = calloc_profile3D ( atmos->n_caoth,
						    tocalloc,
						    atmos->Nz,
						    atmos->Nx,
						    atmos->Ny,
						    atmos->threed );
      if (atmos->ksca3D [ids][iris]==NULL)
	return mem_err_out ( "atmos->ksca3D[][]", ERROR_POSITION );
    }

  /* HG parameters are needed for HG 3d caoth only */
  for (isp=0; isp<=atmos->n_caoth; isp++)
    tocalloc[isp]=0;
  for (isp=1; isp<=atmos->n_caoth; isp++)
    if (atmos->nthreed[isp]>=1 && !(atmos->scatter_type[isp]==MCSCAT_PFT) )
      tocalloc[isp]=1;

  atmos->g1_3D = calloc_profile3D ( atmos->n_caoth,
				    tocalloc,
				    atmos->Nz,
				    atmos->Nx,
				    atmos->Ny,
				    atmos->threed );
  if (atmos->g1_3D==NULL)
    return mem_err_out ( "atmos->g1_3D", ERROR_POSITION );
  atmos->g2_3D = calloc_profile3D ( atmos->n_caoth,
				    tocalloc,
				    atmos->Nz,
				    atmos->Nx,
				    atmos->Ny,
				    atmos->threed );
  if (atmos->g2_3D==NULL)
    return mem_err_out ( "atmos->g2_3D", ERROR_POSITION );
  atmos->ff_3D = calloc_profile3D ( atmos->n_caoth,
				    tocalloc,
				    atmos->Nz,
				    atmos->Nx,
				    atmos->Ny,
				    atmos->threed );
  if (atmos->ff_3D==NULL)
    return mem_err_out ( "atmos->ff_3D", ERROR_POSITION );

  /* r_eff are needed for pft 3d caoth only */
  for (isp=0; isp<=atmos->n_caoth; isp++)
    tocalloc[isp]=0;
  for (isp=1; isp<=atmos->n_caoth; isp++)
    if (atmos->nthreed[isp]>=1 && atmos->scatter_type[isp]==MCSCAT_PFT)
      tocalloc[isp]=1;

  atmos->reff_3D = calloc_profile3D ( atmos->n_caoth,
				      tocalloc,
				      atmos->Nz,
				      atmos->Nx,
				      atmos->Ny,
				      atmos->threed );
  if (atmos->reff_3D==NULL)
    return mem_err_out ( "atmos->reff_3D", ERROR_POSITION );

#ifdef CLDPRP
	
  if (sample_cldprp){
    atmos->sample_cldprp = 1;
    atmos->dxlwc_3D = calloc_profile3D ( atmos->n_caoth, tocalloc, atmos->Nz, atmos->Nx, atmos->Ny, atmos->threed );
    atmos->dylwc_3D = calloc_profile3D ( atmos->n_caoth, tocalloc, atmos->Nz, atmos->Nx, atmos->Ny, atmos->threed );
    atmos->dzlwc_3D = calloc_profile3D ( atmos->n_caoth, tocalloc, atmos->Nz, atmos->Nx, atmos->Ny, atmos->threed );
  }
  
#endif

  /* copy raytracing properties */
  
  found_raytracing_prop = 0;
  for (isp=MCCAOTH_FIR; isp<=atmos->n_caoth; isp++) {
    ispo=isp-MCCAOTH_FIR;

    if ( strncasecmp(caoth3d[ispo].name, "ic", 2 ) == 0 && caoth3d[ispo].n_raytracing_prop > 0) {
      found_raytracing_prop = 1;
      break;
    }
  }
      
  if (found_raytracing_prop>0) {
    atmos->raytracing_prop = calloc (caoth3d[ispo].n_raytracing_prop, sizeof (crystal_prop_struct));
    atmos->n_raytracing_prop = caoth3d[ispo].n_raytracing_prop;
    
    for (i=0; i<atmos->n_raytracing_prop; i++) {

      fprintf (stderr, "CAOTH = %s\n", caoth3d[ispo].raytracing_prop[i].name);

      
      (atmos->raytracing_prop)[i].name = calloc (strlen(caoth3d[ispo].raytracing_prop[i].name)+1, sizeof(char));
      strcpy (atmos->raytracing_prop[i].name, caoth3d[ispo].raytracing_prop[i].name);
      
      atmos->raytracing_prop[i].fraction = caoth3d[ispo].raytracing_prop[i].fraction;
      atmos->raytracing_prop[i].oriented_fraction = caoth3d[ispo].raytracing_prop[i].oriented_fraction;
      atmos->raytracing_prop[i].angdist_width = caoth3d[ispo].raytracing_prop[i].angdist_width;
      atmos->raytracing_prop[i].orientation_dof = caoth3d[ispo].raytracing_prop[i].orientation_dof;
    }
  }

  
  /**********************/
  /* now fill 3D arrays */
  /**********************/


  /* copy 1D properties to 3D profiles                                   */ 
  /* liquid, ice, and total (=liquid + ice + molecules + aerosol + ...)  */
  /* The newest version only copies what really is needed.               */

  status = cp_1D_to_3D (atmos->kabs3D,  atmos->kabs);
  if (status!=0)
    return err_out ("Error %d returned by cp_1D_to_3D(kabs)\n", status);

  for (ids=0; ids<atmos->nscaDS; ids++)
    for (iris=0; iris<atmos->nscaRIS; iris++) {
      status = cp_1D_to_3D (atmos->ksca3D[ids][iris],  atmos->ksca[ids][iris]);
      if (status!=0)
	return err_out ("Error %d returned by cp_1D_to_3D(ksca)\n", status);
    }

  for (ids=0; ids<atmos->nscaDS; ids++)
    for (iris=0; iris<atmos->nscaRIS; iris++)
      for (ivs=0; ivs<atmos->nscaVIS; ivs++) {
	status = cp_1D_to_3D (atmos->kext3D[ids][iris][ivs],  atmos->kext[ids][iris][ivs]);
	if (status!=0)
	  return err_out ("Error %d returned by cp_1D_to_3D(kext)\n", status);
      }

  status = cp_1D_to_3D (atmos->g1_3D, atmos->g1);
  if (status!=0)
    return err_out ("Error %d returned by cp_1D_to_3D(g1)\n", status);
  status = cp_1D_to_3D (atmos->g2_3D, atmos->g2);
  if (status!=0)
    return err_out ("Error %d returned by cp_1D_to_3D(g2)\n", status);
  status = cp_1D_to_3D (atmos->ff_3D, atmos->ff);
  if (status!=0)
    return err_out ("Error %d returned by cp_1D_to_3D(ff)\n", status);


  status = cp_1D_to_3D (atmos->reff_3D, atmos->reff);
  if (status!=0)
    return err_out ("Error %d returned by cp_1D_to_3D(reff)\n", status);

  /* fill 3D optical properties */
  if (atmos->nthreed[MCCAOTH_TOT]>0) {
    for (isp=MCCAOTH_FIR;isp<=atmos->n_caoth;isp++) {
      ispo=isp-MCCAOTH_FIR;
      if (atmos->nthreed[isp]>=1) {
	for (kc=0; kc<atmos->Nz; kc++) {
	  if (atmos->threed[isp][kc]>=1) {
	    for (ic=0; ic<atmos->Nx; ic++) {
	      for (jc=0; jc<atmos->Ny; jc++) {

		/* for now no IS */
		iris = MCRIS_MODE_NORMAL;

		atmos->ksca3D [MCSC_MODE_NORMAL       ][iris]->prof [isp][kc][ic][jc]
		  = caoth3d[ispo].ext[kc][ic][jc] * caoth3d[ispo].ssa[kc][ic][jc];

		if (delta_scaling>-1)
		  atmos->ksca3D [MCSC_MODE_DELTA_SCALE][iris]->prof [isp][kc][ic][jc]
		    = caoth3d[ispo].ext[kc][ic][jc] * caoth3d[ispo].ssa[kc][ic][jc]
		    * (1.0 - caoth3d[ispo].dscale[kc][ic][jc]);

		if (atmos->scatter_type[isp]==MCSCAT_PFT)
		  atmos->reff_3D ->prof [isp][kc][ic][jc] = caoth3d[ispo].reff[kc][ic][jc];
		else {
		  atmos->g1_3D->prof [isp][kc][ic][jc] = caoth3d[ispo].g1[kc][ic][jc];
		  atmos->g2_3D->prof [isp][kc][ic][jc] = caoth3d[ispo].g2[kc][ic][jc];
		  atmos->ff_3D->prof [isp][kc][ic][jc] = caoth3d[ispo].ff[kc][ic][jc];
		}

#ifdef CLDPRP
                
		if (sample_cldprp) {
                  
                  /* fill 3D absorption properties, when sampling cloud optical properties */
                  
                  if (atmos->scatter_type[isp]==MCSCAT_PFT){
                    atmos->kabs3D ->prof [isp][kc][ic][jc] += caoth3d[ispo].ext[kc][ic][jc]
                    * (1.0 - caoth3d[ispo].ssa[kc][ic][jc]);
                  }
                  
                  
		  /* calculate dx/dy/dz gradients of the lwc/iwc */
                  
		  if (ic<atmos->Nx-1) {
                    atmos->dxlwc_3D->prof [isp][kc][ic][jc] = 
                    (caoth3d[ispo].lwc[kc][ic+1][jc]-caoth3d[ispo].lwc[kc][ic][jc]) / caoth3d[ispo].delX;
                  }
		  if (jc<atmos->Ny-1) {
                    atmos->dylwc_3D->prof [isp][kc][ic][jc] = 
                    (caoth3d[ispo].lwc[kc][ic][jc+1]-caoth3d[ispo].lwc[kc][ic][jc]) / caoth3d[ispo].delY;
                  }
		  if (kc<atmos->Nz-1 && atmos->threed[isp][kc+1]>=1) {
                    atmos->dzlwc_3D->prof [isp][kc][ic][jc] = 
                    (caoth3d[ispo].lwc[kc+1][ic][jc]-caoth3d[ispo].lwc[kc][ic][jc]) / (caoth3d[ispo].zd[kc+1]-caoth3d[ispo].zd[kc]);
                  }
                }
                
#endif

		/* add 3D caoths to total scattering and absorption coefficients */
		for (ids=0; ids<atmos->nscaDS; ids++)
		  atmos->ksca3D [ids][iris]->prof [MCCAOTH_TOT][kc][ic][jc] +=
		    atmos->ksca3D [ids][iris]->prof [isp][kc][ic][jc];
		
		atmos->kabs3D ->prof [MCCAOTH_TOT][kc][ic][jc] += caoth3d[ispo].ext[kc][ic][jc]
		  * (1.0 - caoth3d[ispo].ssa[kc][ic][jc]);
		
		/* if ( strncasecmp(caoth3d[ispo].name, "molecular_3d", 3) && ic == 1 && jc==1  ){ */
		/*   fprintf(stderr, "kc %d Z %.2f ext %.4e \n", kc, atmos->Z[kc], caoth3d[ispo].ext[kc][ic][jc]); */
		/* } */
		
		/* copy into IS modes */
		for (ids=0; ids<atmos->nscaDS; ids++)
		  for (iris=1; iris<atmos->nscaRIS; iris++) {
		    atmos->ksca3D [ids][iris]->prof [isp][kc][ic][jc] =
		      atmos->ksca3D [ids][MCRIS_MODE_NORMAL]->prof [isp][kc][ic][jc];
		    atmos->ksca3D [ids][iris]->prof [MCCAOTH_TOT][kc][ic][jc] +=
		      atmos->ksca3D [ids][MCRIS_MODE_NORMAL]->prof [isp][kc][ic][jc];
		  }

	      } /* end for jc */
	    } /* end for ic */
	  } /* end if atmos->threed[isp][kc] */
	} /* end for kc */
      } /* end if atmos->nthreed[isp] */
    } /* end for isp */
  } /* end if atmos->nthreed[0] */

  /* setup extinction coefficient profile; Virtual IS may not be set up before here! */
  /* the memory consumption here could be reduced if we point instead of callocating */
  if (atmos->nthreed[MCCAOTH_TOT]>=1)
    for (kc=0; kc<atmos->Nz; kc++)
      if (atmos->threed[MCCAOTH_TOT][kc]>=1)
	for (ids=0; ids<atmos->nscaDS; ids++)
	  for (iris=0; iris<atmos->nscaRIS; iris++)
	    for (ivs=0; ivs<atmos->nscaVIS; ivs++)
	      for (ic=0; ic<atmos->Nx; ic++)
		for (jc=0; jc<atmos->Ny; jc++)
		  atmos->kext3D [ids][iris][ivs]->prof [MCCAOTH_TOT][kc][ic][jc]
		    = atmos->ksca3D [ids][iris]->prof [MCCAOTH_TOT][kc][ic][jc];

  /* finally, create X[], Y[] arrays of grid coordinates */
  if (spherical3D) {
#ifdef HAVE_SPHER
    status = setup_spherical3D ( "atmos",
				 atmos->Nx,
				 atmos->Ny,
				 spherical3D_scene,
				 spherical3D_scene_lon_min,
				 spherical3D_scene_lon_max,
				 spherical3D_scene_lat_min,
				 spherical3D_scene_lat_max,
				 &(atmos->delX),
				 &(atmos->delY),
				 &(atmos->xmax),
				 &(atmos->ymax),
				 &(atmos->X),
				 &(atmos->Y),
				 &(atmos->Nyg),
				 &(atmos->Yg),
				 &(atmos->realY),
				 &(atmos->sin2_lat),
				 &(atmos->tan2_lat),
				 &(atmos->jc_eq),
				 quiet );
    if (status)
      return fct_err_out (status, "setup_spherical3D", ERROR_POSITION);
#else
    fprintf(stderr,"Error! you are not allowed to use spherical 3D!\n");
    return -1;
#endif
  }
  else {
    atmos->X = calloc ((size_t) (atmos->Nx+1), sizeof(double));
    if (atmos->X==NULL)
      return mem_err_out ("atmos->X", ERROR_POSITION );
  
    atmos->Y = calloc ((size_t) (atmos->Ny+1), sizeof(double));
    if (atmos->Y==NULL)
      return mem_err_out ("atmos->Y", ERROR_POSITION );

    for (ic=0; ic<=atmos->Nx; ic++)
      atmos->X[ic] = (double) ic * atmos->delX;
  
    for (jc=0; jc<=atmos->Ny; jc++)
      atmos->Y[jc] = (double) jc * atmos->delY;

    /* horizontal boundaries */
    atmos->xmax = (double) atmos->Nx * atmos->delX;
    atmos->ymax = (double) atmos->Ny * atmos->delY;

  }

  /* set Ris Factor to 1 */
  atmos->ris_factor = 1.;

  #if HAVE_OPENGL
  if (visualize)
    GLmystic_write_caoth (atmos);
  #endif

  free(tocalloc);
  return 0;
}


/***********************************************************************************/
/* Function: setup_albedo2D                                               @62_30i@ */
/* Description:                                                                    */
/*  Initialize the 2D albedo_struct with the albedo data read from                 */
/*  filename.                                                                      */ 
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int setup_albedo2D (char *albedo_filename, char *albedo_spectralfilename, char *rpv_filename, char *rossli_filename, int isAmbralsFile,
                           atmosphere_struct *atmos,
                           albedo_struct *albedo,
                           double alb, 
			   float *alb_type,
                           float *rpv_rho0, float *rpv_k, float *rpv_theta, 
                           float *rpv_scale, float *rpv_sigma,
			   float *rpv_t1, float *rpv_t2, 
                           char **rpv_labels, int rpv_nlabels,
			   float *hapke_h,float *hapke_b0, float *hapke_w,
			   float *rossli_iso, float *rossli_vol, float *rossli_geo, 
			   int rossli_hotspot,
                           float u10, float pcl, float xsal, float uphi, int solar_wind,
			   float bpdf_u10, 
                           int spectral, float wavelength, int polarisation,
			   int spherical3D, int spherical3D_scene,
			   double spherical3D_scene_lon_min,
			   double spherical3D_scene_lon_max,
			   double spherical3D_scene_lat_min,
			   double spherical3D_scene_lat_max,
			   int nlambda_abs,
			   double *alis_albedo,
			   double **alis_alb_type,
			   int ixmin,
			   int ixmax,
			   int iymin,
			   int iymax,
                           int quiet)
{

  int il=0, status=0, nfi=0, therewereCaMs=0;
  char tmp_filename   [FILENAME_MAX], tmpchar[6];
  int iv=0;

  if (spectral) {
    if (!quiet)
      fprintf(stderr,"Reading spectrally dependent albedo file...\n");
    nfi=strlen(albedo_filename);
    strcpy (tmp_filename, "");
    strncpy (tmp_filename, albedo_filename, nfi-4);
    tmp_filename[nfi-4]=0;
    strcat (tmp_filename, "_");
    sprintf (tmpchar, "%5d", 10000 + (int) wavelength);
    strcat (tmp_filename, tmpchar+1);
    strcat (tmp_filename, albedo_filename+nfi-4);
  }
  else{
    strcpy (tmp_filename, albedo_filename);
  }

  switch (albedo->method) {
  case MCALB_LAM:
    albedo->albedo = alb; 
    if(nlambda_abs > 0){
      albedo->spectral_albedo=calloc(nlambda_abs, sizeof(double));
      for (iv=0; iv<nlambda_abs; iv++)
	albedo->spectral_albedo[iv] = alis_albedo[iv];
    }
    break;
    
  case MCALB_LAM2D:

    if (!quiet)
      fprintf (stderr, " ... reading 2D albedo data\n");
    
    /* read 2D albedo data */
    status = read_2D_albedo (tmp_filename,
			     &(albedo->Nx), &(albedo->Ny),
			     &(albedo->delX), &(albedo->delY),
			     &(albedo->albedo2D),
			     ixmin, ixmax, iymin, iymax,
			     quiet);
    
    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n", status, tmp_filename);
      return status;
    }
    
    if (!spherical3D) {
      if (!equal ((double) (albedo->Nx) * albedo->delX, (double) atmos->Nx * atmos->delX)) {
	fprintf (stderr, "Error, x-size of albedo area (%g) does not equal x-size of cloud area (%g)\n",
		 (double) (albedo->Nx) * albedo->delX, (double) atmos->Nx * atmos->delX);
	return -1;
      }
    
      if (!equal ((double) (albedo->Ny) * albedo->delY, (double) atmos->Ny * atmos->delY)) {
	fprintf (stderr, "Error, y-size of albedo area (%g) does not equal y-size of cloud area (%g)\n",
		 (double) (albedo->Ny) * albedo->delY, (double) atmos->Ny * atmos->delY);
	return -1;
      }
    }

    break;

  case MCALB_LAM2D_SPECTRAL:

    if (!quiet)
      fprintf (stderr, " ... reading 2D albedo data from %s\n", albedo_spectralfilename);


    /* 2D surface albedo - we "abuse" the RPV structure - should rename it */
    status = read_2D_surface_labels (albedo_spectralfilename, 
				     rpv_labels, rpv_nlabels,
				     &(albedo->Nx), &(albedo->Ny),
				     &(albedo->delX), &(albedo->delY),
				     &(albedo->rpv_index), &(albedo->rpv_count),
				     &therewereCaMs, quiet);
    
    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n", status, albedo_spectralfilename);
      return status;
    }
    
    if (!spherical3D) {
      if (!equal ((double) (albedo->Nx) * albedo->delX, (double) atmos->Nx * atmos->delX)) {
	fprintf (stderr, "Error, x-size of 2D albedo area (%g) does not equal x-size of cloud area (%g)\n",
		 (double) (albedo->Nx) * albedo->delX, (double) atmos->Nx * atmos->delX);
	return -1;
      }
    
      if (!equal ((double) (albedo->Ny) * albedo->delY, (double) atmos->Ny * atmos->delY)) {
	fprintf (stderr, "Error, y-size of 2D albedo area (%g) does not equal y-size of cloud area (%g)\n",
		 (double) (albedo->Ny) * albedo->delY, (double) atmos->Ny * atmos->delY);
	return -1;
      }
    }

    if (!quiet)
      for (il=0; il<rpv_nlabels; il++)
        fprintf (stderr, " ... found %d pixels with %s\n", albedo->rpv_count[il], rpv_labels[il]);

    albedo->alb_type = calloc ((size_t) rpv_nlabels, sizeof(double));
        
    for (il=0; il<rpv_nlabels; il++) {
      if (albedo->rpv_count[il]) {
        albedo->alb_type[il] = alb_type[il];
	fprintf (stderr, "alb_type[%d] = %f\n", il, alb_type[il]);
      }
    }      

    /* spectral albedo type for ALIS */
    if (nlambda_abs > 0) {
      albedo->spectral_alb_type = calloc(rpv_nlabels, sizeof(double *));
      for (il=0; il<rpv_nlabels; il++) {
	albedo->spectral_alb_type[il] =	calloc(nlambda_abs, sizeof(double));

	for (iv=0; iv<nlambda_abs; iv++)
	  albedo->spectral_alb_type[il][iv] = alis_alb_type[il][iv];
      }
    }
    
    break;
    
  case MCALB_RPV:

    if (!albedo->reflectalways) {
      fprintf (stderr, "Error, BRDF only implemented together with reflectalways!\n");
      return -1;
    }

    albedo->rpv = calloc ((size_t) 1, sizeof(rpv_struct));

    /* this is no bug, for RPV 1d, il=0 is correct */
    albedo->rpv[0].rho0  = rpv_rho0[0];
    albedo->rpv[0].k     = rpv_k[0];
    albedo->rpv[0].theta = rpv_theta[0];
    albedo->rpv[0].scale = rpv_scale[0];
    albedo->rpv[0].sigma = rpv_sigma[0];
    albedo->rpv[0].t1    = rpv_t1[0];
    albedo->rpv[0].t2    = rpv_t2[0];
    break;

  case MCALB_RPV2D_SPECTRAL:

    if (!albedo->reflectalways) {
      fprintf (stderr, "Error, BRDF only implemented together with reflectalways!\n");
      return -1;
    }

    if (!quiet)
      fprintf (stderr, " ... reading 2D RPV data from %s\n", rpv_filename);

    status = read_2D_surface_labels (rpv_filename, 
				     rpv_labels, rpv_nlabels,
				     &(albedo->Nx), &(albedo->Ny),
				     &(albedo->delX), &(albedo->delY),
				     &(albedo->rpv_index), &(albedo->rpv_count),
				     &therewereCaMs, quiet);
    
    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n", status, rpv_filename);
      return status;
    }
    
    if (therewereCaMs) {
      if (u10<=0 && !polarisation) {
	fprintf (stderr, "Error, need to specify a wind speed for Cox and Munk!\n");
	return -1;
      }
      if (bpdf_u10<=0 && polarisation) {
	fprintf (stderr, "Error, need to specify a wind speed for Tsang!\n");
	return -1;
      }

      if (!polarisation)
	albedo->u10      = u10;
      else
	albedo->u10      = bpdf_u10;
      albedo->pcl        = pcl;
      albedo->xsal       = xsal;
      albedo->uphi       = uphi;
      albedo->solar_wind = solar_wind;

    }
    
    if (!spherical3D) {
      if (!equal ((double) (albedo->Nx) * albedo->delX, (double) atmos->Nx * atmos->delX)) {
	fprintf (stderr, "Error, x-size of rpv2D area (%g) does not equal x-size of cloud area (%g)\n",
		 (double) (albedo->Nx) * albedo->delX, (double) atmos->Nx * atmos->delX);
	return -1;
      }
    
      if (!equal ((double) (albedo->Ny) * albedo->delY, (double) atmos->Ny * atmos->delY)) {
	fprintf (stderr, "Error, y-size of rpv2D area (%g) does not equal y-size of cloud area (%g)\n",
		 (double) (albedo->Ny) * albedo->delY, (double) atmos->Ny * atmos->delY);
	return -1;
      }
    }

    if (!quiet)
      for (il=0; il<rpv_nlabels; il++)
        fprintf (stderr, " ... found %d pixels with %s\n", albedo->rpv_count[il], rpv_labels[il]);

    albedo->rpv = calloc ((size_t) rpv_nlabels, sizeof(rpv_struct));
        
    for (il=1; il<rpv_nlabels; il++) {  /* first element is reserved for Cox and Munk */
      if (albedo->rpv_count[il]) {
        albedo->rpv[il].rho0  = rpv_rho0[il];
        albedo->rpv[il].k     = rpv_k[il];
        albedo->rpv[il].theta = rpv_theta[il];
	albedo->rpv[il].scale = rpv_scale[il];
	albedo->rpv[il].sigma = rpv_sigma[il];
	albedo->rpv[il].t1    = rpv_t1[il];
	albedo->rpv[il].t2    = rpv_t2[il];
      }
    }      

    break;

  case MCALB_COXANDMUNK:

    if (!albedo->reflectalways) {
      fprintf (stderr, "Error, BRDF only implemented together with reflectalways!\n");
      return -1;
    }

    albedo->u10        = u10;
    albedo->pcl        = pcl;
    albedo->xsal       = xsal;
    albedo->uphi       = uphi;
    albedo->solar_wind = solar_wind;

    break;

  case MCALB_TSANG:
    if (!(albedo->reflectalways && polarisation)) {
      fprintf (stderr, "Error, BPRF only implemented together with reflectalways and polarisation!\n");
      return -1;
    }
    albedo->u10  = bpdf_u10;
    
    break; 
    
  case MCALB_HAPKE:

    if (!albedo->reflectalways) {
      fprintf (stderr, "Error, BRDF only implemented together with reflectalways!\n");
      return -1;
    }

    albedo->hapke = calloc ((size_t) 1, sizeof(hapke_brdf_spec));

    albedo->hapke[0].h  = hapke_h [0];
    albedo->hapke[0].b0 = hapke_b0[0];
    albedo->hapke[0].w  = hapke_w [0];

    break;

  case MCALB_ROSSLI:

    if (!albedo->reflectalways) {
      fprintf (stderr, "Error, BRDF only implemented together with reflectalways!\n");
      return -1;
    }

    albedo->rossli = calloc ((size_t) 1, sizeof(rossli_brdf_spec));

    albedo->rossli[0].iso = rossli_iso[0];
    albedo->rossli[0].vol = rossli_vol[0];
    albedo->rossli[0].geo = rossli_geo[0];
    albedo->rossli[0].hotspot = rossli_hotspot;

    break;

  case MCALB_ROSSLI2D:

    if (!albedo->reflectalways) {
      fprintf (stderr, "Error, BRDF only implemented together with reflectalways!\n");
      return -1;
    }

    if (!quiet)
      fprintf (stderr, " ... reading 2D Ross-Li data from %s\n", rossli_filename);

    status = read_2D_rossli (rossli_filename,
			     &(albedo->Nx), &(albedo->Ny),
			     &(albedo->delX), &(albedo->delY),
			     &(albedo->rossli2D), &therewereCaMs,
			     rossli_hotspot, isAmbralsFile,
			     ixmin, ixmax, iymin, iymax,
			     quiet);

    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n", status, rossli_filename);
      return status;
    }
    
    if (therewereCaMs) {
      if (u10<=0) {
	fprintf (stderr, "Error, need to specify a wind speed for Cox and Munk!\n");
	return -1;
      }

      albedo->u10        = u10;
      albedo->pcl        = pcl;
      albedo->xsal       = xsal;
      albedo->uphi       = uphi;
      albedo->solar_wind = solar_wind;

    }
    
    if (!spherical3D) {
      if (!equal ((double) (albedo->Nx) * albedo->delX, (double) atmos->Nx * atmos->delX)) {
	fprintf (stderr, "Error, x-size of RossLi2D area (%g) does not equal x-size of cloud area (%g)\n",
		 (double) (albedo->Nx) * albedo->delX, (double) atmos->Nx * atmos->delX);
	return -1;
      }
    
      if (!equal ((double) (albedo->Ny) * albedo->delY, (double) atmos->Ny * atmos->delY)) {
	fprintf (stderr, "Error, y-size of RossLi2D area (%g) does not equal y-size of cloud area (%g)\n",
		 (double) (albedo->Ny) * albedo->delY, (double) atmos->Ny * atmos->delY);
	return -1;
      }
    }

    break;

  default:
    fprintf (stderr, "Error, unknown albedo->method %d\n", albedo->method);
    return -1;
  }

#ifdef HAVE_SPHER
  if (spherical3D)
    switch (albedo->method) {
    case MCALB_LAM2D:
    case MCALB_LAM2D_SPECTRAL:
    case MCALB_RPV2D_SPECTRAL:
    case MCALB_ROSSLI2D:
      status = setup_spherical3D ("albedo",
				  albedo->Nx,
				  albedo->Ny,
				  spherical3D_scene,
				  spherical3D_scene_lon_min,
				  spherical3D_scene_lon_max,
				  spherical3D_scene_lat_min,
				  spherical3D_scene_lat_max,
				  &(albedo->delX),
				  &(albedo->delY),
				  &(albedo->xmax),
				  &(albedo->ymax),
				  &(albedo->X),
				  &(albedo->Y),
				  NULL,
				  NULL,
				  NULL,
				  NULL,
				  NULL,
				  NULL,
				  quiet);
      if (status!=0)
	return err_out ("Error %d returned by setup_spherical3D()\n", status);
    case MCALB_LAM:
    case MCALB_RPV:
    case MCALB_COXANDMUNK: 
    case MCALB_HAPKE:
    case MCALB_ROSSLI:
    case MCALB_TSANG:
      break;
    default:
      fprintf (stderr, "Fatal error, unknown albedo type %d\n", albedo->method);
      return -1;
    }
#endif

  return 0;
}


/***********************************************************************************/
/* Function: setup_umu2D                                                  @62_30i@ */
/* Description:                                                                    */
/*  Initialize the 2D umu_struct with the umu data read from                       */
/*  filename.                                                                      */ 
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int setup_umu2D (char *umu_filename, sample_struct *sample, 
			atmosphere_struct *atmos, int quiet)
{
  int status=0;
  int Nx=0, Ny=0;
  double delX=0.0, delY=0.0;

  if (!quiet)
    fprintf (stderr, " ... reading 2D umu data\n");
    
  /* read 2D umu data, identical to albedo! */
  status = read_2D_umu (umu_filename,
			&Nx, &Ny,
			&delX, &delY,
			&(sample->umu2D),
			&(sample->phi2D),
			sample->ixmin,
			sample->ixmax,
			sample->iymin,
			sample->iymax,
			quiet);
    
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, umu_filename);
    return status;
  }

  /* tests */
  if ( !equal ((double) delX, (double) atmos->delX) ||
       !equal ((double) delY, (double) atmos->delY) ||
       !equal ( Nx, atmos->Nx ) ||
       !equal ( Ny, atmos->Ny ) ) {

    fprintf (stderr, "Error, size of umu file size of cloud (%d %d %g %g ne %d %d %g %g)\n",
	     Nx , Ny, delX, delY, atmos->Nx , atmos->Ny, atmos->delX, atmos->delY);
    return -1;
  }

  return 0;
}


/* switch from forward tracing to backward tracing;  */
/* basically, exchange observation and sun direction */
/* as well as start and sampling altitudes           */

static int forward2backward (atmosphere_struct *atmos,
                             float **zout, int *nzout,
                             float *sza, float *phi0,
                             sample_struct *sample,
                             result_struct * result,
                             int quiet)
{
  int kc=0;
  int id=0;
  
  if (!quiet) {
    if ((*zout)[0]==-999.0)
      fprintf (stderr, "BACKWARD! before: starting photons at %gm, sampling at the surface\n", sample->zstart);
    else 
      fprintf (stderr, "BACKWARD! before: starting photons at %gm, sampling at %gm\n", sample->zstart, (*zout)[0]*1000.0);
  }
  
  if (*nzout<1) {
    fprintf (stderr, "Error, no zout defined!\n");
    return -1;
  }

  /* save altitude string for backward sampling file; */
  /* to make sure that the output goes to the same    */
  /* location where the forward results would have    */
  /* gone, determine filename from "forward" zout[0]  */
  /* before switching zout and zstart                 */
  strcpy (sample->backward_altstr, "");

  /* need special treatment for surface */

  if ((*zout)[0] != -999.0) {
    for (kc=0; kc<=atmos->Nz; kc++) {
      if (fabs(atmos->Z[kc] - (*zout)[0]*1000.0) < 1E-3) {
        sprintf (sample->backward_altstr, "%d", kc);
        break;
      }
    }
    if (kc==atmos->Nz+1) {
      fprintf (stderr, "\nFATAL error: altitude grid does not contain level %g\n", (*zout)[0]);
      fprintf (stderr, "which has been specified as output altitude in the sample file.\n");
      return -1;
    }
  }


  switch (sample->backward) {
  case MCBACKWARD_EDIR:
  case MCBACKWARD_EGLOB:
  case MCBACKWARD_EDN:
  case MCBACKWARD_EUP:
  case MCBACKWARD_FDIR:
  case MCBACKWARD_FDN:
  case MCBACKWARD_FUP:
  case MCBACKWARD_ABS:
  case MCBACKWARD_EMIS:
  case MCBACKWARD_HEAT:
  case MCBACKWARD_ACT:
  case MCBACKWARD_EXP:
  case MCBACKWARD_EXN:
  case MCBACKWARD_EYP:
  case MCBACKWARD_EYN:
    /* need to allocate memory for one sampling direction */

    if (!quiet)
      fprintf (stderr, " ... backward sampling irradiance / actinic flux!\n");
    
    /* allocate memory for sample radiances */
    sample->Nd  = 1;
    sample->rad = calloc (1+sample->DoLE, sizeof (radang));

    /* copy uvspec radiance to sample array */
    sample->rad[0].theta = sample->backward_vza;
    sample->rad[0].phi   = sample->backward_phi;
    /* sample->rad[0].alpha = 5.0;  default, can be modified using input option mc_rad_alpha */

    if (sample->DoLE) {
      sample->rad[1].theta = 180.0 - sample->backward_vza;
      sample->rad[1].phi   = sample->backward_phi;
      /*  sample->rad[1].alpha = 5.0; default, can be modified using input option mc_rad_alpha */
    }
    
    break;
    
  case MCBACKWARD_RADIANCE:

    if (!quiet) {
      fprintf (stderr, "BACKWARD!         sza = %f, phi0 = %f\n", *sza, *phi0);
      fprintf (stderr, "BACKWARD!         vza = %f, phi  = %f\n", sample->rad[0].theta, sample->rad[0].phi);
    }

    if (sample->Nd < 1) {
      fprintf (stderr, "Error, need at least one radiance direction for MCBACKWARD_RADIANCE.\n");
      fprintf (stderr, "This indicates a programming error!\n");
      return -1;
    }
      
    sample->Nd = 1;  /* we use only the first direction */
    
    /* exchange zenith and azimuth angles */
    *sza  = sample->backward_sza;
    *phi0 = sample->backward_phi0;
    sample->rad[0].theta = sample->backward_vza;
    sample->rad[0].phi   = sample->backward_phi;

    if (sample->DoLE) {
      sample->rad[1].theta = 180.0 - sample->backward_vza;
      sample->rad[1].phi   = sample->backward_phi;
    }
    
    if (!quiet) {
      fprintf (stderr, "BACKWARD! after:  sza = %f, phi0 = %f\n", *sza, *phi0);
      fprintf (stderr, "BACKWARD!         vza = %f, phi  = %f\n", sample->rad[0].theta, sample->rad[0].phi);
    }

    break;
    
  case MCBACKWARD_NONE:
    break;

  default:
    fprintf (stderr, "Error, unknown backward quantity %d\n", sample->backward);
    return -1;
  }
  
  for (id=0; id<1+sample->DoLE; id++) {
    /* need to initialize new radiance direction */
    sample->rad[id].cosalpha = cosd (sample->rad[id].alpha);
    sample->rad[id].omega    = 2.0 * PI * (1.0 - sample->rad[id].cosalpha);
  
    init_direction (-sind(sample->rad[id].theta),
		    -cosd(sample->rad[id].theta), 
		    sind(sample->rad[id].phi), 
		    cosd(sample->rad[id].phi),
		    &(sample->rad[id].dir));
  }

  /* Backup rad.dir */
  cp_direction(&(sample->rad_dir0),&(sample->rad[0].dir));

  /* exchange zstart / zout */
  sample->zstart = sample->backward_zstart;

  /* if no zout's allocated, allocate one */
  if (*nzout<=0)
    *zout = calloc (1, sizeof(float));
  
  (*zout)[0] = sample->backward_zout;
  *nzout = 1;   /* reset to one because we obviously only look at one direction */
  
  /* no surface sampling for backward */
  sample->surface=0;

  /* store first sampling level - we might need it again */
  /* in particular for forward/backward sampling         */
  for (kc=0; kc<=atmos->Nz; kc++)
    if (sample->sample[kc]) {
      sample->sample_forward_kc = kc;
      break;
    }
  
  /* reset sample flag and determine new sample flag */
  for (kc=0; kc<=atmos->Nz; kc++)
    sample->sample[kc]=0;

  if (*nzout>=1) {
    if ((*zout)[0] != -999.0) {
      for (kc=0; kc<=atmos->Nz; kc++) {
        if (fabs(atmos->Z[kc] - (*zout)[0]*1000.0) < 1E-3) {
          sample->sample[kc]=1;
          break;
        }
      }
      if (kc==atmos->Nz+1) {
        fprintf (stderr, "\nFATAL error: altitude grid does not contain level %g\n", (*zout)[0]);
        fprintf (stderr, "which has been specified as output altitude in the sample file.\n");
        return -1;
      }
    }
  }
  
  /* store backward sampling level */
  sample->sample_backward_kc = kc;

  if (!quiet) {
    if (sample->zstart == -999.0)
      fprintf (stderr, "BACKWARD!         starting photons at the surface, sampling at %gm\n", (*zout)[0]*1000.0);
    else
      fprintf (stderr, "BACKWARD!         starting photons at %gm, sampling at %gm\n", sample->zstart, (*zout)[0]*1000.0);
  }

  /* if no user-defined sample range, */
  /* sample all pixels                */
  if (sample->backward_islower==-999 &&
      sample->backward_jslower==-999 &&
      sample->backward_isupper==-999 &&
      sample->backward_jsupper==-999) {
    sample->backward_islower=0;
    sample->backward_jslower=0;
    sample->backward_isupper=sample->Nx-1;
    sample->backward_jsupper=sample->Ny-1;
  }

  return 0;
}


/***********************************************************************************/
/* Function: setup_sample2D                                               @62_30i@ */
/* Description:                                                                    */
/*  Initialize sample_struct with data read from samplefilename and passed as      */
/*  function arguments.                                                            */ 
/*                                                                                 */
/* Parameters:                                                                     */
/*  float *zout:  array of altitudes to be sampled; overwrites data from           */
/*                samplefilename if nzout>0                                        */
/*  int   nzout:  number of altitudes to be sampled; if nzout<=0, data from        */
/*                samplefilename are used.                                         */
/*                                                                                 */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int setup_sample2D (atmosphere_struct *atmos,
                           float *zout, int nzout,
                           float *sza, float *phi0,
                           float *wavelength,
                           sample_struct *sample,
                           int quiet)
{
  int status=0;
  int id=0;

  /********************************************************/
  /* set start and sampling altitude for forward/backward */
  /********************************************************/

  /* default: forward photons start at TOA */
  sample->forward_zstart = atmos->Z[atmos->Nz];
  sample->zstart         = sample->forward_zstart;
  
  /* backward photons start at the first user-defined zout */
  /* or, if no zout is defined, use -999                   */
  if (nzout>0 && zout[0]>=0)
    sample->backward_zstart=zout[0]*1000;  /* convert from km to m */
  else 
    sample->backward_zstart=-999;

  if (nzout<=0)
    sample->forward_zout = -999;
  else 
    sample->forward_zout = zout[0];

  sample->backward_zout = sample->forward_zstart/1000.0;

  /* quick dirty fix for blitz... BCA: clean up generate_photon,
     setup_sample_grid, forward2backward, and setup_sample2d*/

  if (sample->blitz_position != NULL) {
    *sza  = sample->rad[0].theta;
    *phi0 = 180 + sample->rad[0].phi;
  }

  /**********************************************************************/
  /* save original directions - they might be reversed in backward mode */
  /**********************************************************************/

  sample->forward_sza  = *sza;
  sample->forward_phi0 = *phi0;

  sample->backward_vza = sample->forward_sza;
  sample->backward_phi = sample->forward_phi0;

  if (sample->Nd > 0) {
    sample->forward_vza   = sample->rad[0].theta;
    sample->forward_phi   = sample->rad[0].phi;

    sample->backward_sza  = sample->forward_vza;
    sample->backward_phi0 = sample->forward_phi;
  }

  /* CE probably not needed, CHECK!!! */ 
  /* #ifdef NEW_REFLECT */
  /*   if (sample->backward) */
  /*     sample->phi_target=sample->forward_phi0; */
  /*   else */
  /*     sample->phi_target=sample->forward_phi; */
  /* #endif */

  /* now switch forward to backward */ 
  if (sample->backward) {
    status = forward2backward (atmos,
                               &zout, &nzout,
                               sza, phi0,
                               sample, NULL,
                               quiet);
    if (status!=0)
      return err_out ("Error %d returned by forward2backward()\n", status);
  }

  /* initialize radiance directions, after the first one has */
  /* possibly been modified in the above backward section;   */
  for (id=0; id<sample->Nd; id++) {
    sample->rad[id].cosalpha = cosd (sample->rad[id].alpha);
    sample->rad[id].omega    = 2.0 * PI * (1.0 - sample->rad[id].cosalpha);
    
    init_direction (-sind(sample->rad[id].theta),
                    -cosd(sample->rad[id].theta), 
                    sind(sample->rad[id].phi), 
                    cosd(sample->rad[id].phi),
                    &(sample->rad[id].dir));

    if (sample->rad[id].dir.dx[2] == 0.0 && !sample->spherical ) {
      fprintf (stderr,"Error! Horizontal escape radiances will never end for non-spherical geometry!!! Exiting...\n");
      return -1;
    }
  }

#if HAVE_LIDAR
  /*WARNING: Accuracy of wavelength as double is only 6 digits*/
  sample->wavelength = float_zeros(*wavelength,6);
  sample->wavenumber = (2.e9*PI)/sample->wavelength;
#endif


  return 0;
}


/***********************************************************************************/
/* Function: setup_mc_result                                              @62_30i@ */
/* Description:                                                                    */
/*  Allocate memory for the result structures.                                     */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int setup_mc_result (atmosphere_struct *atmos, sample_struct *sample,
                            result_struct *result, int absorption, int thermal_heating_method,
			    int quiet)
{
  int is=0, js=0, kc=0, ip=0, ic=0,il=0;
  
  /* surface quantities */
  result->surf = calloc_radiation_field(sample->Nx, sample->Ny, sample->Nd, sample->Nr, sample->Nt,
                                        sample->nstokes, sample->std, sample->cldprp);

  
  
  /* altitude profiles */
  result->alt = calloc ((size_t) (atmos->Nz+1), sizeof(radiation_field *));

  if (sample->spectral_is || sample->concentration_is){
    result->surf_t =  calloc_radiation_field_t(atmos->Nc, sample->Nx, sample->Ny,  sample->nstokes, atmos->nlambda_abs);
    result->rad_t = calloc ((size_t) (atmos->Nz+1), sizeof(radiation_field_t *));
  }
  
#if HAVE_LIDAR
  /* lidar quantities */
  if (sample->LidarLocEst)
    result->lidar = calloc_lidar_result_field ( sample->Nli,
						sample->LLE_Nt,
						atmos->n_caoth,
						1+(sample->LLE_channels>0),
						sample->LLE_scatterout+1,
						sample->nstokes,
						sample->LLE_No,
						sample->LLE_Na,
						sample->std,
						sample->LLE_jacobian,
						sample->Nz_jac );
  if (sample->abs_jacobian)
    result->jacobian = calloc_jacobian_result_field (sample->std, sample->Nx, sample->Ny, sample->Nz_jac, atmos->n_caoth);
#endif

  for (kc=0; kc<=atmos->Nz; kc++) {
    if (sample->sample[kc]){
      result->alt[kc] = calloc_radiation_field (sample->Nx, sample->Ny, sample->Nd, sample->Nr, sample->Nt, 
                                                sample->nstokes, sample->std, sample->cldprp);
    
      if (sample->spectral_is || sample->concentration_is )
        result->rad_t[kc]=calloc_radiation_field_t(atmos->Nc, sample->Nx, sample->Ny,  sample->nstokes, atmos->nlambda_abs);
    }
  }
  
  

  /* backward radiance / irradiance */
  if (sample->backward) {

    if (!quiet) 
      fprintf (stderr, " ... allocating %d x %d pixels for backward fields\n", sample->Nx, sample->Ny);

    result->back     = calloc (sample->Nx, sizeof (double **));
    result->back2    = calloc (sample->Nx, sizeof (double **));
    result->backemis = calloc (sample->Nx, sizeof (double *));
    if (sample->backward == MCBACKWARD_HEAT)/* **CK  2013.09.23*/
      if (thermal_heating_method == MCBACKWARD_HEAT_DENET || thermal_heating_method == MCBACKWARD_HEAT_HYBRID){ /* **CK  2013.08.27*/
	result->back_dEnet  = calloc (sample->Nx, sizeof (double **)); 
 	result->back_dEnet2  = calloc (sample->Nx, sizeof (double **));     
      }
    for (is=0; is< sample->Nx; is++) {
      result->back[is]     = calloc (sample->Ny, sizeof (double *));
      result->back2[is]    = calloc (sample->Ny, sizeof (double *));
      result->backemis[is] = calloc (sample->Ny, sizeof (double));
      if (sample->backward == MCBACKWARD_HEAT)/* **CK  2013.09.23*/
	if (thermal_heating_method == MCBACKWARD_HEAT_DENET || thermal_heating_method == MCBACKWARD_HEAT_HYBRID){ /* **CK  2013.08.27*/
	  result->back_dEnet[is]  = calloc (sample->Ny, sizeof (double *)); 
	  result->back_dEnet2[is]  = calloc (sample->Ny, sizeof (double *)); 
	}
      for (js=0; js< sample->Ny; js++) {
        result->back[is][js]     = calloc (sample->nstokes, sizeof (double));
        result->back2[is][js]    = calloc (sample->nstokes, sizeof (double));
	if (sample->backward == MCBACKWARD_HEAT)  /* **CK  2013.09.23*/
	  if (thermal_heating_method == MCBACKWARD_HEAT_DENET || thermal_heating_method == MCBACKWARD_HEAT_HYBRID){  /* **CK  2013.08.27*/
	    result->back_dEnet[is][js]  = calloc (12, sizeof (double));  
	    result->back_dEnet2[is][js]  = calloc (12, sizeof (double));
	  }
      }
    }


    if (sample->spectral_is || sample->concentration_is){
      result->back_t  = calloc( atmos->Nc, sizeof (double ****));
      for (ic=0; ic<atmos->Nc; ic++){
        result->back_t[ic]     = calloc (sample->Nx, sizeof (double ***));
	for (is=0; is< sample->Nx; is++){
	  result->back_t[ic][is]     = calloc (sample->Ny, sizeof (double **));
	  for (js=0; js< sample->Ny; js++) {
	    result->back_t[ic][is][js]   = calloc (sample->nstokes, sizeof (double *));
	    for (ip=0; ip<sample->nstokes; ip++)
	      result->back_t[ic][is][js][ip]    = calloc (atmos->nlambda_abs, sizeof (double));
	  }
        }
      }
    }
  }

  if (absorption) {
    result->absorption3D = calloc_ddprofile3D (atmos->Nz, atmos->Nx, atmos->Ny,
					       atmos->threed[MCCAOTH_TOT]);
    if (result->absorption3D==NULL) {
      fprintf (stderr, "Error allocating memory for absorption3D\n");
      return -1;
    }

    if (sample->std) {
      result->absorption3D2 = calloc_ddprofile3D (atmos->Nz, atmos->Nx, atmos->Ny,
						  atmos->threed[MCCAOTH_TOT]);
      if (result->absorption3D2==NULL) {
        fprintf (stderr, "Error allocating memory for absorption3D\n");
        return -1;
      }
    }      
  }    
  
  if(sample->boxairmass)
    result->pathlength_per_layer_tot=calloc((size_t) atmos->Nz, sizeof(double));
  

  if (sample->ncirc) {
    /* QUICK FIX*/
    result->circcontr = calloc (10, sizeof (double *));
    for (il=0;il<10;il++)
      result->circcontr[il] = calloc (sample->ncirc, sizeof (double));
  }

  result->mish = calloc(1,sizeof(mishchenko_cb_struct));
#if HAVE_LIDAR
  result->lidcb = calloc(1,sizeof(lidar_cb_struct));
#endif

  return 0;
}


/***********************************************************************************/
/* Function: incone                                                       @62_30i@ */
/* Description:                                                                    */
/*  Check if a direction dx lies within the cone define in rad.                    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int incone (double *dx, radang rad) 
{
  
  if (dx[0]*rad.dir.dx[0] + dx[1]*rad.dir.dx[1] + dx[2]*rad.dir.dx[2] >= rad.cosalpha) 
    return 1;
  else 
    return 0;
}


/***********************************************************************************/
/* Function: count_photon                                                 @62_30i@ */
/* Description:                                                                    */
/*  Register a photon; that is, add its contribution to the irradiance,            */
/*  actinic flux, and radiance arrays in result_struct *res.                       */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */ 
/***********************************************************************************/

static void count_photon (radiation_field *res, 
                          photon_struct *p, 
                          sample_struct *sample, 
                          double *totweight, 
                          int surfaceparallel, int incoming)
{
  
  int id=0;
  int ir=0, is=0, it=0, js=0, ip=0;
  
  /* determine sample coordinates */
  sample_coord (p, sample, &is, &js);

  if (p->direct) {   /* direct beam source */
    /* direct irradiance and actinic flux */
    res->ndir[is][js] ++;
    res->edir[is][js] += totweight[0];    res->fdir[is][js] += totweight[0] / p->dir.cotheta;
    
    if (sample->std) {
      p->wtree = addtree_stddev (p->wtree, &(res->edir2[is][js]), totweight[0]);
      p->wtree = addtree_stddev (p->wtree, &(res->fdir2[is][js]), totweight[0]/p->dir.cotheta);
    }
    
    /* direct radiance */
    if (!sample->panorama_forward) {
      for (id=0; id<sample->Nd; id++) {
	/* Stokes components loop not necessary, for randomly oriented particles, direct  */
	/* radiation is not polarized. THIS IS NOT TRUE FOR LIDAR!!! BUG!!! CE, please fix... CHECK!!! */ 
	if (incone(p->dir.dx, sample->rad[id])) {
	  res->raddir[id][is][js][0] += totweight[0] / p->dir.cotheta;
	  
	  if (sample->std)
	    p->wtree = addtree_stddev (p->wtree, &(res->raddir2[id][is][js][0]), totweight[0]/p->dir.cotheta);
	}
      }
    }
  }
  else {          /* diffuse */
    
    if ((surfaceparallel && incoming) || (!surfaceparallel && p->dir.dx[2]<0)) {  /* down-welling */
      res->ndn[is][js] ++;
      res->edn[is][js] += totweight[0];
      if (sample->std) 
        p->wtree = addtree_stddev (p->wtree, &(res->edn2[is][js]), totweight[0]);
        

      /* ??? Attention, this is only to reduce noise ??? CHECK!!! */ 
      /* ??? but it actually introduces systematic   ??? */
      /* ??? error! Tests show that we loose up to   ??? */
      /* ??? 1 - 2% of the diffuse components when   ??? */
      /* ??? ACTINIC_CUTOFF is set to 0.01           ??? */
      /* ???                                         ??? */
      /* ??? We correct that partly by multiplying   ??? */
      /* ??? the result by 1/(1-ACTINIC_CUTOFF),     ??? */
      /* ??? assuming that the actinic flux in the   ??? */
      /* ??? cutoff range equals the average actinic ??? */
      /* ??? flux in the hemisphere                  ??? */

      if (fabs(p->dir.cotheta) >= ACTINIC_CUTOFF) {
        res->fdn[is][js] += totweight[0] / p->dir.cotheta / (1.0-ACTINIC_CUTOFF);
        
        if (sample->std)
          p->wtree = addtree_stddev (p->wtree, &(res->fdn2[is][js]), 
                                     totweight[0]/p->dir.cotheta/(1.0-ACTINIC_CUTOFF));
      }
    }
    
    if ((surfaceparallel && !incoming) || (!surfaceparallel && p->dir.dx[2]>0)) {  /* up-welling */
      res->nup[is][js] ++;
      res->eup[is][js] += totweight[0];
      if (sample->std) 
        p->wtree = addtree_stddev (p->wtree, &(res->eup2[is][js]), totweight[0]);

      /* ??? Attention, this is only to reduce noise ??? CHECK!!! */ 
      /* ??? but it actually introduces systematic   ??? */
      /* ??? error! Tests show that we loose up to   ??? */
      /* ??? 1 - 2% of the diffuse components when   ??? */
      /* ??? ACTINIC_CUTOFF is set to 0.01           ??? */
      /* ???                                         ??? */
      /* ??? We correct that partly by multiplying   ??? */
      /* ??? the result by 1/(1-ACTINIC_CUTOFF),     ??? */
      /* ??? assuming that the actinic flux in the   ??? */
      /* ??? cutoff range equals the average actinic ??? */
      /* ??? flux in the hemisphere                  ??? */
      
      if (fabs(p->dir.cotheta) >= ACTINIC_CUTOFF) {
        res->fup[is][js] += totweight[0] / p->dir.cotheta / (1.0-ACTINIC_CUTOFF);
        if (sample->std)
          p->wtree = addtree_stddev (p->wtree, &(res->fup2[is][js]), 
                                     totweight[0]/p->dir.cotheta/(1.0-ACTINIC_CUTOFF));
      }
    }
    
    /* diffuse radiance */
    if (sample->panorama_forward) {
      double theta = acosd(p->dir.dx[2]);

      double phi = 0;
      if (p->dir.dx[0] == 0 && p->dir.dx[1] ==0)
	phi=0;
      else 
	phi=atan2 (p->dir.dx[0],p->dir.dx[1])/M_PI*180.0;

      while (phi<0)
	phi = 360.0 + phi;
      is = (int) ((phi)   / 360.0 * sample->Nx);
      js = (int) ((theta) / 180.0 * sample->Ny);

      /* collect edges of forward panorama in last bin (LF 2016/08/19) */
      if (is==sample->Nx)
          is-=is;
      if (js==sample->Ny)
          js-=js;

      for (ip=0; ip<sample->nstokes; ip++){
	res->raddif[0][is][js][ip] += totweight[ip] / p->dir.cotheta;
	if (sample->std)
	  p->wtree = addtree_stddev (p->wtree, &(res->raddif2[0][is][js][ip]),
				     totweight[ip]/p->dir.cotheta);
      }
    }
    else {
      for (id=0; id<sample->Nd; id++) {
	if (incone(p->dir.dx, sample->rad[id])) {
	  for (ip=0; ip<sample->nstokes; ip++){
	    res->raddif[id][is][js][ip] += totweight[ip] / p->dir.cotheta;
	    if (sample->std)
	      p->wtree = addtree_stddev (p->wtree, &(res->raddif2[id][is][js][ip]), 
					 totweight[ip]/p->dir.cotheta);
	  }
	}
      }
    }

    /* radial pathlength radiance */
    if (sample->Nr>0 && sample->Nt>0) 
      if (sample_radius (p, sample, &ir)) 
        if (sample_time (p, sample, &it)) 
          for (id=0; id<sample->Nd; id++)
            if (incone(p->dir.dx, sample->rad[id])) {
              for (ip=0; ip<sample->nstokes; ip++){
                res->radpat[id][ir][it][ip] += totweight[ip] / p->dir.cotheta;
                if (sample->std)
                  p->wtree = addtree_stddev (p->wtree, &(res->radpat2[id][ir][it][ip]), 
                                             totweight[ip]/p->dir.cotheta);
              }
            }
  }
}


/*TZ bt ****************************************************************************/
/* Function: count_thermal_backward_photon                                @62_30i@ */
/* Description:                                                                    */
/*  Register a thermal backward photon; that is, add its contribution to the       */
/*  irradiance, actinic flux, and radiance arrays in result_struct *res.           */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static void count_thermal_backward_photon (result_struct *result,
                                           photon_struct *p, 
                                           atmosphere_struct *atmos,
                                           surftemp_struct *surftemp,
                                           sample_struct *sample, 
                                           double btemp_plkavg,
                                           double albedo,
                                           int quiet, int toa_surf_flag)
{
  double *add=NULL;
  double *totweight=NULL;

  double *totweight_emis=NULL;  /* **CK  2013.08.27 Include totweight_emis*/
  double *add_emis=NULL; /* **CK  2013.08.27 Include add_emis*/
  
int it=0, jt=0, ip=0;

#if HAVE_ALIS
  int iv=0, ic=0;
  double* totweight_spectral=NULL; 
  double* totweight_concentration=NULL; 
  double factor=0.0;
  double planck_ref=0.0;
  double planck_iv=0.0;
  double planck_ratio=1.0;
#endif

  totweight=calloc(sample->nstokes, sizeof(double));
  add=calloc(sample->nstokes, sizeof(double));

  add_emis=calloc(sample->nstokes, sizeof(double));/* **CK  2013.08.27 Allocate add_emis*/
  totweight_emis=calloc(sample->nstokes, sizeof(double));/* **CK  2013.08.27 Allocate totweight_emis*/
  

#if HAVE_ALIS
  if (sample->spectral_is||sample->concentration_is){
    totweight_spectral=calloc(atmos->nlambda_abs, sizeof(double));
    totweight_concentration=calloc(atmos->Nc, sizeof(double));
    totweight_spectral[0]=1.0;
    totweight_concentration[0]=1.0;
    if (sample->spectral_is)
      spectral_is_weight(&totweight_spectral, p, p, atmos); /* RpA???  CHECK!!! */
    if (sample->concentration_is)
      concentration_is_weight(&totweight_concentration, p, p, atmos);
  }
#endif


  for (ip=0; ip< sample->nstokes; ip++){

/* **CK  2013.08.27 */
    if (sample->backward != MCBACKWARD_HEAT) 
      totweight[ip] = p->phamat[ip][0] * p->weight * exp(-p->tauabs.tot);
    else {
      switch (sample->heat_flag[sample->backward_is][sample->backward_js]) {
      case MCBACKWARD_HEAT_EMABS:
      case MCBACKWARD_HEAT_EMABSOPT:
	totweight[ip] = p->phamat[ip][0] * p->weight * exp(-p->tauabs.tot);
	totweight_emis[ip] = p->phamat[ip][0] * p->weight_emis;  /* **CK 2013.08.27 p->weight_emis only set for new photon distribution for thermal backward heating rates EMABS_Method*/
	break;
      case MCBACKWARD_HEAT_DENET:
	totweight[ip] = p->phamat[ip][0] * p->weight * exp(-p->tauabs.tot);
	break;
      } 
    }


    switch (toa_surf_flag) {
    case MCSTATUS_ABSORB:
      /*linear interpolation of Bplanck for comparison with DISORT*/
      totweight[ip] *= ( atmos->Bplanck->prof  [MCCAOTH_TOT][p->kc][p->ic][p->jc] + 
                         ( ( ( p->x[2] - atmos->Z[p->kc] ) / ( atmos->Z[p->kc+1] - atmos->Z[p->kc] ) ) * 
                           ( atmos->Bplanck->prof  [MCCAOTH_TOT][p->kc+1][p->ic][p->jc]
                             - atmos->Bplanck->prof  [MCCAOTH_TOT][p->kc][p->ic][p->jc] ) ) );
      break;
      
    case MCSTATUS_BOUNDARY_UPPER:
      /* TOA exit => zero contribution and not Bplanck(top level) contribution */
      totweight[ip] = 0.0;
      break;
      
    case MCSTATUS_BOUNDARY_LOWER:
      /* ??? BM:don't think so. We don't need to weight with the emissivity CHECK!!! */ 
      /* ??? because this is already done by reflecting the photon or not    */
      /* ??? reflecting it. When multiplied by (1.0-albedo) the upward       */
      /* ??? irradiance at the surface was divided by 4 if the albedo was    */
      /* ??? set to 0.5!                                                     */

      /* surface absorption => Bplanck(surface temp)*/
      /* totweight *= btemp_plkavg*(1.0-albedo);    */
      
      if (surftemp->surf2D) {
        /* 2D surface temperature */
        surftemp_coord (p, surftemp, &it, &jt);
        totweight[ip] *= surftemp->plkavg2D[it][jt];
      }
      else
        totweight[ip] *= btemp_plkavg;
      
      break;
      
    case MCSTATUS_SURFACE:
      if (surftemp->surf2D) {
        /* 2D surface temperature */
        surftemp_coord (p, surftemp, &it, &jt);
        totweight[ip] *= surftemp->plkavg2D[it][jt];
      }
      else {
        /* if surface temperature not explicitely defined by user      */
	/* interpolate atmospheric planck function to actual location  */
        /* assuming that the surface temperature equals the air        */
        /* temperature at this height.                                 */
	
	if (surftemp->user_defined_btemp)
          totweight[ip] *= btemp_plkavg;
        
        else 
          totweight[ip] *= atmos->Bplanck->prof [MCCAOTH_TOT][p->kc][p->ic][p->jc] + 
            ((( p->x[2] - atmos->Z[p->kc] ) / ( atmos->Z[p->kc+1] - atmos->Z[p->kc] )) * 
             ( atmos->Bplanck->prof [MCCAOTH_TOT][p->kc+1][p->ic][p->jc]
               - atmos->Bplanck->prof [MCCAOTH_TOT][p->kc][p->ic][p->jc]));
      }
      
      break;
      
    default:
      fprintf (stderr, "Error, unknown toa_surf_flag %d\n", toa_surf_flag);
      totweight[0]=0.0;
    }
    
    switch (sample->backward) {
    case MCBACKWARD_EGLOB:
    case MCBACKWARD_EDN:
    case MCBACKWARD_EUP:
    case MCBACKWARD_EXP:
    case MCBACKWARD_EXN:
    case MCBACKWARD_EYP:
    case MCBACKWARD_EYN:
    case MCBACKWARD_FDN:
    case MCBACKWARD_FUP:
    case MCBACKWARD_ABS:
    case MCBACKWARD_ACT:
    case MCBACKWARD_HEAT:   /* emission is counted separately in sample->backemis, therefore no difference to absorption */
      if (sample->backward != MCBACKWARD_HEAT) 
	add[ip] = totweight[ip] * 4.0 * PI; 
      else {
	switch (sample->heat_flag[sample->backward_is][sample->backward_js]) {  /* **CK  2013.08.27 */
	case MCBACKWARD_HEAT_EMABS:
	case MCBACKWARD_HEAT_EMABSOPT:
	  add[ip] = totweight[ip] * 4.0 * PI; 
	  add_emis[ip]=totweight_emis[ip] * 4.0 * PI;
	  break;
	case MCBACKWARD_HEAT_DENET:
	  add[ip] = totweight[ip] * 4.0 * PI; 
	  break;
	}
      }
    
      break;
      
      
    case MCBACKWARD_RADIANCE:
      add[ip]=totweight[ip];
      break;

    case MCBACKWARD_EMIS:
      fprintf (stderr, "\n");
      fprintf (stderr, "*** Warning, with MCBACKWARD_EMIS we are not supposed to end\n");
      fprintf (stderr, "*** in count_thermal_backward_photon()\n");
      fprintf (stderr, "\n");
      free(add);
      free(totweight);
      free(add_emis);  /* **CK  2013.08.27 */
      free(totweight_emis);  /* **CK  2013.08.27 */
      return;
      break;
      
    default:
      fprintf (stderr, "Error, unknown sample->backward %d in count_thermal_backward_photon()\n", sample->backward);
      add[ip]=0;
    }
    
    /* If backward, then relate the escape result to a certain xy index */
    if (sample->backward != MCBACKWARD_HEAT) 
      result->back[sample->backward_is][sample->backward_js][ip] += add[ip];
    else {
      /* **CK  2013.08.27 */
      switch (sample->heat_flag[sample->backward_is][sample->backward_js]) {
      case MCBACKWARD_HEAT_EMABS:
      case MCBACKWARD_HEAT_EMABSOPT:
	result->backemis[sample->backward_is][sample->backward_js] += add_emis[ip];
	result->back[sample->backward_is][sample->backward_js][ip] += add[ip];
	break;
      case MCBACKWARD_HEAT_DENET:
	result->back[sample->backward_is][sample->backward_js][ip] += add[ip];
	break;
      }
    }
    
#if HAVE_ALIS
      if (sample->spectral_is || sample->concentration_is){
        factor=( p->x[2] - atmos->Z[p->kc] ) / ( atmos->Z[p->kc+1] - atmos->Z[p->kc] );
        planck_ref=atmos->Bplanck->prof  [MCCAOTH_TOT][p->kc][p->ic][p->jc] + 
          factor * ( atmos->Bplanck->prof  [MCCAOTH_TOT][p->kc+1][p->ic][p->jc]
                     - atmos->Bplanck->prof  [MCCAOTH_TOT][p->kc][p->ic][p->jc] );
        
        for (iv=0; iv<atmos->nlambda_abs; iv++){
          planck_iv= atmos->Bplanck_spectral[iv][p->kc] + 
            factor * (atmos->Bplanck_spectral[iv][p->kc+1] 
                      - atmos->Bplanck_spectral[iv][p->kc]) ;

          if (sample->spectral_is)
            planck_ratio=planck_iv / planck_ref;
          
          for (ic=0; ic<atmos->Nc; ic++){
            result->back_t[ic][sample->backward_is][sample->backward_js][ip][iv] +=
              add[ip]* totweight_spectral[iv] * totweight_concentration[ic] 
              * planck_ratio;
            /*   fprintf(stderr, "FIXCE ic %d  specweight %g concweight %g back_t  %g add[ip] %g planck %g\n",ic, totweight_spectral[iv], totweight_concentration[ic],  result->back_t[ic][iv][sample->backward_is][sample->backward_js][ip], add[ip], planck_ratio ); CHECK!!! */
          }
        }
      }
#endif
      
 
      if (sample->std){  /* **CK  2013.08.27*/
	if (sample->backward != MCBACKWARD_HEAT)
	  p->wtree = addtree_stddev (p->wtree, &(result->back2[sample->backward_is][sample->backward_js][ip]), add[ip]);  
	else {
	  switch (sample->heat_flag[sample->backward_is][sample->backward_js]) {
	  case MCBACKWARD_HEAT_EMABS:
	  case MCBACKWARD_HEAT_EMABSOPT:
	    p->wtree = addtree_stddev (p->wtree, &(result->back2[sample->backward_is][sample->backward_js][ip]), add[ip]+add_emis[ip]);
	    break;
	  case MCBACKWARD_HEAT_DENET:
	    result->back_dEnet[sample->backward_is][sample->backward_js][p->dEnet_component] += add[ip]/sqrt(p->n_dEnet);  
	    result->back_dEnet2[sample->backward_is][sample->backward_js][p->dEnet_component] += add[ip]*add[ip]; 
	    break;
	  }
	}

      }
 
    }
  //}

#if HAVE_ALIS
  free(totweight_spectral);
  free(totweight_concentration); 
#endif

  free(add);
  free(totweight);
  free(add_emis);  /* **CK  2013.08.27 */
  free(totweight_emis);  /* **CK  2013.08.27*/

}


/***********************************************************************************/
/* Function: intersection1D                                               @62_30i@ */
/* Description:                                                                    */
/*  Calculate intersection with a 1D layer boundary.                               */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void intersection1D (photon_struct *p,
		     atmosphere_struct *atmos,
		     double tau, double tausca,
		     double *step, int source)
{
  double alpha=0;
  
  if (p->dir.hit[2]>=0)   /* upward or downward */
    alpha = ( atmos->Z[p->kc+p->dir.hit[2]] - p->x[2] ) / p->dir.dx[2];
  else {                   /* horizontal         */
    /* RPB: the extra MC_EPSILON is to prevent the photon    */
    /*      from being trapped. The photon overshoots and is */
    /*      corrected to the right position in travel_tau    */
    if (source != MCSRC_THERMAL_BACKWARD)
      alpha = ( tau - tausca + MC_EPSILON ) / get_kscaIS (atmos, p, MCCAOTH_TOT);
    else
      alpha = ( tau - tausca + MC_EPSILON ) / ( get_kscaIS (atmos, p, MCCAOTH_TOT)
					      + get_kabs_tot (atmos, p) );
  }

  *step = alpha;
}


/***********************************************************************************/
/* Function: intersection1D_spherical                                     @62_30i@ */
/* Description:                                                                    */
/*  Calculate intersection with a 1D spherical layer boundary.                      */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Claudia Emde                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int intersection1D_spherical (photon_struct *p,
			      atmosphere_struct *atmos,
			      int refraction, 
			      float *refind, double *step)
{
  int i=0;
  double alpha=-1.;
  double x=0., y=0., z=0., r=0., r_sq=0., r1=0., r2=0.;
  double costheta=0., sintheta=0.; 
  double rdr=0.;
  double r2diff=1.0;
  double r1diff=1.0;
  double N=1.0, b1=0.0, a2=0.0, b2=0.0, fac=0.0;
  double e_r[3], e_t[3];
  
  x=p->x[0]-atmos->xmax/2.;
  y=p->x[1]-atmos->ymax/2.;
  z=p->x[2]+atmos->r_earth;

  /* r of layers below and above */
  r2=atmos->Z[p->kc+1]+atmos->r_earth;
  r1=atmos->Z[p->kc]+atmos->r_earth;

  r_sq=x*x+y*y+z*z;
  r=sqrt(r_sq);

  if (z<atmos->r_earth*1e-5) {
    fprintf (stderr, "\n");
    fprintf(stderr, "*** Warning: Photon travelled around the Earth and is absorbed.\n");
    fprintf (stderr, "\n");
    
    p->weight = 0.0; 
    alpha = 0.;
    return -1;
  }
  
  /* Normalize direction vector */
  e_r[0]=x/r;
  e_r[1]=y/r;
  e_r[2]=z/r;
  
  rdr=0.0;
  for (i=0; i<3; i++)
    rdr+=e_r[i]*p->dir.dx[i];
  
  /* Angle between dr and r */
  costheta=rdr;
  sintheta=sqrt(1.-costheta*costheta);
  
  /* fprintf(stderr, "=====================================\n"); */
/*   fprintf(stderr, "refind kc+1 %g  kc %g \n", refind[p->kc+1], refind[p->kc]); */
/*   fprintf(stderr, "before refrac rdr %g dx %g dy %g dz %g costheta %g sintheta %g \n", */
/*           rdr, p->dir.dx[0], p->dir.dx[1], p->dir.dx[2], costheta, sintheta); */
  
/*   fprintf(stderr, "p %g \n", (1+refind[p->kc])*r*sintheta); */
  
  if (refraction){
  
    /* if direction not parallel to radial or tangential direction */
    if (rdr!=1.0 && sintheta < 0.999 ){
      
      if (costheta >= 0.0)/*upward*/
        fac=(1+refind[p->kc])/(1+refind[p->kc+1]);
      else{/* downward and limb*/
        fac=(1+refind[p->kc+1])/(1+refind[p->kc]);
        for (i=0; i<3; i++)
          e_r[i]*=-1;
        rdr*=-1;
      }
      
      /* Calculate tangential vector */
      N=0; 
      for (i=0; i<3; i++){
        e_t[i]=rdr*e_r[i]-p->dir.dx[i];
        N+=e_t[i]*e_t[i]; 
      }
      
      for (i=0; i<3; i++)
        e_t[i]/=sqrt(N); 
      
      /*    fprintf(stderr, "et %g %g %g \n", e_t[0], e_t[1], e_t[2]);  */
      /*       fprintf(stderr, "er %g %g %g \n", e_r[0], e_r[1], e_r[2]); */
      
      b1=0; 
      for (i=0; i<3; i++)
        b1+=p->dir.dx[i]*e_t[i];
      
      
      b2=fac*b1; 
      a2=sqrt(1.-b2*b2);
      
      /*       fprintf(stderr, "fac %g b1 %g a2 %g b2 %g \n",fac, b1, a2, b2);  */
      for (i=0; i<3; i++)
        p->dir.dx[i]=a2*e_r[i]+b2*e_t[i]; 
      
      /* Recalculate rdr  */
      rdr=0.0;
      for (i=0; i<3; i++)
        rdr+=e_r[i]*p->dir.dx[i];
      
      if (!(costheta >= 0.0)){
        for (i=0;i<3; i++)
          e_r[i]*=-1;
        rdr*=-1;
      }
      
      /* Special case limb ??? CHECK!!! */ 
      if (costheta/rdr <0.0)
        rdr*=-1;
      
      /* Angle between dr and r */
      costheta=rdr;
      sintheta=sqrt(1.-costheta*costheta);
    }
/*       fprintf(stderr, "after refrac rdr %g dx %g dy %g dz %g costheta %g sintheta %g \n", */
/*               rdr, p->dir.dx[0], p->dir.dx[1], p->dir.dx[2], costheta, sintheta); */
/*       fprintf(stderr, "p %g \n", (1+refind[p->kc])*r*sintheta); */
      
      /* Check for normalisation */
      /* N=sqrt(p->dir.dx[0]*p->dir.dx[0] + p->dir.dx[1]*p->dir.dx[1] + p->dir.dx[2]*p->dir.dx[2]); */
/*       fprintf(stderr, "normalisation %g \n \n", N); */
    
  } /* end refraction */
  
  rdr*=r;

  /* upward and limb */
  if (costheta >= 0. || sintheta >= r1/r){
    /* Numerical problems likely to occur here, especially if atmosphere */
    /* has very thin layers, r_sq and r1*r1 are very close, so this threshold has been included. */
    /* 0.01 is a very small value compared to r_sq, which is something like 6370e3*6370e3.*/
    r2diff=r_sq-r2*r2;
    if (fabs(r2diff) > 0.01)
      alpha = r2diff/(-rdr - sqrt( rdr*rdr - r2diff));
    /* special case, intersection with the same layer */
    else
      alpha = -2.*r*costheta;
    
    /* Direction forced to upward for "limb" cases.*/
    p->dir.hit[2]=1;
  }
  
  /* downward 
     else if (costheta < 0. && sintheta <r1/r)*/
  else{
    
    /* Numerical problems likely to occur here, see above */
    r1diff=r_sq-r1*r1;
    if (fabs(r2diff) > 0.01){
      alpha = r1diff/(-rdr + sqrt( rdr*rdr - r1diff));
      p->dir.hit[2]=0;
    }
    
    else{
      /* Photon is absorberd. This can only happen due to numerical 
         problems. */ 
      p->weight = 0.0; 
      alpha = 0.0;
      
      /*Some debugging output:*/
      fprintf( stderr, "Photon absorbed by surface!!\n");  
      fprintf( stderr, "direction: %g \n", acosd(p->dir.dx[2])); 
      fprintf( stderr, "position:  z_p  %g z_i %g z_i+1 %g, r %g \n", p->x[2],atmos->Z[p->kc], atmos->Z[p->kc+1], r-atmos->r_earth);
      fprintf( stderr, "x %f, y %f, z %f, p->x[0] %f, p->x[1], %f \n", x, y, z, p->x[0],  p->x[1] );
      fprintf( stderr, " r1/r %g costheta %g  sintheta %g \n \n" , r1/r, costheta, sintheta); 
      return -1; 
    }
  }
  *step = alpha;
  
  if (alpha < 0.0  ||  alpha != alpha){
    /*Some debugging output:*/
    fprintf( stderr, "direction: %f \n", acosd(p->dir.dx[2])); 
    fprintf( stderr, "position:  z_p  %f z_i %f z_i+1 %f, r %f \n", p->x[2],atmos->Z[p->kc], atmos->Z[p->kc+1], r-atmos->r_earth);
    fprintf( stderr, "x %f, y %f, z %f, p->x[0] %f, p->x[1], %f \n", x, y, z, p->x[0],  p->x[1] );
    fprintf( stderr, " r1/r %g costheta %g  sintheta %g \n" , r1/r, costheta, sintheta) ; 
    fprintf( stderr, " alpha %g  rdr*rdr %g r1diff %g r2diff %g r_sq %g r2 %g \n", alpha,  rdr*rdr, r1diff, r2diff, r_sq, r2 ); 
    fprintf(stderr, "UUPS!!!: Steplength less or equal to zero: %.6e. If this does not happen\n", alpha);
    fprintf(stderr, "         too often, you can probably still trust your result.\n");
    p->weight = 0.0;
    alpha = 0.0;
    return -1;  
  }
  return 0; 
}


/***********************************************************************************/
/* Function: step1D                                                       @62_30i@ */
/* Description:                                                                    */
/*  Move photon to the next layer boundary or to the next scattering point and     */
/*  and update layer number kc.                                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int step1D (photon_struct *p, atmosphere_struct *atmos, double step, 
	    int bcond, int photonpath, int visualize)
{  
  int periodic=0;

  #if HAVE_OPENGL
  double x=0, y=0, z=0;

  /* save start position and direction */
  if (visualize) {
    x = p->x[0];
    y = p->x[1];
    z = p->x[2];
  }
  #endif
  
  /* new position */
  if (!p->ipa) {
    p->x[0] += step * p->dir.dx[0];
    p->x[1] += step * p->dir.dx[1];
  }
  
  p->x[2] += step * p->dir.dx[2];
  
  switch (bcond) {
  case MCBCOND_PERIODIC:
    periodic = perform_periodic_bc (p, atmos);

#if HAVE_OPENGL
    if (visualize)
      if (periodic && photonpath)
        GLmystic_add_periodic_nodes_to_photon_path (p, atmos, x, y, z, step);
#endif

    break;
    
  case MCBCOND_ABSORB:
    if (p->x[0]<0 || p->x[0]>atmos->xmax || p->x[1]<0 || p->x[1]>atmos->ymax) {
      p->weight = 0.0;
      p->photon_status = MCSTATUS_OUTOFDOMAIN;
      return 0;
    }
    break;
    
  case MCBCOND_MIRROR: 
    fprintf (stderr, "Error, mirroring boundary conditions not implemented!\n");
    return -1;
    
  default:
    fprintf (stderr, "Error, unknown boundary conditions %d\n", bcond);
    return -1;
  }
  
  /* determine indices of current box */
  atmos_coord (p, atmos, &(p->ic), &(p->jc));
  
  if (p->photon_status == MCSTATUS_TRAVEL) {
    if (p->dir.hit[2]==0)   /* downward */          
      (p->kc)--;
    
    if (p->dir.hit[2]==1)   /* upward */
      (p->kc)++;
  }

  return periodic;
}


/***********************************************************************************/
/* Function: perform_periodic_bc                                          @62_30i@ */
/* Description:                                                                    */
/*  Check whether photon is outside of grid, and apply periodic bc if necessary.   */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static inline int perform_periodic_bc (photon_struct *p, atmosphere_struct *atmos)
{
  int periodic=0;

  if (p->dir.hit[0]) {
    /* photon moving in +x direction */
    while(p->x[0]<0) {
      p->x[0] += atmos->xmax;
      periodic=1;
    }
    while (p->x[0] >= atmos->xmax) {
      p->x[0] -= atmos->xmax;
      periodic=1;
    }
  }
  else {
    /* photon moving in -x direction */
    while(p->x[0]<=0) {
      p->x[0] += atmos->xmax;
      periodic=1;
    }
    while (p->x[0] > atmos->xmax) {
      p->x[0] -= atmos->xmax;
      periodic=1;
    }
  }

  if (p->dir.hit[1]) {
    /* photon moving in +y direction */
    while (p->x[1] < 0) {
      p->x[1] += atmos->ymax;
      periodic=1;
    }
    while (p->x[1] >= atmos->ymax) {
      p->x[1] -= atmos->ymax;
      periodic=1;
    }
  }
  else {
    /* photon moving in -y direction */
    while (p->x[1] <= 0) {
      p->x[1] += atmos->ymax;
      periodic=1;
    }
    while (p->x[1] > atmos->ymax) {
      p->x[1] -= atmos->ymax;
      periodic=1;
    }
  }

  return periodic;
}


/***********************************************************************************/
/* Function: calloc_photon                                                @62_30i@ */
/* Description:                                                                    */
/*  Allocate memory for struct photon_struct.                                      */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/
/* static */
photon_struct *calloc_photon ( sample_struct *sample,
			       int            Nz,
			       int            source,
			       int            nlambda_abs,
                               int            Nc, 
			       int            n_caoth )
{
  int i=0;

  photon_struct *p = calloc (1, sizeof(photon_struct));
    
  /* allocate memory for radiance probabilities */
  if (sample->Nd>0)
    p->pdir = calloc ((size_t) sample->Nd, sizeof(double));

  /* allocate memory for weight vector */
  p->stokes0 = calloc (4, sizeof(double));
  p->stokes = calloc (4, sizeof(double));
  
#if HAVE_LIDAR
  calloc_photon_lidar_part (p, sample, sample->Nz_jac, source, n_caoth);
#endif
  
  p->phamat=calloc( 4, sizeof(double *));
  for (i=0; i<4; i++)
    p->phamat[i]=calloc( 4, sizeof(double));

  /* spectral absorption */
  if(sample->spectral_is||sample->concentration_is){
    /* p->nlambda=nlambda_abs; */
    /* p->dtauabs_spectral = calloc (nlambda_abs, sizeof(double)); */
    p->Nz_alis = Nz;
    p->pathlength_per_layer = calloc (Nz, sizeof(double));
    p->q_spectral= calloc (nlambda_abs, sizeof(double));
    p->q2_spectral= calloc (nlambda_abs, sizeof(double));
    p->q_albedo_spectral= calloc (nlambda_abs, sizeof(double));
    p->q_concentration = calloc (Nc, sizeof(double));
    p->q2_concentration = calloc (Nc, sizeof(double));
  }
  
  if(sample->boxairmass){
    p->pathlength_per_layer = calloc (Nz, sizeof(double));
    p->Nz_alis = Nz;
  }

  if (sample->ncirc)
    p->tocirco=calloc(100, sizeof(double *));
  
  return p;
}


/***********************************************************************************/
/* Function: gen_default_sunshape                                         @62_30i@ */
/* Description:                                                                    */
/*  Generate a wavelength dependend sunshape according to Köpke 2001               */
/*                                                                                 */
/* Parameters: wvnmlo, wvnmhi: Min and max wavenumbers          		   */
/*             N: number of requested supporting points for sunshape		   */
/*             pd: (output) sunhape						   */
/*             alpha: (output) relative position on sundisc (0-1)		   */
/* Return value									   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Bernhard Reinhardt                                                      */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void gen_default_sunshape(float wvnmlo, float wvnmhi, int N, float *pd, double *alpha)
{
  /* Definition of constants:*/ 
  double planck=6.63e-34;
  double speed_of_light=2.998e8;
  double boltzmann=1.38e-23;
  double T_S=5740.0;
  
  /* Wavelength */
  double lambda=2.0/((wvnmlo+wvnmhi)*100.0);

  /* Limb darkening coefficient*/
  double beta=3.0*planck*speed_of_light*sqrt(sqrt(2.0))/
    (8.0*boltzmann* lambda * T_S); 

  int i;
  
  /*pd = (float*) calloc(N+1,sizeof(float)); */
  /*alpha = (double*) calloc(N+1,sizeof(double)); */

  alpha[0]=0.;
  pd[0]=1.;
  
  for (i=1; i<N; i++){
      alpha[i] = alpha[i-1]+1.0/(N-1);
      pd[i] = (1.0+beta*sqrt(fabs(1.0-alpha[i]*alpha[i])))/(1+beta);
  }     
}


/***********************************************************************************/
/* Function: limb_dark_lest_dir                                           @62_30i@ */
/* Description:                                                                    */
/*  Generate a photon escape direction according to a limb darkening function      */
/*  of the sun. I.e. choose a direction pointing to a spot in the sun disk         */
/*                                                                                 */
/* Parameters: "dir_sun_center" is modified according to sunshape-information      */
/*             "sample" contains information about sunshape                        */
/*             "phi" and "alpha" are filled with direction information about       */
/*             local estimate direction in regard to sun center                    */
/* Return value: 0 if ok, <=-1 if error                                            */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Bernhard Reinhardt                                                      */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int limb_dark_lest_dir (direction *dir_sun_center, sample_struct *sample,
			double *phi, double *alpha)
{

  int j = 0;
  double m = 0.0; /*  slope of linear interpolation */
  double random;
 
  /* Sloppy check of sun shape  */
  if (sample->sample_backward_sunshape_n_lines < 3) {
    fprintf ( stderr, "%s %s %s",
	      "sun shape file must contain more than one line and must be in correct format:",
	      " 1. Col: probability, 2. Col.: Normalized angular distance to center of sun:",
	      " 0-> center of sun, 1->limb of sun\n");
    return -1;
  }

  random = uvspec_random();

  j = locate ( sample->sample_backward_sunshape_F, sample->sample_backward_sunshape_n_lines, random );

  m = ( sample->sample_backward_sunshape_alpha[j+1] - sample->sample_backward_sunshape_alpha[j] ) /
    ( sample->sample_backward_sunshape_F[j+1] - sample->sample_backward_sunshape_F[j] );

  *alpha = sample->sun_radius * ( sample->sample_backward_sunshape_alpha[j]
				  + m * ( random - sample->sample_backward_sunshape_F[j] ) );

  *phi = sc_Isotropic_phi();

  new_direction( cosd(*alpha), *phi, dir_sun_center, 0.);

  return 0;
}


/***********************************************************************************/
/* Function: generate_photon                                              @62_30i@ */
/* Description:                                                                    */
/*  Generate a photon and inititalize its location and direction according to      */
/*  the source.                                                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static photon_struct *generate_photon ( int                source,
					int                photoncounter,
					atmosphere_struct *atmos,
					albedo_struct     *albedo,
					elevation_struct  *elev,
					sample_struct     *sample,
					float              wvnmlo,
					float              wvnmhi,
					float              wavelength,
					double             sza,
					double             phi,
					int                ipa,
					long int           nphotons)   /* **CK 2013.08.27 Include nphotons */
{

/* ============================================================================== */
/* ==================== Declarations and Preparations =========================== */
/* ============================================================================== */
  double B=0, emis=0, phase=0.0, tmpalb=0.0;
  int ia=0, ja=0;
  double* norm = NULL;
  int i=0, j=0, k=0, ip=0, iv=0;
  int hit=0;
  
#if HAVE_LIDAR
  int status=0;
#endif
  
  static int counter1=0; /* **CK 2013.08.27 */
 
  norm = calloc(3,sizeof(double));
  double ran, apar, aperp, dd;

  photon_struct *p = calloc_photon ( sample,
				     atmos->Nz,
				     source,
				     atmos->nlambda_abs,
                                     atmos->Nc,
				     atmos->n_caoth );

#ifdef MUCHOUT
  p->muchoutcounter=photoncounter;
  if (p->muchoutcounter==MUCHOUTPHOTON)
    fprintf(stderr,"Muchouting this photon: %d\n",p->muchoutcounter);
#endif

/* ============================================================================== */
/* ==================== Prepare Solar Eclipse =================================== */
/* ============================================================================== */
  const int sofi = 0;  /* Internal Switch? */
  static int first=1;
  static double *pd = NULL;
  if (sofi&&first) {  /* solar eclipse */
    if (sample->backward) {
      fprintf (stderr, "Error, sofi only works in forward mode!\n");
      return NULL;
    }

    first=0;
    pd = (double*) calloc (2000,sizeof(double));
    sample_photons_sofi (wvnmlo, wvnmhi, pd);
  }
/* ============================================================================== */
/* ==================== Set Basic Stuff ========================================= */
/* ============================================================================== */
  p->photoncounter=photoncounter;
  p->backward_is=sample->backward_is;
  p->backward_js=sample->backward_js;
  p->scattercounter=0;     /* number of scattering events */
  p->spikewarningcounter=0;     /* number of scattering events */
  p->rayleighcounter=0;    /* number of rayleigh scattering events */
  p->q_isoene=1.0;         /* weight of rayleigh scattering events */
  p->direct=1;             /* indicates direct of diffuse */
  p->pathlength=0;         /* photon pathlength           */
  p->weight=1.0;           /* photon weight for e.g. importance sampling */
  p->weight_emis=1.0;      /*  **CK 2013.08.27 photon weight for thermal backward emission - set to 1 here, changed later!  */
  p->phi0=phi;
  if (sample->backward) {
    p->fw_phi0=sample->forward_phi0;
    p->fw_phi=phi;
  }
  else {
    p->fw_phi0=phi;
    p->fw_phi=sample->forward_phi;
  }
  p->special_weight=1.0;   /* photon weight for sponti-clone; has nothing to do with spectral IS! */
  p->vis_beta=0.0;
  p->reflectcounter=0;
  p->tauext_tot=0.0;
  p->photon_status=MCSTATUS_TRAVEL;
  p->maxpathlength = 1e33; /* if photon travels further, kill it */

  /* Initially the polarization of the photon is defined in (x,y,z) coordinates, 
     later with respect to the respective scattering frames */
  p->dir_ort.dx[0]=0.0;
  p->dir_ort.dx[1]=0.0;
  p->dir_ort.dx[2]=1.0;

  /* Assign unity matrix to phase matrix */
  p->phamat[0][0]=1.0;
  p->phamat[1][1]=1.0;
  p->phamat[2][2]=1.0;
  p->phamat[3][3]=1.0;

  /* reset pathlengths */
  p->tauabs.mol = 0;
  p->tauabs.aer = 0;
  p->tauabs.cld = 0;
  p->tauabs.ice = 0;
  p->tauabs.tot = 0;
  p->tauris     = 0;
  p->tauext_tot = 0;
  
  /* reset reallyhit */
  p->reallyhit[0]=0;
  p->reallyhit[1]=0;
  p->reallyhit[2]=0;

  p->ncircos=0;
  p->isclone=0;
  p->clonescattercounter=0;
  p->escapescattercounter=0;
  p->update_atmos_coord=0;
  p->doling=0;

/* ============================================================================== */
/* ==================== Stokes Vector =========================================== */
/* ============================================================================== */
  p->stokes0[0]=1.0;       /* starting photon stokes vector intensity */
  for (i=1; i<4; i++)      /* starting photon weight vector required for polarization */ 
    p->stokes0[i]=0.0; 
  switch (sample->polarisation_state) {
    case 0:                /* In case of natural solar/thermal radiation w=(1,0,0,0) */
      break;
    case 1:
      p->stokes0[1]=1.0;   /* +Q polarisation */
      break;
    case 2:
      p->stokes0[2]=1.0;   /* +U polarisation */
      break;
    case 3:
      p->stokes0[3]=1.0;   /* +V polarisation */
      break;
    case -1:
      p->stokes0[1]=-1.0;  /* -Q polarisation */
      break;
    case -2:
      p->stokes0[2]=-1.0;  /* -U polarisation */
      break;
    case -3:
      p->stokes0[3]=-1.0;  /* -V polarisation */
      break;
    case 4:                /* randomized Stokes Vector, fulfills I^2 = Q^2 + U^2 + V^2 */
      ran = 2*PI*uvspec_random();
      apar = fabs(cos(ran));
      aperp = fabs(sin(ran));
      dd = 2*PI*uvspec_random();
      p->stokes0[1] = apar*apar - aperp*aperp;
      p->stokes0[2] = 2*apar*aperp*cos(dd);
      p->stokes0[3] = -2*apar*aperp*sin(dd);
      break;
    case 5:                /* Manually set initial Stokes Vector here! */
      p->stokes0[1] = sqrt(0.5); /* example values */
      p->stokes0[2] = 0.5;
      p->stokes0[3] = 0.5;
      break;
    default:
      fprintf (stderr,"ERROR, no valid polarisation mode in generate_photon()!\n");
      return NULL;
  }

#ifdef CLDPRP
/* ============================================================================== */
/* ==================== Cloud Prop ============================================== */
/* ============================================================================== */
  if (sample->cldprp){
    p->cldprp.reff_wc=0.0;
    p->cldprp.reff_ic=0.0;
    p->cldprp.rhit_wc=0.0;
    p->cldprp.rhit_ic=0.0;
    p->cldprp.tau_wc=0.0;
    p->cldprp.tau_ic=0.0;
    p->cldprp.dxlwc=0.0;
    p->cldprp.dylwc=0.0;
    p->cldprp.dzlwc=0.0;
    p->cldprp.dxiwc=0.0;
    p->cldprp.dyiwc=0.0;
    p->cldprp.dziwc=0.0;
  }
#endif
	
/* ============================================================================== */
/* ==================== Importance Sampling ===================================== */
/* ============================================================================== */
  p->RIS_mode=MCRIS_MODE_NORMAL; /* standard mode for (Real) Importance Sampling */
  p->VIS_mode=MCVIS_MODE_NORMAL; /* standard mode for Virtual Importance Sampling */
  if (sample->delta_scaling==0) { 
    p->SC_mode=MCSC_MODE_DELTA_SCALE; /* delta scaling mode for Scattering */
    p->DDIS_SC_mode=MCSC_MODE_DELTA_SCALE; /* delta scaling mode for DDIS */
  }
  else {
    p->SC_mode=MCSC_MODE_NORMAL; /* standard mode for Scattering */
    p->DDIS_SC_mode=MCSC_MODE_NORMAL; /* standard mode for DDIS */
  }

  p->iv_alis=-1;
  if (sample->spectral_is){
    p->iv_alis = atmos->ilambda_ref[0];
    p->nlambda = atmos->nlambda_abs; 
    for (iv=0; iv<atmos->nlambda_abs; iv++){
      p->q_spectral[iv] = 1.0;  /* weight for single scattering albedo */
      p->q2_spectral[iv] = 1.0; /* weight for phase function */
      p->q_albedo_spectral[iv] = 1.0; /* weigth for spectral albedo */
      /*  p->dtauabs_spectral[iv]=0.0; old method */
    }
    for (k=0; k<atmos->Nz; k++)
      p->pathlength_per_layer[k] = 0.0;
  }
  
  if (sample->concentration_is){
    p->Nc = atmos->Nc; 
    for (i=0; i<atmos->Nc; i++){
      p->q_concentration[i] = 1.0;  /* weight for single scattering albedo */
      p->q2_concentration[i] = 1.0; /* weight for phase function */
    }
    for (k=0; k<atmos->Nz; k++)
      p->pathlength_per_layer[k] = 0.0;
  }

/* ============================================================================== */
/* ==================== Miscellaneous =========================================== */
/* ============================================================================== */
  if(sample->boxairmass)
   for (k=0; k<atmos->Nz; k++)
     p->pathlength_per_layer[k] = 0.0; 

  /* independent pixel */
  if (sample->tipadir==3) /* should be TIPA_DIR3D */
    p->ipa = 0;
  else
    p->ipa = ipa;

  if (sample->LLE_moon) /* BCA clean up! */
    source = MCSRC_SOLAR;

    /* horizontal start position */
    /* within the pixel: x = i*dx .. (i+1)*dx; y = j*dy .. (j+1)*dy*/

/* ============================================================================== */
/* ========== Set Starting Position and Direction Depending on Source =========== */
/* ============================================================================== */
  switch (source) {
/* ============================================================================== */
/* ========== SOLAR, BLITZ, THERMAL_BACKWARD: Starting Position and Direction === */
/* ============================================================================== */
  case MCSRC_SOLAR:
  case MCSRC_BLITZ:
  case MCSRC_THERMAL_BACKWARD:
  
    if ((p = generate_photon_solar_thermal_backward (sample, 
                                                     atmos, 
                                                     p, 
                                                     elev,
                                                     &hit, 
                                                     nphotons, 
                                                     &counter1,
                                                     sza,
                                                     phi,
                                                     phase,
                                                     source,
                                                     sofi,
                                                     norm,
                                                     pd)) == NULL) {
      fprintf (stderr,"Error in generate_photon_solar_thermal_backward()!\n");
      return NULL;
    }
    break; 

/* ============================================================================== */
/* ========== THERMAL_SURFACE: Starting Position and Direction ================== */
/* ============================================================================== */
  case MCSRC_THERMAL_SURFACE:

    if (elev->elev2D && !sample->backward) {
      fprintf (stderr, "Error, topography and thermal emission only available with backward MYSTIC!\n");
      return NULL;
    }

    if (sample->backward) {
      /* ??? don't forget to consider sensorposition and sensordirection when implementing backward CHECK!!! */ 
      fprintf (stderr, "Error, backward not implemented for MCSRC_THERMAL_SURFACE\n");
      return NULL;
    }

    p->direct=0;   /* no direct beam source */
    
    /* horizontal start position:                */
    /* a rectangle: x = 0 .. xmax; y = 0 .. ymax */
    
    /* loop until photon is emitted */
    while (1) {
      p->x[0] = atmos->xmax * uvspec_random();
      p->x[1] = atmos->ymax * uvspec_random();

      /* albedo coordinates */
      if (albedo->method==MCALB_LAM2D || albedo->method==MCALB_LAM2D_SPECTRAL) {
	if (sample->spherical3D) {
#ifdef HAVE_SPHER
	  coord_spherical3D (p,
			     albedo->Nx, albedo->X,
			     albedo->Ny, albedo->Y,
			     0, NULL,
			     0, NULL,
			     0,
			     &ia, &ja, NULL, NULL,
			     0);
#else
	  fprintf(stderr,"Error! you are not allowed to use spherical 3D!\n");
	  return NULL;
#endif
	}
	else
	  albedo_coord (p, albedo, &ia, &ja);
      }


      switch (albedo->method) {
      case MCALB_LAM:
	tmpalb = albedo->albedo;
	break;
	
      case MCALB_LAM2D:
	tmpalb = albedo->albedo2D[ia][ja];
	break;

      case MCALB_LAM2D_SPECTRAL:
	tmpalb = albedo->alb_type[(int) albedo->rpv_index[ia][ja]];
	break;

      default:
	fprintf (stderr, "Error, albedo->method %d not yet implemented!\n", albedo->method);
      }
      
      if (emission(tmpalb))  
        break;
    }

    atmos_coord (p, atmos, &(p->ic), &(p->jc));
    
    /* initialize direction: isotropic upward */
    norm[0]=0.0;
    norm[1]=0.0;
    norm[2]=1.0;
    random_Lambertian_normal (&(p->dir), norm);

    /* vertical start position: surface */
    /* Make always sure to call set_photon_z() AFTER the new photon     */
    /* direction has been assigned because the start index and position */
    /* might depend on direction                                        */
    set_photon_z (atmos->Z[0], atmos, p);

    break;
    
/* ============================================================================== */
/* ========== THERMAL_ATMOSPHERE: Starting Position and Direction =============== */
/* ============================================================================== */
  case MCSRC_THERMAL_ATMOSPHERE:
    
    if (elev->elev2D && !sample->backward) {
      fprintf (stderr, "Error, topography and thermal emission only available with backward MYSTIC!\n");
      return NULL;
    }

    if (sample->backward) {
      /* ??? don't forget to consider sensorposition and sensordirection when implementing backward CHECK!!! */ 
      fprintf (stderr, "Error, backward not implemented for MCSRC_THERMAL_ATMOSPHERE\n");
      return NULL;
    }

    p->direct=0;   /* no direct beam source */
    
    /* horizontal start position:                */
    /* a rectangle: x = 0 .. xmax; y = 0 .. ymax */
      
    /* initialize direction: isotropic */
    random_direction (&(p->dir));

    /* loop until photon is emitted */
    while (1) {
      p->x[0] = atmos->xmax * uvspec_random();
      p->x[1] = atmos->ymax * uvspec_random();

      /* Make always sure to call set_photon_z() AFTER the new photon     */
      /* direction has been assigned because the start index and position */
      /* might depend on direction                                        */
      set_photon_z ( atmos->Z[0] + uvspec_random() * (atmos->Z[atmos->Nz] - atmos->Z[0]),
		     atmos, p );
      
      atmos_coord (p, atmos, &(p->ic), &(p->jc));

      /* interpolate Bplanck linearely between levels */
      B = atmos->Bplanck->prof [MCCAOTH_TOT][p->kc][p->ic][p->jc] + ( p->x[2] - atmos->Z[p->kc] )
	/ ( atmos->Z[p->kc+1] - atmos->Z[p->kc] ) * 
        ( atmos->Bplanck->prof [MCCAOTH_TOT][p->kc+1][p->ic][p->jc]
	  - atmos->Bplanck->prof [MCCAOTH_TOT][p->kc][p->ic][p->jc] );

      if (atmos->kabs3D->threed[MCCAOTH_TOT][p->kc]>=1)
        emis = B * atmos->kabs3D->prof [MCCAOTH_TOT][p->kc][p->ic][p->jc];
      else 
        emis = B * atmos->kabs->prof [MCCAOTH_TOT][p->kc];
      
      if (uvspec_random() <= emis)
        break;
    }

    break;
      
/* ============================================================================== */
/* ========== LIDAR: Set everything Lidar-related (Function in lidar.c) ========= */
/* ============================================================================== */
  case MCSRC_LIDAR:
#if HAVE_LIDAR
    status = generate_photon_lidar (p, sample, atmos);
    if (status!=0) {
      fprintf (stderr,"Error %d returned by generate_photon_lidar_part()\n", status);
      return NULL;
    }
#endif

    break;
  default:
    fprintf (stderr, "Error, unknown or unimplemented MC source type\n");
    return NULL;
  }

/* ============================================================================== */
/* ========== More Miscellaneous Stuff ========================================== */
/* ============================================================================== */
#if HAVE_LIDAR
  if (sample->LLE_moon)
    status = generate_photon_lidar (p, sample, atmos);
#endif

  /* this is a quick and dirty solution: photons starting horizontally in a non-absorbing, non-scattering */
  /* layer might cause an infinite loop - therefore we set the photon weight of these photons to 0.       */
  /* Should find a better solution - at least check the extinction coefficient and stop only if the       */
  /* latter is 0.                                                                                         */
    
  if (p->dir.dx[2]==0 && sample->bcond != MCBCOND_ABSORB) {
    fprintf (stderr, "Photon starting in horizontal direction - this might cause trouble\n");
    fprintf (stderr, "and the photon weight is set to 0 for this reason. If you see this\n");
    fprintf (stderr, "message just occasionally everything is ok. If it appears many many\n");
    fprintf (stderr, "times, you should not trust the result! Photon coordinates:\n");
    fprintf (stderr, "(%d, %d, %d), (%g, %g, %g), (%g, %g, %g)\n", p->ic, p->jc, p->kc, p->x[0], p->x[1], p->x[2], p->dir.dx[0], p->dir.dx[1], p->dir.dx[2]);

    p->weight=0;
  }

  for(ip=0; ip<4; ip++) {
    p->stokes[ip] = 0.0;
    for(j=0; j<4; j++)
      p->stokes[ip] += p->phamat[ip][j]*p->stokes0[j];
  }

  if (sample->use_p_norm) {
    double mu_norm=0.0;
    if (sample->escape)
      v_mult_mu ( p->dir.dx, sample->rad[0].dir.dx, &mu_norm );
    if (sample->LidarLocEst)
      v_mult ( sample->lidar[sample->ili].dir.dx, p->dir.dx, &mu_norm );
    p->p_norm = get_phase_max ( sample->phase_max, sample->n_phase_max, mu_norm, 0 );
  }
  else
    p->p_norm = 0.0;
  
#ifdef PRINT_PHOTON_PATH
  /* add photon to photon path */
  add_to_photonpath (&(p->path), p->x);
#endif

#if HAVE_LIDAR
  /* CB - add coordinate to cohebasca */
  if (sample->coherent_backscatter && sample->LidarLocEst) {
    add_coord_to_cohebasca ( p );
    p->pss.ic_s = p->ic;
    p->pss.jc_s = p->jc;
    p->pss.kc_s = p->kc;
  }
#endif
  fflush(stderr);
  free(norm);

  return p;
}


/***********************************************************************************/
/* Function: generate_photon_solar_thermal_backward                       @62_30i@ */
/* Description:                                                                    */
/*                                                                                 */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/
photon_struct* generate_photon_solar_thermal_backward (sample_struct* sample, 
                                                       atmosphere_struct* atmos, 
                                                       photon_struct* p, 
                                                       elevation_struct* elev,
                                                       int* hit, 
                                                       long int nphotons, 
                                                       int* counter1,
                                                       double sza,
                                                       double phi,
                                                       double phase,
                                                       int source,
                                                       int sofi,
                                                       double* norm,
                                                       double* pd)

{
  double local_elevation=0.0;
  int ie=0, je=0;
  double alpha=0;
#ifdef HAVE_SPHER
  int status=0;
#endif
/* ------------------------------------------------------------------------------ */
/* ---------- 1. SOLAR BLITZ THERMAL_BW: Various Checks ------------------------- */
/* ------------------------------------------------------------------------------ */
    /* check consistency of settings, exit if error */
    if (source==MCSRC_THERMAL_BACKWARD) {
      if (!sample->backward) {
	fprintf (stderr, "Fatal error, we shouldn't end up in MCSRC_THERMAL_BACKWARD\n");
	fprintf (stderr, "in a forward calculation - programming error!\n");
	return NULL;
      }
      if (sofi) {
	fprintf (stderr, "Fatal error, sofi and MCSRC_THERMAL_BACKWARD\n");
	fprintf (stderr, "does not make sense - programming error!\n");
	return NULL;
      }
      if (sample->backward == MCBACKWARD_EDIR || sample->backward == MCBACKWARD_FDIR ) {
	fprintf (stderr, "Fatal error, direct radiation and MCSRC_THERMAL_BACKWARD\n");
	fprintf (stderr, "does not make sense - programming error!\n");
	return NULL;
      }

      if (sample->spherical) {
	fprintf (stderr, "Warning, spherical and MCSRC_THERMAL_BACKWARD\n");
	fprintf (stderr, "has never been tested, might be bullshit!\n");
      }
    }

    if (source==MCSRC_SOLAR || source==MCSRC_BLITZ) {
      switch(sample->backward) {
      case MCBACKWARD_ABS:
      case MCBACKWARD_HEAT:
      case MCBACKWARD_EMIS:
	//      case MCBACKWARD_ACT:
	fprintf (stderr, "Fatal error, MCSRC_SOLAR does not make sense\n");
	fprintf (stderr, "with backward mode %d - programming error!\n",sample->backward);
	return NULL;
      }
    }

/* ------------------------------------------------------------------------------ */
/* ---------- 2. SOLAR BLITZ THERMAL_BW: Set Horizontal Starting Position ------- */
/* ------------------------------------------------------------------------------ */
    if (source==MCSRC_THERMAL_BACKWARD)
      p->direct=0;   /* no direct beam source */
    else
      p->direct=1;   /* direct beam source */

    /* horizontal start position */
    /* within the pixel: x = i*dx .. (i+1)*dx; y = j*dy .. (j+1)*dy*/

    if (!sample->backward) {  /* forward */
      if (source==MCSRC_BLITZ) {
	/* simulate blitz */
	alpha = uvspec_random();
	p->x[0] = ( sample->blitz_position[3] - sample->blitz_position[0] ) * alpha
	  + sample->blitz_position[0];
	p->x[1] = ( sample->blitz_position[4] - sample->blitz_position[1] ) * alpha
	  + sample->blitz_position[1];
	p->x[2] = ( sample->blitz_position[5] - sample->blitz_position[2] ) * alpha
	  + sample->blitz_position[2];
      }
      else {
	p->x[0] = atmos->xmax * uvspec_random();
	p->x[1] = atmos->ymax * uvspec_random();

	/* radial/pathlength calculation */
	if (sample->Nr>0 && sample->dt>0) {
	  p->x[0] = atmos->xmax / 2.0;
	  p->x[1] = atmos->ymax / 2.0;
	}
      }
    }
    else {  /* backward */
      if (!sample->spherical){   
	if (!sample->sensorposition) {  /* random       */
	  p->x[0] = ((double) sample->backward_is + uvspec_random()) * sample->delX;
	  p->x[1] = ((double) sample->backward_js + uvspec_random()) * sample->delY;
	  
	  /* for horizontal irradiances we start on the side face */
	  switch (sample->backward) {
	  case MCBACKWARD_EXP:
	  case MCBACKWARD_EXN:
	    p->x[0] = ((double) sample->backward_is) * sample->delX;
	    break;

	  case MCBACKWARD_EYP:
	  case MCBACKWARD_EYN:
	    p->x[1] = ((double) sample->backward_js) * sample->delY;
	    break;

	  default:
	    break;
	  }
	}
	else {                          /* user-defined */
	  p->x[0] = sample->sensorposition_x[0];
	  p->x[1] = sample->sensorposition_x[1];
	}
      }
      else {  /* spherical, start position is center of domain */
	p->x[0] = atmos->xmax / 2.0;
	p->x[1] = atmos->ymax / 2.0;
      }
    }

    if (!sample->spherical3D)
      atmos_coord (p, atmos, &(p->ic), &(p->jc));

/* ------------------------------------------------------------------------------ */
/* ---------- 3. SOLAR BLITZ THERMAL_BW: Set Solar Eclipse Starting Position ---- */
/* ------------------------------------------------------------------------------ */

    if (sofi) 
      generate_photon_sofi (atmos, p, sza, phi, pd);
    else {

/* ------------------------------------------------------------------------------ */
/* ---------- 4. SOLAR BLITZ THERMAL_BW: Set Vertical Starting Position 1 ------- */
/* ------------------------------------------------------------------------------ */
      /*********************************************/
      /* first define third photon position (x[2]) */
      /*********************************************/

      if (sample->spherical3D)
	/* set_photon_z is called later for spherical3D; the sensor
	   may still be outside the atmosphere at this point */
	p->x[2] = sample->sensorposition_x[2];
      else {
	if (source!=MCSRC_BLITZ) {
	  if (!sample->sensorposition) {
	    if (sample->zstart != -999.0) /* XXX??? or should it be (sample->zstart>0) for TH_BW? CHECK!!! */
	      p->x[2] = sample->zstart;  
	    else
	      p->x[2] = atmos->Z[0]; 
	  }
	  else
	    p->x[2] = sample->sensorposition_x[2];
	}
	/* else already set for blitz */
      }


/* ------------------------------------------------------------------------------ */
/* ---------- 5. SOLAR BLITZ THERMAL_BW: Initialize Photon Direction ------------ */
/* ------------------------------------------------------------------------------ */
      /******************************************************/
      /* second calculate normal vector relative to surface */
      /******************************************************/

      calc_normal_vector ( p, elev, sample, atmos,
			   sample->surfaceparallel && sample->zstart==-999.0,
			   sample->sensordirection && sample->zstart==-999.0,
			   norm );

      /*************************************/
      /* third initialize photon direction */
      /*************************************/

      if (generate_photon_backward_photon_direction (sample, 
                                                     atmos, 
                                                     p, 
                                                     elev,
                                                     hit, 
                                                     nphotons, 
                                                     counter1,
                                                     sza,
                                                     phi,
                                                     phase,
                                                     source,
                                                     norm) == NULL) {
        fprintf (stderr,"Error in generate_photon_backward_photon_direction()!\n");
        return NULL;
      }

/* ------------------------------------------------------------------------------ */
/* ---------- 6. SOLAR BLITZ THERMAL_BW: Set Vertical Starting Position 2 ------- */
/* ------------------------------------------------------------------------------ */
      /**********************************************************/
      /* fourth determine vertical start position of the photon */
      /**********************************************************/

      /* flat surface */
      
      /* Make always sure to call set_photon_z() AFTER the new photon     */
      /* direction has been assigned because the start index and position */
      /* might depend on direction                                        */

      if (sample->spherical3D) {
#ifdef HAVE_SPHER
	/* if photon is outside of atmosphere, move it to border */
	status = move_satellite_photon_to_toa ( p, atmos );
	if (status<0) {
	  fprintf (stderr,"Error %d returned by move_satellite_photon_to_toa()\n", status);
	  return NULL;
	}
	/* return photon with weight 0 if status>0 */
	if (status>0)
	  return p;

	coord_spherical3D (p,
			   atmos->Nx,  atmos->X,
			   atmos->Ny,  atmos->Y,
			   atmos->Nz,  atmos->Z,
			   atmos->Nyg, atmos->Yg,
			   atmos->r_earth,
			   &(p->ic), &(p->jc), &(p->kc), &(p->jcg),
			   0);

#endif
      } /* spherical3D end */
      else
	set_photon_z(p->x[2], atmos, p);

      if (generate_photon_backward_vertical_position (sample, 
                                                      atmos, 
                                                      p, 
                                                      hit, 
                                                      nphotons, 
                                                      counter1)) {
        fprintf (stderr,"Error in generate_photon_backward_vertical_position()!\n");
        return NULL;
      }

/* ------------------------------------------------------------------------------ */
/* ---------- 7. SOLAR BLITZ THERMAL_BW: Set Topography ------------------------- */
/* ------------------------------------------------------------------------------ */
      /* topography */
      if (elev->elev2D) {

	switch (sample->backward) {
	case MCBACKWARD_NONE:
	  /* nothing to do */
	  break;
	default:

	  /* elevation coordinates */
	  elev_coord (p, elev, &ie, &je);

	  /* elevation at the given coordinates */
	  local_elevation = elevation (elev->surf[ie][je], 
				       p->x[0] - (double) ie * elev->delX, 
				       p->x[1] - (double) je * elev->delY);

	  /* if surface sampling then set start altitude to elevation at start location */
	  if (!sample->sensorposition) {
	    if (sample->zstart == -999.0) {

	      /* in any case, go away from the surface a little bit */
	      /* as the photon might get a direction directly into  */
	      /* the surface and we want to avoid rounding problems */

	      /* Make always sure to call set_photon_z() AFTER the new photon     */
	      /* direction has been assigned because the start index and position */
	      /* might depend on direction                                        */
	      set_photon_z (local_elevation*(1.0+MC_EPSILON), atmos, p);
	    }
	  }     
	  else /* user-defined sensor position */
	    set_photon_z (sample->sensorposition_x[2], atmos, p);

	  /* Sampling altitude below the ground; this happens when    */
	  /* the sampling altitude zout[0] is lower than the highest  */
	  /* point of the model area; this is reasonable, but we need */
	  /* to set the photon weight to 0 and start the photon above */
	  /* ground; it will not be traced anyway and will therefore  */
	  /* not contribute to the result.                            */
	  if (p->x[2] < local_elevation) {
	    set_photon_z (local_elevation*(1.0+MC_EPSILON), atmos, p);

	    /* weight is not consistently handled - need to set both weight and tauabs.tot to NAN */
	    p->tauabs.tot = NOT_A_NUMBER;  /* NAN if photons starts below the ground */
	    p->weight = NOT_A_NUMBER;      /* NAN if photons starts below the ground */
	  }
	} /* endif sample->backward */
      } /* endif elev2D */

/* ------------------------------------------------------------------------------ */
/* ---------- 8. SOLAR BLITZ THERMAL_BW: Reference to NN ------------------------ */
/* ------------------------------------------------------------------------------ */
	/* sample pixel is referenced to NN according to photon direction */
      if (sample->reference_to_NN) {
	if (p->dir.dx[2]==0) {
	  fprintf (stderr, "Fatal error, referencing to NN does not work with\n");
	  fprintf (stderr, "horizontally moving photons! Exiting...\n");
	  return NULL;
	}

	/*	fprintf(stderr,"position before referencing %e %e %d %d\n",p->x[0],p->x[1],p->ic,p->jc); */
	/*	fprintf(stderr,"direction %e %e %e\n",p->dir.dx[0],p->dir.dx[1],p->dir.dx[2]);           */

	/* move photon along photon direction from NN to starting altitude */
	p->x[0] += p->x[2]*p->dir.dx[0]/p->dir.dx[2];
	p->x[1] += p->x[2]*p->dir.dx[1]/p->dir.dx[2];

	/* photon might have left domain, use periodic boundary conditions to move it back into it */
	/* maybe replace the following with
	   status = perform_periodic_bc (p, atmos); */
	while (p->x[0] < 0)
	  p->x[0] += atmos->xmax;
	while (p->x[0] >= atmos->xmax)
	  p->x[0] -= atmos->xmax;
	while (p->x[1] < 0)
	  p->x[1] += atmos->ymax;
	while (p->x[1] >= atmos->ymax)
	  p->x[1] -= atmos->ymax;

	/* update indices of photon */
	atmos_coord (p, atmos, &(p->ic), &(p->jc));

	/*	fprintf(stderr,"position after referencing %e %e %d %d \n",p->x[0],p->x[1],p->ic,p->jc); */

      } /* endif reference_to_NN */

      /* need to weight with the area of the actual surface for backward surface-parallel */
      if (elev->elev2D && sample->surfaceparallel && sample->backward!=MCBACKWARD_NONE) {
	if (sample->sensorposition || sample->zstart==-999.0) {  /* only for photons starting at the surface */
	  /* BMSURFACEPARALLEL */
	  p->weight        /= fabs(norm[2]);
	}
      }
    } /* !!!!!!!!!!!!!!!!!!!! endelse (sofi) !!!!!!!!!!!!!!!!!!!! */
  return p;
}

/***********************************************************************************/
/* Function: generate_photon_backward_photon_direction                    @62_30i@ */
/* Description:                                                                    */
/*                                                                                 */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/
photon_struct* generate_photon_backward_photon_direction (sample_struct* sample, 
                                                          atmosphere_struct* atmos, 
                                                          photon_struct* p, 
                                                          elevation_struct* elev,
                                                          int* hit, 
                                                          long int nphotons, 
                                                          int* counter1,
                                                          double sza,
                                                          double phi,
                                                          double phase,
                                                          int source,
                                                          double* norm)
{
                                  
  double sinsza=0,cossza=0,sinphi=0,cosphi=0;
  int status = 0;
  int i=0;
#if HAVE_VROOM
  double mu=0.0, phidummy=0.0, mu_horz=0.0, mu_ddis=0.0, P_norm=0.0, stheta_ddis=0.0, cphi_ddis=0.0;
  int DDISsing=MCDDIS_NONE, FODDISsing=0, out_of_cone=0;
  locest_struct lest;
  scadis_struct scadis;
#endif

  switch (sample->backward) {
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~ 5.1 SOLAR BLITZ THERMAL_BW / NONE RAD: Photon Direction  ~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    case MCBACKWARD_NONE:
    case MCBACKWARD_RADIANCE:
    
      /* initialize direction using the specified solar zenith and azimuth angles which,    */
      /*  backward mode is active, contain the actual viewing angles (see setup_sample2D) */
            
      /* sine and cosine of solar zenith */
      cossza = cosd(sza);
      sinsza = sind(sza);
            
      /* sine and cosine of solar azimuth */
      cosphi = cosd(phi);
      sinphi = sind(phi);

      /* use different umu for each pixel */
      if (sample->umu2D != NULL && sample->backward==MCBACKWARD_RADIANCE) {
            
        /* sine and cosine of solar zenith */
        cossza = sample->umu2D[p->ic][p->jc];
        sinsza = sqrt(1.0 - cossza*cossza);
            
        /* sine and cosine of solar azimuth */
        cosphi = cosd(sample->phi2D[p->ic][p->jc]);
        sinphi = sind(sample->phi2D[p->ic][p->jc]);
      }
    
      if (!sample->panorama) {
/* .............................................................................. */
/* ..................... 5.1.1 BLITZ / NONE RAD / Photon Direction .............. */
/* .............................................................................. */
        if (source==MCSRC_BLITZ) {
          if ( !sample->escape_eps_ddis_upf )
            random_direction (&(p->dir));
          else {
            /* First start the photon directly into the direction of
         the detector; for blitz, cos(sza)==umu, phi==phi0 */
            init_direction (sinsza, -cossza, sinphi, cosphi, &(p->dir));

            /* Choose a random starting mu according to phase
         function phase_max for the photon by abusing the
         scattering routine sc_mu */
            cossza = sc_mu ( sample->phase_max[0], 1, p->SC_mode, 0.0, 0.0 );
    
            /* Now adjust photon weight for the "unphysical"
         choice of direction: Divide photon-weight by
         probability from phase_max for sending it into this
         direction */
            status = get_phase_matrix_pft ( sample->phase_max[0], cossza, 0, 1, &phase );
            if (status) {
        fct_err_out (status, "get_phase_matrix_pft", ERROR_POSITION);
        return NULL;
            }
            p->weight /=  phase;
    
            /* Then modify direction according to give mu and phi */
            phi = sc_Isotropic_phi();
            new_direction (cossza, phi-90., &(p->dir), 0.);
          }
        }
/* .............................................................................. */
/* .............. 5.1.2 SOLAR THERMAL_BW / NONE RAD / Photon Direction .......... */
/* .............................................................................. */
        else {
                init_direction (sinsza, cossza, sinphi, cosphi, &(p->dir));
        }
                
      }  
/* .............................................................................. */
/* ... 5.1.3 SOLAR BLITZ THERMAL_BW / NONE RAD / Panorama / Photon Direction .... */
/* .............................................................................. */
      else { /* PANORAMA */

              generate_photon_backward_photon_direction_panorama ( sample, p, phi, sza, 
                                                                   cossza, sinsza, 
                                                                   cosphi, sinphi );
      }
      cp_direction (&(p->dir0), &(p->dir));
    
            /* CE: Probably not needed, CHECK!!! */ 
            /* #ifdef NEW_REFLECT */
            /*   p->phi0=calc_phi_horz ( p->dir.dx, NULL); */
            /* #endif */
    
      break;
         
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~ 5.2 SOLAR BLITZ THERMAL_BW / EDIR EGLOB EDN / Photon Direction ~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    case MCBACKWARD_EDIR:
    case MCBACKWARD_EGLOB:
    case MCBACKWARD_EDN:
      if ( sample->escape_eps_ddis_upf!=0.0 ) {
#if HAVE_VROOM
        status = mc_vroom_prep_DDIS (sample, p, atmos, &DDISsing, &FODDISsing, &out_of_cone, &lest, &scadis);
        if (status!=0) {
          fprintf(stderr,"Error %d returned by mc_vroom_prep_DDIS ()\n", status);
          return NULL;
        }
    
        random_reflection_special (sample, atmos, elev, p,
                 sample->phase_max, sample->n_phase_max,
                 &mu, &phidummy, norm,
                 &scadis, lest, 0, DDISsing, FODDISsing, out_of_cone); 
    
        status = mc_vroom_set_mus_and_phis (sample, p, DDISsing,
                    lest, scadis,
                    &mu, &mu_ddis, &stheta_ddis, phi, &cphi_ddis);
        if (status!=0) {
          fprintf(stderr,"Error %d returned by mc_vroom_set_mus_and_phis ()\n", status);
          return NULL;
        }
    
        v_mult_mu ( norm, p->dir.dx, &mu_horz);
        P_norm = get_Lambertian_phase(mu_horz);
    
        if (p->weight==0) {
          p->weight=0;
          return p;
        }
    
        status = mc_vroom_DDIS_weight_and_prep_stuff (sample, atmos, mu_ddis, stheta_ddis, cphi_ddis,
                  DDISsing, out_of_cone,
                  P_norm, P_norm, P_norm, lest, scadis, p);
        if (status!=0) {
          fprintf(stderr,"Error %d returned by mc_vroom_DDIS_weight_and_prep_stuff ()\n", status);
          return NULL;
        }
#endif
      }
      else
        random_Lambertian_normal (&(p->dir), norm);
      break;
        
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~ 5.3 SOLAR BLITZ THERMAL_BW / EUP / Photon Direction ~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    case MCBACKWARD_EUP: 
  /* reverse direction because we need to go downward */
      for (i=0; i<3; i++)
        norm[i] = - norm[i];
      random_Lambertian_normal (&(p->dir), norm);
      break;

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~ 5.4 SOLAR BLITZ THERMAL_BW / FDIR / Photon Direction ~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    case MCBACKWARD_FDIR:
      if (sample->spherical3D) {
        fprintf(stderr,"Error, backward fdir and spherical3d bugs!\n");
        return NULL;
      }
    
      random_Isotropic_normal (&(p->dir), norm); 
      p->weight        = 1.0 / fabs(cosd(sza));
      break;

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~ 5.5 SOLAR BLITZ THERMAL_BW / FDN / Photon Direction ~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    case MCBACKWARD_FDN:
      random_Isotropic_normal (&(p->dir), norm); 
      p->weight        = 2.0;
      break;

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~ 5.6 SOLAR BLITZ THERMAL_BW / FUP / Photon Direction ~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    case MCBACKWARD_FUP:
  /* reverse direction because we need to go inward */
      for (i=0; i<3; i++)
        norm[i] = - norm[i];
      random_Isotropic_normal (&(p->dir), norm); 
      p->weight        = 2.0;
      break;


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~ 5.7 SOLAR BLITZ THERMAL_BW / ABS EMIS ACT / Photon Direction ~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    case MCBACKWARD_ABS:
    case MCBACKWARD_EMIS:
    case MCBACKWARD_ACT:  /* **CK 2013.08.27 remove MCBACKWARD_HEAT from here and add in following lines the specific calculations */

      random_direction(&(p->dir));
      break;
    
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~ 5.8 SOLAR BLITZ THERMAL_BW / HEAT / Photon Direction ~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    case MCBACKWARD_HEAT: /* **CK 2013.08.27 Apply different thermal backward heating rate calculations */
      switch (sample->heat_flag[sample->backward_is][sample->backward_js]) {
      case MCBACKWARD_HEAT_DENET:
        hunt_modified (atmos->Z, atmos->Nz, p->x[2], &(p->kc), hit);  /* **CK 2014.01.13 added hunt_modified here! */
        status = denet_set_direction (p, atmos, sample, nphotons, counter1, norm); /* **CK 2013.08.27 */
        if (status!=0) {
          fprintf (stderr, "Error %d returned by denet_set_direction()\n", status);
          return NULL;
        }
        
        break;
    
      case MCBACKWARD_HEAT_EMABS:
      case MCBACKWARD_HEAT_EMABSOPT:
        random_direction(&(p->dir));
        break;
    
      default: 
        fprintf (stderr, "Error, sample->heat_flag shouldn't be %d here!\n", sample->heat_flag[sample->backward_is][sample->backward_js]);
        return NULL;
      }
      break;
    
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~ 5.9 SOLAR BLITZ THERMAL_BW / EXP / Photon Direction ~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    case MCBACKWARD_EXP: 
      norm[0]=-1; norm[1]=0; norm[2]=0; 
      random_Lambertian_normal (&(p->dir), norm);
      break;

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~ 5.10 SOLAR BLITZ THERMAL_BW / EXN / Photon Direction ~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    case MCBACKWARD_EXN: 
      norm[0]=1; norm[1]=0; norm[2]=0; 
      random_Lambertian_normal (&(p->dir), norm);
      break;

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~ 5.11 SOLAR BLITZ THERMAL_BW / EYP / Photon Direction ~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    case MCBACKWARD_EYP: 
      norm[0]=0; norm[1]=-1; norm[2]=0; 
      random_Lambertian_normal (&(p->dir), norm);
      break;

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~ 5.12 SOLAR BLITZ THERMAL_BW / EYN / Photon Direction ~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    case MCBACKWARD_EYN: 
      norm[0]=0; norm[1]=+1; norm[2]=0; 
      random_Lambertian_normal (&(p->dir), norm);
      break;

    default:
      fprintf (stderr, "Error, unknown backward quantity %d\n",
	       sample->backward);
    return NULL;
  } 


  return p;
}

/***********************************************************************************/
/* Function: generate_photon_backward_photon_direction_panorama           @62_30i@ */
/* Description:                                                                    */
/*                                                                                 */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/
void generate_photon_backward_photon_direction_panorama ( sample_struct* sample, 
                                                          photon_struct* p,
                                                          double phi, 
                                                          double sza, 
                                                          double cossza,
                                                          double sinsza, 
                                                          double cosphi, 
                                                          double sinphi )

{
  double cossza_forward=0.0, sinsza_forward=0.0, cosphi0_forward=0.0, sinphi0_forward=0.0;

  /* initialize direction                             */
  /* and copy direction to "initial direction" struct */
  if (sample->pan_distr_photons_over_pixel) {
    phi = random_variation_above ( phi, sample->pan_dphi );
    /* Bug fix, the following two lines were missing. RPB */
    cosphi = cosd(phi);
    sinphi = sind(phi);

    if( ( sample->pan_circumsolar_var_red && sample->escape_eps_ddis_upf ) ){
       /* Distribute photon directions according to
       phase_max. 
       BE AWARE THAT dist_photon_direction_by_phasemax alters 
       sza, cossza and sinsza*/
       dist_photon_direction_by_phasemax(p, sample, 
           &sza, &cossza, &sinsza);
    }
    else {
      /* Distribute photon directions geometrically */
      sza = random_theta_above ( sza, sample->pan_dtheta );
      cossza = cosd(sza);
      sinsza = sind(sza);
    }
  } /* end if pan_distr_photons_over_pixel */

  if (sample->pan_weight_with_cos)
    p->weight1 = cossza;

  if (sample->pan_alignment!=MCPAN_ALIGNMENT_NONE && sample->align_theta != 180.0 ) {
    /* Get the direction that points from the instrument to
       the sun/alignment */

    sinsza_forward  = sind(sample->align_theta);
    cossza_forward  = cosd(sample->align_theta);
    sinphi0_forward = sind(sample->align_phi);
    cosphi0_forward = cosd(sample->align_phi);

    /* First start the photon directly into the direction of
       alignment */
    init_direction (sinsza_forward, cossza_forward,
        sinphi0_forward, cosphi0_forward, &(p->dir));
    cp_direction (&(p->dir00), &(p->dir));

    /* Then modify direction according to give mu and phi */
    new_direction (-cossza, 90.-phi+sample->add_phi_when_aligning, &(p->dir), 0.);
  }
  else {
    /* this is the standard case!!! */
          init_direction (sinsza, cossza, sinphi, cosphi, &(p->dir));
    init_direction (0,-1,0,0,&(p->dir00));
  }

  if (sample->sample_backward_sunshape) {
    /* p->dir is pointing to the center of the sun modify
       direction slightly to point to somewhere on the sun
       disc directions are chosen according sun shape which
       acts a pdf for the directions (Bernhard Reinhardt) */
    cp_direction (&(sample->rad[0].dir),&(sample->rad_dir0));
    limb_dark_lest_dir ( &(sample->rad[0].dir), sample,
       &(p->lest_phi), &(p->lest_alpha) );
  }
}

/***********************************************************************************/
/* Function: generate_photon_backward_vertical_position                   @62_30i@ */
/* Description:                                                                    */
/*                                                                                 */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int generate_photon_backward_vertical_position (sample_struct* sample, 
                                                atmosphere_struct* atmos, 
                                                photon_struct* p, 
                                                int* hit, 
                                                long int nphotons, 
                                                int* counter1)
{
  int kc_save=0;
  int oopscounter=0;

  switch (sample->backward) {

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~ 6.1 SOLAR BLITZ THERMAL_BW / ABS / Vertical Position ~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    case MCBACKWARD_ABS:
      if (sample->sensorposition) {
        fprintf (stderr, "Error, don't know how to calculate box-averaged absorption/emission/heating\n");
        fprintf (stderr, "for a user-defined sensor position!\n");
        return -1;
      }
            
      /* for absorption/emission/heating/actinic calculations, the backward photon starts in the layer, not at the level */
      
      /* Make always sure to call set_photon_z() AFTER the new photon     */
      /* direction has been assigned because the start index and position */
      /* might depend on direction                                        */
      
      /* the call to hunt_modified() might seem a bit weird - we need it, however, */
      /* as we need the index of the layer above p->x[2] irrespective of the       */
      /* photon direction! This basically undoes the "p->kc--;" in set_photon_z()  */
      hunt_modified (atmos->Z, atmos->Nz, p->x[2], &(p->kc), hit);
    
      /* Use absorption coefficient as photon weight */
      if (atmos->kabs3D->threed[MCCAOTH_TOT][p->kc]>=1)
        p->weight = atmos->kabs3D->prof [MCCAOTH_TOT][p->kc][p->ic][p->jc]
          * ( atmos->Z[p->kc+1] - atmos->Z[p->kc] );
      else
        p->weight = atmos->kabs->prof [MCCAOTH_TOT][p->kc]
          * ( atmos->Z[p->kc+1] - atmos->Z[p->kc] );
      
      /* finally, set the real photon position to a random location in the layer */
      /* set_photon_z (p->x[2] + uvspec_random()*(atmos->Z[p->kc+1]-atmos->Z[p->kc]), atmos, p); */
      set_photon_z (p->x[2] + uvspec_random()*(atmos->Z[p->kc+1]-p->x[2]), atmos, p);
            
      
      break;
    
      
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~ 6.2 SOLAR BLITZ THERMAL_BW / EMIS / Vertical Position ~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    case MCBACKWARD_EMIS:
      
      if (sample->sensorposition) {
        fprintf (stderr, "Error, don't know how to calculate box-averaged absorption/emission/heating\n");
        fprintf (stderr, "for a user-defined sensor position!\n");
        return -1;
      }
            
      /* for absorption/emission/heating/actinic calculations, the backward photon starts in the layer, not at the level */
      
      /* Make always sure to call set_photon_z() AFTER the new photon     */
      /* direction has been assigned because the start index and position */
      /* might depend on direction                                        */
      
      /* the call to hunt_modified() might seem a bit weird - we need it, however, */
      /* as we need the index of the layer above p->x[2] irrespective of the       */
      /* photon direction! This basically undoes the "p->kc--;" in set_photon_z()  */
      hunt_modified (atmos->Z, atmos->Nz, p->x[2], &(p->kc), hit);
    
      /* Use absorption coefficient as photon weight */
      if (atmos->kabs3D->threed[MCCAOTH_TOT][p->kc]>=1)
        p->weight = atmos->kabs3D->prof [MCCAOTH_TOT][p->kc][p->ic][p->jc]
          * ( atmos->Z[p->kc+1] - atmos->Z[p->kc] );
      else
        p->weight = atmos->kabs->prof [MCCAOTH_TOT][p->kc]
          * ( atmos->Z[p->kc+1] - atmos->Z[p->kc] );
      
      /* finally, set the real photon position to a random location in the layer */
      /* set_photon_z (p->x[2] + uvspec_random()*(atmos->Z[p->kc+1]-atmos->Z[p->kc]), atmos, p); */
      set_photon_z (p->x[2] + uvspec_random()*(atmos->Z[p->kc+1]-p->x[2]), atmos, p);
            
      
      break;
    
    
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~ 6.3 SOLAR BLITZ THERMAL_BW / HEAT / Vertical Position ~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    case MCBACKWARD_HEAT:   /* **CK  2013.08.27 */
      if (generate_photon_backward_vertical_position_heat (sample, atmos, p, 
							   hit, nphotons, 
							   counter1 )) {
        fprintf (stderr, "Error in generate_photon_backward_vertical_position_heat()!\n");
        return -1;
      }

      break;
    
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~ 6.4 SOLAR BLITZ THERMAL_BW / ACT / Vertical Position ~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    case MCBACKWARD_ACT: 
      
      if (sample->sensorposition) {
        fprintf (stderr, "Error, don't know how to calculate box-averaged actinic flux\n");
        fprintf (stderr, "for a user-defined sensor position!\n");
        return -1;
      }
      
      /* for absorption/emission/heating/actinic calculations, the backward photon starts in the layer, not at the level */  

      /* Make always sure to call set_photon_z() AFTER the new photon     */
      /* direction has been assigned because the start index and position */
      /* might depend on direction                                        */

      /* the call to hunt_modified() might seem a bit weird - we need it, however, */
      /* as we need the index of the layer above p->x[2] irrespective of the       */
      /* photon direction! This basically undoes the "p->kc--;" in set_photon_z()  */
      hunt_modified (atmos->Z, atmos->Nz, p->x[2], &(p->kc), hit);

      kc_save = p->kc;

      /* for safety reasons we stay away from the upper and lower boundaries */
      set_photon_z (p->x[2] + (MC_EPSILON+(1.0-2.0*MC_EPSILON)*uvspec_random())*(atmos->Z[p->kc+1]-atmos->Z[p->kc]), atmos, p);
      
      oopscounter=0;
      do {
	if (oopscounter>0)
	  fprintf (stderr, "Oops %d, need to throw the dice more than once because %.10f was too close to boundary\n", oopscounter, p->x[2]); 
	
	  /* for safety reasons we stay away from the upper and lower boundaries */
	set_photon_z (atmos->Z[kc_save]*(1+MC_EPSILON) + (atmos->Z[kc_save+1]*(1-MC_EPSILON) - atmos->Z[kc_save]*(1+MC_EPSILON))*uvspec_random(), atmos, p);
	
	  oopscounter++;
	  
	  /* and we roll the dice until we are really away from the boundaries: careful: z coordinate = float! */
      } while (kc_save != p->kc);
    
      break;
            
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* -~~~~ 6.5 SOLAR BLITZ THERMAL_BW / horizontal fluxes / Vertical Position ~~~~~ */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    case MCBACKWARD_EXP: 
    case MCBACKWARD_EXN: 
    case MCBACKWARD_EYP: 
    case MCBACKWARD_EYN: 

      if (p->kc==atmos->Nz-1) {
	fprintf (stderr, "Error, it doesn't make sense to calculate horizontal fluxes above TOA!\n");
	return -1;
      }

      /* the call to hunt_modified() might seem a bit weird - we need it, however, */
      /* as we need the index of the layer above p->x[2] irrespective of the       */
      /* photon direction! This basically undoes the "p->kc--;" in set_photon_z()  */
      hunt_modified (atmos->Z, atmos->Nz, p->x[2], &(p->kc), hit);

      /* random location on vertical grid box face for horizontal irradiances */
      if (!sample->sensorposition) {
	kc_save = p->kc;

	oopscounter=0;
	do {
	  if (oopscounter>0)
	    fprintf (stderr, "Oops %d, need to throw the dice more than once because %.10f was too close to boundary\n", oopscounter, p->x[2]); 

	  /* for safety reasons we stay away from the upper and lower boundaries */
	  set_photon_z (atmos->Z[kc_save]*(1+MC_EPSILON) + (atmos->Z[kc_save+1]*(1-MC_EPSILON) - atmos->Z[kc_save]*(1+MC_EPSILON))*uvspec_random(), atmos, p);

	  oopscounter++;

	  /* and we roll the dice until we are really away from the boundaries: careful: z coordinate = float! */
	} while (kc_save != p->kc);
      }
      
      break;

    default: 
      break;
  } /* end switch sample->backward */
  return 0;
}


/***********************************************************************************/
/* Function: generate_photon_backward_vertical_position_heat              @62_30i@ */
/* Description:                                                                    */
/*                                                                                 */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int generate_photon_backward_vertical_position_heat (sample_struct* sample, 
                                                     atmosphere_struct* atmos, 
                                                     photon_struct* p, 
                                                     int* hit, 
                                                     long int nphotons, 
                                                     int* counter1 )
{
  int status = 0;
  int kc_save = 0;
  int oopscounter = 0;
  
  switch (sample->heat_flag[sample->backward_is][sample->backward_js]) {

/* .............................................................................. */
/* ........ 6.3.1 SOLAR BLITZ THERMAL_BW / HEAT / DENET / Vertical Position ..... */
/* .............................................................................. */
  case MCBACKWARD_HEAT_DENET:
      
    if (sample->sensorposition) {
      fprintf (stderr, "Error, don't know how to calculate box-averaged absorption/emission/heating\n");
      fprintf (stderr, "for a user-defined sensor position!\n");
      return -1;
    }
    hunt_modified (atmos->Z, atmos->Nz, p->x[2], &(p->kc), hit);
    p->x[0] = ((double) sample->backward_is + uvspec_random()) * sample->delX;
    p->x[1] = ((double) sample->backward_js + uvspec_random()) * sample->delY;

    kc_save=p->kc;
    do {
      if (oopscounter>0)
	fprintf (stderr, "Oops %d, need to throw the dice more than once because %.10f was too close to boundary\n", oopscounter, p->x[2]);
      
      /* for safety reasons we stay away from the upper and lower boundaries */
      set_photon_z (atmos->Z[kc_save]*(1+MC_EPSILON) + (atmos->Z[kc_save+1]*(1-MC_EPSILON) - atmos->Z[kc_save]*(1+MC_EPSILON))*uvspec_random(), atmos, p);
      
      oopscounter++;
      
      /* and we roll the dice until we are really away from the boundaries: careful: z coordinate = float! */
    } while (kc_save != p->kc);
    
    status = denet_weight(p, atmos, sample, nphotons, counter1);
    if (status!=0) {
      fprintf (stderr, "Error %d returned by denet_weight()\n", status);
      return -1;
    }
    break;
    
/* .............................................................................. */
/* ......... 6.3.2 SOLAR BLITZ THERMAL_BW / HEAT / EMABS / Vertical Position .... */
/* .............................................................................. */
  case MCBACKWARD_HEAT_EMABS:
    if (sample->sensorposition) {
      fprintf (stderr, "Error, don't know how to calculate box-averaged absorption/emission/heating\n");
      fprintf (stderr, "for a user-defined sensor position!\n");
      return -1;
    }
        
    /* for absorption/emission/heating/actinic calculations, the backward photon starts in the layer, not at the level */
  
    /* Make always sure to call set_photon_z() AFTER the new photon     */
    /* direction has been assigned because the start index and position */
    /* might depend on direction                                        */
  
    /* the call to hunt_modified() might seem a bit weird - we need it, however, */
    /* as we need the index of the layer above p->x[2] irrespective of the       */
    /* photon direction! This basically undoes the "p->kc--;" in set_photon_z()  */
    hunt_modified (atmos->Z, atmos->Nz, p->x[2], &(p->kc), hit);

    if (atmos->kabs3D->threed[MCCAOTH_TOT][p->kc]>=1){   /* if 3D */
      /* calculate weights */
      p->weight      = atmos->kabs3D->prof [MCCAOTH_TOT][p->kc][p->ic][p->jc] * ( atmos->Z[p->kc+1] - atmos->Z[p->kc] );
      p->weight_emis = - atmos->kabs3D->prof [MCCAOTH_TOT][p->kc][p->ic][p->jc] * ( atmos->Z[p->kc+1] - atmos->Z[p->kc] );  
    }  /* end 3D */
    else{   /* if 1D */
      p->weight      =  atmos->kabs->prof [MCCAOTH_TOT][p->kc] * ( atmos->Z[p->kc+1] - atmos->Z[p->kc] );
      p->weight_emis = - atmos->kabs->prof [MCCAOTH_TOT][p->kc] * ( atmos->Z[p->kc+1] - atmos->Z[p->kc] );
    }  /* end 1D */


    /* finally, set the real photon position to a random location in the layer */
    /* set_photon_z (p->x[2] + uvspec_random()*(atmos->Z[p->kc+1]-atmos->Z[p->kc]), atmos, p); */
    set_photon_z (p->x[2] + uvspec_random()*(atmos->Z[p->kc+1]-p->x[2]), atmos, p);

    /* Calculate Emission at the emission-location of the photon - not analytically! */        
    p->weight_emis *= ( atmos->Bplanck->prof  [MCCAOTH_TOT][p->kc][p->ic][p->jc] +
            ( ( ( p->x[2] - atmos->Z[p->kc] ) / ( atmos->Z[p->kc+1] - atmos->Z[p->kc] ) ) *
        ( atmos->Bplanck->prof  [MCCAOTH_TOT][p->kc+1][p->ic][p->jc]
          - atmos->Bplanck->prof  [MCCAOTH_TOT][p->kc][p->ic][p->jc] ) ) );
    
    break;

/* .............................................................................. */
/* ....... 6.3.3 SOLAR BLITZ THERMAL_BW / HEAT / EMABSOPT / Vertical Position ... */
/* .............................................................................. */
  case MCBACKWARD_HEAT_EMABSOPT:
    if (sample->sensorposition) {
      fprintf (stderr, "Error, don't know how to calculate box-averaged absorption/emission/heating\n");
      fprintf (stderr, "for a user-defined sensor position!\n");
      return -1;
    }

    /* for absorption/emission/heating/actinic calculations, the backward photon starts in the layer, not at the level */
  
    /* Make always sure to call set_photon_z() AFTER the new photon     */
    /* direction has been assigned because the start index and position */
    /* might depend on direction                                        */
  
    /* the call to hunt_modified() might seem a bit weird - we need it, however, */
    /* as we need the index of the layer above p->x[2] irrespective of the       */
    /* photon direction! This basically undoes the "p->kc--;" in set_photon_z()  */
    hunt_modified (atmos->Z, atmos->Nz, p->x[2], &(p->kc), hit);

  
    if (atmos->kabs3D->threed[MCCAOTH_TOT][p->kc]>=1){   /* if 3D */
      /* calculate weights */
      p->weight      = atmos->kabs3D->prof [MCCAOTH_TOT][p->kc][p->ic][p->jc] * ( atmos->Z[p->kc+1] - atmos->Z[p->kc] );
      p->weight_emis = - atmos->kabs3D->prof [MCCAOTH_TOT][p->kc][p->ic][p->jc] * ( atmos->Z[p->kc+1] - atmos->Z[p->kc] );  
    }  /* end 3D */
    else{   /* if 1D */
      p->weight      =  atmos->kabs->prof [MCCAOTH_TOT][p->kc] * ( atmos->Z[p->kc+1] - atmos->Z[p->kc] );
      p->weight_emis = - atmos->kabs->prof [MCCAOTH_TOT][p->kc] * ( atmos->Z[p->kc+1] - atmos->Z[p->kc] );
    }  /* end 1D */
    

    status = emabsopt_location (p, atmos, sample, nphotons);
    if (status!=0) {
      fprintf (stderr, "Error %d returned by emabsopt_location()\n", status);
      return -1;
    }
    
    p->weight_emis *= ( atmos->Bplanck->prof  [MCCAOTH_TOT][p->kc][p->ic][p->jc] +
            ( ( ( p->x[2] - atmos->Z[p->kc] ) / ( atmos->Z[p->kc+1] - atmos->Z[p->kc] ) ) *
        ( atmos->Bplanck->prof  [MCCAOTH_TOT][p->kc+1][p->ic][p->jc]
          - atmos->Bplanck->prof  [MCCAOTH_TOT][p->kc][p->ic][p->jc] ) ) );

    break;

      
  default:
    fprintf (stderr, "No Method used for setup of thermal backward heating rates!!!\n");
    return -1;

  }
  return 0;
}


/***********************************************************************************/
/* Function: destroy_photon                                               @62_30i@ */
/* Description:                                                                    */
/*  Free memory of a photon_struct.                                                */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void destroy_photon (photon_struct *p, int n_caoth)
{
  int i=0;

  if (p->pdir != NULL)
    free (p->pdir);

  if (p->stokes != NULL) {
    free(p->stokes0);
    free(p->stokes);
  }

#if HAVE_LIDAR
  destroy_photon_lidar_part(p, n_caoth);
#endif

  if (p->phamat != NULL){
    for(i=0; i<4; i++)
      free(p->phamat[i]);
    free(p->phamat);
  }
  
  /* ALIS or concentration importance sampling */
  if (p->pathlength_per_layer!=NULL){ 
    free(p->pathlength_per_layer);
    free(p->q_spectral);
    free(p->q2_spectral);
    free(p->q_albedo_spectral);
    free(p->q_concentration);
    free(p->q2_concentration);
    /*     free(p->dtauabs_spectral); */
  } 

  free(p->tocirco);
  
  free(p);
}


/***********************************************************************************/
/* Function: slt2hrz                                                      @62_30i@ */
/* Description:                                                                    */
/*  Calculate the conversion factor for irradiance on a slant surface to           */
/*  irradiance on a horizonal surface.                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static double slt2hrz (elevation_struct *elev, photon_struct *p, int surfaceparallel, int incoming)
{

  double slant2horz = 1;
  int ie=0, je=0;
  double cophi=0, coalp=0;
  double norm[3] = {0,0,0};
  
  if (elev->elev2D) {
    /* phi:   incidence angle on slant surface                  */
    /* theta:   incidence angle on horizontal surface           */
    /* alpha: angle between slant and horizontal surface        */
              
    /* elevation upward normal vector */
    elev_coord_normal_vector (p, elev, norm);

    cophi = (norm[0]*p->dir.dx[0] + norm[1]*p->dir.dx[1] + norm[2]*p->dir.dx[2]);

    if (incoming)
      cophi = -cophi;

    /* and finally, the angle between the slant surface and the horizontal */
    coalp = norm[2];
    
              
    if (cophi<=0) {
      elev_coord (p, elev, &ie, &je);
      fprintf (stderr, "OUCH!\n");
      fprintf (stderr, "%d %d %g %g %g %g\n", ie, je, cophi, p->x[0], p->x[1], p->x[2]);
      slant2horz = 0;   /* don't count photon */
      return -1;
    }
    else
      slant2horz = coalp*p->dir.cotheta/cophi;
    
    /* much simpler if surface-parallel instead of horizontal */

    if (surfaceparallel)
      /* BMSURFACEPARALLEL */
      slant2horz = 1.0;
      /* slant2horz = coalp; */
  }
  
  
  return slant2horz;
}


/***********************************************************************************/
/* Function: area_average                                                 @62_30i@ */
/* Description:                                                                    */
/*  Calculate area averages of irradiance and actinic flux.                        */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

 static void area_average (radiation_field *res, double ***back,
                           int Nx, int Ny,
			   int islower, int isupper, int jslower, int jsupper, 
			   int   *ndir, int   *ndn, int   *nup,
			   float *edir, float *edn, float *eup,
			   float *fdir, float *fdn, float *fup,
			   int backward, double **surface_area, int **surface_area_counter,
			   int elev2D, int surfaceparallel)
{
  int is=0, js=0, counter=0;
  double edir_d=0, edn_d=0, eup_d=0;
  double fdir_d=0, fdn_d=0, fup_d=0;
  double normalization=0;
  double factor=0;

  if (islower<0)
    islower=0;

  if (isupper<0)
    isupper=Nx-1;

  if (jslower<0)
    jslower=0;

  if (jsupper<0)
    jsupper=Ny-1;


  if (elev2D && surfaceparallel) {
    for (is=islower; is<=isupper; is++)
      for (js=jslower; js<=jsupper; js++) {
	factor+=surface_area[is][js];
	counter++;
      }
    
    factor/=counter;
  }
  else 
    factor = 1.0;
  

  if (backward)
    normalization = (double) ((isupper-islower+1)*(jsupper-jslower+1)) * factor;
  else 
    normalization = (double) (Nx * Ny) * factor;
  
  /* reset result */
  *ndir = 0; *ndn  = 0; *nup  = 0;
  *edir = 0; *edn  = 0; *eup  = 0;
  *fdir = 0; *fdn  = 0; *fup  = 0;
  
  /* initialize counter to NAN if backward;                         */
  /* in the forward case we are sure that we will get contributions */
  /* to all quantities while backward we only get one quantitiy     */
  /* while the others should be NAN                                 */

  if (backward!=MCBACKWARD_NONE) {
    edir_d = 0.0 / 0.0; 
    edn_d  = 0.0 / 0.0;
    eup_d  = 0.0 / 0.0;
    fdir_d = 0.0 / 0.0; 
    fdn_d  = 0.0 / 0.0; 
    fup_d  = 0.0 / 0.0;
  }

  switch (backward) {
  case MCBACKWARD_NONE:
    break;

  case MCBACKWARD_EDIR:
    edir_d = 0.0;
    break;

  case MCBACKWARD_EDN:
  case MCBACKWARD_EXN:
  case MCBACKWARD_EYN:
    edn_d = 0.0;
    break;

  case MCBACKWARD_EUP:
  case MCBACKWARD_EXP:
  case MCBACKWARD_EYP:
    eup_d = 0.0;
    break;

  case MCBACKWARD_FDIR:
    fdir_d = 0.0;
    break;

  case MCBACKWARD_FDN:
    fdn_d = 0.0;
    break;

  case MCBACKWARD_FUP:
    fup_d = 0.0;
    break;

  case MCBACKWARD_EGLOB: 
  case MCBACKWARD_ABS:
  case MCBACKWARD_ACT:
  case MCBACKWARD_EMIS:
  case MCBACKWARD_HEAT:
    /* nothing to do */
    break;



  case MCBACKWARD_RADIANCE:
    break;

  default:
    fprintf (stderr, "Error, unknown backward quantity %d in area_average()\n", backward);
    return;
  }

  
  /* need to calculate the sum in double precision because */
  /* the results may cover a high dynamic range            */
  for (is=0; is<Nx; is++) {
    for (js=0; js<Ny; js++) {
      
      *ndir   += res->ndir[is][js];
      *ndn    += res->ndn [is][js];
      *nup    += res->nup [is][js];

      /* weight with surface area of each pixel */
      if (elev2D && surfaceparallel) {                                                                   
	if (is>=islower && is<=isupper && js>=jslower && js<=jsupper)                                    
	  factor = surface_area[is][js];                                                                   
	else                                                                                               
	  factor = 1.0;                                                                                    
      }                                                                                                  
      else 
	factor = 1.0;
      
      if (!backward) {
        edir_d += res->edir[is][js]*factor;
        edn_d  += res->edn [is][js]*factor;
        eup_d  += res->eup [is][js]*factor;
        
        fdir_d += res->fdir[is][js]*factor;
        fdn_d  += res->fdn [is][js]*factor;
        fup_d  += res->fup [is][js]*factor;
      }
      else {
        if (is>=islower && is<=isupper && js>=jslower && js<=jsupper) {
	  
          switch (backward) {
          case MCBACKWARD_NONE:
          case MCBACKWARD_RADIANCE:
            break;
            
          case MCBACKWARD_EDIR:
            edir_d += back [is][js][0]*factor;
            break;
            
          case MCBACKWARD_EDN:
          case MCBACKWARD_EXN:
          case MCBACKWARD_EYN:
            edn_d  += back [is][js][0]*factor;
            break;
            
          case MCBACKWARD_EUP:
          case MCBACKWARD_EXP:
          case MCBACKWARD_EYP:
            eup_d  += back [is][js][0]*factor;
            break;
            
          case MCBACKWARD_FDIR:
            fdir_d += back [is][js][0]*factor;
            break;
            
          case MCBACKWARD_FDN:
            fdn_d  += back [is][js][0]*factor;
            break;
            
          case MCBACKWARD_FUP:
            fup_d  += back [is][js][0]*factor;
            break;
            
          case MCBACKWARD_EGLOB: 
          case MCBACKWARD_ABS: 
          case MCBACKWARD_ACT:  
          case MCBACKWARD_HEAT:  
          case MCBACKWARD_EMIS:  
            /* nothing to do */  
            break;

          default:
            fprintf (stderr, "Error, unknown backward quantity %d in area_average()\n", backward);
            return;
          }
        }
      }
    }
  }
    
  /* UVSPEC wants averaged intensity, not actinic flux */
  fdir_d /= (4.0 * PI);
  fdn_d  /= (4.0 * PI);
  fup_d  /= (4.0 * PI);


  *edir = (float) (edir_d / normalization);
  *edn  = (float) (edn_d  / normalization);
  *eup  = (float) (eup_d / normalization);
  *fdir = (float) (fdir_d / normalization);
  *fdn  = (float) (fdn_d  / normalization);
  *fup  = (float) (fup_d  / normalization);
}


/***********************************************************************************/
/* Function: summarize_result                                             @62_30i@ */
/* Description:                                                                    */
/*  Summarize result, create result arrays, write ASCII output files.              */ 
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int summarize_result (result_struct *result, sample_struct *sample,
                             atmosphere_struct *atmos, albedo_struct *albedo,
                             int source, long int mcphotons, int mcsimulations, double sza, 
                             int absorption, int write_files, int write_output_as_netcdf,
                             int tenstream,
                             char *basename,
                             float *rfldir, float *rfldn,  float *flup,
                             float *uavgso, float *uavgdn, float *uavgup,
                             float ***rfldir3d, float ***rfldn3d, float ***flup3d, 
                             float ***uavgso3d, float ***uavgdn3d, float ***uavgup3d, 
                             float ****radiance3d, float ***absback3d,
                             float ***abs3d, float ******radiance3d_is,
                             float ***rfldir3d_var, float ***rfldn3d_var, float ***flup3d_var, 
                             float ***uavgso3d_var, float ***uavgdn3d_var, float ***uavgup3d_var, 
                             float ****radiance3d_var, float ***absback3d_var,
                             float ***abs3d_var, 
			     int escape, int elev2D, double wavelength, int quiet)
{
  int lu=0, nu=0;
  double incident=0;
  int *ndir=NULL, *ndn=NULL, *nup=NULL;

  int isp=0;        /* caoth index   */
  int kc=0;         /* z coordinate  */

  /* XXXumu_2d??? CHECK!!! */ 

  int islower=0, jslower=0, isupper=0, jsupper=0;
  int lu3D=0;

  char specradfilename [FILENAME_MAX] ="";

  double factor=0.;
  double dmcphotons = (double) mcphotons;

  char flxext[FILENAME_MAX] = "", flxext2[FILENAME_MAX] = "";
  char radext[FILENAME_MAX] = "", radext2[FILENAME_MAX] = "";
  char rplext[FILENAME_MAX] = "", rplext2[FILENAME_MAX] = "";
  char absext[FILENAME_MAX] = "", absext2[FILENAME_MAX] = "";
  char bacext[FILENAME_MAX] = "", bacext2[FILENAME_MAX] = "";
  char bacjacext[FILENAME_MAX] = "", bacjacext2[FILENAME_MAX] = "";
  char picext[FILENAME_MAX] = "";
  char extension[FILENAME_MAX] = "", extension2[FILENAME_MAX] = "";

  if(tenstream) 
    return 0;

/* ======================= IMPORTANT ======================== */
/* If you want to implement a new output method please put it */
/* into a new function and use this function to call it!      */
/* Do not directly add your method into this function!        */
/* ========================================================== */
 
/* ====== 1. Lidar Output, special, needs nothing of the rest! */
  if (sample->LidarLocEst) {
#if HAVE_LIDAR
    int status=0;
    factor = (double) mcsimulations / dmcphotons;
    if ((status = summarize_result_lidar (sample, atmos, result, basename, write_files, 
                                          write_output_as_netcdf, factor)) ) {
      fprintf (stderr, "Error, summarize_result_lidar_part returned status %d\n", status);
      return -1;
    }
#endif
  }
  /* not lidar, else brace goes until end of function! */
  else {
 
/* ====== 2. Preparing everything */
    /* allocate memory for photon counters */
    ndir = calloc (atmos->Nz+1, sizeof(int));
    ndn  = calloc (atmos->Nz+1, sizeof(int));
    nup  = calloc (atmos->Nz+1, sizeof(int));
  
    /* extension for absorption/heating/actinic etc. file */
    switch (absorption) {
    case MCFORWARD_ABS_ACTINIC:
      strcpy (extension,  ".act");
      strcpy (extension2, ".act.std");
      break;
  
    case MCFORWARD_ABS_ABSORPTION:
    case MCFORWARD_ABS_EMISSION:
    case MCFORWARD_ABS_HEATING:
    case MCFORWARD_ABS_NONE:
      strcpy (extension,  ".abs");
      strcpy (extension2, ".abs.std");
      break;
  
    default:
      fprintf (stderr, "Error, unknown absorption type %d\n", absorption);
      return -1;
    }
  
    switch (source) {
    case MCSRC_THERMAL_ATMOSPHERE:
      strcpy (radext, ".atm.rad");
      strcpy (rplext, ".atm.rpl");
      strcpy (flxext, ".atm.flx");
  
      strcpy (absext, ".atm");
      strcat (absext, extension);
  
      strcpy (radext2, ".atm.rad.std");
      strcpy (rplext2, ".atm.rpl.std");
      strcpy (flxext2, ".atm.flx.std");
  
      strcpy (absext2, ".atm");
      strcat (absext2, extension2);
      break;
  
    case MCSRC_THERMAL_SURFACE:
      strcpy (radext, ".sur.rad");
      strcpy (rplext, ".sur.rpl");
      strcpy (flxext, ".sur.flx");
  
      strcpy (absext, ".sur");
      strcat (absext, extension);
  
      strcpy (radext2, ".sur.rad.std");
      strcpy (rplext2, ".sur.rpl.std");
      strcpy (flxext2, ".sur.flx.std");
  
      strcpy (absext2, ".sur");
      strcat (absext2, extension2);
      break;
    
    case MCSRC_SOLAR:
    case MCSRC_BLITZ:
    case MCSRC_THERMAL_BACKWARD:/*TZ bt ...*/
      strcpy (radext, ".rad");
      strcpy (rplext, ".rpl");
      strcpy (flxext, ".flx");
      
      strcpy (absext, extension);
  
      strcpy (radext2, ".rad.std");
      strcpy (rplext2, ".rpl.std");
      strcpy (flxext2, ".flx.std");
  
      strcpy (absext2, extension2);
        
      /* special treatment for backward: here we want the TOA radiance   */
      /* in a special file *.bac which contains the contribution of each */
      /* TOA pixel to the irradiance at zout                             */
  
      switch (sample->backward) {
      case MCBACKWARD_NONE:
        break;
        
      case MCBACKWARD_EDIR:
      case MCBACKWARD_EGLOB: 
      case MCBACKWARD_EDN:
      case MCBACKWARD_EUP:
      case MCBACKWARD_EXP:
      case MCBACKWARD_EXN:
      case MCBACKWARD_EYP:
      case MCBACKWARD_EYN:
      case MCBACKWARD_FDIR:
      case MCBACKWARD_FDN:
      case MCBACKWARD_FUP:
      case MCBACKWARD_ABS: 
      case MCBACKWARD_EMIS:
      case MCBACKWARD_HEAT:
      case MCBACKWARD_ACT:   
        strcpy (bacext,  ".flx");
        strcpy (bacext2, ".flx.std");
        
        strcpy (radext,  ".bac");
        strcpy (radext2, ".bac.std");
        break;
        
      case MCBACKWARD_RADIANCE:
        strcpy (bacext,  ".rad");
        strcpy (bacext2, ".rad.std");
        strcpy (picext,  ".ppm");
      
        strcpy (radext,  ".bac");
        strcpy (radext2, ".bac.std");
  
        if (sample->abs_jacobian) {
  	strcpy (bacjacext,  ".jac");
  	strcpy (bacjacext2, ".jac.std");
        }
  
        break;
        
      default:
        fprintf (stderr, "Error, unknown backward quantity %d\n", sample->backward);
        return -1;
      }
      break;
  
    case MCSRC_LIDAR:
      break;
  
    default:
      fprintf (stderr, "Error, unknown or unimplemented MC source type\n");
      return -1;
    }
  
    if (sample->spectral_is || sample->concentration_is){
      strcpy(specradfilename, basename);
      strcat(specradfilename, ".is.spc");
    }
  
    /* "absolute calibration" */
    switch (source) {
    case MCSRC_THERMAL_ATMOSPHERE:
      /* normalize to average column emission */
      if (absorption==MCFORWARD_ABS_EMISSION)
        incident = 1.0 / (4.0 * PI); /* integrate radiance over 4 pi */
      else 
        incident = 1.0 / atmos->Watm;
  
      break;
  
    case MCSRC_THERMAL_BACKWARD:/*TZ bt ...*/
      if (absorption==MCFORWARD_ABS_EMISSION)
        incident = 1.0 / (4.0 * PI); /* integrate radiance over 4 pi */
      else
        incident = 1.0;
   
      break;/*...TZ bt*/
      
    case MCSRC_THERMAL_SURFACE:
  
      /* normalize to average surface emission */
      incident = 1.0 / albedo->Wsurf;
      break;
        
    case MCSRC_SOLAR:
    case MCSRC_BLITZ:
  
      /* correct for slant incidence at TOA or BOA;                */
      /* use fabs() because cosine might be negative for backward; */
      /* the factor is the same for forward and backward because   */
      /* this is the cosine of the solar zenith angle in the       */
      /* definition of the transmittance outside mystic()          */                                     
      
      incident = fabs (1.0 / cosd(sample->forward_sza));
      break;
  
    case MCSRC_LIDAR:
      break;
  
    default:
      fprintf (stderr, "Error, unknown or unimplemented MC source type\n");
      return -1;
    }
  
    if (sample->backward || source == MCSRC_BLITZ) 
      factor = 1.0 / incident; 
    else
      factor = (double) sample->Nx * (double) sample->Ny / incident; 

    
    /* 3. Air Mass Factors! */
    if (sample->boxairmass){
      if (summarize_result_boxairmass (sample, atmos, result, basename,
				       islower, isupper, jslower, jsupper)) {
	fprintf (stderr, "Error in summarize_result_boxairmass().\n");
          return -1;
      }
    }
  
/* ======================================== Forward ===================================================== */
/* ====== 4. Forward Surface Irradiance, Radiance, RPL! */
/* (RPL = Radial Path Length) */
    if (!sample->backward) {
      if (summarize_result_forward_surface(sample, atmos, result, basename, dmcphotons, 
                                           factor, incident, write_files, elev2D, 
                                           flxext, flxext2, radext, radext2, rplext, 
                                           rplext2)) {
        fprintf (stderr, "Error in summarize_result_forward.\n");
        return -1;
      }
    }
  
/* ======================================== Backward ==================================================== */
    else {
      if (sample->Nd < 1) {
        fprintf (stderr, "Fatal error, need at least one radiance direction in backward mode!\n");
        fprintf (stderr, "Please inform the programmer!\n");
        return -1;
      }
      
      if (sample->backward_writeallpixels) {
        islower = 0;
        isupper = sample->Nx-1;
        jslower = 0;
        jsupper = sample->Ny-1;
      }
      else {
        islower = sample->backward_islower;
        isupper = sample->backward_isupper;
        jslower = sample->backward_jslower;
        jsupper = sample->backward_jsupper;
      }
  
      if(sample->backward != MCBACKWARD_RADIANCE) {
/* ====== 5. Backward Irradiance! */
        if (summarize_result_backward_irradiance(sample, atmos, result, basename, specradfilename,
                                                 bacext, bacext2, dmcphotons, factor,
                                                 wavelength, mcsimulations, write_files,
                                                 elev2D, islower, isupper, jslower, jsupper)) {
          fprintf (stderr, "Error in summarize_result_backward_irradiance.\n");
          return -1;
        }
      }
      else {
/* ====== 6. Cloudprop (backward)! */
  #ifdef CLDPRP
        if (sample->cldprp) {
          summarize_result_cldprp (sample, atmos, result, basename, elev2D, islower, isupper, 
                                   jslower, jsupper);
        }
  #endif
  
/* ====== 7. Backward Radiance! */
        if (summarize_result_backward_radiance(sample, atmos, result, basename, specradfilename,
                                               bacext, bacext2, picext, bacjacext, bacjacext2,
                                               dmcphotons, factor, sza, source, mcsimulations,
                                               write_files, elev2D, islower, isupper, jslower,
                                               jsupper)) {
            fprintf (stderr, "Error in summarize_result_backward_radiance.\n");
            return -1;
        }
      }
    }
  
    /* photons at altitude levels */
    for (kc=0; kc<=atmos->Nz; kc++) {
      if (sample->sample[kc]) {
        /* this loop is needed in forward mode; in backward mode it is only */
        /* needed if the user wants the .bac files (mc_writeback)           */
	/* ======================================== Forward ===================================================== */
	/* ====== 8. Forward Altitude Irradiance! */
        if (!sample->backward) {
          if (summarize_result_forward_altitude_irradiance (sample,  atmos,  result,  basename,
                                                            dmcphotons, factor, write_files, kc,
                                                            flxext, flxext2)) {
            fprintf (stderr, "Error in summarize_result_forward_altitude_irradiance.\n");
            return -1;
          }
        }

	/* ======================================== Forward (& Backward) ======================================== */
	/* ====== 9. Forward Altitude Radiance (Backward: bac files)! */
        if (!sample->backward || (sample->backward && sample->backward_writeback)) {
          if (summarize_result_forward_altitude_radiance(sample, atmos, result, 
                                                         basename, specradfilename,
                                                         radext, radext2, dmcphotons, 
                                                         factor, write_files,kc)) {
            fprintf (stderr, "Error in summarize_result_fw_alt_rad().\n");
            return -1;
          }
        }

	/* ======================================== Forward & Backward ========================================== */
	/* ====== 10. Altitude RPL! */
        if (summarize_result_altitude_rpl(sample, atmos, result, basename, rplext, rplext2, 
                                          dmcphotons, incident, write_files, kc)) {
            fprintf (stderr, "Error in summarize_result_alt_rpl().\n");
            return -1;
        }
      }
    }

    /* ======================================== Forward & Backward ========================================== */
    /* ====== 11. Passback3D! */
    lu3D=0;     /* counter for user levels */
    if (summarize_result_passback3D(sample, atmos, result, rfldir3d, rfldn3d, flup3d, 
                                    uavgso3d, uavgdn3d, uavgup3d, radiance3d, absback3d,
                                    rfldir3d_var, rfldn3d_var, flup3d_var, uavgso3d_var, 
                                    uavgdn3d_var, uavgup3d_var, radiance3d_var, absback3d_var,
                                    abs3d_var, radiance3d_is, escape, islower, isupper,
                                    jslower, jsupper, &lu3D)) {
      fprintf (stderr, "Error in summarize_result_passback3D.\n");
      return -1;
    }
    
   
  
/* ====== 12. Absorbed Energy! */
    if (absorption!=0) {
      if (summarize_result_absorption(sample,  atmos, result, abs3d,  abs3d_var, 
                                      basename, absext, absext2, dmcphotons,
                                      incident, write_files, absorption)) {
        fprintf (stderr, "Error in summarize_result_absorption().\n");
        return -1;
      }
    }
  
/* ====== 13. Information for the user! */
    if (!quiet) {
      fprintf (stderr, "\n\n");
      fprintf (stderr, "lower boundary     = %g m\n", atmos->Z[0]);
      fprintf (stderr, "upper boundary     = %g m\n", atmos->Z[atmos->Nz]);
      fprintf (stderr, "number of photons  = %ld\n",  mcphotons);
      fprintf (stderr, "solar zenith angle = %g\n",   sample->forward_sza);
      fprintf (stderr, "                                  \n");
      for (isp=1; isp<=atmos->n_caoth; isp++)
        fprintf (stderr, "total 1D %10s scattering       OD = %g\n",
  	       atmos->caoth_name[isp], atmos->tscatot[isp]);
      fprintf (stderr, "total scattering                OD = %g\n", atmos->tscatot[MCCAOTH_TOT]);
      fprintf (stderr, "                                  \n");
      for (isp=1; isp<=atmos->n_caoth; isp++)
        fprintf (stderr, "total 1D %10s absorption       OD = %g\n",
  	       atmos->caoth_name[isp], atmos->tabstot[isp]);
      fprintf (stderr, "total absorption                OD = %g\n", atmos->tabstot[MCCAOTH_TOT]);
    }
  
    lu=0;
      
    if (sample->surface) {
      area_average (result->surf, result->back, sample->Nx, sample->Ny,
                    sample->backward_islower, sample->backward_isupper, sample->backward_jslower, sample->backward_jsupper, 
                    &(ndir[0]),   &(ndn[0]),    &(nup[0]),
                    &(rfldir[0]), &(rfldn[0]),  &(flup[0]),
                    &(uavgso[0]), &(uavgdn[0]), &(uavgup[0]),
                    sample->backward, sample->surface_area, sample->surface_area_counter, elev2D, sample->surfaceparallel);
      lu++;
    }
  
  
    for (kc=0; kc<=atmos->Nz; kc++) { 
      if (sample->sample[kc]) {
        area_average (result->alt[kc], result->back, sample->Nx, sample->Ny,
                      sample->backward_islower, sample->backward_isupper, sample->backward_jslower, sample->backward_jsupper, 
                      &(ndir[lu]),   &(ndn[lu]),    &(nup[lu]),
                      &(rfldir[lu]), &(rfldn[lu]),  &(flup[lu]),
                      &(uavgso[lu]), &(uavgdn[lu]), &(uavgup[lu]),
                      sample->backward, sample->surface_area, sample->surface_area_counter, elev2D, sample->surfaceparallel);
        lu++;
      }
    }
    nu=lu;
    
    
    if (!quiet) {
      fprintf (stderr, "\ntotal number of (unweighted) photons:\n");
      fprintf (stderr, " user level,      direct, diffuse down,   diffuse up\n");
      for (lu=0; lu<nu; lu++)
        fprintf (stderr, "        %3d %12d  %12d  %12d\n", lu, ndir[lu], ndn[lu], nup[lu]);
      fprintf (stderr, "\n");
    }
  
    free (ndir); free (ndn); free (nup);
  
    if (!quiet)
      fprintf (stderr, " ... passing 3D fields for %d user levels back to uvspec()\n", lu3D);
  
  } /* End brace else-Lidar-loop */
  return 0;
}

/***********************************************************************************/
/* Function: print_sr_error_message                                       @62_30i@ */
/* Description: Prints out certain error messages for summarize_result             */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void print_sr_error_message(char* filename) {
  perror(NULL);
  fprintf (stderr, "errno = %d\n", errno);
  fprintf (stderr, "Error opening %s for writing\n", filename);
}

/***********************************************************************************/
/* Function: summarize_result_forward_surface                             @62_30i@ */
/* Description: Outputs forward surface quantities                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int summarize_result_forward_surface(sample_struct* sample, 
                                     atmosphere_struct* atmos,
                                     result_struct* result,
                                     char* basename,
                                     double dmcphotons,
                                     double factor,
                                     double incident,
                                     int write_files,
                                     int elev2D,
                                     char* flxext, char* flxext2,
                                     char* radext, char* radext2,
                                     char* rplext, char* rplext2) 
{
  int is=0,js=0; /* sampling grid coordinates */
  int id=0;      /* radiance index */
  int ip=0;      /* stokes component index */
  int iv=0;      /* wavelength index for spectral calculations with importance sampling */
  int ir=0,it=0;
  int ic=0; /* atmospheric grid coordinates - why no jc in most functions? */
  double var, rplfactor, areafactor, avg;
  double normomega=0;
  FILE *histfile   =NULL, *histfile2   =NULL;
  FILE *histradfile=NULL, *histradfile2=NULL;
  FILE *histrplfile=NULL, *histrplfile2=NULL;
  char histfilename    [FILENAME_MAX] ="", histfilename2    [FILENAME_MAX] ="";
  char histradfilename [FILENAME_MAX] ="", histradfilename2 [FILENAME_MAX] ="";
  char histrplfilename [FILENAME_MAX] ="", histrplfilename2 [FILENAME_MAX] ="";

  strcpy (histfilename,     basename);
  strcpy (histfilename2,    basename);
  strcpy (histradfilename,  basename);
  strcpy (histradfilename2, basename);
  strcpy (histrplfilename,  basename);
  strcpy (histrplfilename2, basename);

  strcat (histfilename,     flxext);
  strcat (histfilename2,    flxext2);
  strcat (histradfilename,  radext);
  strcat (histradfilename2, radext2);
  strcat (histrplfilename,  rplext);
  strcat (histrplfilename2, rplext2);
    
  if (write_files) {
        /* open histogram file */
    if ((histfile = fopen (histfilename, "w")) == NULL) {
            print_sr_error_message(histfilename);
      return -1;
    }
    /* open radiance histogram file */
    if (sample->Nd>0) {
      if ((histradfile = fopen (histradfilename, "w")) == NULL) {
              print_sr_error_message(histradfilename);
        return -1;
      }
    }
    /* open radiance/pathlength histogram file */
    if (sample->Nr>0 && sample->Nt>0 && sample->Nd>0) {
      if ((histrplfile = fopen (histrplfilename, "w")) == NULL) {
              print_sr_error_message(histrplfilename);
        return -1;
      }
    }
    if (sample->std) {  /* standard deviations */
      
        /* open histogram file */
      if ((histfile2 = fopen (histfilename2, "w")) == NULL) {
              print_sr_error_message(histfilename2);
        return -1;
      }
      /* open radiance histogram file */
      if (sample->Nd>0) {
        if ((histradfile2 = fopen (histradfilename2, "w")) == NULL) {
                print_sr_error_message(histradfilename2);
          return -1;
        }
      }
      /* open radiance/pathlength histogram file */
      if (sample->Nr>0 && sample->Nt>0 && sample->Nd>0) {
        if ((histrplfile2 = fopen (histrplfilename2, "w")) == NULL) {
                print_sr_error_message(histrplfilename2);
          return -1;
        }
      }
    }
  }

// =========================== Surface Irradiance
  for (is=0; is<sample->Nx; is++) {
    for (js=0; js<sample->Ny; js++) {
      if (elev2D && sample->surfaceparallel)
        areafactor = 1.0 / sample->surface_area[is][js];
      else 
        areafactor = 1.0;

      /* calculate average */
      result->surf->edir[is][js] /= dmcphotons;
      result->surf->edn [is][js] /= dmcphotons;
      result->surf->eup [is][js] /= dmcphotons;
      result->surf->fdir[is][js] /= dmcphotons;
      result->surf->fdn [is][js] /= dmcphotons;
      result->surf->fup [is][js] /= dmcphotons;
      
      /* calculate standard deviation from second moment */
      if (sample->std) {
        
        avg = result->surf->edir[is][js];
        var = result->surf->edir2[is][js] / dmcphotons;
        result->surf->edir2[is][js] = std_noisy (var, avg, dmcphotons);

        avg = result->surf->edn[is][js];
        var = result->surf->edn2[is][js] / dmcphotons;
        result->surf->edn2[is][js] = std_noisy (var, avg, dmcphotons);

        avg = result->surf->eup[is][js];
        var = result->surf->eup2[is][js] / dmcphotons;
        result->surf->eup2[is][js] = std_noisy (var, avg, dmcphotons);

        avg = result->surf->fdir[is][js];
        var = result->surf->fdir2[is][js] / dmcphotons;
        result->surf->fdir2[is][js] = std_noisy (var, avg, dmcphotons);

        avg = result->surf->fdn[is][js];
        var = result->surf->fdn2[is][js] / dmcphotons;
        result->surf->fdn2[is][js] = std_noisy (var, avg, dmcphotons);

        avg = result->surf->fup[is][js];
        var = result->surf->fup2[is][js] / dmcphotons;
        result->surf->fup2[is][js] = std_noisy (var, avg, dmcphotons);
 
        result->surf->edir2[is][js] *= factor * areafactor;
        result->surf->edn2 [is][js] *= factor * areafactor;
        result->surf->eup2 [is][js] *= factor * areafactor;
        result->surf->fdir2[is][js] *= factor * areafactor;
        result->surf->fdn2 [is][js] *= factor * areafactor;
        result->surf->fup2 [is][js] *= factor * areafactor;
      }
      
      result->surf->edir[is][js] *= factor * areafactor;
      result->surf->edn [is][js] *= factor * areafactor;
      result->surf->eup [is][js] *= factor * areafactor;
      result->surf->fdir[is][js] *= factor * areafactor;
      result->surf->fdn [is][js] *= factor * areafactor;
      result->surf->fup [is][js] *= factor * areafactor;
      
      
      if (write_files) {
        fprintf (histfile, "%g %g  %.6e %.6e %.6e %.6e %.6e %.6e\n",
           ((double) is + 0.5) * sample->delX,
           ((double) js + 0.5) * sample->delY,
           result->surf->edir[is][js], result->surf->edn [is][js], result->surf->eup [is][js],
           result->surf->fdir[is][js], result->surf->fdn [is][js], result->surf->fup [is][js]);
  
        if (sample->std) 
          fprintf (histfile2, "%g %g  %.6e %.6e %.6e %.6e %.6e %.6e\n",
             ((double) is + 0.5) * sample->delX,
             ((double) js + 0.5) * sample->delY,
             result->surf->edir2[is][js], result->surf->edn2 [is][js], result->surf->eup2 [is][js],
             result->surf->fdir2[is][js], result->surf->fdn2 [is][js], result->surf->fup2 [is][js]);
      }

// =========================== Surface Radiance
      for (id=0; id<sample->Nd; id++) {
        for (ip=0; ip<sample->nstokes; ip++) {
          
          /* calculate average */
          result->surf->raddir[id][is][js][ip] /= dmcphotons;
          result->surf->raddif[id][is][js][ip] /= dmcphotons;
          result->surf->radesc[id][is][js][ip] /= dmcphotons;
          
          
          /* calculate standard deviation from second moment */
          if (sample->std) {
            avg = result->surf->raddir[id][is][js][ip];
            var = result->surf->raddir2[id][is][js][ip] / dmcphotons;
            result->surf->raddir2[id][is][js][ip] = std_noisy (var, avg, dmcphotons);

            avg = result->surf->raddif[id][is][js][ip];
            var = result->surf->raddif2[id][is][js][ip] /dmcphotons;
            result->surf->raddif2[id][is][js][ip] = std_noisy (var, avg, dmcphotons);

            avg = result->surf->radesc[id][is][js][ip];
            var = result->surf->radesc2[id][is][js][ip] / dmcphotons;
            result->surf->radesc2[id][is][js][ip] = std_noisy (var, avg, dmcphotons);
 
            /* divide by solid angle to get radiance */
	    if (!sample->panorama_forward) {
	      result->surf->raddir2[id][is][js][ip] /= sample->rad[id].omega;
	      result->surf->raddif2[id][is][js][ip] /= sample->rad[id].omega;
	    }
	    else {
	      /* do forward panorama normalization here; divide by individual pixel solid angle */
	      /* and multiply with number of sample pixels sample->Nx * sample->Ny              */
	      normomega = (cos ((double) (js) / (double) sample->Ny * M_PI) - cos ((double) (js+1) / (double) sample->Ny * M_PI)) / sample->Nx * 2 * M_PI;
	      normomega *= sample->Nx * sample->Ny;

	      result->surf->raddir2[id][is][js][ip] /= normomega;
	      result->surf->raddif2[id][is][js][ip] /= normomega;
	    }

            result->surf->radesc2[id][is][js][ip] /= (4.0 * PI);
            result->surf->raddir2[id][is][js][ip] *= factor * areafactor;
            result->surf->raddif2[id][is][js][ip] *= factor * areafactor;
            result->surf->radesc2[id][is][js][ip] *= factor * areafactor;
          }

          /* divide by solid angle to get radiance */
	  if (!sample->panorama_forward) {
	    result->surf->raddir[id][is][js][ip] /= sample->rad[id].omega;
	    result->surf->raddif[id][is][js][ip] /= sample->rad[id].omega;
	  }
	  else {
	      /* do forward panorama normalization here; divide by individual pixel solid angle */
	      normomega = (cos ((double) (js) / (double) sample->Ny * M_PI) - cos ((double) (js+1) / (double) sample->Ny * M_PI)) / sample->Nx * 2 * M_PI;
	      normomega *= sample->Nx * sample->Ny;

	      result->surf->raddir[id][is][js][ip] /= normomega;
	      result->surf->raddif[id][is][js][ip] /= normomega;
	  }
	  result->surf->radesc[id][is][js][ip] /= (4.0 * PI);
          result->surf->raddir[id][is][js][ip] *= factor * areafactor;
          result->surf->raddif[id][is][js][ip] *= factor * areafactor;
          result->surf->radesc[id][is][js][ip] *= factor * areafactor;

          if (sample->spectral_is || sample->concentration_is){
            for (ic=0; ic<atmos->Nc; ic++){
              for (iv=0; iv<atmos->nlambda_abs; iv++){
                result->surf_t->rad_t[ic][is][js][ip][iv] /= dmcphotons;
                result->surf_t->rad_t[ic][is][js][ip][iv] /= (4.0 * PI);
                result->surf_t->rad_t[ic][is][js][ip][iv] *= factor * areafactor;
              }
            }
          }

          if (write_files) {
	    fprintf (histradfile, "%g %g  %.6e %.6e %.6e %.6e %.6e %.6e\n",
                     ((double) is + 0.5) * sample->delX,
                     ((double) js + 0.5) * sample->delY,
                     sample->rad[id].theta,
                     sample->rad[id].externalphi,
                     sample->rad[id].omega,
                     result->surf->raddir[id][is][js][ip],
                     result->surf->raddif[id][is][js][ip],
                     result->surf->radesc[id][is][js][ip]);
            
            if (sample->std)
              fprintf (histradfile2, "%g %g  %.6e %.6e %.6e %.6e %.6e %.6e\n",
                       ((double) is + 0.5) * sample->delX,
                       ((double) js + 0.5) * sample->delY,
                       sample->rad[id].theta,
                       sample->rad[id].externalphi,
                       sample->rad[id].omega,
                       result->surf->raddir2[id][is][js][ip],
                       result->surf->raddif2[id][is][js][ip],
                       result->surf->radesc2[id][is][js][ip]);
            
            /* CE - spectral radiance at surface, write to file ??? CHECK!!!
               for (iv=0; iv<atmos->nlambda_abs; iv++){
               fprintf(specradfile, "%d, %.6e",  iv, result->surf_t->rad_t[iv][ip]);
               }
            */
          }
        }
      }
    }
  }
  
  // =========================== Surface RPL
  for (ir=0; ir<sample->Nr; ir++) {
    
    /* normalize to area of circular ring and pathlength interval */
    rplfactor = 1.0 / (PI * sample->dr * sample->dr * (double) (2*ir+1)) / incident / (3.0E8 * sample->dt); 

    for (it=0; it<sample->Nt; it++) {
      for (id=0; id<sample->Nd; id++) {
        for (ip=0; ip<sample->nstokes; ip++){
          
          /* calculate average */
          result->surf->radpat[id][ir][it][ip] /= dmcphotons;
          result->surf->radpes[id][ir][it][ip] /= dmcphotons;
          
          /* calculate standard deviation from second moment */
          if (sample->std) {

            avg = result->surf->radpat[id][ir][it][ip];
            var = result->surf->radpat2[id][ir][it][ip] / dmcphotons;
            result->surf->radpat2[id][ir][it][ip] = std_noisy (var, avg, dmcphotons);

            avg = result->surf->radpes[id][ir][it][ip];
            var = result->surf->radpes2[id][ir][it][ip] / dmcphotons;
            result->surf->radpes2[id][ir][it][ip] = std_noisy (var, avg, dmcphotons);
 
            /* divide by solid angle to get radiance */
            result->surf->radpat2[id][ir][it][ip] /= sample->rad[id].omega;
            result->surf->radpat2[id][ir][it][ip] *= rplfactor;
            
            result->surf->radpes2[id][ir][it][ip] /= (4.0 * PI);
            result->surf->radpes2[id][ir][it][ip] *= rplfactor;
          }
          
          /* divide by solid angle to get radiance */
          result->surf->radpat[id][ir][it][ip] /= sample->rad[id].omega;
          result->surf->radpat[id][ir][it][ip] *= rplfactor;
          
          result->surf->radpes[id][ir][it][ip] /= (4.0 * PI);
          result->surf->radpes[id][ir][it][ip] *= rplfactor;
          
          if (write_files) {
            fprintf (histrplfile, "%4d %4d  %.6e %.6e %.6e %.6e %.6e %.6e\n",
                     ir, it,
                     sample->rad[id].theta,
                     sample->rad[id].externalphi,
                     sample->rad[id].omega,
                     0.0,
                     result->surf->radpat[id][ir][it][ip],
                     result->surf->radpes[id][ir][it][ip]);

            if (sample->std) {
              fprintf (histrplfile2, "%4d %4d  %.6e %.6e %.6e %.6e %.6e %.6e\n",
                       ir, it,
                       sample->rad[id].theta,
                       sample->rad[id].externalphi,
                       sample->rad[id].omega,
                       0.0,
                       result->surf->radpat2[id][ir][it][ip],
                       result->surf->radpes2[id][ir][it][ip]);
            }
          }
        }
      }
    }
  }

  /* close files */
  if (write_files) {
    (void) fclose (histfile);

    if (sample->Nd>0)
      (void) fclose (histradfile);
    
    if (sample->Nr>0 && sample->Nt>0 && sample->Nd>0)
      (void) fclose (histrplfile);
  
    if (sample->std) {
      (void) fclose (histfile2);
      
      if (sample->Nd>0)
        (void) fclose (histradfile2);
      
      if (sample->Nr>0 && sample->Nt>0 && sample->Nd>0)
        (void) fclose (histrplfile2);

    }
  }

  return 0;
}

/***********************************************************************************/
/* Function: summarize_result_forward_altitude_irradiance                 @62_30i@ */
/* Description: Outputs forward quantities                                         */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int summarize_result_forward_altitude_irradiance (sample_struct* sample, 
                                                  atmosphere_struct* atmos,
                                                  result_struct* result,
                                                  char* basename,
                                                  double dmcphotons,
                                                  double factor,
                                                  int write_files,
                                                  int kc,
                                                  char* flxext, char* flxext2)
{
  int is=0,js=0; /* sampling grid coordinates */
  double var, avg;
  FILE *altfile    =NULL, *altfile2    =NULL;
  char altfilename     [FILENAME_MAX] ="", altfilename2     [FILENAME_MAX] ="";

// =========================== Altitude Irradiance
  if (write_files) {

    sprintf (altfilename,    "%s%d%s", basename, kc, flxext);
    sprintf (altfilename2,    "%s%d%s", basename, kc, flxext2);

    if ((altfile = fopen(altfilename, "w")) == NULL) {
      print_sr_error_message(altfilename);
      return -1;
    }

    if (sample->std) {
      if ((altfile2 = fopen(altfilename2, "w")) == NULL) {
        print_sr_error_message(altfilename2);
        return -1;
      }
    }
  } 
  for (is=0; is<sample->Nx; is++) {
    for (js=0; js<sample->Ny; js++) {

      /* calculate average */
      result->alt[kc]->edir[is][js] /= dmcphotons;
      result->alt[kc]->edn [is][js] /= dmcphotons;
      result->alt[kc]->eup [is][js] /= dmcphotons;

      result->alt[kc]->fdir[is][js] /= dmcphotons;
      result->alt[kc]->fdn [is][js] /= dmcphotons;
      result->alt[kc]->fup [is][js] /= dmcphotons;

      /* calculate standard deviation from second moment */
      if (sample->std) {

        avg = result->alt[kc]->edir[is][js];
        var = result->alt[kc]->edir2[is][js] / dmcphotons;
        result->alt[kc]->edir2[is][js] = std_noisy (var, avg, dmcphotons);

        avg = result->alt[kc]->edn[is][js];
        var = result->alt[kc]->edn2[is][js] / dmcphotons;
        result->alt[kc]->edn2[is][js] = std_noisy (var, avg, dmcphotons);

        avg = result->alt[kc]->eup[is][js];
        var = result->alt[kc]->eup2[is][js] / dmcphotons;
        result->alt[kc]->eup2[is][js] = std_noisy (var, avg, dmcphotons);

        avg = result->alt[kc]->fdir[is][js];
        var = result->alt[kc]->fdir2[is][js] / dmcphotons;
        result->alt[kc]->fdir2[is][js] = std_noisy (var, avg, dmcphotons);

        avg = result->alt[kc]->fdn[is][js];
        var = result->alt[kc]->fdn2[is][js] / dmcphotons;
        result->alt[kc]->fdn2[is][js] = std_noisy (var, avg, dmcphotons);

        avg = result->alt[kc]->fup[is][js];
        var = result->alt[kc]->fup2[is][js] / dmcphotons;
        result->alt[kc]->fup2[is][js] = std_noisy (var, avg, dmcphotons);

        result->alt[kc]->edir2[is][js] *= factor;
        result->alt[kc]->edn2 [is][js] *= factor;
        result->alt[kc]->eup2 [is][js] *= factor;

        result->alt[kc]->fdir2[is][js] *= factor;
        result->alt[kc]->fdn2 [is][js] *= factor;
        result->alt[kc]->fup2 [is][js] *= factor;
      }

      result->alt[kc]->edir[is][js] *= factor;
      result->alt[kc]->edn [is][js] *= factor;
      result->alt[kc]->eup [is][js] *= factor;

      result->alt[kc]->fdir[is][js] *= factor;
      result->alt[kc]->fdn [is][js] *= factor;
      result->alt[kc]->fup [is][js] *= factor;

      if (write_files) {
        fprintf (altfile, "%g %g  %.6e %.6e %.6e %.6e %.6e %.6e\n",
            ((double) is + 0.5) * sample->delX, ((double) js + 0.5) * sample->delY,
            result->alt[kc]->edir[is][js], result->alt[kc]->edn [is][js],
            result->alt[kc]->eup [is][js], result->alt[kc]->fdir[is][js],
            result->alt[kc]->fdn [is][js], result->alt[kc]->fup [is][js]);

        if (sample->std) 
          fprintf (altfile2, "%g %g  %.6e %.6e %.6e %.6e %.6e %.6e\n",
              ((double) is + 0.5) * sample->delX, ((double) js + 0.5) * sample->delY,
              result->alt[kc]->edir2[is][js], result->alt[kc]->edn2 [is][js],
              result->alt[kc]->eup2 [is][js], result->alt[kc]->fdir2[is][js],
              result->alt[kc]->fdn2 [is][js], result->alt[kc]->fup2 [is][js]);
      }
    }
  }

  if (write_files) {
    (void) fclose (altfile);

    if (sample->std)
      (void) fclose (altfile2);
  }
  return 0;
}

/***********************************************************************************/
/* Function: write_backward_irradiance                                    @62_30i@ */
/* Description: Writes output for backward irradiance in summarize_result          */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void write_backward_irradiance(sample_struct* sample, double result,
                               FILE* backfile, int is, int js, int nn)
{
  int i;
  fprintf (backfile, "%g %g ", ((double) is + 0.5) * sample->delX, 
                               ((double) js + 0.5) * sample->delY);
  for(i=0;i<6;i++) {
    if (i == nn)
      fprintf(backfile, " %.6e",result);
    else
      fprintf(backfile, " %.6e",NOT_A_NUMBER);
  }
  fprintf(backfile, "\n");
}


/***********************************************************************************/
/* Function: summarize_result_backward_irradiance                         @62_30i@ */
/* Description: Outputs backward irradiance, spectral IS, Circ Rad, HEAT           */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int summarize_result_backward_irradiance(sample_struct* sample, 
                                         atmosphere_struct* atmos,
                                         result_struct* result,
                                         char* basename,
                                         char* specradfilename,
                                         char* bacext,
                                         char* bacext2,
                                         double dmcphotons,
                                         double factor,
                                         double wavelength,
                                         int mcsimulations,
                                         int write_files,
                                         int elev2D,
                                         int islower,
                                         int isupper,
                                         int jslower,
                                         int jsupper)
{
  int is=0,js=0; /* sampling grid coordinates */
  int ip=0;      /* stokes component index */
  int iv=0;      /* wavelength index for spectral calculations with importance sampling */
  int ic=0,il=0,i=0;
  double avg, var, areafactor;
  double norm = 1.;
  FILE *backfile   =NULL, *backfile2  =NULL;
  FILE *specradfile=NULL;
  FILE *circradfile=NULL;
  char backfilename    [FILENAME_MAX] ="", backfilename2    [FILENAME_MAX] ="";
  char circradfilename [FILENAME_MAX] ="";
  char circstr         [FILENAME_MAX] ="";

  if (sample->backward == MCBACKWARD_EDIR || 
      sample->backward == MCBACKWARD_EGLOB || 
      sample->backward == MCBACKWARD_EDN  || 
      sample->backward == MCBACKWARD_EUP  ||   
      sample->backward == MCBACKWARD_EXP  ||   
      sample->backward == MCBACKWARD_EXN  ||   
      sample->backward == MCBACKWARD_EYP  ||   
      sample->backward == MCBACKWARD_EYN  ||   
      sample->backward == MCBACKWARD_FDIR || 
      sample->backward == MCBACKWARD_FDN  || 
      sample->backward == MCBACKWARD_FUP)
    norm=4.0;
 
  sprintf (backfilename,  "%s%s%s", basename, sample->backward_altstr, bacext);
  sprintf (backfilename2, "%s%s%s", basename, sample->backward_altstr, bacext2);

  if (write_files) {
    if ((backfile = fopen (backfilename, "w")) == NULL) {
      print_sr_error_message(backfilename);
      return -1;
    }
          
    if (sample->std) {
      if ((backfile2 = fopen (backfilename2, "w")) == NULL) {
        print_sr_error_message(backfilename2);
        return -1;
      }
    }
  }
  
/* =========== Irradiance Factors */ 
  for (is=islower; is<=isupper; is++) {
    for (js=jslower; js<=jsupper; js++) {

      if (elev2D && sample->surfaceparallel)
        areafactor = 1.0 / sample->surface_area[is][js];
      else 
        areafactor = 1.0;

      for (ip=0; ip<sample->nstokes; ip++) {
  
        result->back[is][js][ip] /= (dmcphotons / (double) mcsimulations); /* **CK 2013.08.27 Add dmcphotons/mcsimulations*/
      
        if (sample->backward == MCBACKWARD_HEAT)
          if (sample->heat_flag[is][js] == MCBACKWARD_HEAT_EMABS || sample->heat_flag[is][js] == MCBACKWARD_HEAT_EMABSOPT)
            result->backemis[is][js] /= (dmcphotons / (double) mcsimulations);  /* **CK 2013.08.27 Add dmcphotons/mcsimulations*/
      
      
        if (sample->std) {
          if (sample->backward != MCBACKWARD_HEAT) { /* **CK 2013.09.27 */
            avg = result->back[is][js][ip];
            var = result->back2[is][js][ip] / (dmcphotons / (double) mcsimulations); 
          }
          else {
            switch (sample->heat_flag[is][js]) {
            case MCBACKWARD_HEAT_EMABS:
            case MCBACKWARD_HEAT_EMABSOPT:
              avg = result->back [is][js][ip] + result->backemis [is][js];
              var = result->back2[is][js][ip] / (dmcphotons / (double) mcsimulations); 
              break;
            case MCBACKWARD_HEAT_DENET:
              avg = 0;
              var = 0;
              
              if (atmos->threed[MCCAOTH_TOT][sample->backemis_kc]>=1) { /* **CK 1014.01.17*/
                for (i=0;i<12;i++) {
                  avg += result->back_dEnet[is][js][i]*result->back_dEnet[is][js][i];
                  var += result->back_dEnet2[is][js][i];
                }
              }
              else {
                for (i=0;i<4;i++){ 
                  avg += result->back_dEnet[is][js][i]*result->back_dEnet[is][js][i];
                  var += result->back_dEnet2[is][js][i];
                }
              }
              avg = sqrt(avg); 
              break;
          
            default:
              fprintf (stderr, "Error, unknown sample->heat_flag %d\n", sample->heat_flag[is][js]);
              return -1;
            }
          }
        
          if (sample->backward == MCBACKWARD_HEAT && sample->heat_flag[is][js]==MCBACKWARD_HEAT_DENET)  /* **CK 2014.01.17 */
            result->back2[is][js][ip] = std_noisy (var, avg, (dmcphotons / (double) mcsimulations)*(dmcphotons / (double) mcsimulations));
          else
            result->back2[is][js][ip] = std_noisy (var, avg, dmcphotons / (double) mcsimulations);  /* **CK 2013.09.27 Add dmcphotons/mcsimulations*/   
        
          result->back2[is][js][ip] *= factor / norm * areafactor;
            
        } 
       
        if (sample->backward !=  MCBACKWARD_HEAT){ /* **CK 2013.09.27 */
          result->backemis[is][js] *= (4.0 * PI);
          result->back[is][js][ip] *= factor/ norm * areafactor;
        }
        else {
          switch (sample->heat_flag[is][js]) {
          case MCBACKWARD_HEAT_EMABS:
          case MCBACKWARD_HEAT_EMABSOPT:
            result->backemis[is][js] *= factor/ norm * areafactor; 
            result->back[is][js][ip] *= factor/ norm * areafactor;
          break;
          case MCBACKWARD_HEAT_DENET:
            result->back[is][js][ip] *= factor/ norm * areafactor;
          break;
          default: 
            fprintf (stderr, "Error, unknown sample->heat_flag %d\n", sample->heat_flag[is][js]);
            return -1;
          }
        }
        
        //if(sample->spectral_is || sample->concentration_is) {
          //for (ic=0; ic<atmos->Nc; ic++) {
            //for (iv=0; iv<atmos->nlambda_abs; iv++) {
              //result->back_t[ic][is][js][ip][iv] /= (dmcphotons / (double) mcsimulations); /* BM 2013.10.27 - mcsimulations was missing */
              //result->back_t[ic][is][js][ip][iv] *= factor/ norm * areafactor;
            //} 
          //}
        //}   
      } /* ip */
    }  /* js */
  }  /* is */
	
  if (write_files) {

/* =========== Write Irradiance Files */ 
    for (is=islower; is<=isupper; is++) {
      for (js=jslower; js<=jsupper; js++) {
        switch (sample->backward) {
        case MCBACKWARD_EDIR:
                write_backward_irradiance(sample,result->back[is][js][0],backfile,is,js,0);
          break;
  
        case MCBACKWARD_EGLOB: 
        case MCBACKWARD_EDN:
        case MCBACKWARD_EXN:
        case MCBACKWARD_EYN:
                write_backward_irradiance(sample,result->back[is][js][0],backfile,is,js,1);
          break;
  
        case MCBACKWARD_EUP:
        case MCBACKWARD_EXP:
        case MCBACKWARD_EYP:
                write_backward_irradiance(sample,result->back[is][js][0],backfile,is,js,2);
          break;
                
        case MCBACKWARD_FDIR:
                write_backward_irradiance(sample,result->back[is][js][0],backfile,is,js,3);
          break;
                
        case MCBACKWARD_FDN:
                write_backward_irradiance(sample,result->back[is][js][0],backfile,is,js,4);
          break;
                
        case MCBACKWARD_FUP:
                write_backward_irradiance(sample,result->back[is][js][0],backfile,is,js,5);
          break;
                
        case MCBACKWARD_ABS:
        case MCBACKWARD_ACT:
        case MCBACKWARD_EMIS:
        case MCBACKWARD_HEAT:
          fprintf (stderr, "MCBACKWARD_ABS, MCBACKWARD_ACT etc not written to mc.flx\n");
          break;
  
        case MCBACKWARD_NONE:
          fprintf (stderr, "Error, we are not supposed to be here in a forward simulation.\n");
          return -1;
          break;
  
        default:
          fprintf (stderr, "Error, unknown backward quantity %d\n", sample->backward);
          return -1;
        }
          
        if (sample->std) {
          switch (sample->backward) {
          case MCBACKWARD_EDIR:
                  write_backward_irradiance(sample,result->back2[is][js][0],backfile2,is,js,0);
            break;
                  
          case MCBACKWARD_EGLOB: 
          case MCBACKWARD_EDN:
          case MCBACKWARD_EXN:
          case MCBACKWARD_EYN:
                  write_backward_irradiance(sample,result->back2[is][js][0],backfile2,is,js,1);
            break;
                  
          case MCBACKWARD_EUP:
          case MCBACKWARD_EXP:
          case MCBACKWARD_EYP:
                  write_backward_irradiance(sample,result->back2[is][js][0],backfile2,is,js,2);
            break;
                  
          case MCBACKWARD_FDIR:
                  write_backward_irradiance(sample,result->back2[is][js][0],backfile2,is,js,3);
            break;
                  
          case MCBACKWARD_FDN:
                  write_backward_irradiance(sample,result->back2[is][js][0],backfile2,is,js,4);
            break;
                  
          case MCBACKWARD_FUP:
                  write_backward_irradiance(sample,result->back2[is][js][0],backfile2,is,js,5);
            break;
                  
          case MCBACKWARD_ABS:
          case MCBACKWARD_ACT:
          case MCBACKWARD_EMIS:
          case MCBACKWARD_HEAT:
            fprintf (stderr, "MCBACKWARD_ABS, MCBACKWARD_ACT etc not written to mc.flx\n");
            break;
                  
          case MCBACKWARD_NONE:
            fprintf (stderr, "Error, we are not supposed to be here in a forward simulation.\n");
            return -1;
            break;
                  
          default:
            fprintf (stderr, "Error, unknown backward quantity %d\n", sample->backward);
            return -1;
          }
        }
      }
    }
    fclose (backfile);
    if (sample->std)
      fclose (backfile2);

/* =========== Spectral Importance Sampling */ 
    if (sample->spectral_is || sample->concentration_is) {
      if ((specradfile = fopen (specradfilename, "w")) == NULL) {
        print_sr_error_message(specradfilename);
        return -1;
      }
      for (ic=0; ic<atmos->Nc; ic++) {
        for (iv=0; iv<atmos->nlambda_abs; iv++) {
          for (is=islower; is<=isupper; is++) {
            for (js=jslower; js<=jsupper; js++) {
              if (elev2D && sample->surfaceparallel)
                areafactor = 1.0 / sample->surface_area[is][js];
              else 
                areafactor = 1.0;
              for (ip=0; ip<sample->nstokes; ip++) {
                result->back_t[ic][is][js][ip][iv] /= (dmcphotons / (double) mcsimulations); /* BM 2013.10.27 - mcsimulations was missing */
                result->back_t[ic][is][js][ip][iv] *= factor/ norm * areafactor;
              }
              fprintf(specradfile, "%.8f %d %d %.6e \n",  atmos->lambda[iv], is, js, 
                      result->back_t[ic][is][js][0][iv]);
            }
          }
        }
      }
      fclose(specradfile); 
    }
  }

/* =========== Circ Rad */
  if (sample->ncirc) {
    if (elev2D && sample->surfaceparallel)
      areafactor = 1.0 / sample->surface_area[0][0];
    else 
      areafactor = 1.0;

    for (il=0;il<10;il++) {
      for(ic=0;ic<sample->ncirc;ic++) {
        result->circcontr[il][ic] *= factor/norm * areafactor / dmcphotons;
      }
    }

    strcpy(circradfilename, basename);
    sprintf (circstr, ".%5.1f.circ", wavelength);
    strcat(circradfilename, circstr);
    if ((circradfile = fopen (circradfilename, "w")) == NULL) {
      print_sr_error_message(circradfilename);
      return -1;
    }
    for(ic=0;ic<sample->ncirc;ic++) {
      fprintf(circradfile, "%9.3f %12.3f", wavelength,sample->dcirc*((double)(ic+1)));
      for(il=0;il<10;il++) {
        fprintf(circradfile, " %12.4e",result->circcontr[il][ic]);
        fprintf(circradfile, "\n");
      }
    }
    fclose(circradfile); 
  }
      
  return 0;
}

/***********************************************************************************/
/* Function: summarize_result_cldprp                                      @62_30i@ */
/* Description: summarize_result output for cloud properties                       */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

#ifdef CLDPRP
int summarize_result_cldprp (sample_struct* sample, atmosphere_struct* atmos,
                              result_struct* result, char* basename, int elev2D,
                              int islower, int isupper, int jslower, int jsupper)
                              
{ 
  int kc=0, is=0, js=0;
  FILE *cldprp_file=NULL;
  char cldprp_filename[FILENAME_MAX] ="";
  /* double areafactor; */

  sprintf (cldprp_filename, "%s_%f", basename, result->wavelength);
  strcat(cldprp_filename, ".cldprp");

  for (kc=0; kc<=atmos->Nz; kc++) { 
    if (sample->sample[kc]) {
              
      if ((cldprp_file = fopen (cldprp_filename, "w")) == NULL) {
              print_sr_error_message(cldprp_filename);
        return -1;
      }
              
      for (is=islower; is<=isupper; is++) {
        for (js=jslower; js<=jsupper; js++) {  
                  
          /* normalize estimated effective radius with summed-up tau */
          if (result->alt[kc]->cldprp.tau_wc[is][js]>0.0)
            result->alt[kc]->cldprp.reff_wc[is][js] /= (result->alt[kc]->cldprp.tau_wc[is][js]);
          if (result->alt[kc]->cldprp.tau_ic[is][js]>0.0)
            result->alt[kc]->cldprp.reff_ic[is][js] /= (result->alt[kc]->cldprp.tau_ic[is][js]);
                        
          /* normalize with number of traced photons */
          if (result->alt[kc]->cldprp.totweights[is][js] > 0.0){
            result->alt[kc]->cldprp.tau_wc[is][js] /= (result->alt[kc]->cldprp.totweights[is][js]);
            result->alt[kc]->cldprp.tau_ic[is][js] /= (result->alt[kc]->cldprp.totweights[is][js]);              
          }
                        
          /* normalize with number of traced photons which flew through clouds */
          if (result->alt[kc]->cldprp.wc_weights[is][js] > 0.0){
            result->alt[kc]->cldprp.rhit_wc[is][js] /= (result->alt[kc]->cldprp.wc_weights[is][js]);
            result->alt[kc]->cldprp.dylwc[is][js] /= (result->alt[kc]->cldprp.wc_weights[is][js]);
            result->alt[kc]->cldprp.dylwc[is][js] /= (result->alt[kc]->cldprp.wc_weights[is][js]);              
            result->alt[kc]->cldprp.dzlwc[is][js] /= (result->alt[kc]->cldprp.wc_weights[is][js]);  
          }
                        
          if (result->alt[kc]->cldprp.ic_weights[is][js] > 0.0){
            result->alt[kc]->cldprp.rhit_ic[is][js] /= (result->alt[kc]->cldprp.ic_weights[is][js]);
            result->alt[kc]->cldprp.dxiwc[is][js] /= (result->alt[kc]->cldprp.ic_weights[is][js]);
            result->alt[kc]->cldprp.dyiwc[is][js] /= (result->alt[kc]->cldprp.ic_weights[is][js]);              
            result->alt[kc]->cldprp.dziwc[is][js] /= (result->alt[kc]->cldprp.ic_weights[is][js]);  
          }

/* why is this even here? */
/*          if (elev2D && sample->surfaceparallel) */ 
/*            areafactor = 1.0 / sample->surface_area[is][js]; */
/*          else */
/*            areafactor = 1.0;  */
                        
          fprintf (cldprp_file, "%4d %4d %7.3f %7.3f %7.3f %7.3f %8.3f %8.3f %+7.5f %+7.5f %+7.5f %+7.5f %+7.5f %+7.5f\n",
             is,
             js,
             result->alt[kc]->cldprp.reff_wc[is][js],
             result->alt[kc]->cldprp.reff_ic[is][js],
             result->alt[kc]->cldprp.rhit_wc[is][js],
             result->alt[kc]->cldprp.rhit_ic[is][js],
             result->alt[kc]->cldprp.tau_wc[is][js],
             result->alt[kc]->cldprp.tau_ic[is][js],
             result->alt[kc]->cldprp.dylwc[is][js],
             result->alt[kc]->cldprp.dylwc[is][js],
             result->alt[kc]->cldprp.dzlwc[is][js],
             result->alt[kc]->cldprp.dxiwc[is][js],
             result->alt[kc]->cldprp.dyiwc[is][js],
             result->alt[kc]->cldprp.dziwc[is][js]);
        }
      }
      fclose(cldprp_file);
    }
  }
  return 0;
}
#endif

/***********************************************************************************/
/* Function: summarize_result_backward_radiance                           @62_30i@ */
/* Description: Outputs backward radiance                                          */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int summarize_result_backward_radiance(sample_struct* sample, 
                                       atmosphere_struct* atmos,
                                       result_struct* result,
                                       char* basename,
                                       char* specradfilename,
                                       char* bacext,
                                       char* bacext2,
                                       char* picext,
                                       char* bacjacext,
                                       char* bacjacext2,
                                       double dmcphotons,
                                       double factor,
                                       double sza,
                                       int source,
                                       int mcsimulations,
                                       int write_files,
                                       int elev2D,
                                       int islower,
                                       int isupper,
                                       int jslower,
                                       int jsupper)
{

  int is=0,js=0,ic=0,ip=0,iv=0,isp=0,ijac=0,kc=0;
  double avg, var, areafactor;
  FILE *backfile   =NULL, *backfile2   =NULL;
  FILE *backjacfile=NULL, *backjacfile2=NULL;
  FILE *picturefile=NULL, *specradfile =NULL;;
  char backfilename    [FILENAME_MAX] ="", backfilename2    [FILENAME_MAX] ="";
  char backjacfilename [FILENAME_MAX] ="", backjacfilename2 [FILENAME_MAX] ="";
  char picturefilename [FILENAME_MAX] =""; 
  double picfac=0.;
  float xval=0.0,yval=0.0;
 
  if (write_files) {
    sprintf (backfilename,  "%s%s%s", basename, sample->backward_altstr, bacext);
    sprintf (backfilename2, "%s%s%s", basename, sample->backward_altstr, bacext2);
          
    if ((backfile = fopen (backfilename, "w")) == NULL) {
      print_sr_error_message(backfilename);
      return -1;
    }
          
    if (sample->std) {
      if ((backfile2 = fopen (backfilename2, "w")) == NULL) {
        print_sr_error_message(backfilename2);
        return -1;
      }
    }
  
    if (sample->abs_jacobian) {
      sprintf (backjacfilename,  "%s%s%s", basename, sample->backward_altstr, bacjacext);
      sprintf (backjacfilename2, "%s%s%s", basename, sample->backward_altstr, bacjacext2);
          
      if ((backjacfile = fopen (backjacfilename, "w")) == NULL) {
        print_sr_error_message(backjacfilename);
        return -1;
      }
          
      if (sample->std) {
        if ((backjacfile2 = fopen (backjacfilename2, "w")) == NULL) {
           print_sr_error_message(backjacfilename2);
          return -1;
        }
      }
    }
          
    if (sample->spectral_is || sample->concentration_is) {
      if ((specradfile = fopen (specradfilename, "w")) == NULL) {
        print_sr_error_message(specradfilename);
        return -1;
      }
    }
  }
/* ============ Backward Radiance */ 
  for (is=islower; is<=isupper; is++) {
    for (js=jslower; js<=jsupper; js++) {
  
      if (elev2D && sample->surfaceparallel)
        areafactor = 1.0 / sample->surface_area[is][js];
      else 
        areafactor = 1.0;
  
      for (ip=0; ip<sample->nstokes; ip++){
        result->back[is][js][ip] /= (dmcphotons / (double) mcsimulations); /* BM 2013.10.27 - mcsimulations was missing */
                
        if (sample->std) {
          avg = result->back [is][js][ip];
          var = result->back2[is][js][ip] / (dmcphotons / (double) mcsimulations); /* BM 2013.10.27 - mcsimulations was missing */
          result->back2[is][js][ip] = std_noisy (var, avg, dmcphotons / (double) mcsimulations);  /* BM 2013.10.27 - mcsimulations was missing */
          if (source==MCSRC_THERMAL_BACKWARD)  /*TZ bt*/
            result->back2[is][js][ip] *= factor * areafactor;/*TZ bt */
          else/*TZ bt*/
            result->back2[is][js][ip] *= factor / 4.0 / PI * areafactor;
        } 
                
        if (source==MCSRC_THERMAL_BACKWARD){  /*TZ bt*/
          result->back[is][js][ip] *= factor * areafactor;/*TZ bt*/
        }
        else{/*TZ bt*/
          result->back[is][js][ip] *= factor / 4.0 / PI * areafactor;
        }
      }
      if (sample->abs_jacobian) {
        for (isp=1; isp<=atmos->n_caoth; isp++) {
          for (ijac=0; ijac<2; ijac++) { /* scattering and absorption */
            for (kc=0; kc<atmos->Nz; kc++) {
    
              result->jacobian->jacobian_t[is][js][isp][ijac][kc] /= (dmcphotons / (double) mcsimulations); /* BM 2013.10.27 - mcsimulations was missing */
                    
              if (sample->std) {
                avg = result->jacobian->jacobian_t [is][js][isp][ijac][kc];
                var = result->jacobian->jacobian_t2[is][js][isp][ijac][kc] / (dmcphotons / (double) mcsimulations); /* BM 2013.10.27 - mcsimulations was missing */
                result->jacobian->jacobian_t2[is][js][isp][ijac][kc] = std_noisy (var, avg, dmcphotons / (double) mcsimulations); /* BM 2013.10.27 - mcsimulations was missing */
                if (source==MCSRC_THERMAL_BACKWARD)  /*TZ bt*/
                  result->jacobian->jacobian_t2[is][js][isp][ijac][kc] *= factor * areafactor;/*TZ bt */
                else/*TZ bt*/
                  result->jacobian->jacobian_t2[is][js][isp][ijac][kc] *= factor / 4.0 / PI * areafactor;
              }
        
              if (source==MCSRC_THERMAL_BACKWARD)  /*TZ bt*/
                result->jacobian->jacobian_t[is][js][isp][ijac][kc] *= factor * areafactor;/*TZ bt*/
              else/*TZ bt*/
                result->jacobian->jacobian_t[is][js][isp][ijac][kc] *= factor / 4.0 / PI * areafactor;
            } 
          }
        }
      }
    }
  }

/* ============ Write Files */
  if (write_files) {
    for (is=islower; is<=isupper; is++) {
      for (js=jslower; js<=jsupper; js++) {
        if (sample->panorama) {
          sample->backward_is=is;
          sample->backward_js=js;
          set_panorama_dir(sample,&yval,&xval);
          xval-=180.0;
          if (sample->pan_distr_photons_over_pixel) {
            xval+=0.5*sample->pan_dphi;
            yval+=0.5*sample->pan_dtheta;
          }
        }
        else {
          xval = ((double) is + 0.5) * sample->delX;
          yval = ((double) js + 0.5) * sample->delY;
        }
  
        for(ip=0; ip<sample->nstokes; ip++)
          fprintf (backfile, "%.8e %.8e  %.6e %.6e %.6e %.6e %.6e %.6e\n",
                   xval,
                   yval,
                   sza,
                   sample->rad[0].externalphi,
                   0.0,
                   0.0,
                   0.0,
                   result->back[is][js][ip]);
          
        if (sample->std)
          for(ip=0; ip<sample->nstokes; ip++)
            fprintf (backfile2, "%.8e %.8e  %.6e %.6e %.6e %.6e %.6e %.6e\n",
                     xval,
                     yval,
                     sample->forward_vza,
                     sample->forward_phi,
                     0.0,
                     0.0,
                     0.0,
                     result->back2[is][js][ip]);
          
        if (sample->abs_jacobian) {
          for (kc=0; kc<atmos->Nz; kc++) {
            fprintf (backjacfile, "%.8e %.8e %d %.6e %.6e %.6e",
                     xval,
                     yval,
                     kc,
                     atmos->Z[kc],
                     sza,
                     sample->rad[0].externalphi);
            for (isp=1; isp<=result->jacobian->n_caoth; isp++)
              fprintf (backjacfile, " %.6e %.6e",
              result->jacobian->jacobian_t[is][js][isp][0][kc],
              result->jacobian->jacobian_t[is][js][isp][1][kc]);
              fprintf (backjacfile, "\n");
          }
  
          if (sample->std) {
            for (kc=0; kc<atmos->Nz; kc++) {
              fprintf (backjacfile2, "%.8e %.8e %d %.6e %.6e %.6e",
                       xval,
                       yval,
                       kc,
                       atmos->Z[kc],
                       sza,
                       sample->rad[0].externalphi);
              for (isp=1; isp<=result->jacobian->n_caoth; isp++)
                fprintf (backjacfile2, " %.6e %.6e",
                result->jacobian->jacobian_t2[is][js][isp][0][kc],
                result->jacobian->jacobian_t2[is][js][isp][1][kc]);
              fprintf (backjacfile2, "\n");
            }
          }
        }
      }
    }
    fclose (backfile);
    if (sample->std)
      fclose (backfile2);
    if (sample->abs_jacobian) {
      fclose (backjacfile);
      if (sample->std)
        fclose (backjacfile2);
    }

/* ============= PanPicture */
    if (sample->pan_picture) {
      sprintf (picturefilename, "%s%s%s", basename, sample->backward_altstr, picext);
      if ((picturefile = fopen (picturefilename, "w")) == NULL) {
        print_sr_error_message(picturefilename);
        return -1;
      }
      fprintf(picturefile, "P3\n# mystic\n%d %d\n255\n",jsupper-jslower+1,isupper-islower+1);
      for (is=islower; is<=isupper; is++) {
        for (js=jslower; js<=jsupper; js++) {
          picfac = (picfac>result->back[is][js][0]?picfac:result->back[is][js][0]);
        }
      }
      picfac=255.0/picfac;
      for (is=islower; is<=isupper; is++) {
        for (js=jslower; js<=jsupper; js++) {
          fprintf ( picturefile, "%d %d %d\n",
                   (int) (result->back[is][js][0]*picfac),
                   (int) (result->back[is][js][0]*picfac),
                   (int) (result->back[is][js][0]*picfac) );
        }
      }
      fclose(picturefile);
    }
  }

/* ============ Spectral IS */  
  if (sample->spectral_is || sample->concentration_is) {
    for (is=islower; is<=isupper; is++) {
      for (js=jslower; js<=jsupper; js++) {
        if (elev2D && sample->surfaceparallel)
          areafactor = 1.0 / sample->surface_area[is][js];
        else 
          areafactor = 1.0;
        for (ip=0; ip<sample->nstokes; ip++){
          if(sample->spectral_is || sample->concentration_is){
            for (ic=0; ic<atmos->Nc; ic++){
              for (iv=0; iv<atmos->nlambda_abs; iv++){
                result->back_t[ic][is][js][ip][iv] /= (dmcphotons / (double) mcsimulations); /* BM 2013.10.27 - mcsimulations was missing */
              }
            }
          }
          if (source==MCSRC_THERMAL_BACKWARD){  /*TZ bt*/
            for (ic=0; ic<atmos->Nc; ic++) {
              for (iv=0; iv<atmos->nlambda_abs; iv++) {
                result->back_t[ic][iv][is][js][ip] *= factor * areafactor;
              }
            }
          }
          else{/*TZ bt*/
            for (ic=0; ic<atmos->Nc; ic++) {
              for (iv=0; iv<atmos->nlambda_abs; iv++) {
                result->back_t[ic][is][js][ip][iv] *= factor / 4.0 / PI * areafactor;
              }         
            }
          }
        }
      }
    }  
    if (write_files) {
      for (ic=0; ic<atmos->Nc; ic++){
        for (iv=0; iv<atmos->nlambda_abs; iv++){
          for (is=islower; is<=isupper; is++){
            for (js=jslower; js<=jsupper; js++){
              if(sample->nstokes==1) {
                fprintf(specradfile, "%.8f %d %d %.6e \n",  atmos->lambda[iv], is, js, 
                        result->back_t[ic][is][js][0][iv]);
              }
              else {
                fprintf(specradfile, "%.8f %d %d %.6e %.6e %.6e %.6e\n",    
                        atmos->lambda[iv], is, js, result->back_t[ic][is][js][0][iv], 
                        result->back_t[ic][is][js][1][iv], result->back_t[ic][is][js][2][iv], 
                        result->back_t[ic][is][js][3][iv]);
              }
            }
          }
        }
      }
        
    fclose(specradfile); 
    }
  }
  return 0;
}


/***********************************************************************************/
/* Function: summarize_result_forward_altitude_radiance                   @62_30i@ */
/* Description: Forward radiances at altitudes, backward bac files                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int summarize_result_forward_altitude_radiance(sample_struct* sample, 
                                               atmosphere_struct* atmos,
                                               result_struct* result,
                                               char* basename,
                                               char* specradfilename,
                                               char* radext,
                                               char* radext2,
                                               double dmcphotons,
                                               double factor,
                                               int write_files,
                                               int kc)
{

  int is=0,js=0,id=0,ip=0,ic=0,iv=0;
  double avg, var;
  FILE *altradfile=NULL, *altradfile2=NULL, *specradfile=NULL;
  char altradfilename[FILENAME_MAX] ="";
  char altradfilename2[FILENAME_MAX] ="";

  double normomega=0;
  
// =========================== Forward Altitude Radiance  -> radext = rad
// =========================== Backward Altitude Radiance -> radext = bac
  if (sample->Nd>0) {
    sprintf (altradfilename, "%s%d%s", basename, kc, radext);
    sprintf (altradfilename2, "%s%d%s", basename, kc, radext2);

    if (write_files) {
      if ((altradfile = fopen(altradfilename, "w")) == NULL) {
        print_sr_error_message(altradfilename);
        return -1;
      }

      if (sample->std) {
        if ((altradfile2 = fopen(altradfilename2, "w")) == NULL) {
          print_sr_error_message(altradfilename2);
          return -1;
        }
      }
    }
  }

  for (is=0; is<sample->Nx; is++) {
    for (js=0; js<sample->Ny; js++) {
      for (id=0; id<sample->Nd; id++) {
        for (ip=0; ip<sample->nstokes; ip++){

          /* calculate average */ 
          result->alt[kc]->raddir[id][is][js][ip] /= dmcphotons;
          result->alt[kc]->raddif[id][is][js][ip] /= dmcphotons;
          result->alt[kc]->radesc[id][is][js][ip] /= dmcphotons;

          /* calculate standard deviation from second moment */
          if (sample->std) {

            avg = result->alt[kc]->raddir[id][is][js][ip];
            var = result->alt[kc]->raddir2[id][is][js][ip] / dmcphotons;
            result->alt[kc]->raddir2[id][is][js][ip] = std_noisy (var, avg, dmcphotons);

            avg = result->alt[kc]->raddif[id][is][js][ip];
            var = result->alt[kc]->raddif2[id][is][js][ip] / dmcphotons;
            result->alt[kc]->raddif2[id][is][js][ip] = std_noisy (var, avg, dmcphotons);

            avg = result->alt[kc]->radesc[id][is][js][ip];
            var = result->alt[kc]->radesc2[id][is][js][ip] / dmcphotons;
            result->alt[kc]->radesc2[id][is][js][ip] = std_noisy (var, avg, dmcphotons);

            /* divide by solid angle to get radiance */
	    if (sample->panorama_forward) {
	      result->alt[kc]->raddir2[id][is][js][ip] /= sample->rad[id].omega;
	      result->alt[kc]->raddif2[id][is][js][ip] /= sample->rad[id].omega;
	    }
	    else {
	      /* do forward panorama normalization here; divide by individual pixel solid angle */
	      normomega = (cos ((double) (js) / (double) sample->Ny * M_PI) - cos ((double) (js+1) / (double) sample->Ny * M_PI)) / sample->Nx * 2 * M_PI;
	      normomega *= sample->Nx * sample->Ny;

	      result->alt[kc]->raddir2[id][is][js][ip] /= normomega;
	      result->alt[kc]->raddif2[id][is][js][ip] /= normomega;
	    }

	    result->alt[kc]->radesc2[id][is][js][ip] /= (4.0 * PI);

            result->alt[kc]->raddir2[id][is][js][ip] *= factor;
            result->alt[kc]->raddif2[id][is][js][ip] *= factor;
            result->alt[kc]->radesc2[id][is][js][ip] *= factor;

            switch (sample->backward) {
              case MCBACKWARD_EGLOB: 
              case MCBACKWARD_EDIR:
              case MCBACKWARD_EDN:
              case MCBACKWARD_EXP:
              case MCBACKWARD_EXN:
              case MCBACKWARD_EYP:
              case MCBACKWARD_EYN:
              case MCBACKWARD_EUP:
              case MCBACKWARD_FDIR:
              case MCBACKWARD_FDN:
              case MCBACKWARD_FUP:
              case MCBACKWARD_ABS:
              case MCBACKWARD_ACT:
                result->alt[kc]->raddir2[id][is][js][ip] *= PI;
                result->alt[kc]->raddif2[id][is][js][ip] *= PI;
                result->alt[kc]->radesc2[id][is][js][ip] *= PI;
                break;

              case MCBACKWARD_NONE:
              case MCBACKWARD_RADIANCE:
                break;

              default:
                fprintf (stderr, "Error, unknown backward quantity %d\n", sample->backward);
                return -1;
            }
          }

	  if (!sample->panorama_forward) {
	    result->alt[kc]->raddir[id][is][js][ip] /= sample->rad[id].omega;
	    result->alt[kc]->raddif[id][is][js][ip] /= sample->rad[id].omega;
	  }
	  else {
	      /* do forward panorama normalization here; divide by individual pixel solid angle */
	      normomega = (cos ((double) (js) / (double) sample->Ny * M_PI) - cos ((double) (js+1) / (double) sample->Ny * M_PI)) / sample->Nx * 2 * M_PI;
	      normomega *= sample->Nx * sample->Ny;

	      result->alt[kc]->raddir[id][is][js][ip] /= normomega;
	      result->alt[kc]->raddif[id][is][js][ip] /= normomega;
	  }
	  
	  result->alt[kc]->radesc[id][is][js][ip] /= (4.0 * PI);
	  
          result->alt[kc]->raddir[id][is][js][ip] *= factor;
          result->alt[kc]->raddif[id][is][js][ip] *= factor;
          result->alt[kc]->radesc[id][is][js][ip] *= factor;

          switch (sample->backward) {
            case MCBACKWARD_EGLOB: 
            case MCBACKWARD_EDIR:
            case MCBACKWARD_EDN:
            case MCBACKWARD_EXP:
            case MCBACKWARD_EXN:
            case MCBACKWARD_EYP:
            case MCBACKWARD_EYN:
            case MCBACKWARD_EUP:
            case MCBACKWARD_FDIR:
            case MCBACKWARD_FDN:
            case MCBACKWARD_FUP:
            case MCBACKWARD_ABS:
            case MCBACKWARD_ACT:  
              result->alt[kc]->raddir[id][is][js][ip] *= PI;
              result->alt[kc]->raddif[id][is][js][ip] *= PI;
              result->alt[kc]->radesc[id][is][js][ip] *= PI;
              break;

            case MCBACKWARD_NONE:
            case MCBACKWARD_RADIANCE:
              break;

            default:
              fprintf (stderr, "Error, unknown backward quantity %d\n", sample->backward);
              return -1;
          }

          if (write_files) {
            fprintf (altradfile, "%g %g  %.6e %.6e %.6e %.6e %.6e %.6e\n",
		     ((double) is + 0.5) * sample->delX, ((double) js + 0.5) * sample->delY,
		     sample->rad[id].theta, sample->rad[id].externalphi, sample->rad[id].omega, 
		     result->alt[kc]->raddir[id][is][js][ip],
		     result->alt[kc]->raddif[id][is][js][ip], result->alt[kc]->radesc[id][is][js][ip]);
	    
            if (sample->std)
              fprintf (altradfile2, "%g %g  %.6e %.6e %.6e %.6e %.6e %.6e\n",
		       ((double) is + 0.5) * sample->delX, ((double) js + 0.5) * sample->delY,
		       sample->rad[id].theta, sample->rad[id].externalphi, sample->rad[id].omega, 
		       result->alt[kc]->raddir2[id][is][js][ip],
		       result->alt[kc]->raddif2[id][is][js][ip], result->alt[kc]->radesc2[id][is][js][ip]);
          }
	} // end ip
      } // end id
    } // end js
  } // end is
  if (write_files) {
    if (sample->Nd>0)
      (void) fclose (altradfile);

    if (sample->std) {
      if (sample->Nd>0)
        (void) fclose (altradfile2);
    }
  }


/* ============ Spectral IS */  
  if (sample->spectral_is || sample->concentration_is){
    for (is=0; is<sample->Nx; is++) {
      for (js=0; js<sample->Ny; js++) {
        for (ip=0; ip<sample->nstokes; ip++){
          for (ic=0; ic<atmos->Nc; ic++){
            for (iv=0; iv<atmos->nlambda_abs; iv++){
              result->rad_t[kc]->rad_t[ic][is][js][ip][iv] /= (4.0 * PI); 
              result->rad_t[kc]->rad_t[ic][is][js][ip][iv] *= factor;
              result->rad_t[kc]->rad_t[ic][is][js][ip][iv] /= dmcphotons;
            }
          }
        }
      }
    }
    if (write_files) {
      if ((specradfile = fopen(specradfilename, "w")) == NULL) {
        print_sr_error_message(specradfilename);
        return -1;
      }
      for (is=0; is<sample->Nx; is++) {
        for (js=0; js<sample->Ny; js++) {
          for (ic=0; ic<atmos->Nc; ic++){
            for (iv=0; iv<atmos->nlambda_abs; iv++){
              if(sample->nstokes==1){
                fprintf(specradfile, "%.8f %d %d %.6e \n",  atmos->lambda[iv], 
                    is, js, 
                    result->rad_t[kc]->rad_t[ic][is][js][0][iv]);
              }
              else {
                fprintf(specradfile, "%.8f %d %d %.6e %.6e %.6e %.6e\n",    
                    atmos->lambda[iv], is, js, 
                    result->rad_t[kc]->rad_t[ic][is][js][0][iv], 
                    result->rad_t[kc]->rad_t[ic][is][js][1][iv],
                    result->rad_t[kc]->rad_t[ic][is][js][2][iv],
                    result->rad_t[kc]->rad_t[ic][is][js][3][iv]);
              }
            }
          }
        }
      }
      if (sample->Nd>0)
        (void) fclose (specradfile);
    }
  }

  return 0;
}

/***********************************************************************************/
/* Function: summarize_result_passback3D                                  @62_30i@ */
/* Description: Outputs 3D fields                                                  */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int summarize_result_passback3D(sample_struct* sample, 
                                atmosphere_struct* atmos,
                                result_struct* result,
                                float ***rfldir3d, float ***rfldn3d, float ***flup3d, 
                                float ***uavgso3d, float ***uavgdn3d, float ***uavgup3d, 
                                float ****radiance3d, float ***absback3d,
                                float ***rfldir3d_var, float ***rfldn3d_var, float ***flup3d_var, 
                                float ***uavgso3d_var, float ***uavgdn3d_var, float ***uavgup3d_var, 
                                float ****radiance3d_var, float ***absback3d_var,
                                float ***abs3d, float ******radiance3d_is,
                                int escape,
                                int islower,
                                int isupper,
                                int jslower,
                                int jsupper,
                                int* lu3D)
{
  int is=0,js=0,ip=0,kc=0,ic=0,iv=0;

/* ============ Forward */
  if (!sample->backward) {
    if (sample->passback3D && sample->surface) {
      for (is=0; is<sample->Nx; is++) {
        for (js=0; js<sample->Ny; js++) {
          rfldir3d [(*lu3D)][is][js] = (float) result->surf->edir [is][js];
          rfldn3d  [(*lu3D)][is][js] = (float) result->surf->edn  [is][js];
          flup3d   [(*lu3D)][is][js] = (float) result->surf->eup  [is][js];
          uavgso3d [(*lu3D)][is][js] = (float) result->surf->fdir [is][js];
          uavgdn3d [(*lu3D)][is][js] = (float) result->surf->fdn  [is][js];
          uavgup3d [(*lu3D)][is][js] = (float) result->surf->fup  [is][js];
          
          if (sample->Nd > 0) {
            for (ip=0; ip<sample->nstokes; ip++){
              if (escape) {
                radiance3d [(*lu3D)][is][js][ip] = (float) result->surf->radesc[0][is][js][ip];

		if (sample->std)
		  radiance3d_var [(*lu3D)][is][js][ip] = (float) result->surf->radesc2[0][is][js][ip];
	      }
              else {                                                       
                radiance3d [(*lu3D)][is][js][ip] = (float) result->surf->raddir[0][is][js][ip] +
                  (float) result->surf->raddif[0][is][js][ip];

		if (sample->std)
		  radiance3d_var [(*lu3D)][is][js][ip] = (float) result->surf->raddir2[0][is][js][ip] +
		    (float) result->surf->raddif2[0][is][js][ip];
	      }
	    }
          }

	  /* variance */
	  if (sample->std) {
	    rfldir3d_var [(*lu3D)][is][js] = (float) (result->surf->edir2 [is][js] * result->surf->edir2 [is][js]);
	    rfldn3d_var  [(*lu3D)][is][js] = (float) (result->surf->edn2  [is][js] * result->surf->edn2  [is][js]);
	    flup3d_var   [(*lu3D)][is][js] = (float) (result->surf->eup2  [is][js] * result->surf->eup2  [is][js]);
	    uavgso3d_var [(*lu3D)][is][js] = (float) (result->surf->fdir2 [is][js] * result->surf->fdir2 [is][js]);
	    uavgdn3d_var [(*lu3D)][is][js] = (float) (result->surf->fdn2  [is][js] * result->surf->fdn2  [is][js]);
	    uavgup3d_var [(*lu3D)][is][js] = (float) (result->surf->fup2  [is][js] * result->surf->fup2  [is][js]);


	    if (sample->Nd > 0) {
	      for (ip=0; ip<sample->nstokes; ip++){

		if (escape)
		  radiance3d_var [(*lu3D)][is][js][ip] = (float) (result->surf->radesc2[0][is][js][ip] * 
							       result->surf->radesc2[0][is][js][ip]);
		else
		  radiance3d_var [(*lu3D)][is][js][ip] = (float) ((result->surf->raddir2[0][is][js][ip] + result->surf->raddif2[0][is][js][ip]) * 
							       (result->surf->raddir2[0][is][js][ip] + result->surf->raddif2[0][is][js][ip]));
	      }
	    }
	  }
        } /* for (js=0; is<sample->Ny; ... */
      } /* for (is=0; is<sample->Nx; ... */

      (*lu3D)++;
    }
  }
/* ============ Backward */
  else { 
    if (sample->passback3D) {
      if (sample->Nd < 1) {
        fprintf (stderr, "Fatal error, need at least one radiance direction in backward mode!\n");
        fprintf (stderr, "Please inform the programmer!\n");
        return -1;
      }
      
      for (is=islower; is<=isupper; is++) {
        for (js=jslower; js<=jsupper; js++) {
          
          /* initialize with NAN */
          rfldir3d  [(*lu3D)][is][js] = 0.0 / 0.0;
          rfldn3d   [(*lu3D)][is][js] = 0.0 / 0.0;
          flup3d    [(*lu3D)][is][js] = 0.0 / 0.0;
          uavgso3d  [(*lu3D)][is][js] = 0.0 / 0.0;
          uavgdn3d  [(*lu3D)][is][js] = 0.0 / 0.0;
          uavgup3d  [(*lu3D)][is][js] = 0.0 / 0.0;

          if (sample->backward == MCBACKWARD_ABS || sample->backward == MCBACKWARD_ACT)
            absback3d [(*lu3D)][is][js] = 0.0 / 0.0;

	  /* variances */
	  if (sample->std) {
	    rfldir3d_var [(*lu3D)][is][js] = 0.0 / 0.0;
	    rfldn3d_var  [(*lu3D)][is][js] = 0.0 / 0.0;
	    flup3d_var   [(*lu3D)][is][js] = 0.0 / 0.0;
	    uavgso3d_var [(*lu3D)][is][js] = 0.0 / 0.0;
	    uavgdn3d_var [(*lu3D)][is][js] = 0.0 / 0.0;
	    uavgup3d_var [(*lu3D)][is][js] = 0.0 / 0.0;
	  }

          /* we want only numbers if something was calculated; otherwise NAN */
          if (is>=sample->backward_islower && is<=sample->backward_isupper &&
              js>=sample->backward_jslower && js<=sample->backward_jsupper) {
            switch (sample->backward) {
	    case MCBACKWARD_EGLOB: 
            case MCBACKWARD_EDIR:
              rfldir3d [(*lu3D)][is][js]     = (float) result->back[is][js][0];

	      if (sample->std)
		rfldir3d_var [(*lu3D)][is][js] = (float) (result->back2[is][js][0] * result->back2[is][js][0]);
	      if (sample->backward==MCBACKWARD_EDIR)
		break;
              
            case MCBACKWARD_EDN:
            case MCBACKWARD_EXN:
            case MCBACKWARD_EYN:
              rfldn3d  [(*lu3D)][is][js] = (float) result->back[is][js][0];

	      if (sample->std)
		rfldn3d_var [(*lu3D)][is][js] = (float) (result->back2[is][js][0] * result->back2[is][js][0]);
              break;
              
            case MCBACKWARD_EUP:
            case MCBACKWARD_EXP:
            case MCBACKWARD_EYP:
              flup3d   [(*lu3D)][is][js] = (float) result->back[is][js][0];

	      if (sample->std)
		flup3d_var [(*lu3D)][is][js] = (float) (result->back2[is][js][0] * result->back2[is][js][0]);
              break;
              
            case MCBACKWARD_FDIR:
              uavgso3d [(*lu3D)][is][js] = (float) result->back[is][js][0];

	      if (sample->std)
		uavgso3d_var [(*lu3D)][is][js] = (float) (result->back2[is][js][0] * result->back2[is][js][0]);
              break;
              
            case MCBACKWARD_FDN:
              uavgdn3d [(*lu3D)][is][js] = (float) result->back[is][js][0];

	      if (sample->std)
		uavgdn3d_var [(*lu3D)][is][js] = (float) (result->back2[is][js][0] * result->back2[is][js][0]);
              break;
              
            case MCBACKWARD_FUP:
              uavgup3d [(*lu3D)][is][js] = (float) result->back[is][js][0];

	      if (sample->std)
		uavgup3d_var [(*lu3D)][is][js] = (float) (result->back2[is][js][0] * result->back2[is][js][0]);
              break;
              
            case MCBACKWARD_ABS:
            case MCBACKWARD_ACT:
              absback3d  [(*lu3D)][is][js]  = (float) result->back[is][js][0];

	      if (sample->std)
		absback3d_var [(*lu3D)][is][js] = (float) (result->back2[is][js][0] * result->back2[is][js][0]);
              break;
              
            case MCBACKWARD_EMIS:
              absback3d  [(*lu3D)][is][js]  = (float) result->backemis[is][js];

	      if (sample->std)
		absback3d_var [(*lu3D)][is][js] = (float) 0.0;   /* standard deviation = 0 because analytical */
              break;


            case MCBACKWARD_HEAT: /* **CK 2013.09.27 */
	      switch (sample->heat_flag[is][js]) {
	      case MCBACKWARD_HEAT_EMABS:
	      case MCBACKWARD_HEAT_EMABSOPT:
		absback3d  [(*lu3D)][is][js]  = (float) result->backemis[is][js] + (float) result->back[is][js][0]; 
		break;
	      case MCBACKWARD_HEAT_DENET:
		absback3d  [(*lu3D)][is][js]  = (float) result->back[is][js][0]; 
		break;
	      }
	      
	      if (sample->std)
		absback3d_var [(*lu3D)][is][js] =(float) (result->back2[is][js][0] * result->back2[is][js][0]);   /* 27.02.2013 **CK **BM: Calculation for thermal heating rates std */
              break;

              
            case MCBACKWARD_RADIANCE:
              for (ip=0; ip<sample->nstokes; ip++) {
                radiance3d [(*lu3D)][is][js][ip] = (float) result->back[is][js][ip];
                
                if (sample->spectral_is || sample->concentration_is){
                  for (ic=0; ic<atmos->Nc; ic++){
                    for (iv=0; iv<atmos->nlambda_abs; iv++){
                      radiance3d_is [(*lu3D)][ic][is][js][ip][iv] = (float) result->back_t[ic][is][js][ip][iv];
                    }
                  }
                }
                
                if (sample->std)
		  radiance3d_var [(*lu3D)][is][js][ip] = (float) (result->back2[is][js][ip] * result->back2[is][js][ip]);
	      }

              break;
              
            default:
              fprintf (stderr, "Error, unknown backward quantity %d\n", sample->backward);
              return -1;
            }
          }
        }  /* for (js=0; js<sample->Ny; ... */
      }  /* for (is=0; is<sample->Nx; ... */

      (*lu3D)++;
    }
  }


  /* ============ Write Results */
  for (kc=0; kc<=atmos->Nz; kc++) {
    if (sample->sample[kc]) {
      if (!sample->backward) {
        
        /* no need to consider backward here - there is only one backward */
        /* level which has been considered above.                         */
        if (sample->passback3D) {
          for (is=0; is<sample->Nx; is++) {
            for (js=0; js<sample->Ny; js++) {
              rfldir3d [(*lu3D)][is][js] = (float) result->alt[kc]->edir[is][js];
              rfldn3d  [(*lu3D)][is][js] = (float) result->alt[kc]->edn [is][js];
              flup3d   [(*lu3D)][is][js] = (float) result->alt[kc]->eup [is][js];
              uavgso3d [(*lu3D)][is][js] = (float) result->alt[kc]->fdir[is][js];
              uavgdn3d [(*lu3D)][is][js] = (float) result->alt[kc]->fdn [is][js];
              uavgup3d [(*lu3D)][is][js] = (float) result->alt[kc]->fup [is][js];
              

              for (ip=0; ip<sample->nstokes; ip++){
                if (sample->Nd > 0) {
                  if (escape)
                    radiance3d [(*lu3D)][is][js][ip] = (float) result->alt[kc]->radesc[0][is][js][ip];
                  else                                                  
                    radiance3d [(*lu3D)][is][js][ip] = (float) result->alt[kc]->raddir[0][is][js][ip] +
                      (float) result->alt[kc]->raddif[0][is][js][ip];
                  if (sample->spectral_is || sample->concentration_is){
                    for (ic=0; ic<atmos->Nc; ic++){
                      for (iv=0; iv<atmos->nlambda_abs; iv++){
                        radiance3d_is [(*lu3D)][ic][is][js][ip][iv] = 
                          result->rad_t[kc]->rad_t[ic][is][js][ip][iv];
                      }
                    } 
                  }
                }
              }

	      /* variances */
	      if (sample->std) {
		rfldir3d_var [(*lu3D)][is][js] = (float) (result->alt[kc]->edir2[is][js] * result->alt[kc]->edir2[is][js]);
		rfldn3d_var  [(*lu3D)][is][js] = (float) (result->alt[kc]->edn2 [is][js] * result->alt[kc]->edn2 [is][js]);
		flup3d_var   [(*lu3D)][is][js] = (float) (result->alt[kc]->eup2 [is][js] * result->alt[kc]->eup2 [is][js]);
		uavgso3d_var [(*lu3D)][is][js] = (float) (result->alt[kc]->fdir2[is][js] * result->alt[kc]->fdir2[is][js]);
		uavgdn3d_var [(*lu3D)][is][js] = (float) (result->alt[kc]->fdn2 [is][js] * result->alt[kc]->fdn2 [is][js]);
		uavgup3d_var [(*lu3D)][is][js] = (float) (result->alt[kc]->fup2 [is][js] * result->alt[kc]->fup2 [is][js]);
		
		for (ip=0; ip<sample->nstokes; ip++){
		  if (sample->Nd > 0) {
		    if (escape)
		      radiance3d_var [(*lu3D)][is][js][ip] = (float) (result->alt[kc]->radesc2[0][is][js][ip] * result->alt[kc]->radesc2[0][is][js][ip]);
		    else                                                  
		      radiance3d_var [(*lu3D)][is][js][ip] = (float) ((result->alt[kc]->raddir2[0][is][js][ip] + result->alt[kc]->raddif2[0][is][js][ip]) * 
								   (result->alt[kc]->raddir2[0][is][js][ip] + result->alt[kc]->raddif2[0][is][js][ip]));
		  }
		}
	      }

            }  /* for (js=0; js<sample->Ny; ... */
          }  /* for (is=0; is<sample->Nx; ... */
          
          (*lu3D)++;
        }
      }
    }  
  } // end Photons at height levels!
  return 0; 
}

/***********************************************************************************/
/* Function: summarize_result_alt_rpl                                     @62_30i@ */
/* Description: Altitude RPL                                                       */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int summarize_result_altitude_rpl(sample_struct* sample, 
                                  atmosphere_struct* atmos,
                                  result_struct* result,
                                  char* basename,
                                  char* rplext,
                                  char* rplext2,
                                  double dmcphotons,
                                  double incident,
                                  int write_files,
                                  int kc)
{
  int ir=0,it=0,id=0,ip=0;
  double rplfactor=0.;
  char altrplfilename  [FILENAME_MAX] ="", altrplfilename2  [FILENAME_MAX] ="";
  FILE *altrplfile =NULL, *altrplfile2 =NULL;
  double avg=0.,var=0.;

      if (sample->Nr>0 && sample->Nt>0) {
        
        if (kc!=2)
          fprintf (stderr, "KC = %d\n", sample->sample[kc]);
        
        if (write_files) {
          sprintf (altrplfilename, "%s%d%s", basename, kc, rplext);
          sprintf (altrplfilename2, "%s%d%s", basename, kc, rplext2);
          
          if (sample->Nd>0)
            if ((altrplfile = fopen(altrplfilename, "w")) == NULL) {
              print_sr_error_message(altrplfilename);
              return -1;
            }
        
          if (sample->std)
            if (sample->Nd>0)
              if ((altrplfile2 = fopen(altrplfilename2, "w")) == NULL) {
                print_sr_error_message(altrplfilename2);
                return -1;
              }
        }

        for (ir=0; ir<sample->Nr; ir++) {
          
          /* area of circular ring */
          rplfactor = 1.0 / (PI * sample->dr * sample->dr * (double) (2*ir+1)) / incident / (3.0E8 * sample->dt); 
          
          for (it=0; it<sample->Nt; it++) {
            for (id=0; id<sample->Nd; id++) {
              
              for (ip=0; ip<sample->nstokes; ip++){
                
                /* calculate average */ 
                result->alt[kc]->radpat[id][ir][it][ip] /= dmcphotons;
                result->alt[kc]->radpes[id][ir][it][ip] /= dmcphotons;
              
                /* calculate standard deviation from second moment */
                if (sample->std) {
                  
                  avg = result->alt[kc]->radpat[id][ir][it][ip];
                  var = result->alt[kc]->radpat2[id][ir][it][ip] / dmcphotons;
                  result->alt[kc]->radpat2[id][ir][it][ip] = std_noisy (var, avg, dmcphotons);
                  
                  avg = result->alt[kc]->radpes[id][ir][it][ip];
                  var = result->alt[kc]->radpes2[id][ir][it][ip] / dmcphotons;
                  result->alt[kc]->radpes2[id][ir][it][ip] = std_noisy (var, avg, dmcphotons);
                  
                  /* divide by solid angle to get radiance */
                  result->alt[kc]->radpat2[id][ir][it][ip] /= sample->rad[id].omega;
                  result->alt[kc]->radpat2[id][ir][it][ip] *= rplfactor;

                  result->alt[kc]->radpes2[id][ir][it][ip] /= (4.0 * PI);
                  result->alt[kc]->radpes2[id][ir][it][ip] *= rplfactor;
                }
                
                /* divide by solid angle to get radiance */
                result->alt[kc]->radpat[id][ir][it][ip] /= sample->rad[id].omega;
                result->alt[kc]->radpat[id][ir][it][ip] *= rplfactor;
                
                result->alt[kc]->radpes[id][ir][it][ip] /= (4.0 * PI);
                result->alt[kc]->radpes[id][ir][it][ip] *= rplfactor;
                
                if (write_files) {
                  fprintf (altrplfile, "%4d %4d  %.6e %.6e %.6e %.6e %.6e %.6e\n",
                           ir, it, sample->rad[id].theta,
                           sample->rad[id].externalphi, sample->rad[id].omega, 0.0,
                           result->alt[kc]->radpat[id][ir][it][ip], result->alt[kc]->radpes[id][ir][it][ip]);
                
                  if (sample->std)
                    fprintf (altrplfile2, "%4d %4d  %.6e %.6e %.6e %.6e %.6e %.6e\n",
                             ir, it, sample->rad[id].theta,
                             sample->rad[id].externalphi, sample->rad[id].omega, 0.0,
                             result->alt[kc]->radpat2[id][ir][it][ip], result->alt[kc]->radpes2[id][ir][it][ip]);
                }
              }
            }
          }
        }
        
        if (write_files) {
          if (sample->Nd>0)
            (void) fclose (altrplfile);
        
          if (sample->Nd>0) {
            if (sample->std)
              (void) fclose (altrplfile2);
          }
        }
      }
  return 0;
}

/***********************************************************************************/
/* Function: summarize_result_absorption                                  @62_30i@ */
/* Description: Absorption                                                         */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int summarize_result_absorption(sample_struct* sample, 
                                atmosphere_struct* atmos,
                                result_struct* result,
                                float ***abs3d, 
                                float ***abs3d_var, 
                                char* basename,
                                char* absext,
                                char* absext2,
                                double dmcphotons,
                                double incident,
                                int write_files,
                                int absorption)
{
  int kc=0,ic=0,jc=0;
  double tempabs=0.;
  double var=0.,avg=0.;
  FILE *altfile    =NULL, *altfile2    =NULL;
  char altfilename     [FILENAME_MAX] ="", altfilename2     [FILENAME_MAX] ="";


    /* multiplication factor to consider the area of a single pixel */
    /* attention: absorption is defined on the 3D caoth grid, not   */
    /* the sample grid!                                             */

    if (absorption == MCFORWARD_ABS_EMISSION)
      tempabs = 1;
    else 
      tempabs = (double) (atmos->Nx * atmos->Ny);

    /* photons at altitude levels */
    for (kc=0; kc<atmos->Nz; kc++) { 
      if (atmos->threed[MCCAOTH_TOT][kc]>=1) {
          
        if (write_files) {
          sprintf (altfilename, "%s%d%s", basename, kc, absext);
          sprintf (altfilename2, "%s%d%s", basename, kc, absext2);
          
          if ((altfile = fopen(altfilename, "w")) == NULL) {
            print_sr_error_message(altfilename);
            return -1;
          }
          
          
          if (sample->std) {
            if ((altfile2 = fopen(altfilename2, "w")) == NULL) {
              print_sr_error_message(altfilename2);
              return -1;
            }
          }
        }

        for (ic=0; ic<atmos->Nx; ic++) {
          for (jc=0; jc<atmos->Ny; jc++) {
            
            if (absorption!=MCFORWARD_ABS_EMISSION) {
                
              /* calculate average */
              result->absorption3D->tot[kc][ic][jc] /= dmcphotons;
            
              /* calculate standard deviation from second moment */
              if (sample->std) {
                avg = result->absorption3D->tot[kc][ic][jc];
                var = result->absorption3D2->tot[kc][ic][jc] / dmcphotons;
                result->absorption3D2->tot[kc][ic][jc] = std_noisy (var, avg, dmcphotons);
                

                result->absorption3D2->tot[kc][ic][jc] /= incident;
              }
            }

            result->absorption3D->tot[kc][ic][jc] /= incident;
            
            if (absorption == MCFORWARD_ABS_ACTINIC) {  /* dann sagt Rolle: */
              result->absorption3D->tot[kc][ic][jc] /=
		( atmos->kabs3D->prof [MCCAOTH_TOT][kc][ic][jc] *
		  ( atmos->Z[kc+1] - atmos->Z[kc] ) );

              if (sample->std)
		result->absorption3D2->tot[kc][ic][jc] /=
		  ( atmos->kabs3D->prof [MCCAOTH_TOT][kc][ic][jc] *
		    ( atmos->Z[kc+1]-atmos->Z[kc] ) );
            }
            
            if (write_files) {
              fprintf (altfile, "%g %g %g\n",
                       ((double) ic + 0.5) * atmos->delX, ((double) jc + 0.5) * atmos->delY,
                       result->absorption3D->tot[kc][ic][jc]*tempabs);
              
              if (sample->std)
                fprintf (altfile2, "%g %g %g\n",
                         ((double) ic + 0.5) * atmos->delX, ((double) jc + 0.5) * atmos->delY,
                         result->absorption3D2->tot[kc][ic][jc]*tempabs);
            }         
          }
        }

        if (write_files) {
          (void) fclose (altfile);
          
          if (sample->std)
            (void) fclose (altfile2);
        }
      }
    }

    /* copy to result array */
    if (sample->passback3D) {
      for (kc=0; kc<atmos->Nz; kc++) {
        for (ic=0; ic<atmos->Nx; ic++) {
          for (jc=0; jc<atmos->Ny; jc++) {
            if (atmos->threed[MCCAOTH_TOT][kc]>=1) { /* **CK added bracket */
              abs3d [kc][ic][jc] = result->absorption3D->tot[kc][ic][jc]*tempabs;
              if (sample->std) { /* **CK added for forward mc_std */
		abs3d_var [kc][ic][jc] = result->absorption3D2->tot[kc][ic][jc]*tempabs*result->absorption3D2->tot[kc][ic][jc]*tempabs;
              }
	    }
          }
        }
      }
    }
  return 0;
}

/***********************************************************************************/
/* Function: summarize_result_boxairmass                                  @62_30i@ */
/* Description: Outputs box airmass factors                                        */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int summarize_result_boxairmass (sample_struct* sample, 
                                 atmosphere_struct* atmos,
                                 result_struct* result,
                                 char* basename,
				 int islower,
                                 int isupper,
                                 int jslower,
                                 int jsupper)
{
  int is=0,js=0,ip=0,kc=0,id=0;
  FILE *airmassfile=NULL;
  char airmassfilename    [FILENAME_MAX] ="";  
  strcpy(airmassfilename, basename);
  strcat(airmassfilename, ".amf");  
    
  if ((airmassfile = fopen (airmassfilename, "w")) == NULL) {
    perror(NULL);
    fprintf (stderr, "errno = %d\n", errno);
    fprintf (stderr, "Error opening %s for writing\n", airmassfilename);
    return -1; 
  }
  if (sample->backward) {
    for (is=islower; is<=isupper; is++) {
      for (js=jslower; js<=jsupper; js++) {
        for (ip=0; ip<sample->nstokes; ip++){
          for (kc=0; kc<atmos->Nz; kc++){
            result->pathlength_per_layer_tot[kc]/=
              (result->back[is][js][ip]*(atmos->Z[kc+1]-atmos->Z[kc]));
          }
        }
      }
    }
  }
  else {
    for (is=0; is<sample->Nx; is++) {
      for (js=0; js<sample->Ny; js++) {
	for (id=0; id<sample->Nd; id++) {
	  for (ip=0; ip<sample->nstokes; ip++){
	    for (kc=0; kc<atmos->Nz; kc++){
	      result->pathlength_per_layer_tot[kc]/=
		(result->alt[0]->radesc[id][is][js][ip]*(atmos->Z[kc+1]-atmos->Z[kc]));
	    }
          }
        }
      }
    }
  }
  for (kc=0; kc<atmos->Nz; kc++){
    fprintf(airmassfile, "%10.3f %12.6f\n", atmos->Z[kc]/1000., 
            result->pathlength_per_layer_tot[kc]);
  }
  
  fclose(airmassfile); 
  return 0;
}


/***********************************************************************************/
/* Function: reflect                                                      @62_30i@ */
/* Description:                                                                    */
/*  Surface reflection. If albedo->reflectalways is set, each photons is reflected */
/*  isotropically and the weight is reduced according to the BRDF. Otherwise,      */
/*  the photon is randomly reflected of not (works only with Lambertian albedo)    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int reflect (photon_struct *p, albedo_struct *albedo, elevation_struct *elev,
                    atmosphere_struct *atmos, sample_struct *sample, 
                    result_struct *result,
                    float wvnmlo, float wvnmhi, float *refind)
{
  int il=0, id=0, isp=0, status=0, ic=0, dolamb=0;
  double mu1=0, mu2=0, phi1=0, phi2=0, paw=0;

  double norm[3] = {0,0,1};
  double cotheta=0;

  int ia=0, ja=0;
  photon_struct *p_escape=NULL;
  int ie=0, je=0;
  double nr=0.0, ni=0.0, zsurf=0.0;

  /* Rotation matrix for transformation of surface in spherical geometry */
  
#if HAVE_VROOM
  double mu=0.0, mu_ddis=0.0, phi=0.0, mu_horz=0.0, P_norm=0.0, stheta_ddis=0.0, cphi_ddis=0.0;
  int DDISsing=MCDDIS_NONE, FODDISsing=0, out_of_cone=0;
  locest_struct lest;
  scadis_struct scadis;
#endif

  int i=0;
  double pdir_inc[3];
  
#ifdef MUCHOUT
  double mua=0.0;
#endif

#if HAVE_LIDAR
  int it=0;
#endif

#ifdef MUCHOUT
  if (p->muchoutcounter==MUCHOUTPHOTON) {
    if (p->scattercounter==0)
      fprintf(stderr,"counters %d %d %d %d --- reflect; satdir: %e %e %e\n",
	      p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
	      sample->rad[0].dir.dx[0],sample->rad[0].dir.dx[1],sample->rad[0].dir.dx[2]);
    fprintf(stderr,"reflect x %e %e %e i %d %d %d dir %e %e %e spec %e\n",p->x[0],p->x[1],p->x[2],p->ic,p->jc,p->kc,p->dir.dx[0],p->dir.dx[1],p->dir.dx[2],p->x[2]-400.0);
  }
#endif

  if(sample->polarisation){
    /*Copy incoming direction*/
    for (i=0; i<3; i++)
      pdir_inc[i]=p->dir.dx[i];
  }
  
  /* determine coordinates in albedo/BRDF grid */
    switch (albedo->method) {
    case MCALB_LAM2D:
    case MCALB_LAM2D_SPECTRAL:
    case MCALB_RPV2D_SPECTRAL:
    case MCALB_ROSSLI2D:
      if (sample->spherical3D) {
#ifdef HAVE_SPHER
	coord_spherical3D (p,
			   albedo->Nx, albedo->X,
			   albedo->Ny, albedo->Y,
			   0, NULL,
			   0, NULL,
			   0,
			   &ia, &ja, NULL, NULL,
			   0);
#else
	fprintf(stderr,"Error! you are not allowed to use spherical 3D!\n");
	return -1;
#endif
      }
      else
	albedo_coord (p, albedo, &ia, &ja);
    case MCALB_LAM:
    case MCALB_RPV:
    case MCALB_COXANDMUNK: 
    case MCALB_HAPKE:
    case MCALB_ROSSLI:
    case MCALB_TSANG:
      break;
    default:
      fprintf (stderr, "Fatal error, unknown albedo type %d\n", albedo->method);
      return -1;
    }

  /* get surface type and use corresponding table */
  if (albedo->method == MCALB_RPV2D_SPECTRAL || albedo->method == MCALB_LAM2D_SPECTRAL)
    il = albedo->rpv_index[ia][ja];

  if(sample->ncirc) {
    ic = (int) (sqrt ( ( sample->sensorposition_x[0] - p->x[0] )
		       *( sample->sensorposition_x[0] - p->x[0] )
		       + ( sample->sensorposition_x[1] - p->x[1] )
		       *( sample->sensorposition_x[1] - p->x[1] ) )
		/ sample->dcirc);
    if (il>9)
      fprintf(stderr,"Increase number of indices lcirc\n");
    if (ic<sample->ncirc && il!=-1) {
      if (p->ncircos>99) {
	fprintf(stderr,"Increase number of ncircos\n");
	exit(1);
      }
      p->tocirco[p->ncircos]=&(result->circcontr[il][ic]);
      p->ncircos++;
    }
  }

  /* explicit treatment of reflection: photons are */
  /* randomly reflected according to the albedo;   */
  /* if a photon is not reflected, return 0 and    */
  /* kill the photon immediately afterwards.       */

  if (!albedo->reflectalways) { 
    switch (albedo->method) {
    case MCALB_LAM2D:
      if (!reflection(albedo->albedo2D[ia][ja]))
        return MCSTATUS_ABSORB;

      break;
      
    case MCALB_LAM2D_SPECTRAL:
      if (!reflection(albedo->alb_type[(int) albedo->rpv_index[ia][ja]]))    
        return MCSTATUS_ABSORB;

      break;
      
    case MCALB_LAM:
      if (!reflection(albedo->albedo))
        return MCSTATUS_ABSORB;
      break;
        
    case MCALB_RPV:
    case MCALB_RPV2D_SPECTRAL:
    case MCALB_COXANDMUNK: 
    case MCALB_TSANG:
    case MCALB_HAPKE:
    case MCALB_ROSSLI:
    case MCALB_ROSSLI2D:
        
      fprintf (stderr, "Error, BRDF only implemented together with reflectalways!\n");
      return -1;
    break;

    default:
      fprintf (stderr, "Error, albedo->method %d not yet implemented!\n", albedo->method);
      return -1;
    }
  }
  
  /* Normal surface vector */
  calc_normal_vector ( p, elev, sample, atmos,
		       1, 0, /* in any case, use elev if present, never use sensordir */
		       norm );

  /* initial photon direction (only used for implicit albedo method) */
  v_mult_mu (norm, p->dir.dx, &mu1);
  mu1 = -mu1;
  phi1 = calc_phi_horz(p->dir.dx, norm);

  /* relative azimuth between incoming photon and wind */
  paw = phi1 - albedo->uphi;
  /* if solar wind, set paw to zero, is old version */
  if (albedo->solar_wind)
    paw = 0.0;
  
  if (sample->escape && sample->DoLE &&
      p->scattercounter+1 >= sample->minscatters &&
      p->escapescattercounter <= p->scattercounter+1) {

  
    if (!p->doling) {
      if (il!=0) { /* double reflection from ocean impossible! */
	/* copied from escape stuff */
	id = sample->Nd;
      
	v_mult_mu (norm, sample->rad[id].dir.dx, &cotheta);

#ifdef MUCHOUT
	if(p->muchoutcounter==MUCHOUTPHOTON) {
	  fprintf(stderr,"doling part one (reflect): cotheta %e \n",cotheta);
	  fprintf(stderr,"                  norm %e %e %e\n",norm[0],norm[1],norm[2]);
	  fprintf(stderr,"                  raddir %e %e %e\n",sample->rad[id].dir.dx[0],sample->rad[id].dir.dx[1],sample->rad[id].dir.dx[2]);
	}
#endif

	if (cotheta > 1e-5) { /* quick fix, too small values lead to problems! */

	  if (sample->refraction){
	    fprintf(stderr,"Error! Refraction and reflect do not work together!\n");
	    return -1;
	  }

	  p_escape = create_escape_photon ( p, atmos, sample, id );
	  /* fprintf(stderr,"direction dole 1: %e %e %e\n",p_escape->dir.dx[0],p_escape->dir.dx[1],p_escape->dir.dx[2]); */
	  /* very dirty trick: */
	  p_escape->wtree = p->wtree;

	  phi2 = calc_phi_horz(sample->rad[id].dir.dx, norm);

	  /* no polarization*/
	  if(!sample->polarisation ){
	    p_escape->pdir[0] = 4.0 * cotheta *
	      reflection_probability_tot (albedo, p_escape,
					  cotheta, mu1, phi2 - phi1,
					  ia, ja, il, sample->DoLE,
					  wvnmlo, wvnmhi, paw,
					  sample->spherical3D,
					  sample->spectral_is, 
					  atmos->nlambda_abs,
					  &status );
	    if (status)
	      return fct_err_out (status, "reflection_probability_tot", ERROR_POSITION);
	    
	  }
	  else{
	    /* Calculate reflection matrix */
	    p_escape->pdir[0] = 4.0 * cotheta *
	      reflection_polarized (albedo, p_escape, 
				    sample->nstokes, sample->backward,
				    p->dir.dx, sample->rad[id].dir.dx, norm,
				    mu1, phi1*PI/180,
				    cotheta,
				    phi2*PI/180,
				    sample->phi_target,
				    ia, ja, il,
				    wvnmlo, wvnmhi, 1,
				    sample->spherical3D,
				    sample->spectral_is, atmos->nlambda_abs, atmos->lambda,
				    &status);
	    
	    if (status!=0)
	      return err_out ("Error, reflection_polarized() returned status %d\n", status);
	  }

	  /* jacobian for surface does not make sense */
	  if (sample->abs_jacobian || sample->LLE_jacobian)
	    for (isp=0; isp<atmos->n_caoth; isp++)
	      p_escape->lest.pdir_sct [isp][0] = 0.0; 
        

	  /* do not calculate escape probability */
	  p_escape->doling=0; /* skip count escape */
	  status = escape_probability (p_escape, atmos,
				       sample, result,
				       elev, 0, il, -1,
				       1, refind);
      
	  if (status!=0)
	    return err_out ("Error, escape_probability() returned status %d\n", status);

	  /* Call reflect for Fresnel */
	  /* from photon_journey */
	  switch(p_escape->photon_status) {
	  case MCSTATUS_SURFACE:
	    elev_coord (p_escape, elev, &ie, &je);
	    zsurf = elevation (elev->surf[ie][je],
			       p_escape->x[0] - (double) ie * elev->delX, 
			       p_escape->x[1] - (double) je * elev->delY);
	    if (fabs(p_escape->x[2] - zsurf) > MC_EPSILON * elev->surfmax) {
	      fprintf (stderr, "FATAL error:\n");
	      fprintf (stderr, "cross_surface() found that photon hit the surface at\n");
	      fprintf (stderr, "x = %g, y = %g, z = %.10f\n", p->x[0], p->x[1], p->x[2]);
	      fprintf (stderr, "while the surface elevation at this xy-location\n");
	      fprintf (stderr, "is actually %.10f.\n", zsurf);
	      return -1;
	    }
	    p_escape->x[2] = zsurf;
	  case MCSTATUS_BOUNDARY_LOWER:
	    status = reflect (p_escape, albedo, elev, atmos, sample, result,
			      wvnmlo, wvnmhi, refind);              
	    break;
	  }

	  /* free memory */
	  destroy_photon (p_escape, atmos->n_caoth);
#ifdef MUCHOUT
	  if(p->muchoutcounter==MUCHOUTPHOTON)
	    fprintf(stderr,"end doling part one (reflect)\n");
#endif
	}
      }
    }
    else { /* second part of dole */
#ifdef MUCHOUT
      if(p->muchoutcounter==MUCHOUTPHOTON)
	fprintf(stderr,"doling?\n");
#endif
      if (il==0) { /* only apply if ocean=CaM */
#ifdef MUCHOUT
	if(p->muchoutcounter==MUCHOUTPHOTON)
	  fprintf(stderr,"doling part two\n");
#endif

	/* copied for escape stuff */
	id = 0;
      
	v_mult_mu (norm, sample->rad[id].dir.dx, &cotheta);

	if (cotheta > 0) {

	  if (sample->refraction){
	    fprintf(stderr,"Error! Refraction and reflect do not work together!\n");
	    return -1;
	  }

	  p_escape = create_escape_photon ( p, atmos, sample, id );
	  /* very dirty trick: */
	  p_escape->wtree = p->wtree;

	  /* change photon direction */
	  cp_direction (&(p_escape->dir), &(sample->rad[0].dir));
	  /* fprintf(stderr,"direction dole 2: %e %e %e\n",p->dir.dx[0],p->dir.dx[1],p->dir.dx[2]); */

	  phi2 = calc_phi_horz(sample->rad[id].dir.dx, norm);

	  /* Multiply locest with Fresnel */
	  index_water ( 0.5*1e4*(1.0/wvnmlo + 1.0/wvnmhi),
			albedo->xsal,
			&nr, &ni );
	  p_escape->pdir[0] *= Fresnel(nr,ni,cotheta);
	  /* fprintf(stderr,"fresnel %e %e \n",cotheta,Fresnel(nr,ni,cotheta)); */

	  /* calculate escape probability */
	  p_escape->doling=0; /* now do count escape */
	  if (p_escape->kc<0)
	    p_escape->kc=0;
	  status = escape_probability (p_escape, atmos,
				       sample, result,
				       elev, id, il, -1,
				       1, refind);

	  if (status!=0)
	    return err_out ("Error, escape_probability() returned status %d\n", status);

	  /* free memory */
	  destroy_photon (p_escape, atmos->n_caoth);
	}
      }

#ifdef MUCHOUT
      if(p->muchoutcounter==MUCHOUTPHOTON)
	fprintf(stderr,"end doling part two\n");
#endif
      return 0;
    }
  }

  /* calculate probabilities for escape radiance */
  if (sample->escape && !(sample->DoLE && il==0) &&
      p->scattercounter >= sample->minscatters &&
      p->escapescattercounter <= p->scattercounter) {

    /* calculate probabilities for the radiance angles;    */
    /* for the case of a Lambertian surface this is very   */
    /* simple because all directions are equally likely;   */
    /* in case of a bi-directional reflection function,    */
    /* the thus introduced angular dependence needs to be  */
    /* taken into account (as it is done in scattering())  */

    for (id=0; id<sample->Nd; id++) {
      
      v_mult_mu (norm, sample->rad[id].dir.dx, &cotheta);

      if (cotheta > 0) {

	if (sample->refraction){
	  fprintf(stderr,"Error! Refraction and reflect do not work together!\n");
	  return -1;
	}

	p_escape = create_escape_photon ( p, atmos, sample, id );

	phi2 = calc_phi_horz(sample->rad[id].dir.dx, norm);

	/* no polarization*/
	if(!sample->polarisation ){
	  p_escape->pdir[id] = 4.0 * cotheta *
	    reflection_probability_tot (albedo, p_escape,  
					cotheta, mu1, phi2 - phi1,
					ia, ja, il, sample->DoLE,
					wvnmlo, wvnmhi, paw,
					sample->spherical3D, 
					sample->spectral_is, atmos->nlambda_abs, 
					&status );
	  if (status)
	    return fct_err_out (status, "reflection_probability_tot", ERROR_POSITION);
	}
	else{
	  /* Calculate reflection matrix */
	  p_escape->pdir[id] = 4.0 * cotheta *
	    reflection_polarized (albedo, p_escape, 
				  sample->nstokes, sample->backward,
				  p->dir.dx, sample->rad[id].dir.dx, norm,
				  mu1, phi1*PI/180,
				  cotheta,
				  phi2*PI/180,
				  sample->phi_target,
				  ia, ja, il,
				  wvnmlo, wvnmhi, 1,
				  sample->spherical3D, 
				  sample->spectral_is, atmos->nlambda_abs, atmos->lambda,
				  &status);
	  
	  if (status!=0)
	    return err_out ("Error, reflection_polarized() returned status %d\n", status);
	}

	/* jacobian for surface does not make sense */
	if (sample->abs_jacobian || sample->LLE_jacobian)
	  for (isp=0; isp<atmos->n_caoth; isp++)
	    p_escape->lest.pdir_sct [isp][0] = 0.0; 
        
	/* calculate escape probability */
	status = escape_probability (p_escape, atmos,
				     sample, result,
				     elev, id, il, -1,
				     1, refind);
      
	if (status!=0)
	  return err_out ("Error, escape_probability() returned status %d\n", status);

	/* free memory */
	destroy_photon (p_escape, atmos->n_caoth);
      }
    }
  }

#if HAVE_LIDAR
  /* Single Local Estimator */
  if ( sample->LidarLocEst && !sample->DoLE &&
       p->scattercounter >= sample->minscatters &&
       p->escapescattercounter <= p->scattercounter) {
    
    status = lidar_detector_kernel ( sample->lidar[sample->ili], p->x, &(p->lest) );
    if (status!=0)
      return err_out ("Error %d returned by lidar_detector_kernel()\n", status);
    
    /* stop immediately if photon passed last range bin */
    if (p->pathlength > p->maxpathlength - p->lest.dist) {
      p->photon_status=MCSTATUS_PURGE;
      return MCSTATUS_PURGE;
    }
    
    /* test whether local estimate of this scatter can contribute to result */
    /* find range bin */
    if (sample->LLE_Nt > 0) {
      it = locate (sample->lidar_t, sample->LLE_Nt, p->pathlength + p->lest.dist);
      if (it >= sample->LLE_Nt || it < 0) {
	it = -1;
	p->lest.cosalpha = -1.0;
      }
    }
    else
      it = -1;

    /* only if event within opening angle of lidar */
    if ( sample->lidar[sample->ili].cosalpha[0] <= p->lest.cosalpha ) {
      /* calculate angle between initial photon direction and lidar direction */

      v_mult_mu (norm, p->lest.dir.dx, &cotheta);

      if (cotheta > 0) {
	p_escape = create_escape_photon ( p, atmos, sample, 0 );

	if (sample->refraction){
	  fprintf(stderr,"Error! Refraction and reflect do not work together!\n");
	  return -1;
	}

	phi2 = calc_phi_horz(p_escape->lest.dir.dx, norm);

        /* probability that the photon scatters into the direction of the lidar detector */
        if (!sample->polarisation){
          p_escape->lest.pdir = 4.0 * cotheta *
            reflection_probability_tot (albedo, p, 
                                        cotheta, mu1, phi2 - phi1,
                                        ia, ja, il, sample->DoLE,
                                        wvnmlo, wvnmhi, paw,
					sample->spherical3D, 
					sample->spectral_is, atmos->nlambda_abs,
					&status );
          
	  if (status)
	    return fct_err_out (status, "reflection_probability_tot", ERROR_POSITION);
        }
        else{
          /* Calculate reflection matrix for water surface (Mishchenko code) */
          p_escape->lest.pdir = 4.0 * cotheta *
            reflection_polarized (albedo, p, 
				  sample->nstokes, sample->backward,
				  p->dir.dx, p_escape->lest.dir.dx, norm,
				  mu1, phi1*PI/180,
				  cotheta,
				  phi2*PI/180,
				  sample->phi_target,
				  ia, ja, il,
				  wvnmlo, wvnmhi, 1,
				  sample->spherical3D, 
				  sample->spectral_is, atmos->nlambda_abs, atmos->lambda,
				  &status);
	  
          if (status!=0)
            return err_out ("Error, reflection_polarized() returned status %d\n", status);
        }


	/* assuming that reflection is isoenergetic */
	p_escape->lest.pdir_iso = p_escape->lest.pdir;

	/* jacobian for surface does not make sense */
	if (sample->abs_jacobian || sample->LLE_jacobian)
	  for (isp=0; isp<atmos->n_caoth; isp++)
	    p_escape->lest.pdir_sct [isp][0] = 0.0; 

        /* propagate virtual photon to the lidar detector and add contribution to local estimate */
        status = escape_probability (p_escape, atmos,
                                     sample, result,
                                     elev, 0, il, it,
                                     0, refind);
        if (status!=0)
          return err_out ("Error, escape_probability() returned status %d\n", status);

	/* free memory */
	destroy_photon (p_escape, atmos->n_caoth);
      }
    }
  }
#endif

  /********************************************/
  /* 3. return in case CP is no longer needed */
  /********************************************/

  if ( p->isclone && p->clonescattercounter+2 == sample->ntupelLE ) { /* n-tupel local estimate */
    p->photon_status = MCSTATUS_PURGE;
    return MCSTATUS_PURGE;
  }


#if HAVE_VROOM

  /*******************/
  /* 4. DDIS stuff 1 */
  /*******************/

  /* NOTE: Lidar and DDIS and reflection has not been joined completely */
  /* what is not yet working: UDA */
  /*                          jacobians */
  /* what is working but will evtl. create spikes: */
  /*            DDIS when Detector slightly above horizon */
  /* Latter could be improved by iterative stuff if necessary (RPB,15.1.37) */

  if ( sample->escape_eps_ddis_upf || sample->LLE_D_DIS ) {
    status = mc_vroom_prep_DDIS (sample, p, atmos, &DDISsing, &FODDISsing, &out_of_cone, &lest, &scadis);
    if (status!=0)
      return err_out ("Error %d returned by mc_vroom_prep_DDIS ()\n", status);

    /* UDA not implemented for reflection */
    if (DDISsing == MCDDIS_UDA)
      DDISsing = MCDDIS_UPF;
  }

#ifdef MUCHOUT
  if (p->muchoutcounter==MUCHOUTPHOTON) {
    fprintf(stderr,"counters %d %d %d %d --- reflect DDIS: muhorz %e mumax %e epsfac %e DDISING %d mumin %e Fmin %e \n",
	    p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
	    mu_horz, scadis.mu_max, scadis.epsfac, DDISsing, scadis.mu_min, scadis.F_min[0]);
  }
#endif

#endif


  /******************************/
  /* 7. calculate new direction */
  /******************************/

#if HAVE_VROOM
  if (il==0 && sample->DoLE) { /* specular reflection */
    p->dir.dx[2]=-p->dir.dx[2];
    hitflag (&(p->dir));
    v_mult_mu ( p->dir.dx, norm, &mu_horz );
    if (mu_horz<1e-5) { /* ocean is not horizontal, kill photon */
      p->photon_status = MCSTATUS_PURGE;
      return MCSTATUS_PURGE;
    }
    /* Multiply photon weight with Fresnel */
    index_water ( 0.5*1e4*(1.0/wvnmlo + 1.0/wvnmhi),
		  albedo->xsal,
		  &nr, &ni );
    p->weight *= Fresnel(nr,ni,p->dir.dx[2]);
#ifdef MUCHOUT
    if(p->muchoutcounter==MUCHOUTPHOTON)
      fprintf(stderr,"specular reflect %e %e %e %e %e\n",p->dir.dx[0],p->dir.dx[1],p->dir.dx[2],Fresnel(nr,ni,p->dir.dx[2]),mu_horz);
#endif
    return MCSTATUS_INCREMENT;
  }
  else {
    /* Lambertian reflection for all except CaM */
    dolamb=1;

    switch (albedo->method) {
    case MCALB_LAM:
    case MCALB_LAM2D:
    case MCALB_LAM2D_SPECTRAL:
    case MCALB_RPV:    /* in this case, il=0, and rpv is saved in il=0! */
    case MCALB_HAPKE:
    case MCALB_ROSSLI:
      break;
    case MCALB_TSANG:
    case MCALB_COXANDMUNK:
      dolamb=0;
      break;
    case MCALB_RPV2D_SPECTRAL:  /* in this case, Cox&Munk is in il=0 */
      if (il==0) {
	if (sample->DoLE)
	  return  err_out ("Error, oceabrdfc 1 should not be called when DoLE is %d\n", sample->DoLE);
	dolamb=0;
      }
      break;
    case MCALB_ROSSLI2D:
      if (albedo->rossli2D[ia][ja].isCaM) {
	if (sample->DoLE)
	  return err_out ("Error, oceabrdfc 2 should not be called when DoLE is %d\n", sample->DoLE);
	dolamb=0;
      }
      break;
    default:
      fprintf (stderr, "Error, albedo->method %d not yet implemented!\n", albedo->method);
      return -1;
    }

    random_reflection_special (sample, atmos, elev, p,
			       sample->phase_max, sample->n_phase_max,
			       &mu, &phi, norm,
			       &scadis, lest, dolamb, DDISsing, FODDISsing, out_of_cone); 
  }

#else
  /* photon is always reflected Lambertian, that is, */
  /* with constant intensity in all directions       */
  random_Lambertian_normal (&(p->dir), norm);
#endif

  /* zero weight -> don't need to trace photon any further */
  if (p->weight==0.) {
    p->photon_status=MCSTATUS_PURGE;
    return MCSTATUS_PURGE;
  }

#ifdef MUCHOUT
  if (p->muchoutcounter==MUCHOUTPHOTON) {
    v_mult_mu ( p->dir.dx, sample->rad[0].dir.dx, &mua );
  
    fprintf(stderr,"counters %d %d %d %d --- intdir: %e %e %e mu %e box %d %d %d\n",
	    p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
	    p->dir.dx[0], p->dir.dx[1], p->dir.dx[2], mua, p->ic, p->jc, p->kc);
  }
#endif

#if HAVE_VROOM
  /*******************/
  /* 8. DDIS stuff 2 */
  /*******************/

  if ( p->RIS_mode != MCRIS_MODE_NORMAL ||   /* RIS             */
       sample->LLE_D_DIS ||                  /* DDIS for LIDAR  */
       sample->escape_eps_ddis_upf ||        /* DDIS for escape */
       sample->LLE_channels ||               /* Raman/HSRL      */
       sample->LLE_jacobian ||               /* jacobian LIDAR  */
       atmos->ris_factor != 1. ) {          /* is RIS-factor used? */

    v_mult_mu ( norm, p->dir.dx, &mu_horz );

    if (sample->LLE_D_DIS || (sample->escape_eps_ddis_upf != 0.0) ) {
      status = mc_vroom_set_mus_and_phis (sample, p, DDISsing,
					  lest, scadis,
					  &mu, &mu_ddis, &stheta_ddis, phi, &cphi_ddis);
      if (status!=0)
	return err_out ("Error %d returned by mc_vroom_set_mus_and_phis ()\n", status);
    }

    /* we do not do jacobians here! this is still missing */

    /* For calculating the weight concerning DDIS, we here assume */
    /* that the reflection was done lambertian. Further below the */
    /* weight is corrected with respect to other albedo types,    */
    /* which is done with and without DDIS                        */
    if (dolamb)
      P_norm = get_Lambertian_phase(mu_horz);
    else
      /* factor 2.0 because of only upper half of sphere */
      P_norm = 2.0;

    status = mc_vroom_DDIS_weight_and_prep_stuff (sample, atmos, mu_ddis, stheta_ddis, cphi_ddis,
						  DDISsing, out_of_cone,
						  P_norm, P_norm, P_norm, lest, scadis, p);
    if (status!=0)
      return err_out ("Error %d returned by mc_vroom_DDIS_weight_and_prep_stuff ()\n", status);

  }

#endif
  
  /* store new photon direction (only for implicit albedo method) */
  /* implicit treatment of reflection: photons are */
  /* always reflected Lambertian; afterwards, the  */
  /* photon weight is reduced according to the     */
  /* albedo/BRDF                                   */
  if (albedo->reflectalways) {
    v_mult_mu ( norm, p->dir.dx, &mu2 );
    phi2 = calc_phi_horz(p->dir.dx, norm);
  
    /* calculate photon weight according to albedo/BRDF (only for implicit albedo method) */
    if (!sample->polarisation){
      p->weight *= reflection_probability_tot (albedo, p,
                                               mu2, mu1, phi2 - phi1,
                                               ia, ja, il, sample->DoLE,
                                               wvnmlo, wvnmhi, paw,
					       sample->spherical3D, 
					       sample->spectral_is, atmos->nlambda_abs,
					       &status );
      
      if (status)
	return fct_err_out (status, "reflection_probability_tot", ERROR_POSITION);
    }
    else{
      /* Calculate reflection matrix */
      p->weight *= reflection_polarized(albedo, p,
					sample->nstokes, sample->backward,
					pdir_inc, p->dir.dx, norm,
					mu1, phi1*PI/180,
                                        mu2, phi2*PI/180,
					sample->phi_target,
                                        ia, ja, il,
                                        wvnmlo, wvnmhi, 0, 
					sample->spherical3D, 
					sample->spectral_is, atmos->nlambda_abs,
					atmos->lambda,
                                        &status);
      if (status!=0)
        return err_out ("Error, reflection_polarized() returned status %d\n", status);
    }

    /* in case of isotropic distribution of reflected photons, need factor */
    if (!dolamb) 
      /* factor 2.0 because of only upper half of sphere */
      p->weight *= get_Lambertian_phase(mu2)/2.0;

    /* zero weight -> don't need to trace photon any further */
    if (p->weight==0) {
      p->photon_status=MCSTATUS_PURGE;
      return MCSTATUS_PURGE;
    }
  }
  
  /* after surface reflection, photon is no longer direct radiation */
  p->direct=0;
 
#ifdef MUCHOUT
  if (p->muchoutcounter==MUCHOUTPHOTON) {
    v_mult_mu ( p->dir.dx, sample->rad[0].dir.dx, &mua );
  
    fprintf(stderr,"counters %d %d %d %d --- newdir: %e %e %e mu %e box %d %d %d weight %e\n",
	    p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
	    p->dir.dx[0], p->dir.dx[1], p->dir.dx[2], mua, p->ic, p->jc, p->kc,p->weight);
  }
#endif

  return MCSTATUS_INCREMENT;
}


/***********************************************************************************/
/* Function: calc_phi_horz                                                @62_30i@ */
/* Description:                                                                    */
/*  Calculate phirad from dir                                                      */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

double calc_phi_horz (double *dx, double *n)
{
  double x=0.0, y=0.0, phi=0.0;

  int i=0, horizontal=0;

  double u[3], v[3];
  double u2inv=0;

  if (n==NULL)
    horizontal=1;
  else if (n[2]==1.0)
    horizontal=1;

  if (horizontal) {
    /* special case, horizontal surface */
    x=dx[0];
    y=dx[1];
  }
  else {
    /********************************/
    /* u,v are vectors of length 1; */
    /* u,v,dx form a left-handed    */
    /* orthogonal system            */
    /*                              */
    /* Part is from new_direction   */
    /********************************/

    u[2] = + sqrt( n[0] * n[0] + n[1] * n[1] );

    u2inv = 1./u[2];

    v[0] = - n[1] * u2inv;
    v[1] = + n[0] * u2inv;
    v[2] =   0.;

    u[0] = - n[2] * v[1];
    u[1] = + n[2] * v[0];

    for (i=0; i<3; i++)
      x += dx[i] * v[i];
    for (i=0; i<3; i++)
      y += dx[i] * u[i];
  }

  /* radiance direction */
  if (x==0 && y==0)  /* special case: vertical incidence */
    phi = 0.0;
  else
    if (y==0) { /* special case, 1/dx[1] is infinity */
      if (x<0)
	phi=-90.;
      else
	phi=90.;
    }
    else
      phi = atand(x/y);

  if (y<0)
    phi += 180.;
  else
    if (x<0)
      phi += 360.;

  return phi;
}


/***********************************************************************************/
/* Function: reflection_probability_tot                                   @62_30i@ */
/* Description:                                                                    */
/*  Calculate the total probability for reflection angle mu = cos(theta).          */
/*  CE: included spectal albedo weight for ALIS                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static double reflection_probability_tot (albedo_struct *albedo, photon_struct *p,
                                          double murad, double mu1, double deltaphirad, 
                                          int ia, int ja, int il, int DoLE,
                                          float wvnmlo, float wvnmhi, double paw,
					  int spherical3D, 
					  int spectral_is, int nlambda,
					  int *status)
{

  *status=0;
  int iv=0;
  
  if (murad <= 0.0)
    return 0.0;

  /* calculate probabilities for the radiance angles;    */
  /* for the case of a Lambertian surface this is very   */
  /* simple because all directions are equally likely;   */
  /* in case of a bi-directional reflection function,    */
  /* the thus introduced angular dependence needs to be  */
  /* taken into account (as it is done in scattering())  */

  if (!(albedo->reflectalways))
    return 1.0;

  switch (albedo->method) {
  case MCALB_LAM:
    if(spectral_is && albedo->spectral_albedo!=NULL){
      if ( albedo->albedo != 0.0){ 
	for (iv=0; iv<nlambda; iv++){
	  p->q_albedo_spectral[iv] *=  albedo->spectral_albedo[iv] /
	    albedo->albedo ;
	}
      }
    }
    return albedo->albedo;
      
  case MCALB_LAM2D:
    return albedo->albedo2D[ia][ja];

  case MCALB_LAM2D_SPECTRAL:
    if(spectral_is){
      if(albedo->alb_type[(int) albedo->rpv_index[ia][ja]] != 0.0){
	for (iv=0; iv<nlambda; iv++){
	  if (spherical3D) {
#ifdef HAVE_SPHER
	    coord_spherical3D (p,
			       albedo->Nx, albedo->X,
			       albedo->Ny, albedo->Y,
			       0, NULL,
			       0, NULL,
			       0,
			       &ia, &ja, NULL, NULL,
			       0);
#else
	    fprintf(stderr,"Error! you are not allowed to use spherical 3D!\n");
	    *status=-1;
	    return NOT_A_NUMBER;
#endif
	  }
	  else
	    albedo_coord (p, albedo, &ia, &ja);	
	 
	  p->q_albedo_spectral[iv] *=  
	    albedo->spectral_alb_type[(int) albedo->rpv_index[ia][ja]][iv]/
	    albedo->alb_type[(int) albedo->rpv_index[ia][ja]] ;
	  }
      }	
    }
    return albedo->alb_type[(int) albedo->rpv_index[ia][ja]];
    
  case MCALB_RPV:    /* in this case, il=0, and rpv is saved in il=0! */
  case MCALB_RPV2D_SPECTRAL:  /* in this case, Cox&Munk is in il=0 */
    if (il==0 && albedo->method == MCALB_RPV2D_SPECTRAL) {
      /* attention: mu1, murad, deltaphirad are specified differently in oceabrdfc, see also src_f/BDREF.f */
      if (DoLE) {
	*status=-1;
	return  err_out ("Error, oceabrdfc 1 should not be called when DoLE is %d\n", DoLE);
      }

      return oceabrdfc (wvnmlo, wvnmhi, murad, mu1, deltaphirad/180.0*PI, 
			albedo->u10, albedo->pcl, albedo->xsal, paw);
    }
    else 
      return rpv_brdf (albedo->rpv[il].rho0, albedo->rpv[il].k, albedo->rpv[il].theta,
		       albedo->rpv[il].scale, albedo->rpv[il].sigma,
		       albedo->rpv[il].t1, albedo->rpv[il].t2,
		       mu1, murad, fabs(deltaphirad));

  case MCALB_HAPKE:
    return c_bidir_reflectivity_hapke ( wvnmlo, wvnmhi, mu1, murad, fabs(deltaphirad),
					albedo->hapke[0].b0, albedo->hapke[0].h, albedo->hapke[0].w );
  case MCALB_ROSSLI:
#if 1
    return c_bidir_reflectivity_rossli ( &(albedo->rossli[il]), mu1, murad, fabs(deltaphirad/180.0*PI) );
#else
    // former version
    return ambrals_brdf (albedo->rossli[il].iso, 
			 albedo->rossli[il].vol,
			 albedo->rossli[il].geo,
			 mu1, murad, fabs(deltaphirad));
#endif

  case MCALB_ROSSLI2D:
    if (albedo->rossli2D[ia][ja].isCaM) {
      /* attention: mu1, murad, deltaphirad are specified differently in oceabrdfc, see also src_f/BDREF.f */
      if (DoLE) {
	*status=-1;
	return  err_out ("Error, oceabrdfc 2 should not be called when DoLE is %d\n", DoLE);
      }

      return oceabrdfc (wvnmlo, wvnmhi, murad, mu1, deltaphirad/180.0*PI, 
			albedo->u10, albedo->pcl, albedo->xsal, paw);
    }
    else {
#if 1
      return c_bidir_reflectivity_rossli ( &(albedo->rossli2D[ia][ja]), mu1, murad, fabs(deltaphirad/180.0*PI) );
#else
      // former version
      return ambrals_brdf (albedo->rossli2D[ia][ja].iso,
			   albedo->rossli2D[ia][ja].vol,
			   albedo->rossli2D[ia][ja].geo,
			   mu1, murad, fabs(deltaphirad));
#endif
    }
  case MCALB_COXANDMUNK:
    /* attention: mu1, murad, deltaphirad are specified differently in oceabrdfc, see also src_f/BDREF.f */
    if (DoLE) {
      *status=-1;
      return  err_out ("Error, oceabrdfc 3 should not be called when DoLE is %d\n", DoLE);
    }

    return oceabrdfc (wvnmlo, wvnmhi, murad, mu1, deltaphirad/180.0*PI, 
		      albedo->u10, albedo->pcl, albedo->xsal, paw);

  case MCALB_TSANG:
    fprintf (stderr, "BPDF by Tsang/Mishchenko currently only implemented with polarisation,\n");
    fprintf (stderr, "Please use Cox and Munk BRDF. \n");
    *status = -1;
    return NOT_A_NUMBER;
              
  default:
    fprintf (stderr, "Error, albedo->method %d not yet implemented!\n", albedo->method);
    *status = -1;
    return NOT_A_NUMBER;
  }

  fprintf (stderr, "Error, something wrong in reflection_probability_tot!\n");
  *status = -1;
  return NOT_A_NUMBER;
}


/***********************************************************************************/
/* Function: scattering                                                   @62_30i@ */
/* Description:                                                                    */
/*  Scattering event (Rayleigh, aerosol, water, ice); pdf defined by the           */
/*  properties of the current box or layer.                                        */
/*                                                                                 */
/*  NOTE: Most steps in this routine are only used in case of LIDAR simulations!   */
/*        All other calculations only need the steps 1a, 2, (3)*, 5, 6, 7, 9.      */
/*        Escape can also use DIS: then it further needs steps 4a, 8a, 8c          */
/*        In case some of the tricks used for Lidar are extended to other methods, */
/*        following steps will be needed:                                          */
/*         - Real Importance Sampling : step 8c                                    */
/*         - Virtual Importance Sampling: step 3                                   */
/*                                                                                 */
/*  * Comment: 3 is also needed for normal mode since random_scatter_type can even */
/*             return VIRTUAL in normal mode! This is due to the fact than         */
/*             1D profiles are double while 3D profiles are float.                 */
/*             Small discrepancies in the sum of the scatter coefficient can thus  */
/*             lead to an undefined scatter type in atmospheric layers with both   */
/*             1D and 3D media. The probability of such is ca. 0.5e-7.             */
/*             Within floating point precision, it is safest to ignore the         */
/*             scattering and let the photon move on.                              */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: ...., last big modification by Robert Buras                             */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int scattering ( atmosphere_struct *atmos,
			photon_struct     *p,
			sample_struct     *sample,
			result_struct     *result,
			elevation_struct  *elev,
			albedo_struct     *albedo,
			float              wvnmlo,
			float              wvnmhi, 
			float             *refind )
{
  double mu=0.0, mu2=0.0, stheta2=0.0, phi=0.0, cphi2=0.0;

  int i=0, j=0, id=0, ip=0, iv=0;
  int ie=0, je=0;
  double zsurf=0.0;

  int status=0;

  int isp=0;
  photon_struct *p_escape=NULL;

#if HAVE_VROOM
  double P_norm=0.0, P_spec=0.0, P_isoene=0.0;
  int DDISsing=MCDDIS_NONE, FODDISsing=0, out_of_cone=0;
  locest_struct lest;
  scadis_struct scadis;
#endif
  
#if HAVE_LIDAR
  int it=0;
  int cb_hit=0;
  int within_angle=0;
  double MKn=0.;
  double **Ztemp = NULL;
  double mu_save;
#endif

#ifdef MUCHOUT
  double mua=0.0;
#endif

  double pdir_inc[3]; 
  double phi_inc=0.0; 
  
  /* six independent elements for randomly oriented particles */
  double phase_matrix[6]= { 0, 0, 0, 0, 0, 0 };
  
  /* Set to 1 if method to calculate total phase matrix for whole atmosphere 
     should be used (works forward and backward). 
     Else the Stokes vector is calculated along the photon
     path (works only in forward mode).*/
  int polmat=1;

  double dir_sensor[3];

  if (sample->polarisation){
    /*Copy incoming direction*/
    for (i=0; i<3; i++)
      pdir_inc[i]=p->dir.dx[i];
  }

  /*************************************************************************/
  /* 1. Estimators: Estimate probability that photon contributes to signal */
  /*************************************************************************/

  /* weight for spectral importance sampling */
  /* BCA, ksca_tot in denominator should probably be get_kscaIS!!! XXX */
  /* CE: this is not correct when 3D clouds are included, quick fix with passback3D */
  if (sample->spectral_is){
    
      for (iv=0; iv<atmos->nlambda_abs; iv++){
	
	if (atmos->threed[MCCAOTH_TOT][p->kc]>=1 )
	  p->q_spectral[iv] *= (get_ksca(atmos,p,MCCAOTH_TOT)-(atmos->ksca[MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->prof [MCCAOTH_TOT])[p->kc]+atmos->ksca_spectral[MCCAOTH_TOT][iv][p->kc])/get_ksca(atmos,p,MCCAOTH_TOT);
	
	else
	  p->q_spectral[iv]*= atmos->ksca_spectral[MCCAOTH_TOT][iv][p->kc] 
	    /(atmos->ksca[MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->prof [MCCAOTH_TOT])[p->kc]; 
	
	
	/* debugging output */
      /* 	if (1==1) //iv==0) */
      /* 	  fprintf(stderr, "kc %d p->q_spectral[iv] %g pathlengths %g %g \n q %   g calc: aer %g mol %g,\n  ivs tot %g mol %g aer %g \n", p->kc, p->q_spectral[iv], */
      /* 		  p->pathlength, p->pathlength_per_layer[p->kc], */
      /* 		  atmos->ksca_spectral[MCCAOTH_TOT][iv][p->kc]/ */
      /* 		  (atmos->ksca[MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->prof [MCCAOTH_TOT])[p->kc], */
      /* 		  (atmos->ksca[MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->prof [MCCAOTH_AER])[p->kc], */
      /* 		  (atmos->ksca[MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->prof [MCCAOTH_MOL])[p->kc], */
      /* 		  atmos->ksca_spectral[MCCAOTH_TOT][iv][p->kc], atmos->ksca_spectral[MCCAOTH_MOL][iv][p->kc], */
      /* 		  atmos->ksca_spectral[MCCAOTH_AER][iv][p->kc]);  */
	
      } 
  }
    
  
  if (sample->concentration_is){
    for (i=0; i<atmos->Nc; i++){
      p->q_concentration[i]*=
        ((atmos->ksca[MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->prof [MCCAOTH_TOT])[p->kc]
         - (atmos->ksca[MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->prof[MCCAOTH_AER])[p->kc]
         + atmos->ksca_scaled[i][p->kc])
        / (atmos->ksca[MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->prof [MCCAOTH_TOT])[p->kc];
      /* debugging output */
      /* if (i==0) */
/*         fprintf(stderr, " p->q_concentration[i] %g kscatot %g  ksca_aer %g ksca_scaled %g \n",  p->q_concentration[i], */
/*                 (atmos->ksca[MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->prof [MCCAOTH_TOT])[p->kc], */
/*                 (atmos->ksca[MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->prof[MCCAOTH_AER])[p->kc], */
/*                 atmos->ksca_scaled[i][p->kc]); */
    }
  }

  /***************************************************/
  /* 1a. Escape Radiances, aka Directional Estimator */
  /***************************************************/

#ifdef MUCHOUT
  if (p->muchoutcounter==MUCHOUTPHOTON) {
    if (p->scattercounter==0 && sample->escape)
      fprintf(stderr,"counters %d %d %d %d --- scatter; satdir: %e %e %e\n",
	      p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
	      sample->rad[0].dir.dx[0],sample->rad[0].dir.dx[1],sample->rad[0].dir.dx[2]);
    if (p->scattercounter==0 && sample->LidarLocEst)
      fprintf(stderr,"counters %d %d %d %d --- scatter; satdir: %e %e %e\n",
	      p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
	      sample->lidar[sample->ili].dir.dx[0],
	      sample->lidar[sample->ili].dir.dx[1],
	      sample->lidar[sample->ili].dir.dx[2]);
  }
#endif

  if (sample->escape && sample->DoLE &&
      p->scattercounter+1 >= sample->minscatters &&
      p->escapescattercounter <= p->scattercounter+1) {

    if (!p->doling) {
#ifdef MUCHOUT
      if(p->muchoutcounter==MUCHOUTPHOTON)
	fprintf(stderr,"doling part one (scatter)\n");
#endif

      /* copied from escape stuff */
      id = sample->Nd;

      p_escape = create_escape_photon ( p, atmos, sample, id );
      /* very dirty trick: */
      p_escape->wtree = p->wtree;
      /* fprintf(stderr,"direction dols 1: %e %e %e\n",p_escape->dir.dx[0],p_escape->dir.dx[1],p_escape->dir.dx[2]); */
      if (sample->refraction){
	/* RPB change! I replaced p with p_escape. This should be more correct */
        correct_escape_dir_refract( dir_sensor, p_escape, atmos, sample, refind );
        v_mult_mu ( p->dir.dx, dir_sensor, &mu2 );
      }
      else
        v_mult_mu ( p->dir.dx, sample->rad[id].dir.dx, &mu2 );

      status = get_phase_matrix_total ( atmos, p_escape, mu2, 
					sample->nphamat, 1, sample->spectral_is, sample->concentration_is,
					atmos->ris_factor,
					0,
					phase_matrix,
					p_escape->lest.pdir_sct, NULL, &(p_escape->weight) );
      if (status)
	return fct_err_out (status, "get_phase_matrix_total", ERROR_POSITION);

      id = 0;
      p_escape->pdir[id]=phase_matrix[0];

      /* Calculate weight vector for the sensor direction */
      if(sample->polarisation) {
	/* Compute weights of scattered stokes vector */
	status = calc_stokes_vector_for_escape ( sample,
						 atmos,
						 phase_matrix,
						 mu2,
						 p->dir,
						 sample->rad[1].phi + 180.,
						 sample->rad[1].dir,
						 polmat,
						 p_escape );
	if (status)
	  return fct_err_out (status, "calc_stokes_vector_for_escape", ERROR_POSITION);
      }
 
      /* do not calculate escape probability */
      p_escape->doling=0; /* skip count escape */
      status = escape_probability (p_escape, atmos,
                                   sample, result,
                                   elev, 0, -1, -1,
                                   0, refind);
      if (status)
        return fct_err_out (status, "escape_probability", ERROR_POSITION);

      /* Call reflect for Fresnel */
      /* from photon_journey */
      switch(p_escape->photon_status) {
      case MCSTATUS_SURFACE:
	elev_coord (p_escape, elev, &ie, &je);
	zsurf = elevation (elev->surf[ie][je],
			   p_escape->x[0] - (double) ie * elev->delX, 
			   p_escape->x[1] - (double) je * elev->delY);
	if (fabs(p_escape->x[2] - zsurf) > MC_EPSILON * elev->surfmax) {
	  fprintf (stderr, "FATAL error:\n");
	  fprintf (stderr, "cross_surface() found that photon hit the surface at\n");
	  fprintf (stderr, "x = %g, y = %g, z = %.10f\n", p_escape->x[0], p_escape->x[1], p_escape->x[2]);
	  fprintf (stderr, "while the surface elevation at this xy-location\n");
	  fprintf (stderr, "is actually %.10f.\n", zsurf);
	  return -1;
	}
	p_escape->x[2] = zsurf;
      case MCSTATUS_BOUNDARY_LOWER:
	status = reflect (p_escape, albedo, elev, atmos, sample, result,
			  wvnmlo, wvnmhi, refind);              
	break;
      }
      
      /* free memory */
      destroy_photon (p_escape, atmos->n_caoth);
#ifdef MUCHOUT
      if(p->muchoutcounter==MUCHOUTPHOTON)
	fprintf(stderr,"end doling part one (scatter)\n");
#endif
    }
  }

  if (sample->escape &&
      p->scattercounter >= sample->minscatters &&
      p->escapescattercounter <= p->scattercounter) {

    /* calculate probabilities for all other radiance angles */
    for (id=0; id<sample->Nd; id++) {
      p_escape = create_escape_photon ( p, atmos, sample, id );

      /* Calculated correct direction in case of refraction */
      if (sample->refraction){
        correct_escape_dir_refract( dir_sensor, p, atmos, sample, refind );
        /*fprintf(stderr, "dir_sensor %g %g %g \n", dir_sensor[0], 
                dir_sensor[1], dir_sensor[2]);
                fprintf(stderr, "raddir %g %g %g \n",   sample->rad[id].dir.dx[0],
                sample->rad[id].dir.dx[1],  sample->rad[id].dir.dx[2]);
        */      
        
        /* calculate angle between photon direction before scattering and radiance direction */
        v_mult_mu ( p->dir.dx, dir_sensor, &mu2 );
        
      }
      /* calculate angle between photon direction before scattering and radiance direction */
      else
        v_mult_mu ( p->dir.dx, sample->rad[id].dir.dx, &mu2 );
      
      status = get_phase_matrix_total ( atmos, p_escape, mu2, 
					sample->nphamat, 1, sample->spectral_is, sample->concentration_is, 
                                        atmos->ris_factor,
					0,
					phase_matrix,
					p_escape->lest.pdir_sct, NULL, &(p_escape->weight) );
      if (status)
	return fct_err_out (status, "get_phase_matrix_total", ERROR_POSITION);

      p_escape->pdir[id]=phase_matrix[0];

      /* Calculate weight vector for the sensor direction */
      if (sample->polarisation) {
	/* Compute weights of scattered stokes vector */
	status = calc_stokes_vector_for_escape ( sample,
						 atmos,
						 phase_matrix,
						 mu2,
						 p->dir,
						 sample->rad[id].phi + 180.,
						 sample->rad[id].dir,
						 polmat,
						 p_escape );
	if (status)
	  return fct_err_out (status, "calc_stokes_vector_for_escape", ERROR_POSITION);
      }
      if (sample->coherent_backscatter) 
        save_phase_matrix_elements_temp ( result->mish, p_escape->phamat, p_escape->scattercounter, phase_matrix );
      /* propagate virtual photon in the radiance direction and add contributions to escape */
      status = escape_probability (p_escape, atmos,
                                   sample, result,
                                   elev, id, -1, -1,
                                   0, refind);
      if (status)
        return fct_err_out (status, "escape_probability", ERROR_POSITION);

      /* free memory */
      destroy_photon (p_escape, atmos->n_caoth);
    }
  }

#if HAVE_LIDAR
  /***********************/
  /* 1b. Local Estimator */
  /***********************/

  if ( sample->LidarLocEst &&
       p->scattercounter >= sample->minscatters &&
       p->escapescattercounter <= p->scattercounter) {

    /* EXPERIMENTAL option, leave planar_switch = 0 in cloud3d.c for normal behaviour */
    /* use this to simulate planar waves with opening angle = 0 */
    /* weight is not correct yet !!! */
    /* 0 "spherical waves" (normal case) */ 
    /* 1 transmitter planar waves */ 
    /* 2 detector planar waves */ 
    /* 3 both planar waves */
    /* for CB either use 0 or 3 */
    if (sample->planar_switch > 1) {
      MKn = plane_wave_locest (p, sample);
    }
    /* choose hit point on detector and derive distance weight */
    else {
      status = lidar_detector_kernel ( sample->lidar[sample->ili], p->x, &(p->lest));

      if (status!=0)
  	    return err_out ("Error %d returned by lidar_detector_kernel()\n", status);
    }
    /* stop immediately if photon passed last range bin */
    if (p->pathlength > p->maxpathlength - p->lest.dist) {
      p->photon_status=MCSTATUS_PURGE;
      return MCSTATUS_PURGE;
    }

    /* test whether local estimate of this scatter can contribute to result */
    /* find range bin */
    if (sample->LLE_Nt > 0) {
      it = locate (sample->lidar_t, sample->LLE_Nt, p->pathlength + p->lest.dist);
      if (it >= sample->LLE_Nt || it < 0) {
	it = -1;
	p->lest.cosalpha = -1.0;
      }
    }
    else
      it = -1;

    /* only if event within opening angle of lidar */
    if ( sample->planar_switch < 2 && sample->lidar[sample->ili].cosalpha[0] <= p->lest.cosalpha )
      within_angle = 1;

    else if ( sample->planar_switch > 1 && sample->lidar[sample->ili].radius >= MKn )
      within_angle = 1;

    if (within_angle == 1) {
      p_escape = create_escape_photon ( p, atmos, sample, 0 );

      if (sample->refraction){
	fprintf(stderr,"Error! Refraction and lidar do not work together!\n");
	return -1;
      }

      /* when using CB, first scatter must be in FOV of the detector */
      /* and last scatter must be in FOV of the transmitter! */
      if (sample->coherent_backscatter && p->scattercounter != 0) {
        cb_hit = check_opan_cb(p,p_escape,sample,atmos);
      }

      /* calculate angle between initial photon direction and lidar direction */
      v_mult_mu ( p->dir.dx, p_escape->lest.dir.dx, &mu2 );

      /* CB - save second-to-last p->weight */
      if (sample->coherent_backscatter)
        p_escape->pss.weight_sl = p_escape->weight;

      /* probability that the photon scatters into the direction of the lidar detector */
      status = get_phase_matrix_total ( atmos, p_escape, mu2,
		                        sample->nphamat, 1, sample->spectral_is, sample->concentration_is, 
		                        atmos->ris_factor, 0,
		                        phase_matrix,
		                        p_escape->lest.pdir_sct, NULL, &(p_escape->weight) );

      if (status)
	return fct_err_out (status, "get_phase_matrix_total", ERROR_POSITION);

      p_escape->lest.pdir=phase_matrix[0];

      if (sample->LLE_jacobian)
	/* for HSRL lidar, evtl. move somewhere else */
	p_escape->lest.pdir_iso = p_escape->lest.pdir - p_escape->lest.pdir_sct[MCCAOTH_MOL][0];

      /* Calculate weight vector for the sensor direction */
      if (sample->polarisation) {
	/* Compute weights of scattered stokes vector */
	status = calc_stokes_vector_for_escape ( sample,
						 atmos,
						 phase_matrix,
						 mu2,
						 p->dir,
						 calc_phi_horz(p_escape->lest.dir.dx, NULL),
						 p_escape->lest.dir, /* evtl +180. */
						 polmat,
						 p_escape );
	if (status)
	  return fct_err_out (status, "calc_stokes_vector_for_escape", ERROR_POSITION);
      }

      /* !!! CB Lidar Main stuff happens here !!! */
      if (sample->coherent_backscatter) {
        add_coord_to_cohebasca ( p_escape );
        /* only continue if at least second order of scattering */
        if ( p_escape->scattercounter > 0 && cb_hit == 1 ) {
	    for (i=0; i<3; i++)
	      p_escape->cb.end[i] = p_escape->lest.hitpoint[i];
          status = calc_stokes_alt ( sample, atmos, p_escape, elev, refind);
	  if (status)
	    return fct_err_out (status, "calc_stokes_alt", ERROR_POSITION);

          /* only continue if both stokes vectors are > 0 */
          if (p_escape->stokes[0] > 0.
           && p_escape->cb.stokes_alt[0] > 0.
           && p_escape->cb.tausca < 709.) { /* quick fix to prevent inf in lidar.c in certain situations */
            if (!sample->nocb_switch) {
              calc_stokes_cb (&(p_escape->cb), p_escape->stokes);
              result->lidcb->cPLD += 1u;
              result->lidcb->meanPLD += p_escape->cb.pld;
            }
            p_escape->cb.cb = 1;
          }
        }
      }
      /* End CB */

      /* propagate virtual photon to the lidar detector and add contribution to local estimate */
      status = escape_probability (p_escape, atmos,
                                   sample, result,
                                   elev, 0, -1, it,
                                   0, refind);
      if (status)
	return fct_err_out (status, "escape_probability", ERROR_POSITION);

      /* free memory */
      destroy_photon (p_escape, atmos->n_caoth);
    }
  }
#endif

  /***********************************/
  /* 2. determine what is scattering */
  /***********************************/

  isp = random_scatter_type (atmos, p);
  
  /*******************************/
  /* 3. if VIRTUALSCATTER return */
  /*******************************/

  if (isp == atmos->n_caoth+1) {
    p->photon_status = MCSTATUS_VIRTUALSCATTER;
    return MCSTATUS_VIRTUALSCATTER;
  }

  /* n-tupel local estimate; clone photon no longer needed */
  if ( p->isclone && p->clonescattercounter+2 == sample->ntupelLE ) {
    p->photon_status = MCSTATUS_PURGE;
    return MCSTATUS_PURGE;
  }

#if HAVE_VROOM
  /*******************/
  /* 4. DDIS stuff 1 */
  /*******************/

  if ( sample->escape_eps_ddis_upf || sample->LLE_D_DIS ) {
    status = mc_vroom_prep_DDIS (sample, p, atmos, &DDISsing, &FODDISsing, &out_of_cone, &lest, &scadis);
    if (status!=0)
      return err_out ("Error %d returned by mc_vroom_prep_DDIS ()\n", status);
  }
#endif

  /**************************************/
  /* 5. Calculate scattering angle: phi */
  /**************************************/

#if HAVE_LIDAR
  if (FODDISsing)
    phi = scadis.d_phi * ( 2.0*uvspec_random() - 1.0 );
  else
    /* standard case */
    phi = sc_Isotropic_phi();
#else
  phi = sc_Isotropic_phi();
#endif

  /*************************************/
  /* 6. Calculate scattering angle: mu */
  /*************************************/

#if HAVE_VROOM
  status = mu_scatter_special (atmos, p,
			       sample->phase_max, sample->n_phase_max, isp, &mu,
			       scadis, DDISsing, out_of_cone);
  if (status!=0)
    return err_out ("Error, mu_scatter_special() returned status %d\n", status);
#else
  /* standard case */
  status = mu_scatter (atmos, p, isp, &mu);
  if (status!=0)
    return err_out ("Error, mu_scatter() returned status %d\n", status);
#endif

#if HAVE_LIDAR
  mu_save = mu;
  if (sample->coherent_backscatter && sample->LidarLocEst) {
    result->lidcb->mu_ave += mu;
    result->lidcb->mu_counter += 1;
  }
#endif


  /********************************************/
  /* 7. calculate new direction               */
  /*    NOTE: definition of phi is different! */
  /********************************************/
  
  /* If polarization is calculated the new direction depends on the incoming phi, */
  /* even if the direction of the photon is exactly parallel to the z-axis. */
  /* This happens usually at the first scattering, if umu = +-1 (backward) or */
  /* phi0=0 (forward)*. If the photon is scattered randomly exactly in z direction, */
  /* the plane of polarization will be lost (this is highly unlikely). */

  /* first scattering */
  if (p->scattercounter==0)
    phi_inc=p->phi0;
  else
    phi_inc=0.0; 
  
#if HAVE_VROOM
  /****************/
  /* 4c. for both */
  /****************/

  /* save photon direction for later use and turn it into the direction of the cone center if DDIS active */
  for (i=0; i<3; i++)
    scadis.dirold_dx[i] = p->dir.dx[i];

  if (DDISsing)
    for (i=0; i<3; i++)
      p->dir.dx[i] = lest.dir.dx[i];

#endif

  if (atmos->i_ic[isp] > 0 && atmos->n_raytracing_prop>0)
    new_direction_raytracing (atmos, mu, phi - 90., &(p->dir), phi_inc, wvnmlo, wvnmhi);
  else 
    new_direction (mu, phi - 90., &(p->dir), phi_inc);

#ifdef MUCHOUT
  if (p->muchoutcounter==MUCHOUTPHOTON) {
    if (sample->escape)
      v_mult_mu ( p->dir.dx, sample->rad[0].dir.dx, &mua );
    if (sample->LidarLocEst)
      v_mult_mu ( p->dir.dx, sample->lidar[sample->ili].dir.dx, &mua );

    fprintf(stderr,"counters %d %d %d %d --- newdir: %e %e %e mu %e box %d %d %d mu %e phi %e\n",
	    p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
	    p->dir.dx[0], p->dir.dx[1], p->dir.dx[2], mua, p->ic, p->jc, p->kc,mu,phi);
  }
#endif

#if HAVE_VROOM
  /*******************/
  /* 8. DDIS stuff 2 */
  /*******************/

  if (sample->LLE_D_DIS || (sample->escape_eps_ddis_upf != 0.0) ) {
    status = mc_vroom_set_mus_and_phis (sample, p, DDISsing,
					lest, scadis,
					&mu, &mu2, &stheta2, phi, &cphi2);
    if (status!=0)
      return err_out ("Error %d returned by mc_vroom_set_mus ()\n", status);
  }

  if ( (p->RIS_mode != MCRIS_MODE_NORMAL) ||  /* RIS                 */
       sample->LLE_D_DIS ||                   /* DDIS for LIDAR      */
       sample->escape_eps_ddis_upf ||         /* DDIS for escape     */
       sample->LLE_channels ||                /* Raman/HSRL          */
       sample->LLE_jacobian ||                /* jacobian LIDAR      */
       sample->abs_jacobian ||                /* jacobian SCIA       */
       atmos->ris_factor != 1. ) {           /* is ris-factor used? */

    status = mc_vroom_scattering_calc_phases_and_jacobians ( sample, p, atmos, mu, 0,
							     &P_norm, &P_spec, &P_isoene );
    if (status!=0)
      return err_out ("Error %d returned by mc_vroom_scattering_calc_phases_and_jacobians ()\n", status);

    status = mc_vroom_DDIS_weight_and_prep_stuff (sample, atmos, mu2, stheta2, cphi2,
						  DDISsing, out_of_cone,
						  P_norm, P_spec, P_isoene, lest, scadis, p);
    if (status!=0)
      return err_out ("Error %d returned by mc_vroom_weight_and_prep_stuff ()\n", status);

  }

#endif

  /******************************************************/
  /* 9. calculate polarization, Stokes vector weight    */
  /******************************************************/
  
  if (sample->polarisation){
    if (sample->use_phase_matrix_total) {
      /* another version, evtl. better statistics, but a bit slower */
      /* the scatter type averaged phase_matrix is used */
      status = get_phase_matrix_total ( atmos, p, mu,  /* correct?? */
                                        NPHAMAT, 0, sample->spectral_is,
                                        sample->concentration_is, atmos->ris_factor, 0,
                                        phase_matrix, NULL, NULL, &(p->weight) );
      if (status)
        return fct_err_out (status, "get_phase_matrix_total", ERROR_POSITION);
    }
    else {
      status = get_phase_matrix_caoth ( atmos, p, isp, mu,
  			              NPHAMAT, phase_matrix );
      if (status!=0)
        return fct_err_out (status, "get_phase_matrix_caoth", ERROR_POSITION);
    }

    /* CB - Chris */
#if HAVE_LIDAR
    if (sample->coherent_backscatter && sample->LidarLocEst) {
      add_coord_to_cohebasca ( p );
      result->lidcb->mu_ave += mu_save;
      result->lidcb->mu_counter += 1;
      if (p->scattercounter == 0 )
        p->cb.phamfirst = phase_matrix[0];
    }
#endif

    if (polmat){
      phase_matrix_mult ( p->phamat, phase_matrix,
			  p->dir.dx, pdir_inc,
                          p->fw_phi0, p->fw_phi,0,
			  p->scattercounter, sample->backward );
      
      /* Z*I_0 this weight is only used if photon is counted after this
         scattering event, when p->phamat is total phase matrix */
      for(ip=0; ip<4; ip++) {
	p->stokes[ip] = 0.0;
	for(j=0; j<4; j++)
	  p->stokes[ip] += p->phamat[ip][j]*p->stokes0[j];
      }

#if HAVE_LIDAR
    /* CB - get the correct phase matrices for the reverse path */
      if (sample->coherent_backscatter && sample->LidarLocEst && p->scattercounter > 0) {

        if (p->cb.Zmiddlealt[0][0] == 0.)
          for (i=0;i<4;i++)
            p->cb.Zmiddlealt[i][i] = 1.;

        Ztemp = calloc(4, sizeof(double *));
        for (i=0; i<4; i++)
          Ztemp[i] = calloc(4, sizeof(double));

        for (i=0; i<4; i++)
          for (j=0; j<4; j++)
            Ztemp[i][j] = p->cb.Zmiddlealt[i][j];

        phase_matrix_mult(Ztemp, phase_matrix,
			  p->dir.dx, pdir_inc,
                          p->fw_phi0, p->fw_phi, 0,
			  p->scattercounter, 1-sample->backward);
        
        for (i=0; i<4; i++)
          for (j=0; j<4; j++)
            p->cb.Zmiddlealt[i][j] = Ztemp[i][j];

        for (i=0; i<4; i++)
          free(Ztemp[i]);
      }
#endif
    }  
    else {
      stokes_vector_sca(p->stokes, phase_matrix,  p->dir.dx, pdir_inc, p->fw_phi0, p->fw_phi); 
    }
  }
  
#if HAVE_LIDAR
  free(Ztemp);
#endif
  p->photon_status=MCSTATUS_INCREMENT;
  return MCSTATUS_INCREMENT;
}  


/***********************************************************************************/
/* Function: random_scatter_type                                          @62_30i@ */
/* Description:                                                                    */
/*  determine what is scattering                                                   */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static inline int random_scatter_type (atmosphere_struct *atmos, photon_struct *p)
{
  double alpha=0.0;
  double kscaact=0.0;
  int isp=MCCAOTH_NONE;

  alpha = uvspec_random() * get_kext (atmos, p);
  for (isp=1; isp <= atmos->n_caoth; isp++) {
    kscaact = get_kscaIS (atmos, p, isp);
    if (kscaact > alpha)
      break;
    alpha -= kscaact;
  }

  return isp;
}


/***********************************************************************************/
/* Function: mu_scatter                                                   @62_30i@ */
/* Description:                                                                    */
/*  Return mu for the scattering depending on the scatter type isp                 */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int mu_scatter ( atmosphere_struct    *atmos,
		 photon_struct        *p,
		 int                   isp,
		 double               *mu )
{

  switch(atmos->scatter_type[isp]) {
  case MCSCAT_MOL:    /* Rayleigh scattering */
    p->rayleighcounter++;
    *mu = sc_Rayleigh_mu_depol ( atmos->rayleigh_depol );
    break;
  case MCSCAT_AER:   /* Aerosol scattering   */
    /* ??? need to include phase matrix CHECK!!! */ 
    *mu = sc_mu ( &(atmos->phase_aer[p->kc]), 1, 0, 0., 0. );
    break;
  case MCSCAT_HG1:
    *mu = sc_HG_mu ( get_g1(atmos, p, isp) );
    break;
  case MCSCAT_HG2:
    *mu = sc_HG2_mu ( get_g1(atmos, p, isp),
		      get_g2(atmos, p, isp),
		      get_ff(atmos, p, isp) );
    break;
  case MCSCAT_PFT:
    if (atmos->phase[isp]->n == 1)
      *mu = sc_mu (atmos->phase[isp]->iphase[0], 1, p->SC_mode, 0., 0.);
    else
      *mu = sc_interp_mu (get_reff(atmos, p, isp), atmos->phase[isp], p->SC_mode);
    break;
  default:    
    fprintf(stderr,"Error, no such type %d of scattering!!!\n", atmos->scatter_type[isp]);
    return -1;
  }

  return 0;
}


/***********************************************************************************/
/* Function: get_phase_matrix_total                                       @62_30i@ */
/* Description:                                                                    */
/*  Calculate the total probability for scattering angle mu = cos(theta).          */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int get_phase_matrix_total ( atmosphere_struct *atmos,
			     photon_struct     *p,
			     double             mu,
			     int                np,
			     int                vismode,
			     int                spectral_is,
			     int                concentration_is,
                             double             ris_factor,
			     int                release,
			     /* Output */
			     double            *phase_matrix,
			     double           **P_caoth,
			     double            *P_spec,
			     double            *weight )
{
  int isp=0, ip=0, ic=0;
  static double *ksca=NULL;
  static double **P_save=NULL;
  double norm=0.0, ksca_special=0.0;
  int rismode=0;
  int status=0;
  /*  int iv=0; */

  /* set rismode */
  if (ris_factor != 1.)
    rismode=1;
  if (p!=NULL)
    if (p->RIS_mode==MCRIS_MODE_MAS)
      rismode=1;

  /* this function is called once at end of mystic in order to free memory */
  if (release==1) {
    if (ksca!=NULL) {
      free(ksca);
      ksca=NULL;
    }
    if (P_save!=NULL) {
      for (isp=0; isp<=atmos->n_caoth; isp++)
	free(P_save [isp]);
      free(P_save);
      P_save=NULL;
    }
    return 0.0;
  }

  /* allocate if first call */
  if (ksca==NULL)
    ksca = calloc((size_t) atmos->n_caoth+1, sizeof(double));
  if (P_save==NULL) {
    P_save = calloc((size_t) atmos->n_caoth+1, sizeof(double *));
    for (isp=0; isp<=atmos->n_caoth; isp++)
      P_save [isp] = calloc(NPHAMAT, sizeof(double)); 
  }

  /* point P_caoth to local variable if not given by function call */
  if (P_caoth==NULL)
    P_caoth=P_save;

  /* get scattering coefficients for each caoth */
  for (isp=1; isp<=atmos->n_caoth; isp++)
    ksca [isp] = get_ksca (atmos, p, isp);

  /* normalization */
  norm = 1. / get_ksca (atmos, p, MCCAOTH_TOT);

  if (vismode == 0)
    ksca_special = get_kscaIS (atmos, p, MCCAOTH_TOT);
  else
    ksca_special = get_kext (atmos, p);

  /* need to add to weight!!! */
  if (weight!=NULL)
    *weight /= norm * ksca_special;

  /* calculate phase matrix/function for all caoths */
  for (isp=1; isp<=atmos->n_caoth; isp++) {
    if ( ksca[isp] > 0.0 ) {
      status= get_phase_matrix_caoth ( atmos, p, isp, mu, np, P_caoth[isp] );
      if (status)
	return fct_err_out (status, "get_phase_matrix_caoth", ERROR_POSITION);
    }
  }

  /* P_spec is needed for VROOM */
  if (P_spec!=NULL)
    if (rismode) {
      *P_spec = 0.0;
      for (isp=1; isp<=atmos->n_caoth; isp++)
	*P_spec += get_kscaIS (atmos, p, isp) * P_caoth[isp][0];
      *P_spec /= ksca_special;
    }

  /* add weight to P_caoth */
  for (isp=1; isp<=atmos->n_caoth; isp++)
    for(ip=0; ip<np; ip++)
      P_caoth [isp][ip] *= ksca [isp] * norm; 

  /* initialize scattering phase matrix */
  for(ip=0; ip<np; ip++)
    phase_matrix [ip] = 0.0;

  /* add all caoths to scattering phase matrix */
  for (isp=1; isp<=atmos->n_caoth; isp++)
    for(ip=0; ip<np; ip++)
      phase_matrix [ip] += P_caoth [isp][ip];

  /* save total scattering phase matrix in P_caoth */
  for(ip=0; ip<np; ip++)
    P_caoth [MCCAOTH_TOT][ip] = phase_matrix [ip];

  /*  Weight to correct phase function for ALIS */
  /*   Uncomment the following if you do not want to neglect it. Caution: CPU time increases by a factor of app. 1.5 !!!! */
  /*   Also it is not yet clear whether it works with vroom */
  /*  if (spectral_is){ */
  /*     for (iv=0; iv<atmos->nlambda_abs; iv++){ */
  /*           p->q2_spectral[iv]*= */
  /*             (paer + pcld + pice + */
  /*              atmos->ksca_spectral[iv][p->kc]*Rayleigh_depol (mu, atmos->rayleigh_depol))/ */
  /*             (get_ksca (atmos, p, MCCAOTH_TOT) - get_ksca (atmos, p, MCCAOTH_MOL) + atmos->ksca_spectral[iv][p->kc])/ */
  /*             (( pcld + pice + pmol + paer ) * norm); */
  /*     } */
  /*   } */
  
  /* Phase function correction is very important for aerosol concentration importance sampling !!! */
  if (concentration_is){
    for (ic=0; ic<atmos->Nc; ic++){
      
      /*  if ((P_caoth [MCCAOTH_TOT][0] != 0.0) && */
      /*           (get_ksca (atmos, p, MCCAOTH_TOT)-ksca [MCCAOTH_AER] + atmos->ksca_scaled[ic][p->kc]) !=0.0) { */
      if (ksca [MCCAOTH_AER]!=0.0){
        p->q2_concentration[ic]*=
          (P_caoth [MCCAOTH_TOT][0] + (-1.0 + atmos->ksca_scaled[ic][p->kc]/ksca [MCCAOTH_AER]) *  P_caoth [MCCAOTH_AER][0] ) /
          P_caoth [MCCAOTH_TOT][0] *
          get_ksca (atmos, p, MCCAOTH_TOT) /(get_ksca (atmos, p, MCCAOTH_TOT)-ksca [MCCAOTH_AER] + atmos->ksca_scaled[ic][p->kc]);
      }
    }
  }
  
  if (P_spec!=NULL)
    /* set P_spec for non-special cases */
    if (!rismode)
      *P_spec=phase_matrix [0];

#ifdef MUCHOUT
  if (p->muchoutcounter==MUCHOUTPHOTON) {
    for (isp=1; isp<=atmos->n_caoth; isp++)
      fprintf(stderr,"counters %d %d %d %d --- caoth %d kscapprof %e norm %e mu %e\n",
	      p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
	      isp,P_caoth[isp][0],norm,mu);
  }
#endif

  return 0;
}


/***********************************************************************************/
/* Function: register_at_start                                            @62_30i@ */
/* Description:                                                                    */
/*  Register a photon at the start of its journey.                                 */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int register_at_start (photon_struct *p, sample_struct *sample,
                              atmosphere_struct *atmos, result_struct *result,
                              int source, int absorption)
{

  double *weight; /* Claudia, please tell me this isn't a bug!!! CHECK!!! */ 
  weight=calloc(4, sizeof(double)); 
  weight[0]=1.0;
  
  switch (source) {
  case MCSRC_SOLAR:

    if (sample->sample[atmos->Nz])
      count_photon (result->alt[atmos->Nz], p, sample, weight, 0, 0);
    
    break;
    
  case MCSRC_THERMAL_SURFACE:
    
    /* count outgoing photon */
    /* CE: is outgoing photon polarized?? if yes we may need to change the weight vector. */
    /* probably it is not polarized CHECK!!! */ 
    count_photon (result->surf, p, sample, weight, 0, 0);
    
    /* additionally add to altitude profile if sample[0]   */
    /* is switched on                                      */
    /* ??? careful: does not work with 2D elevation!!! ??? CHECK!!! */
    
    if (sample->sample[0])
      count_photon (result->alt[0], p, sample, weight, 0, 0);
    
    break;
    
  case MCSRC_THERMAL_BACKWARD:/*TZ bt*/
  case MCSRC_THERMAL_ATMOSPHERE:
    
    /* photon is not added to any flux or radiance      */
    /* counter because it is highly unlikely that       */
    /* the photon started exactly at a layer            */
    /* boundary                                         */

    /* However, if absorption is MCFORWARD_ABS_HEATING, */
    /* we need to subtract the photon from the absorbed */
    /* energy (= heating rate) because an emitted       */
    /* photon acts as cooling                           */
    
    if (absorption == MCFORWARD_ABS_HEATING)
      if (atmos->threed[MCCAOTH_TOT][p->kc]>=1) {
        result->absorption3D->tot[p->kc][p->ic][p->jc] -= 1.0;

        if (sample->std)
          p->wtree = addtree_stddev (p->wtree, &(result->absorption3D2->tot[p->kc][p->ic][p->jc]), -1.0);
      }

    break;

  case MCSRC_LIDAR:
  case MCSRC_BLITZ:
    break;

  default:
    fprintf (stderr, "Error, unknown or unimplemented MC source type\n");
    return -1;
  }

  free(weight); 
  
  return 0;
}


#ifdef CLDPRP

/***********************************************************************************/
/* Function: sample_cldprp                                                @62_30i@ */
/* Description:                                                                    */
/*  Sample cloud properties (reff, tau, lwc gradient) along the photon path        */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Florian Ewald                                                           */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int sample_cldprp (photon_struct *p, atmosphere_struct *atmos, double step)
{
  int isp=0.0;
  double tau=0.0;
  double reff=0.0;
  
  double taubox_wc  = 0.0;
  double taubox_ic  = 0.0;
  double reffbox_wc = 0.0;
  double reffbox_ic = 0.0;
  double rhitbox_wc = 0.0; 
  double rhitbox_ic = 0.0; 
  double dxbox_wc   = 0.0;
  double dxbox_ic   = 0.0;
  double dybox_wc   = 0.0;
  double dybox_ic   = 0.0;
  double dzbox_wc   = 0.0;
  double dzbox_ic   = 0.0;
  
  /* start sampling cloud properties at first cloud contact or continue
   sampling when photon had already cloud contact - no sampling if photon
   was cloned/scattered before cloud contact, to avoid values for clearsky pixels  */
  
  if ((p->isclone == 0 && p->scattercounter == 0) || (p->cldprp.rhit_wc+p->cldprp.rhit_ic) > 0.0){
    
    /* loop over n profiles to sum up cloud properties in this cloud box */
    for (isp=MCCAOTH_FIR; isp<=atmos->n_caoth; isp++){
      if (atmos->i_wc[isp] || atmos->i_ic[isp]){ 
        
        /* reset values for every new profile */
        reff   = 0.0;
        tau    = 0.0;
        
        /* lookup tau to test if in cloud */
        tau  = step * (get_ksca (atmos, p, isp) + get_ksca (atmos, p, isp));
        
        /* lookup effective radius in cloud */
        reff = get_reff (atmos, p, isp);
        
        if (atmos->i_wc[isp]){
          
          /* integrate effective radius and tau along tau */
          reffbox_wc += tau*reff;
          taubox_wc  += tau;
          
          /* lookup effective radius for first cloud contact */
          if (p->cldprp.rhit_wc == 0.0 || rhitbox_wc)
            rhitbox_wc += tau*reff;   
          
          /* lookup lwc gradient if photon has been scattered in the cloud for the first time */
          if (p->photon_status == MCSTATUS_SCATTER){          
            if ((p->cldprp.dxlwc+p->cldprp.dylwc+p->cldprp.dzlwc) == 0 || (dxbox_wc+dybox_wc+dzbox_wc)){
              dxbox_wc = get_dxlwc (atmos, p, isp);
              dybox_wc = get_dylwc (atmos, p, isp);
              dzbox_wc = get_dzlwc (atmos, p, isp);
            }
          }
          
        }
        
        if (atmos->i_ic[isp]){
          
          /* integrate effective radius and tau along tau */
          reffbox_ic += tau*reff;
          taubox_ic  += tau;
          
          /* lookup effective radius for first cloud contact */
          if (p->cldprp.rhit_ic == 0.0 || rhitbox_ic)
            rhitbox_ic += tau*reff;   
          
          /* lookup lwc gradient if photon has been scattered in the cloud for the first time */
          if (p->photon_status == MCSTATUS_SCATTER){          
            if ((p->cldprp.dxiwc+p->cldprp.dyiwc+p->cldprp.dziwc) == 0 || (dxbox_ic+dybox_ic+dzbox_ic)){
              dxbox_ic += tau*get_dxlwc (atmos, p, isp);
              dybox_ic += tau*get_dylwc (atmos, p, isp);
              dzbox_ic += tau*get_dzlwc (atmos, p, isp);
            }
          }
          
        }
        
      }
    }
    
    /* sum up the integrated tau for all profiles in this cloud box */
    p->cldprp.tau_wc  += taubox_wc;
    p->cldprp.tau_ic  += taubox_ic;
    
    /* sum up the weighted reff for all profiles in this cloud box */    
    p->cldprp.reff_wc += reffbox_wc;
    p->cldprp.reff_ic += reffbox_ic;
    
    /* register effective radius for first cloud contact */
    if (taubox_wc)
      p->cldprp.rhit_wc = rhitbox_wc/taubox_wc;
    if (taubox_ic)
      p->cldprp.rhit_ic = rhitbox_ic/taubox_ic;
    
    /* register lwc gradient if photon has been scattered in the cloud for the first time*/
    if (dxbox_wc+dybox_wc+dzbox_wc){
      p->cldprp.dxlwc = dxbox_wc/taubox_wc;
      p->cldprp.dylwc = dybox_wc/taubox_wc;
      p->cldprp.dzlwc = dzbox_wc/taubox_wc;
    }
    
    if (dxbox_ic+dybox_ic+dzbox_ic){
      p->cldprp.dxiwc = dxbox_ic/taubox_ic;
      p->cldprp.dyiwc = dybox_ic/taubox_ic;
      p->cldprp.dziwc = dzbox_ic/taubox_ic;
    }
    
  }
  
  return 0;
}

#endif
  

/***********************************************************************************/
/* Function: travel_tau                                                   @62_30i@ */
/* Description:                                                                    */
/*  Travel optical depth tau, or to the next surface or boundary.                  */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: ...., last big modification by Robert Buras                             */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int travel_tau (photon_struct *p, atmosphere_struct *atmos,
		sample_struct *sample, result_struct *result,
		elevation_struct *elev, double tau, int absorption,
		int source, int photonpath, int visualize, float *refind,
		int quiet)
{
  double tausca=0.0, step=0.0, gamma=0.0;
  double ksca=0.0;
  int status=0;
  double omegal=0.0;
  double totweight[4];

  int is_threed=0;
  int ip=0, iv=0, ic=0;
  int intersection_counter=0;
  int need_VIS=0;
  double kscamin=0.0;

#ifdef HAVE_SPHER
  int first=1;
#endif
#if HAVE_LIDAR
  int periodic=0;
  int moving_towards_detector=0, need_preparation=0, need_action=0;
  struct qidd_struct *qidd;

  qidd = calloc ((size_t) 1, sizeof(qidd_struct));
  qidd->recalc = 1;

  need_preparation = (sample->LLE_RIS_MAS || sample->LLE_VIS_QIDD);
#endif

  /* initialize status; if unchanged at end of travel step, the photon crosses a grid boundary */
  p->photon_status = MCSTATUS_TRAVEL;

#if HAVE_LIDAR
  if ( need_preparation )
    mc_lidar_prepare_travel_tau (sample, p, &moving_towards_detector);
#endif

  
  while (1==1) {  /* until photon reaches optical depth tau */

#ifdef MUCHOUT
    if (p->muchoutcounter==MUCHOUTPHOTON)
      fprintf(stderr,"traveltau x %e %e %e i %d %d %d dir %e %e %e tau %e %e scaco %d weight %e\n",p->x[0],p->x[1],p->x[2],p->ic,p->jc,p->kc,p->dir.dx[0],p->dir.dx[1],p->dir.dx[2],tau,tausca,p->scattercounter,p->weight);
#endif

    is_threed = (atmos->threed[MCCAOTH_TOT][p->kc]>=1);
    
    p->RIS_mode = MCRIS_MODE_NORMAL;
    p->VIS_mode = MCVIS_MODE_NORMAL;

#if HAVE_LIDAR
    periodic=0;
    need_action = ( p->lest.will_hit_det_plane || p->lest.will_hit_cone || sample->LLE_VIS_QIDD );

    if (source != 99 && need_preparation ) {
      status = mc_lidar_prepare_travel_tau_step (sample, p, qidd, atmos, tau - tausca, moving_towards_detector);
      if (status!=0)
	return err_out ("Error %d returned by mc_lidar_prepare_travel_tau_step()\n", status);
    }

    if ( sample->RIS_MS && !p->scattercounter )
      p->RIS_mode = MCRIS_MODE_MAS;

#endif

    /**************************************************************/
    /* 1. find out what happens next with photon, and step length */
    /**************************************************************/

    /* calculate next intersection point with the box boundaries   */
    if (sample->spherical3D) {
#ifdef HAVE_SPHER
      status = intersection3D_spherical (p, atmos,
					 sample->spherical3D_scene,
					 &first, &step);
      if (status)
	return err_out("Error %d returned by intersection3D_spherical\n", status);
#else
      fprintf(stderr,"Error! you are not allowed to use spherical 3D!\n");
      return -1;
#endif
    }
    else {
      if (is_threed==1) { /* 3D layer */
#if HAVE_MYSTIC3D
	status = intersection3D (p, atmos, tau, tausca, &step);
	if (status)
	  return err_out("Error %d returned by intersection3D\n", status);
#else
	fprintf(stderr,"Error! you are not allowed to use mystic 3D!\n");
	return -1;
#endif
      }
      else {
	if (sample->spherical){
	  status = intersection1D_spherical (p, atmos,
					     sample->refraction, refind,  &step);
	  if (status<0) {
	    fprintf (stderr, "\n");
	    fprintf(stderr, "*** Warning: Photon absorbed in intersection1D_spherical. \n");
	    fprintf (stderr, "\n");
	  
	    p->photon_status = MCSTATUS_OUTOFDOMAIN;
	    break;
	  }
	}
	else
	  intersection1D (p, atmos, tau, tausca, &step, source);
      }
    }
    intersection_counter++;

    if (intersection_counter==30000 && !quiet)
      fprintf(stderr,"Warning! Photon passed 30000 intersections, probably it has entered an infinite loop!\n");

#if HAVE_LIDAR
    if ( need_action ) {
      status = mc_lidar_travel_tau_action (p, &step, qidd);
      if (status == 1)
	continue;
    }
#endif

    if (sample->VIS_CS)
      p->vis_beta = 0.0; /* dangerous, might destroy lidar ! */

    /* set extinction coefficient */
    ksca = get_kext (atmos, p);

    /* for thermal backward absorption is included in extinction */
    if (source==MCSRC_THERMAL_BACKWARD)
      ksca += get_kabs_tot (atmos, p);

    /* do some VIS if using VROOM */
    if (sample->VIS_CS) {
      /* check whether dtau VIS needed */
      if (step * ksca < 1.0 && p->p_norm > 10.0 ) {
	/* Check whether box is a spiky one */
	if (atmos->threed[MCCAOTH_TOT][p->kc]>=1)
	  need_VIS = atmos->spiky_box[p->kc][p->ic][p->jc];
	else
	  need_VIS = atmos->spiky_box[p->kc][0][0];
	if (need_VIS) {
	  kscamin = 0.001 * p->p_norm;
	  if (kscamin > 1.0)
	    kscamin = 1.0;
	  kscamin /= step;
	  if (ksca < kscamin) {
	    ksca = kscamin;
	    p->vis_beta = kscamin;
	  }
	}
      }
    }

    /* passed next scattering point on the way to the box boundary */
    if (tausca + step * ksca > tau) {
      step = ( tau - tausca ) / ksca;
      p->photon_status = MCSTATUS_SCATTER;
      
      if (source==MCSRC_THERMAL_BACKWARD) {
	/* decide whether scattering or absorption happened */
	omegal = get_kext (atmos, p) / ksca; /* RPB does not know what RIS or VIS do here! CHECK!!! */
        if (sample->spectral_is){
          for (iv=0; iv<atmos->nlambda_abs; iv++){
            p->q_spectral[iv]=atmos->kabs_spectral[MCCAOTH_TOT][iv][p->kc]/
              get_kabs_tot (atmos, p);
          }
        }
        if (sample->concentration_is){
          for (ic=0; ic<atmos->Nc; ic++){
            p->q_concentration[ic]=( get_kabs_tot (atmos, p) 
                                     - (atmos->kabs->prof [MCCAOTH_AER])[p->kc]+atmos->kabs_scaled[ic][p->kc])/
              get_kabs_tot (atmos, p);
            /*   fprintf(stderr, "FIXCE, ic %d, tot %g, aer %g, scaled %g q %g \n", CHECK!!! */
/*                     ic, get_kabs_tot (atmos, p), (atmos->kabs->prof [MCCAOTH_AER])[p->kc], */
/*                     atmos->kabs_scaled[ic][p->kc], p->q_concentration[ic]); */
          }
        }
        if (uvspec_random() > omegal){
	  p->photon_status = MCSTATUS_ABSORB; /* real absorption reached flag*/
        }
      }
    }

    if (step<-MC_EPSILON) {
      fprintf (stderr, "\n!!FATAL error in travel_tau(), step %g < 0 in %dD box (%d, %d, %d).\n",
               step, (is_threed==1?3:1), p->ic, p->jc, p->kc);
      fprintf (stderr,   "Start at (%g, %g, %g), direction (%g, %g, %g)\n", 
               p->x[0], p->x[1], p->x[2], p->dir.dx[0], p->dir.dx[1], p->dir.dx[2]);
      fprintf (stderr,   "hit (%d, %d, %d), cotheta %f, status %d\n", 
               p->dir.hit[0], p->dir.hit[1], p->dir.hit[2], p->dir.cotheta, p->photon_status);
      fprintf (stderr,   "tau %e, tausca %e, beta_local %e # of intersections %d\n", 
               tau, tausca, ksca, intersection_counter);

      
      return -1;
    }
    
    if (step<0)
      step=0;
      
    /* check if photon crossed the 2D surface */

    if (elev->elev2D)  {
      
      status = cross_surface (p, step,
                              elev,
                              sample->bcond,
                              &gamma);
      
      
      /* if cross_surface ran into some inconsistency */
      if (status<0) {
        fprintf (stderr, "Something terrible has happened! cross_surface() returned\n");
        fprintf (stderr, "status %d\n", status);
        return -1;
      }
      
      if (status>0) {
#ifdef MUCHOUT
	if (p->muchoutcounter==MUCHOUTPHOTON)
	  fprintf(stderr,"crossing surface at \n");
#endif
        step = gamma;
        p->photon_status = MCSTATUS_SURFACE;
      }
    }

    
#ifdef CLDPRP
    
    /*******************************************************************/
    /* 1b. integration for visible cloud property estimation (F.Ewald) */
    /*******************************************************************/
    
    if (sample->cldprp){
      
      status = sample_cldprp (p, atmos, step);      
      
      if (status)
        return err_out ("ERROR! sample_cldprp returned status %d\n",status);
      
    }
    
#endif
    
    
    /*******************************************************************/
    /* 2. apply step length, move photon, and exit if photon path ends */
    /*******************************************************************/

    /* add absorption to absorption counter */
    if (absorption && is_threed==1)  /* 3D layer */
      mc_add_absorption_3D (atmos->kabs3D, atmos->kabs, /* RPB does not know what RIS or VIS do here! CHECK!!! */
                            p->kc, p->ic, p->jc, step, exp(-p->tauabs.tot) * p->weight * p->stokes[0], sample->std,
                            result->absorption3D, result->absorption3D2, &(p->wtree));
    
    /* sum up absorption optical depth */
    /* Please, Tobias and Bernhard, look at the following line:*/
    if (source!=MCSRC_THERMAL_BACKWARD) /*TZ bt*/ /* CAUTION! THIS LINE WAS ADDED BY RPB! CHECK WHETHER IT WAS CORRECT TO DO SO! CHECK!!! */
      mc_add_optical_depth (atmos, p, step, sample->coherent_backscatter);

    /* add step to pathlength */
    p->pathlength += step;

#if HAVE_LIDAR
    tausca += step * ksca_for_tausca (p, &step, &ksca, qidd);
    if ( need_action ||
	 sample->LLE_jacobian ||
	 sample->abs_jacobian ||
	 p->RIS_mode != MCRIS_MODE_NORMAL ||
	 atmos->ris_factor != 1. )
      mc_lidar_travel_tau_addition (atmos, sample, p, &step, qidd, &ksca);
#else
    /* update passed add optical depth */
    tausca += step * ksca;
#endif

    if (sample->spectral_is || sample->concentration_is || sample->boxairmass){ 
      /* original method: */
      /*       dtauabs_spectral_calc(atmos, p, step); */
      p->pathlength_per_layer[p->kc] += step;
    }
    
    /* move photon to box boundary, to next scattering point,  */
    /* or to the 2D surface                                    */
    /* TZ bt: or to absorption point (backward thermal)        */
    
    if (is_threed==1 || sample->spherical3D) { /* 3D layer */
#if HAVE_MYSTIC3D
      status = step3D (p, atmos, step, sample->bcond, photonpath, sample->spherical3D, visualize);
#else
      fprintf(stderr,"Error! you are not allowed to use mystic 3D!\n");
      return -1;

#endif
    }
    else
      status = step1D (p, atmos, step, sample->bcond, photonpath, visualize);
    
    if (status<0) {
      fprintf (stderr, "Error %d returned by step%dD()\n", status,
               (is_threed==1?3:1));
      return status;
    }

   
#if HAVE_LIDAR
    if (status==1)
      periodic=1;

    if ( need_action )
      mc_lidar_travel_special_cases (sample, p, qidd, &ksca, tau - tausca);

    /* travelled through periodic boundary condition */
    if (periodic==1 && (sample->LLE_VIS_FOD || sample->LLE_eps_ddis_uda) ) {
      status = initialize_vis_fod (sample, p);
      if (status)
	return err_out ("ERROR! initialize_vis_fod returned status %d\n",status);
    }
#endif

    /* end photon travel unless it simply crosses the boundary or hits the FOV cone */
    if (p->photon_status != MCSTATUS_TRAVEL)
      break;
    
    /********************************/
    /* 3. count photon if necessary */
    /********************************/

    /* count photon at horizontal boundary */
    /* TZ bt: If not backward thermal mode! Then the photons are only */
    /* counted once - when absorbed - and not when passing by. */
    if (sample->sample[p->kc + (1-p->dir.hit[2])] && source!=MCSRC_THERMAL_BACKWARD) {/*TZ bt*/
      if (p->reallyhit[2] || (is_threed==0) ) {
        
        /* count passing photon */
        for (ip=0; ip<sample->nstokes; ip++)
          totweight[ip] = p->weight * p->stokes[ip] * exp(-p->tauabs.tot);
        
        count_photon (result->alt[p->kc + (1-p->dir.hit[2])], 
                      p, sample, totweight, 0, 0);
      }
    }
    
    /**********************************/
    /* 4. check if photon leaves grid */
    /**********************************/

    /* finally check if photon arrived at the lower surface */
    /* or at the top of the atmosphere                      */
    
    if (p->kc<0) {
      p->kc=0;
      
      /* to avoid roundoff errors */
      if (!sample->spherical && !sample->spherical3D)   /* This does not work in spherical geometry */
        p->x[2] = atmos->Z[p->kc];
      
      p->photon_status = MCSTATUS_BOUNDARY_LOWER;

      if (elev->elev2D) {
        fprintf (stderr, "FATAL error! Photon reached the 1D boundary in a\n");
        fprintf (stderr, "             2D elevation case (%dD layer)!\n",
                 (is_threed==1?3:1));
        return -1;
      }
      
      break;
    }
    
    if (p->kc>atmos->Nz-1) {
      p->kc=atmos->Nz-1;

      /* to avoid roundoff errors */
      if (!sample->spherical && !sample->spherical3D)   /* This does not work in spherical geometry */
        p->x[2] = atmos->Z[p->kc+1]; 
      
      p->photon_status = MCSTATUS_BOUNDARY_UPPER;
      break;
    }
  } /* end while loop */

#if HAVE_LIDAR
  free(qidd);
#endif

  return p->photon_status;
}


/***********************************************************************************/
/* Function: create_escape_photon                                         @62_30i@ */
/* Description:                                                                    */
/* Create photon needed for the escape probability from actual photon.             */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static photon_struct *create_escape_photon ( photon_struct     *photon, 
					     atmosphere_struct *atmos,
					     sample_struct     *sample,
					     int                id )
{
  photon_struct *p = calloc_photon ( sample,
				     atmos->Nz,
				     MCSRC_NONE,
				     atmos->nlambda_abs,
                                     atmos->Nc,
				     atmos->n_caoth );

  /* work on a copy of the photon, not the original one */
  cp_photon_struct (p, photon, sample, atmos->n_caoth);
  
  /* the direction to be traced is the radiance direction */
  if (sample->escape)
    cp_direction (&(p->dir), &(sample->rad[id].dir));
#if HAVE_LIDAR
  /* the direction to be traced is the detector/end of backward beam direction */
  if (sample->LidarLocEst)
    cp_direction (&(p->dir), &(photon->lest.dir));
#endif
  return p;
}


/***********************************************************************************/
/* Function: escape_probability                                           @62_30i@ */
/* Description:                                                                    */
/* Calculate the escape probability to all layers to be sampled.                   */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: ...., last big modification by Robert Buras                             */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int escape_probability ( photon_struct     *p,
				atmosphere_struct *atmos,
				sample_struct     *sample,
				result_struct     *result,
				elevation_struct  *elev,
				int                id,
				int                il,
				int                itl,
				int                surface,
				float             *refind )
{
  double step=0.0, gamma=0.0, slant2horz=0.0;
  int status=0;
  int is_threed=0;

#ifdef HAVE_SPHER
  int first=1;
#endif
  
#if HAVE_ALIS
  double* totweight_spectral;
  double* totweight_concentration;
#endif

#ifdef HAVE_SPHER
  double e_r[3], r=0.0, rlev=0.0;
  double sintheta=0.0, costheta=0.0;
#endif

#if HAVE_VROOM
  int LoWeRR_crit=0; /* Iwabuchi trick */
  double tau=0.0, tauabs=0.0;

  tau = sample->LE_taucrit;
#endif

  if (sample->pan_quicklook && !surface)
    p->pdir[id]=1.0;

  /* for spherical 3D, we should check here whether the earth is in the way */
#ifdef HAVE_SPHER
  if (sample->spherical3D && !sample->refraction && sample->backward==MCBACKWARD_RADIANCE) {
    cart_to_spher ( p->x, p->dir.dx,
		    &r, &costheta, &sintheta, e_r );

    if (costheta < 0) {
      /* escape photon is moving downward, check whether it hits earth */
      rlev = atmos->Z[0] + atmos->r_earth;
      if ( sintheta < rlev/r )
	/* photon will hit earth, do not escape it */
	return 0;
    }
  }
#endif

  /* initialize status; if unchanged at end of travel step, the photon crosses a grid boundary */
  p->photon_status = MCSTATUS_TRAVEL;
  p->RIS_mode = MCRIS_MODE_NORMAL;
  p->VIS_mode = MCVIS_MODE_NORMAL;

  /* ??? need to check this ???      CHECK!!!               */ 
  /* special treatment for photons starting at the surface: */
  if (surface && sample->escape) {
    
    /* count escaping photon at the surface */
    if (!sample->spectral_is && !sample->concentration_is && !sample->boxairmass)
      count_escape (p, sample, atmos, result->surf, id, il,
                    0, NULL, NULL, NULL, NULL, NULL, result->surf_t,
                    p->tauabs.tot, 1.0, NULL, NULL, NULL, result->mish);
    /* ????? special treatment for the layer which coincides with the surface; CHECK!!! */ 
    /* e.g. z=0 if no elevation but for safety reasons also z=elevation ...    */
    /* it may also happen at arbitrary locations in the atmosphere but so      */
    /* rare that it may be neglected                                           */
  } /* endif surface && escape */
  
  /* if vertical component is zero, there will be */
  /* no contribution to any level                 */
  if (p->dir.dx[2]==0 && !sample->spherical && !sample->LidarLocEst)
    return 0;

  while (1==1) {  /* until photon reaches the surface or TOA */
    is_threed = (atmos->threed[MCCAOTH_TOT][p->kc]>=1); /* 3D layer */
    
    /**************************************************************/
    /* 1. find out what happens next with photon, and step length */
    /**************************************************************/

    /* calculate next intersection point with the box boundaries   */
    if (sample->spherical3D) {
#ifdef HAVE_SPHER
      status = intersection3D_spherical (p, atmos,
					 sample->spherical3D_scene,
					 &first, &step);
      if (status)
	return err_out("Error %d returned by intersection3D_spherical\n", status);
#else
      fprintf(stderr,"Error! you are not allowed to use spherical 3D!\n");
      return -1;
#endif
    }
    else {
      if (is_threed==1) { /* 3D layer */
#if HAVE_MYSTIC3D
	status = intersection3D (p, atmos, 10000, 0, &step);
	if (status)
	  return err_out("Error %d returned by intersection3D\n", status);
#else
	fprintf(stderr,"Error! you are not allowed to use mystic 3D!\n");
	return -1;
#endif
      }
      else {
	if (sample->spherical){
	  status = intersection1D_spherical (p, atmos,
					     sample->refraction, refind,  &step);
	  /* Photon is absorbed by the surface.*/
	  if (status<0)
	    break;
	}
	else {
	  intersection1D (p, atmos, 10000, 0, &step, 0);
	}
      }
    }

#ifdef MYSTIC_DEBUG
    if (step<-MC_EPSILON) {
      fprintf (stderr, "\nFATAL error in escape_probability(), step %g < 0 in %dD box (%d, %d, %d).\n",
               step, (is_threed==1?3:1), p->ic, p->jc, p->kc);
      fprintf (stderr,   "Start at (%g, %g, %g), direction (%g, %g, %g)\n", 
               p->x[0], p->x[1], p->x[2], p->dir.dx[0], p->dir.dx[1], p->dir.dx[2]);
      fprintf (stderr,   "hit (%d, %d, %d), cotheta %f\n", 
               p->dir.hit[0], p->dir.hit[1], p->dir.hit[2], p->dir.cotheta);
      
      return -1;
    }
#endif
      
    if (step<0)
      step=0;
    
    /* check if photon crossed the 2D surface */
    if (elev->elev2D)  {
      
      status = cross_surface (p, step,
                              elev,
                              sample->bcond,
                              &gamma);
      
      
      /* if cross_surface ran into some inconsistency */
      if (status<0) {
        fprintf (stderr, "Something terrible has happened! cross_surface() returned\n");
        fprintf (stderr, "status %d\n", status);
        return -1;
      }
      
      if (status>0) {
#ifdef MUCHOUT
	if(p->muchoutcounter==MUCHOUTPHOTON) {
	  fprintf(stderr,"escape photon crossing surface at %e %e %e\n",p->x[0],p->x[1],p->x[2]);
	  //	  return -1;
	}
#endif
        step = gamma;
        p->photon_status = MCSTATUS_SURFACE;
      }
    }
    
#if HAVE_LIDAR
    /* stop if detector / end of backward beam has been reached */
    if (sample->LidarLocEst)
      if (p->pathlength + step >= p->parent_photon->pathlength + p->parent_photon->lest.dist) {
        step = p->parent_photon->pathlength + p->parent_photon->lest.dist - p->pathlength;
        p->photon_status = MCSTATUS_HIT_LIDAR;
      }
#endif

#ifdef CLDPRP
    /*******************************************************************/
    /* 1b. integration for visible cloud property estimation (F.Ewald) */
    /*******************************************************************/

    if (sample->cldprp){
      status = sample_cldprp (p, atmos, step);      
      if (status)
        return err_out ("ERROR! sample_cldprp returned status %d\n",status);
    }
#endif

    /*******************************************************************/
    /* 2. apply step length, move photon, and exit if photon path ends */
    /*******************************************************************/

    /* sum up absorption and scattering optical depths */
    /* step is the escape distance within each box between two box boundaries (or between initital photon position and box boundary) depending on the escape direction */
    p->tauabs.tot += step * ( get_ksca (atmos, p, MCCAOTH_TOT) + get_kabs_tot (atmos, p) );
    p->pathlength += step;
    if (sample->coherent_backscatter && sample->LidarLocEst)
      p->tauext_tot += step * (get_ksca(atmos, p, MCCAOTH_TOT) + get_kabs(atmos,p, MCCAOTH_TOT));
  
    if (sample->spectral_is || sample->concentration_is || sample->boxairmass){
      /*   dtauabs_spectral_calc(atmos, p, step); */
      p->pathlength_per_layer[p->kc] += step;
    }
    
    
#if HAVE_LIDAR
    /* calculate Jacobian weight matrix FIJ */
    if (sample->LLE_jacobian || sample->abs_jacobian)
      p->r_jacobian [p->kc] += step;
#endif

#if HAVE_VROOM
    if (sample->LE_taucrit>=0.0) {
      /* Iwabuchi trick */
      if ( p->tauabs.tot - p->parent_photon->tauabs.tot > tau ) {
	tau += random_tau();
	if ( LoWeRR_crit || p->tauabs.tot - p->parent_photon->tauabs.tot > tau ) {
	  p->photon_status=MCSTATUS_PURGE;
	  return 0;
	}
	LoWeRR_crit=1;
      }
    }
#endif

    /* move photon to box boundary, to next scattering point,  */
    /* or to the 2D surface                                    */
      
    if (is_threed==1 || sample->spherical3D) { /* 3D layer */
#if HAVE_MYSTIC3D
      status = step3D (p, atmos, step, sample->bcond, 0, sample->spherical3D, 0);
#else
      fprintf(stderr,"Error! you are not allowed to use mystic 3D!\n");
      return -1;
#endif
    }
    else
      status = step1D (p, atmos, step, sample->bcond, 0, 0);
  
    if (status<0) {
      fprintf (stderr, "Error %d returned by step%dD()\n", status,
               (is_threed==1?3:1));
      return status;
    }

    /* end photon travel unless it simply crosses the boundary */
    if (p->photon_status != MCSTATUS_TRAVEL )
      break;
      
    /********************************/
    /* 3. count photon if necessary */
    /********************************/
    /* register photon at horizontal boundary */
    if (sample->sample[p->kc + (1-p->dir.hit[2])] && sample->escape) {
      if (p->reallyhit[2] || (is_threed==0) ) {

	if (LoWeRR_crit)
	  tauabs = p->parent_photon->tauabs.tot + sample->LE_taucrit;
	else
	  tauabs = p->tauabs.tot;
        
        if (sample->spectral_is || sample->concentration_is){
#if HAVE_ALIS
          totweight_spectral=calloc(atmos->nlambda_abs, sizeof(double));
          totweight_concentration=calloc(atmos->Nc, sizeof(double));
          
          totweight_spectral[0]=1.0; 
          totweight_concentration[0]=1.0;
          
          if (sample->spectral_is)
            spectral_is_weight(&totweight_spectral, p, p->parent_photon, atmos);
          if (sample->concentration_is)
            concentration_is_weight(&totweight_concentration, p, p->parent_photon, atmos);
          
          count_escape (p, sample, atmos,
                        result->alt[p->kc + (1-p->dir.hit[2])], id, il,
                        sample->backward, result->back, result->back2, result->back_t, 
			result->circcontr,
                        result->jacobian, result->rad_t[p->kc + (1-p->dir.hit[2])],
                        tauabs, 1.0, totweight_spectral, totweight_concentration,
                        result->pathlength_per_layer_tot, result->mish);
        
          free(totweight_spectral);
          free(totweight_concentration);
	  
#endif
	}
        else{
	  if (!p->doling)
            count_escape (p, sample, atmos,
			  result->alt[p->kc + (1-p->dir.hit[2])], id, il,
			  sample->backward, result->back, result->back2, result->back_t, result->circcontr,
			  result->jacobian, NULL,
			  tauabs, 1.0, NULL, NULL, result->pathlength_per_layer_tot, result->mish);
        }
      }
    }

    /**********************************/
    /* 4. check if photon leaves grid */
    /**********************************/

    /* finally check if photon arrived at the lower surface */
    /* or at the top of the atmosphere                      */
    
    if (p->kc<0) {
      p->photon_status = MCSTATUS_BOUNDARY_LOWER;
        
      if (elev->elev2D) {
        fprintf (stderr, "FATAL error! Escape photon reached the 1D boundary in a\n");
        fprintf (stderr, "             2D elevation case (%dD layer)!\n\n",
                 (is_threed==1?3:1));
        return -1;
      }
        
      break;
    }
      
    if (p->kc>atmos->Nz-1) {
      p->photon_status = MCSTATUS_BOUNDARY_UPPER;
      break;
    }
  } 
 /* end step loop */ 

#if HAVE_LIDAR
  /**********************************************************************/
  /* 5. if lidar photon has not yet reached satellite detector, move it */
  /**********************************************************************/
  
  if ( sample->LidarLocEst && p->photon_status == MCSTATUS_BOUNDARY_UPPER ) {
    p->pathlength = p->parent_photon->pathlength + p->parent_photon->lest.dist;
    p->photon_status = MCSTATUS_HIT_LIDAR;
    /* we do not need to change the photon position, it is no longer needed */
  }
#endif

  /******************************************/
  /* 6. count photon at ground and in lidar */
  /******************************************/

  /* register photon at the lower boundary */
  if (sample->escape) {
    if (p->photon_status==MCSTATUS_BOUNDARY_LOWER ||
        p->photon_status==MCSTATUS_SURFACE ) {
    
      /* conversion factor to correct for slant surface */
      slant2horz = slt2hrz (elev, p, sample->surfaceparallel, 1);
      if (slant2horz<0.0)
	return fct_err_out (-1, "slt2hrz", ERROR_POSITION);
    
      if (LoWeRR_crit)
	tauabs = p->parent_photon->tauabs.tot + sample->LE_taucrit;
      else
	tauabs = p->tauabs.tot;
      
      if (!sample->spectral_is && !sample->concentration_is && !sample->boxairmass)
        count_escape (p, sample, atmos, result->surf, id, il,
                      0, NULL, NULL, NULL, NULL, NULL,  result->surf_t,
                      tauabs, slant2horz, NULL, NULL, result->pathlength_per_layer_tot, result->mish);
    }
  }

#if HAVE_LIDAR
  /* register local estimate */
  if (itl!=-1) {
    status = mc_lidar_locest (result, sample, atmos, p, itl, LoWeRR_crit);
    if (status!=0)
      return err_out ("Error %d returned by mc_lidar_locest()\n", status);
  }
#endif

  return 0;
}


/***********************************************************************************/
/* Function: count_escape                                                 @62_30i@ */
/* Description:                                                                    */
/* Count the escape probability at the actual layer.                               */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: put into function by Robert Buras                                       */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static void count_escape (photon_struct *p, sample_struct *sample,
			  atmosphere_struct *atmos,
			  radiation_field *result_alt, int id, int il,
			  int backward, double ***result_back, double ***result_back2,
                          double *****result_back_t, 
			  double **circcontr,
			  jacobian_result_field *jacobian,
                          radiation_field_t *result_rad_t,
			  double tauabs, double slant2horz, 
                          double *totweight_spectral,
                          double *totweight_concentration,
                          double *pathlength_per_layer_tot,
                          mishchenko_cb_struct* mish
                          )
{
  int is=0, js=0, ip=0, ir=0, it=0, iv=0, ic=0;
  int isp=0, kc=0;
  double jac_alpha=0.0, jac=0.0;

  double totweight=0.0, add=0.0;
  
  /* count passing photon */
  sample_coord (p, sample, &is, &js);

  totweight = exp(-tauabs) * p->pdir[id] * p->weight / p->dir.cotheta;

  if ( sample->coherent_backscatter && !sample->LidarLocEst )
    save_phase_matrix_elements ( mish, totweight, p->scattercounter, p->stokes0 );
 
#ifdef MUCHOUT
  if (p->muchoutcounter==MUCHOUTPHOTON)
    fprintf(stderr,"counters %d %d %d %d --- escape: add %e pdir %e pweight %e abs %e stokes %e cos %e\n",
	    p->scattercounter,p->escapescattercounter,p->clonescattercounter,p->SC_mode,
	    totweight*p->stokes[0],
	    p->pdir[id],p->weight,exp(-tauabs),p->stokes[0],p->dir.cotheta);
#endif

#ifdef CLDPRP
  
  if (sample->cldprp) {
    
    if (p->rayleighcounter == 0 || (p->cldprp.rhit_wc+p->cldprp.rhit_ic) > 0.0){
      
      result_alt->cldprp.reff_wc[sample->backward_is][sample->backward_js] += totweight * p->stokes[0] * p->cldprp.reff_wc;
      result_alt->cldprp.reff_ic[sample->backward_is][sample->backward_js] += totweight * p->stokes[0] * p->cldprp.reff_ic;
      result_alt->cldprp.rhit_wc[sample->backward_is][sample->backward_js] += totweight * p->stokes[0] * p->cldprp.rhit_wc;
      result_alt->cldprp.rhit_ic[sample->backward_is][sample->backward_js] += totweight * p->stokes[0] * p->cldprp.rhit_ic;
      
      result_alt->cldprp.tau_wc[sample->backward_is][sample->backward_js]  += totweight * p->stokes[0] * p->cldprp.tau_wc;
      result_alt->cldprp.tau_ic[sample->backward_is][sample->backward_js]  += totweight * p->stokes[0] * p->cldprp.tau_ic;
      
      result_alt->cldprp.dxlwc[sample->backward_is][sample->backward_js]   += totweight * p->stokes[0] * p->cldprp.dxlwc;
      result_alt->cldprp.dylwc[sample->backward_is][sample->backward_js]   += totweight * p->stokes[0] * p->cldprp.dylwc;
      result_alt->cldprp.dzlwc[sample->backward_is][sample->backward_js]   += totweight * p->stokes[0] * p->cldprp.dzlwc;
      result_alt->cldprp.dxiwc[sample->backward_is][sample->backward_js]   += totweight * p->stokes[0] * p->cldprp.dxiwc;
      result_alt->cldprp.dyiwc[sample->backward_is][sample->backward_js]   += totweight * p->stokes[0] * p->cldprp.dyiwc;
      result_alt->cldprp.dziwc[sample->backward_is][sample->backward_js]   += totweight * p->stokes[0] * p->cldprp.dziwc;
      
      result_alt->cldprp.totweights[sample->backward_is][sample->backward_js]   += totweight * p->stokes[0];
      if (p->cldprp.rhit_wc > 0.0) 
        result_alt->cldprp.wc_weights[sample->backward_is][sample->backward_js] += totweight * p->stokes[0];
      if (p->cldprp.rhit_ic > 0.0) 
        result_alt->cldprp.ic_weights[sample->backward_is][sample->backward_js] += totweight * p->stokes[0];
    }
  }
  
#endif
  
  
  for (ip=0; ip<sample->nstokes; ip++){
    
    add = totweight * p->stokes[ip];
    result_alt->radesc[id][is][js][ip] += add; /* RADESC */
    
    if(sample->boxairmass)
      for (kc=0; kc<atmos->Nz; kc++){ 
        pathlength_per_layer_tot[kc]+= add * p->pathlength_per_layer[kc];
      }    
    
    if (sample->spectral_is || sample->concentration_is){
      for (ic=0; ic<atmos->Nc; ic++){
        for (iv=0; iv<atmos->nlambda_abs; iv++){
          result_rad_t->rad_t[ic][is][js][ip][iv] += add * totweight_spectral[iv]
            * totweight_concentration[ic];
        }
      }
    }

    if (sample->std)
      p->parent_photon->wtree = addtree_stddev 
        (p->parent_photon->wtree, &(result_alt->radesc2[id][is][js][ip]), add);
    
    if (backward) {
      result_back[sample->backward_is][sample->backward_js][ip] += add;
      
      if (sample->spectral_is || sample->concentration_is)
        for (ic=0; ic<atmos->Nc; ic++)
          for (iv=0; iv<atmos->nlambda_abs; iv++)
            result_back_t[ic][sample->backward_is][sample->backward_js][ip][iv] += 
              add * totweight_spectral[iv]*totweight_concentration[ic];

      if(sample->ncirc) {
	for (ic=0;ic<p->parent_photon->ncircos;ic++)
	  *(p->parent_photon->tocirco[ic]) += add;
      }
      
      if (sample->std)
        p->parent_photon->wtree = addtree_stddev 
          (p->parent_photon->wtree, &(result_back2[sample->backward_is][sample->backward_js][ip]), add);
    }
    
    if (sample->Nr>0 && sample->Nt>0)
      if (sample_radius (p, sample, &ir))
        if (sample_time (p, sample, &it)) {
          result_alt->radpes[id][ir][it][ip] += add; /* radial/pathlength */
          
          if (sample->std)
            p->parent_photon->wtree = addtree_stddev 
              (p->parent_photon->wtree, &(result_alt->radpes2[id][ir][it][ip]), add);
        }
  }
  

  /* jacobian for non lidar; RESINC */
  /* this is done for molecular only! We are interested in changes of
     signal by varying NO2 concentration, i.e. varying molecular
     absorption, and in changes of signal by varying wavelength
     slightly. Latter mainly implies changes in rayleigh scattering
     and molecular absorption, and only in second instance changes in
     mie scattering. However, in order to calculate changes by changes
     in mie scattering, we would have to take into account not only
     changes in scattering and absorption coefficients, but also in
     phase function, which is not possible. Also, 3D jacobian not
     possible. */

  if (backward && sample->abs_jacobian) {
    for (isp=1; isp<=atmos->n_caoth; isp++) {    
      for (kc=0; kc<atmos->Nz; kc++) {

	/* scattering */
	jac_alpha = atmos->ksca [MCSC_MODE_NORMAL][MCRIS_MODE_NORMAL]->prof [isp][kc];

	/* jacobian for beta_sct */
	if (jac_alpha!=0.0)
	  jac = totweight * p->stokes[0] * 
	    ( - p->r_jacobian [kc] * jac_alpha
	      + p->parent_photon->q_jacobian [0][isp-1][kc]
	      + ((double) (kc == p->parent_photon->kc))
	      * p->parent_photon->lest.pdir_sct[isp][0] / p->pdir[id] );
	else
	  jac = 0.0;

	if (jac != 0.0) {
	  jacobian->jacobian_t [sample->backward_is][sample->backward_js][isp][MC_JAC_VAR_BETA][kc] += jac;

	  if (sample->std)
	    jacobian->jacobian_tact[isp][MC_JAC_VAR_BETA][kc] += jac;
	}

	/* absorption */
	jac_alpha = atmos->kabs->prof [isp][kc];

	/* jacobian for beta_abs */
	if (jac_alpha!=0.0)
	  jac = totweight * p->stokes[0] * ( - p->r_jacobian [kc] * jac_alpha );
	else
	  jac = 0.0;

	if (jac != 0.0) {
	  jacobian->jacobian_t [sample->backward_is][sample->backward_js][isp][MC_JAC_VAR_ABS][kc] += jac;

	  if (sample->std)
	    jacobian->jacobian_tact[isp][MC_JAC_VAR_ABS][kc] += jac;
	}
      } /* endfor kc */
    } /* endfor isp */
  } /* endif jacobian */

}

/***********************************************************************************/
/* Function: rpv_brdf                                                     @62_30i@ */
/* Description:  BRDF calculatation according to Rahman, Verstraete,               */
/*               and Pinty [1993]; attention: a factor 1/PI is added which         */
/*               translates from "biconical reflectance factor" to BRDF; this      */
/*               is different from Frank Evans' rpv_reflection.f which was used as */
/*               a basis for rpv_brdf(); also different from Frank Evans' code     */
/*               is the hot spot occurring at phi=180 degree: Reflection into the  */
/*               hot spot direction implies a 180 degree azimuth change!           */
/*               Tested behaviour (March 1, 2004): radiances correct for both      */
/*               cone sampling and escape radiance.                                */
/*               Lately implemented snow stuff from Meerkötter et al. (RPB)        */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: last change by RPB                                                      */
/*                                                                        @i62_30@ */
/***********************************************************************************/

double rpv_brdf (double rho0, double k, double theta, 
		 double scale, double sigma, double t1, double t2,
                 double mu1, double mu2, double phi)
{

  /* Local variables */
  double capg=0, cosg=0, f=0, h=0, m=0, t=0;
  double alpha=0;
  double result=0;
  
  phi = 180.0 - phi;


  /* hot spot */
  if (phi==0 && mu1==mu2)   /* special case (might cause numerical troubles) */
    return rho0 * ( pow(2.0*mu1*mu1*mu1,k-1.0) *
		    (1.0-theta)/(1.0+theta)/(1.0+theta) * (2.0-rho0)
		    + sigma/mu1 ) * ( t1 * exp(PI * t2) + 1.0 ) * scale;

  m = pow(mu1*mu2*(mu1+mu2), k-1.0);
  
  if (mu1<1.0 && mu2<1.0)
    alpha = sqrt((1.0 - mu1*mu1) * (1.0 - mu2*mu2)) * cosd(phi);
  else 
    alpha = 0;
  
  cosg = mu1*mu2 + alpha;
  
  f = (1.0 - theta * theta) / 
    pow(1.0 + 2.0*theta*cosg + theta*theta, 1.5);
  
  capg = sqrt(1.0/mu1/mu1 + 1.0/mu2/mu2 - 2.0 - 2.0*alpha/mu1/mu2);
  
  h = 1.0 + (1.0-rho0) / (1.0+capg) ;

  t = 1.0 + t1 * exp ( t2 * ( PI - acos(cosg) ) );
  
  result = rho0 * ( m * f * h + sigma / mu1 ) * t * scale;

  if (result < 0)
    result = 0;

  return result;
}


/***********************************************************************************/
/* Function: equal                                                        @62_30i@ */
/* Description:                                                                    */
/*  Two numbers equal?                                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int equal (double x, double y)
{
  if ((float) x != (float) y)
    return 0;

  return 1;
}


/***********************************************************************************/
/* Function: sample_coord                                                 @62_30i@ */
/* Description:                                                                    */
/*  Determine sample coordinates of a photon.                                      */          
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static inline void sample_coord (photon_struct *p, sample_struct *sample, int *is, int *js)
{
  /* dummy for spherical3D and panorama; actually, this function should never be called! Clean up! */
  if (sample->spherical3D || sample->panorama) {
    *is=0;
    *js=0;
    return;
  }

  /* determine sample coordinates */
  *is = (int) (p->x[0] / sample->delX);
  *js = (int) (p->x[1] / sample->delY);
  
  /* check boundaries */
  if (*is<0) {
    fprintf (stderr, "OOPS, set sample coordinate is from %d to %d (you may still trust your result?)\n", *is, 0);
    fprintf (stderr, "      x = %f, y = %f\n", p->x[0], p->x[1]);
    *is=0;
  }
  
  if (*is>=sample->Nx) {
    fprintf (stderr, "OOPS, set sample coordinate is from %d to %d (you may still trust your result?)\n", *is, sample->Nx-1);
    fprintf (stderr, "      x = %f, y = %f\n", p->x[0], p->x[1]);
    *is=sample->Nx-1;
  }
  
  if (*js<0) {
    fprintf (stderr, "OOPS, set sample coordinate js from %d to %d (you may still trust your result?)\n", *js, 0);
    fprintf (stderr, "      x = %f, y = %f\n", p->x[0], p->x[1]);
    *js=0;
  }
  
  if (*js>=sample->Ny) {
    fprintf (stderr, "OOPS, set sample coordinate js from %d to %d (you may still trust your result?)\n", *js, sample->Ny-1);
    fprintf (stderr, "      x = %f, y = %f\n", p->x[0], p->x[1]);
    *js=sample->Ny-1;
  }
}


/***********************************************************************************/
/* Function: sample_radius                                                @62_30i@ */
/* Description:                                                                    */
/*  Determine radial coordinate of a photon.                                       */          
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int sample_radius (photon_struct *p, sample_struct *sample, int *ir)
{
  /* determine sample radial coordinate */
  double dx = p->x[0] - 0.5 * (double) (sample->Nx) * sample->delX;
  double dy = p->x[1] - 0.5 * (double) (sample->Ny) * sample->delY;

  *ir = (int) ( sqrt (dx*dx + dy*dy) / sample->dr );
  
  if (*ir<0 || *ir >= sample->Nr)
    return 0;
  else
    return 1;
}


/***********************************************************************************/
/* Function: sample_time                                                  @62_30i@ */
/* Description:                                                                    */
/*  Determine time/pathlength coordinate of a photon.                              */          
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int sample_time (photon_struct *p, sample_struct *sample, int *it)
{
  
  /* determine sample time coordinates */
  double t = p->pathlength / 3.0E8;  /* speed of light is assumed to be 3E8 m/s */

  *it = (int) ( t / sample->dt );

  if (*it<0 || *it >= sample->Nt)
    return 0;
  else
    return 1;
}


/***********************************************************************************/
/* Function: atmos_coord                                                  @62_30i@ */
/* Description:                                                                    */
/*  Determine atmospheric coordinates of a photon.                                 */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void atmos_coord (photon_struct *p, atmosphere_struct *atmos, int *ic, int *jc)
{
  /* determine atmospheric coordinates */
  *ic = (int) (p->x[0] / atmos->delX);
  *jc = (int) (p->x[1] / atmos->delY);
  
  /* check boundaries */
  if (*ic<0) {
    fprintf (stderr, "OOPS, set atmospheric coordinate ic from %d to %d (you may still trust your result?)\n", *ic, 0);
    fprintf (stderr, "      x = %f, y = %f photon %d \n", p->x[0], p->x[1], p->photoncounter);
    *ic=0;
  }
  
  if (*ic>=atmos->Nx) {
    if ( (double) *ic == p->x[0] / atmos->delX )
      /* photon is at periodic boundary ! */
      (*ic)--;
    else {
      fprintf (stderr, "OOPS, set atmospheric coordinate ic from %d to %d (you may still trust your result?)\n", *ic, atmos->Nx-1);
      fprintf (stderr, "      x = %f, y = %f photon %d \n", p->x[0], p->x[1], p->photoncounter);
      *ic=atmos->Nx-1;
    }
  }
  
  if (*jc<0) {
    fprintf (stderr, "OOPS, set atmospheric coordinate jc from %d to %d (you may still trust your result?)\n", *jc, 0);
    fprintf (stderr, "      x = %f, y = %f photon %d \n", p->x[0], p->x[1], p->photoncounter);
    *jc=0;
  }
  
  if (*jc>=atmos->Ny) {
    if ( (double) *jc == p->x[1] / atmos->delY )
      /* photon is at periodic boundary ! */
      (*jc)--;
    else {
      fprintf (stderr, "OOPS, set atmospheric coordinate jc from %d to %d (you may still trust your result?)\n", *jc, atmos->Ny-1);
      fprintf (stderr, "      x = %f, y = %f\n", p->x[0], p->x[1]);
      *jc=atmos->Ny-1;
    }
  }
}


/***********************************************************************************/
/* Function: elev_coord                                                   @62_30i@ */
/* Description:                                                                    */
/*  Determine elevation coordinates of a photon.                                   */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void elev_coord (photon_struct *p, elevation_struct *elev, int *ie, int *je)
{
  /* determine elevation coordinates */
  *ie = (int) (p->x[0] / elev->delX);
  *je = (int) (p->x[1] / elev->delY);
  
  /* if we are at the right / upper boundary, we want to get the     */
  /* index of the last pixel, rather than one more; this is required */
  /* because elev_coord() is called while the photon steps from      */
  /* box boundary where it eventually hits the right or upper        */
  /* boundaries as well                                              */

  if (*ie==elev->Nx)
    if (p->x[0] == elev->delX * (double) elev->Nx)
      (*ie)--;    

  if (*je==elev->Ny)
    if (p->x[1] == elev->delY * (double) elev->Ny)
      (*je)--;    


  /* check boundaries */
  if (*ie<0) {
    fprintf (stderr, "OOPS, set elevation coordinate ie from %d to %d (you may still trust your result?)\n", *ie, 0);
    fprintf (stderr, "      x = %f, y = %f\n", p->x[0], p->x[1]);
    *ie=0;
  }
  
  if (*ie>=elev->Nx) {
    if ( (double) *ie == p->x[0] / elev->delX )
      /* photon is at periodic boundary ! */
      (*ie)--;
    else {
      fprintf (stderr, "OOPS, set elevation coordinate ie from %d to %d (you may still trust your result?)\n", *ie, elev->Nx-1);
      fprintf (stderr, "      x = %f, y = %f\n", p->x[0], p->x[1]);
      *ie=elev->Nx-1;
    }
  }
  
  if (*je<0) {
    fprintf (stderr, "OOPS, set elevation coordinate je from %d to %d (you may still trust your result?)\n", *je, 0);
    fprintf (stderr, "      x = %f, y = %f\n", p->x[0], p->x[1]);
    *je=0;
  }
  
  if (*je>=elev->Ny) {
    if ( (double) *je < p->x[1] / elev->delY * ( 1 + MC_EPSILON ) )
      /* photon is at periodic boundary ! */
      (*je)--;
    else {
      fprintf (stderr, "OOPS, set elevation coordinate je from %d to %d (you may still trust your result?)\n", *je, elev->Ny-1);
      fprintf (stderr, "      x = %f, y = %f, doubleje %e\n", p->x[0], p->x[1],p->x[1]/elev->delY-(double)*je);
      *je=elev->Ny-1;
    }
  }
}


/***********************************************************************************/
/* Function: calc_normal_vector                                           @62_30i@ */
/* Description:                                                                    */
/*  Determine normal vector on surface.                                            */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: RPB                                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static inline void calc_normal_vector ( photon_struct     *p,
					elevation_struct  *elev,
					sample_struct     *sample,
					atmosphere_struct *atmos,
					int                use_elev,
					int                use_sensordirection,
					double            *norm )
{
  double temp=0.0;

  /* default is horizontal surface */
  norm[0] = 0;
  norm[1] = 0;
  norm[2] = 1;

  if (elev->elev2D && use_elev)
    elev_coord_normal_vector (p, elev, norm);
  else if (use_sensordirection) {     /* user-defined surface normal */
    norm[0] = sample->sensordirection_dx[0];
    norm[1] = sample->sensordirection_dx[1];
    norm[2] = sample->sensordirection_dx[2];
  }
  else if ( sample->spherical && !(p->x[0]-atmos->xmax/2.==0. && p->x[1]-atmos->ymax/2.==0.) ) {
    norm[0] = p->x[0] - atmos->xmax/2.;
    norm[1] = p->x[1] - atmos->ymax/2.;
    norm[2] = p->x[2] + atmos->r_earth;
    /* normalize normal vector */
    temp = sqrt ( norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2] ) ;
    norm[0] /= temp;
    norm[1] /= temp;
    norm[2] /= temp;
  }
#ifdef HAVE_SPHER
  else if ( sample->spherical3D && !(p->x[0]==0.0 && p->x[1]==0.0) )
    cart_to_spher ( p->x, NULL, &temp, NULL, NULL, norm );
#endif

  return;
}

/***********************************************************************************/
/* Function: elev_coord_normal_vector                                     @62_30i@ */
/* Description:                                                                    */
/*  Determine normal vector on elevation coordinates of a photon.                  */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: RPB                                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static inline void elev_coord_normal_vector (photon_struct *p, elevation_struct *elev, double *dx)
{
  int ie=0, je=0, i=0;
  double temp=0.0;

  elev_coord (p, elev, &ie, &je);
  /* upward normal to the slant surface */
  dx[0] = -elev->surf[ie][je].a - elev->surf[ie][je].c * (p->x[1] - (double) je * elev->delY);
  dx[1] = -elev->surf[ie][je].b - elev->surf[ie][je].c * (p->x[0] - (double) ie * elev->delX);
  dx[2] = 1;
                
  /* normalize dx */
  temp = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
  for (i=0; i<3; i++)
    dx[i] /= temp;
}


/***********************************************************************************/
/* Function: albedo_coord                                                 @62_30i@ */
/* Description:                                                                    */
/*  Determine albedo coordinates of a photon.                                      */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static inline void albedo_coord (photon_struct *p, albedo_struct *albedo, int *ia, int *ja)
{
  /* albedo coordinates */
  *ia = (int) (p->x[0] / albedo->delX);
  *ja = (int) (p->x[1] / albedo->delY);

  /* check boundaries */
  if (*ia<0) {
    fprintf (stderr, "OOPS, set albedo coordinate ia from %d to %d (you may still trust your result?)\n", 
             *ia, 0);
    fprintf (stderr, "      x = %f, y = %f\n", p->x[0], p->x[1]);
    *ia=0;
  }
  
  if (*ia>=albedo->Nx) {
    fprintf (stderr, "OOPS, set albedo coordinate ia from %d to %d (you may still trust your result?)\n", *ia, albedo->Nx-1);
    fprintf (stderr, "      x = %f, y = %f\n", p->x[0], p->x[1]);
    *ia=albedo->Nx-1;
  }
  
  if (*ja<0) {
    fprintf (stderr, "OOPS, set albedo coordinate ja from %d to %d (you may still trust your result?)\n", *ja, 0);
    fprintf (stderr, "      x = %f, y = %f\n", p->x[0], p->x[1]);
    *ja=0;
  }
  
  if (*ja>=albedo->Ny) {
    fprintf (stderr, "OOPS, set albedo coordinate ja from %d to %d (you may still trust your result?)\n", *ja, albedo->Ny-1);
    fprintf (stderr, "      x = %f, y = %f\n", p->x[0], p->x[1]);
    *ja=albedo->Ny-1;
  }
}


/***********************************************************************************/
/* Function: surftemp_coord                                               @62_30i@ */
/* Description:                                                                    */
/*  Determine surface temperature coordinates of a photon.                         */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static inline void surftemp_coord (photon_struct *p, surftemp_struct *surftemp, int *it, int *jt)
{
  /* surface temperature coordinates */
  *it = (int) (p->x[0] / surftemp->delX);
  *jt = (int) (p->x[1] / surftemp->delY);
  
  /* check boundaries */
  if (*it<0) {
    fprintf (stderr, "OOPS, set surface temperature coordinate ia from %d to %d (you may still trust your result?)\n", *it, 0);
    fprintf (stderr, "      x = %f, y = %f\n", p->x[0], p->x[1]);
    *it=0;
  }
  
  if (*it>=surftemp->Nx) {
    fprintf (stderr, "OOPS, set surface temperature coordinate ia from %d to %d (you may still trust your result?)\n", *it, surftemp->Nx-1);
    fprintf (stderr, "      x = %f, y = %f\n", p->x[0], p->x[1]);
    *it=surftemp->Nx-1;
  }
  
  if (*jt<0) {
    fprintf (stderr, "OOPS, set surface temperature coordinate ja from %d to %d (you may still trust your result?)\n", *jt, 0);
    fprintf (stderr, "      x = %f, y = %f\n", p->x[0], p->x[1]);
    *jt=0;
  }
  
  if (*jt>=surftemp->Ny) {
    fprintf (stderr, "OOPS, set surface temperature coordinate ja from %d to %d (you may still trust your result?)\n", *jt, surftemp->Ny-1);
    fprintf (stderr, "      x = %f, y = %f\n", p->x[0], p->x[1]);
    *jt=surftemp->Ny-1;
  }
}


/***********************************************************************************/
/* Function: calculate_emission                                           @62_30i@ */
/* Description:                                                                    */
/*                                                                                 */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int calculate_emission (atmosphere_struct *atmos, result_struct *result)
{
  int ic=0, jc=0, kc=0;
  
  for (kc=0; kc<atmos->Nz; kc++)
    for (ic=0; ic<atmos->Nx; ic++)
      for (jc=0; jc<atmos->Ny; jc++) 
        if (atmos->threed[MCCAOTH_TOT][kc]>=1){
          
          /* Average value of Planck function */
          /* emission is negative!            */
          result->absorption3D->tot[kc][ic][jc] = - 0.5 * 
            ( atmos->Bplanck->prof [MCCAOTH_TOT][kc][ic][jc] + atmos->Bplanck->prof [MCCAOTH_TOT][kc+1][ic][jc] ) *
            atmos->kabs3D->prof [MCCAOTH_TOT][kc][ic][jc] * ( atmos->Z[kc+1] - atmos->Z[kc] );


          /* Linear interpolation of temperature */
          /* ?????? CE  ???????    CHECK!!!
           For heating rate calculations the emission term calculated here 
           seems not to be consistent with the emission calculated in 
           thermal backward. 
          */
          /* Interpolation of temperature is also not really correct ...*/
          /*Tavg = 0.5* (atmos->T[kc] + atmos->T[kc+1]);
          F77_FUNC (cplkavg, CPLKAVG) (wvnmlo, wvnmhi, &Tavg, &plkavg);
          result->absorption3D->tot[kc][ic][jc] = plkavg*
            atmos->kabs3D->prof [MCCAOTH_TOT][kc][ic][jc] * (atmos->Z[kc+1] - atmos->Z[kc]);
          */
          
        }
  
  return 0;  /* if o.k. */
}


/***********************************************************************************/
/* Function: calculate_backward_emission                                  @62_30i@ */
/* Description:                                                                    */
/*                                                                                 */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int calculate_backward_emission (atmosphere_struct *atmos, sample_struct *sample, int kc, double *emis)
{  
  int ic=0, jc=0, ic1=0, ic2=0, jc1=0, jc2=0;
  double xs1=0, xs2=0, ys1=0, ys2=0;  /* sample pixel coordinates                */
  double xc1=0, xc2=0, yc1=0, yc2=0;  /* atmospheric (=caothy) pixel coordinates */
  double xl=0, xr=0, yl=0, yu=0;
  double sum=0, fraction=0;

  *emis=0;

  if (atmos->threed[MCCAOTH_TOT][kc]>=1){

    /* first determine index of edges of the sample pixel in caoth coordinates; */
    /* this should in principle be done with atmos_coord, but the latter        */
    /* requires a photon struct as input while we provide simply x/y values     */

    /* sample pixel edges */
    xs1 = (double) (sample->backward_is)   * sample->delX;
    xs2 = (double) (sample->backward_is+1) * sample->delX;
    ys1 = (double) (sample->backward_js)   * sample->delY;
    ys2 = (double) (sample->backward_js+1) * sample->delY;

    /* caoth pixels containing the sample pixel edges */
    ic1 = (int) ( xs1 / atmos->delX);
    ic2 = (int) ( xs2 / atmos->delX);
    jc1 = (int) ( ys1 / atmos->delY);
    jc2 = (int) ( ys2 / atmos->delY);

    /* check boundaries */
    if (ic1<0) ic1=0;
    if (ic2<0) ic2=0;
    if (jc1<0) jc1=0;
    if (jc2<0) jc2=0;

    if (ic1 >= atmos->Nx) ic1 = atmos->Nx-1;
    if (ic2 >= atmos->Nx) ic2 = atmos->Nx-1;
    if (jc1 >= atmos->Ny) jc1 = atmos->Ny-1;
    if (jc2 >= atmos->Ny) jc2 = atmos->Ny-1;


    sum=0;
    for (ic=ic1;ic<=ic2;ic++)
      for (jc=jc1;jc<=jc2;jc++) {
	xc1 = (double) ( ic    * atmos->delX);
	xc2 = xc1 + atmos->delX;
	yc1 = (double) ( jc    * atmos->delY);
	yc2 = yc1 + atmos->delY;

	xl=(xc1>xs1?xc1:xs1);
	xr=(xc2<xs2?xc2:xs2);

	yl=(yc1>ys1?yc1:ys1);
	yu=(yc2<ys2?yc2:ys2);

	fraction = (xr-xl)*(yu-yl)/(sample->delX * sample->delY);
	sum+=fraction;

	*emis = *emis - fraction * 0.5 * 
	  ( atmos->Bplanck->prof [MCCAOTH_TOT][kc][ic][jc] + atmos->Bplanck->prof [MCCAOTH_TOT][kc+1][ic][jc] ) *
	  atmos->kabs3D->prof [MCCAOTH_TOT][kc][ic][jc] * (atmos->Z[kc+1] - atmos->Z[kc]);
      }

    /* fprintf (stderr, "ic = %d - %d, jc = %d - %d, sum=%.10f\n", ic1, ic2, jc1, jc2, sum); */
    /* check if the sum of all fractions is 1 */
    if (fabs(sum-1.0)>MC_EPSILON) {
      fprintf (stderr, "Error in calculate_backward_emission, sum = %.14f\n", sum);
      return -1;
    }

  }
  else {
    /* Average value of Planck function */
    /* emission is negative!            */
    *emis = - 0.5 * 
      ( atmos->Bplanck->prof [MCCAOTH_TOT][kc][0][0] + atmos->Bplanck->prof [MCCAOTH_TOT][kc+1][0][0] ) *
      atmos->kabs->prof [MCCAOTH_TOT][kc] * ( atmos->Z[kc+1] - atmos->Z[kc] );
  }

  return 0;
}



/***********************************************************************************/
/* Function: HGint                                                        @62_30i@ */
/* Description:                                                                    */
/*  Calculate the integrated HG phase function.                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static double HGint (double g, double mu)
{
  return (1.0-g*g)/2.0/g/sqrt(1.0+g*g-2.0*g*mu) - (1.0-g)/2.0/g;
}


/***********************************************************************************/
/* Function: HGint                                                        @62_30i@ */
/* Description:                                                                    */
/*  Calculate the integrated double HG phase function.                             */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static double HG2int (double g1, double g2, double ff, double mu)
{
  return ff * HGint(g1,mu) + (1.0-ff) * HGint(g2,mu);
}
 

/***********************************************************************************/
/* Function: sc_HG2_theta_old                                             @62_30i@ */
/* Description:                                                                    */
/*  Calculate a random polar angle for the double Henyey-Greenstein                */
/*  phase function.                                                                */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

double sc_HG2_theta_old (double g1, double g2, double ff)
{
  double theta = 0.0;
  double P = uvspec_random();

  double F=0, FF=0, mu=0;

  double accur=0.001;


  if (ff==1.0) {
    theta = acos(sc_HG_mu (g1));
  }
  else {

    /* ATTENTION: No error check, but very rude correction */
    
    if (g1<-1.0) {
      fprintf (stderr, " ... g1 = %g out of range, setting to -1\n", g1);
      g1=-1;
    }
    if (g1>1.0) {
      fprintf (stderr, " ... g1 = %g out of range, setting to 1\n", g1);
      g1=1;
    }

    
    if (g2<-1.0) {
      fprintf (stderr, " ... g2 = %g out of range, setting to -1\n", g2);
      g2=-1;
    }
    if (g2>1.0) {
      fprintf (stderr, " ... g2 = %g out of range, setting to 1\n", g2);
      g2=1;
    }
    
    
    if (ff<-1.0) {
      fprintf (stderr, " ... ff = %g out of range, setting to -1\n", ff);
      ff=-1;
    }
    if (ff>1.0) {
      fprintf (stderr, " ... ff = %g out of range, setting to 1\n", ff);
      ff=1;
    }
    
    
    /* start value slightly different from 0, otherwise */
    /* we get a local maximum at 0.                     */
    mu =- accur + 2.0*accur*uvspec_random();
    
    F  = HG2int (g1, g2, ff, mu);
    FF = HG2    (g1, g2, ff, mu); 
    
    /* Newton iteration */
    while (fabs(P-F)>accur) {
      F  = HG2int (g1, g2, ff, mu);
      FF = HG2    (g1, g2, ff, mu); 
      
      mu += (P-F)/ FF;

      if (mu<-1)
        mu=-1;
      if (mu>1)
        mu=1;
    }

    theta = acos(mu);
  }  


  return theta;
}


/***********************************************************************************/
/* Function: sc_HG2_mu                                                    @62_30i@ */
/* Description:                                                                    */
/*   Calculate a random polar angle for the double Henyey-Greenstein               */
/*   phase function.                                                               */
/*   Returns mu=cos(theta), for more speed!                                        */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: ...., last big modification by Robert Buras                             */
/*                                                                        @i62_30@ */
/***********************************************************************************/

double sc_HG2_mu (double g1, double g2, double ff)
{
  double mu = 0.0;
  double P = uvspec_random();
  
  /* ATTENTION: No error check, but very rude correction */
  if (ff<-1.0) {
    fprintf (stderr, " ... ff = %g out of range, setting to -1\n", ff);
    ff=-1;
  }

  if (ff>1.0) {
    fprintf (stderr, " ... ff = %g out of range, setting to 1\n", ff);
    ff=1;
  }

  if (P<=ff) /* forward,  with a probability of ff   */
    mu = sc_HG_mu (g1);
  else       /* backward, with a probability of 1-ff */
    mu = sc_HG_mu (g2);

  return mu;
}


/***********************************************************************************/
/* Function: init_direction                                                        */
/* Description:                                                                    */
/*  Initialize direction with solar zenith sza and azimuth phi.                    */
/*  Convention:    phi=0,   South                                                  */
/*                 phi=90,  West                                                   */
/*                 phi=180, North                                                  */
/*                 phi=270, East                                                   */
/*  These definition refers to the position of the sun; e.g. phi=0 means that      */
/*  the sun is in the South and the photons are headed North.                      */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                                 */
/***********************************************************************************/

void init_direction (double sinsza, double cossza, double sinphi, double cosphi, 
                     direction *dir)
{
  dir->dx[0]=  sinsza*sinphi;
  dir->dx[1]=  sinsza*cosphi;;
  dir->dx[2]= -cossza;

  /* cosine of solar zenith angle for radiance calculation */
  dir->cotheta = fabs(cossza);

  hitflag (dir);
}


/***********************************************************************************/
/* Function: hitflag                                                               */
/* Description:                                                                    */
/*  Calculate 'hitflag' which is used to determine the general direction of        */
/*  a photon:                                                                      */
/*   hit[i] =  1   positive                                                        */
/*   hit[i] =  0   negative                                                        */
/*   hit[i] = -1   no component                                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                                 */
/***********************************************************************************/

void hitflag (direction *dir)
{
  int i=0;

  for (i=0; i<3; i++) {
    dir->hit[i] = -1;
    
    if (dir->dx[i]>0)
      dir->hit[i]=1;
    else {
      if (dir->dx[i]<0)
        dir->hit[i] = 0;
    }
  }
}


/***********************************************************************************/
/* Function: addtree_stddev                                                        */
/* Description:                                                                    */
/*  Add a new photon weight to the weights tree.                                   */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                                 */
/***********************************************************************************/

struct tnode *addtree_stddev (struct tnode *p, double *loc, double weight)
{ 
  if (p == NULL){ /* a new word has arrived */
    p = (struct tnode *) calloc(1, sizeof(struct tnode));   /* create a new node */
    p->loc = loc;
    p->weight = weight;
    p->left = p->right = NULL;
  }
  else {
    if (loc == p->loc)
      p->weight+=weight;   /* add photon weight */
    else {
      if (loc < p->loc)    /* smaller goes into left subtree */
        p->left = addtree_stddev (p->left, loc, weight);
      else                 /* larger goes into right subtree */
        p->right = addtree_stddev (p->right, loc, weight);
    }
  }
  return(p);
}


/***********************************************************************************/
/* Function: treeprint_stddev                                                      */
/* Description:                                                                    */
/*  Store the squared photon weights at their destinations.                        */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                                 */
/***********************************************************************************/

void treeprint_stddev (struct tnode *p)
{
  /* ??? check if freeing is ok ??? CHECK!!! */ 
  if (p != NULL){
    treeprint_stddev (p->left);
    free (p->left);
    *(p->loc) += p->weight*p->weight; /* write weight to final destination */
    treeprint_stddev (p->right);
    free (p->right);
  }
}


/***********************************************************************************/
/* Function: add_to_photonpath                                                     */
/* Description:                                                                    */
/*  Add a photon location to the photon path                                       */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                                 */
/***********************************************************************************/

#ifdef PRINT_PHOTON_PATH
static void add_to_photonpath (struct photon_path **p, double *x)
{
  int i=0;

  struct photon_path *new = calloc (1, sizeof (struct photon_path));

  for (i=0; i<3; i++)
    new->x[i] = x[i];

  new->next = *p;
  *p = new;
}
#endif


/***********************************************************************************/
/* Function: print_and_clear_photonpath                                            */
/* Description:                                                                    */
/*  Print all photon locations to stdout and free memory                           */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                                 */
/***********************************************************************************/

#ifdef PRINT_PHOTON_PATH
static void print_and_clear_photonpath (struct photon_path *p)
{
  if (p!=NULL) {
    print_and_clear_photonpath (p->next);
    fprintf (stdout, "%.6e %.6e %.6e\n", p->x[0], p->x[1], p->x[2]);
    free(p);
  }
}
#endif


/***********************************************************************************/
/* Function: setup_profiles1D                                             @62_30i@ */
/* Description:                                                                    */
/*  Initialize the 1D part of the atmosphere_struct with the 1D profiles passed    */
/*  as function arguments.                                                         */ 
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int setup_profiles1D ( int                n_caoth,
			      float            **dt_s,
			      float            **om_s,
			      float            **g1_s,
			      float            **g2_s,
			      float            **f_s,
			      float            **ds_s,
                              float            **re_s,
			      float             *zprof,
			      int                nlyr,
			      sample_struct     *sample, 
			      atmosphere_struct *atmos,
                              alis_struct       *alis 
                              )
{
  int lc=0, status=0;
  double dz_inv=0.0;
  int ids=0, iris=0, ivs=0, isp=0, iv=0;
  int ic=0;
  
  /* allocate memory for 1D profiles in atmosphere_struct *atmos */
  status = calloc_1D_atmos (atmos, nlyr, n_caoth, alis->nlambda_abs, alis->nlambda_ref);
  if (status!=0)
    return err_out ("Error %d returned by calloc_1D_atmos()\n", status);
  
  /* !!! attention !!!  z[0 ... atmos->Nz], (dsca->air)[0 ... atmos->Nz-1], ...     */
  /* distinguish between layer (0..atmos->Nz-1) and level (0..atmos->Nz) properties */
  
  atmos->Nc=alis->Nc;
  
  /* convert altitude km -> m      */
  for (lc=0; lc<atmos->Nz+1; lc++)
    atmos->Z[lc] = zprof[lc]*1000.0;
  
  /* for now only treat standard scattering and non-importance-sampled scattering */
  ids = MCSC_MODE_NORMAL;
  iris = MCRIS_MODE_NORMAL;

  if (sample->spectral_is){
    for (iv=0; iv<alis->nlambda_abs; iv++)
      atmos->lambda[iv]=alis->lambda[iv];
    for (iv=0; iv<alis->nlambda_ref; iv++)
      atmos->ilambda_ref[iv]=alis->ilambda_ref[iv];
  }
  
  if (alis->Nc>1){
    sample->concentration_is=1; 
    ASCII_calloc_double(&(atmos->kabs_scaled), atmos->Nc, atmos->Nz);
    ASCII_calloc_double(&(atmos->ksca_scaled), atmos->Nc, atmos->Nz);
  }
  
  for (lc=0; lc<atmos->Nz; lc++) {
    for (isp=1; isp<=n_caoth; isp++) { /* this necessitates that MCCAOTH_MOL=CAOTH_MOL+1;
                                          MCCAOTH_AER=CAOTH_AER+1*/
      atmos->ksca [ids][iris]->prof [isp][lc] = dt_s[isp-1][lc] * om_s[isp-1][lc];
      atmos->kabs            ->prof [isp][lc] = dt_s[isp-1][lc] * ( 1.0 - om_s[isp-1][lc] );
    }
    
    /* spectral molecular absorption */
    if (sample->spectral_is){
      /* Absorption is considered by importance sampling, for the
	 calculation itself it can be set to 0 */
      /* atmos->kabs->prof [MCCAOTH_MOL][lc]=0.0; */
        
      /* moderate scaling scattering coefficient can improve
	 efficiency, especially if only Rayleigh scattering*/
      /* is considered without surface reflections */
      /* atmos->ksca [ids][iris]->prof [MCCAOTH_MOL][lc]*=2; */
      
      /* FIXCE CHECK!!! */ 
      /* fprintf(stderr, "caoth %d\n", n_caoth); */
      for (iv=0; iv<alis->nlambda_abs; iv++){
	for (isp=1; isp<=n_caoth; isp++) {

        /* layers are counted different in uvspec and mystic*/
	  atmos->ksca_spectral[isp][iv][lc] = alis->dt[iv][isp-1][lc] * alis->om[iv][isp-1][lc];
	  atmos->kabs_spectral[isp][iv][lc] = alis->dt[iv][isp-1][lc] * ( 1.0 - alis->om[iv][isp-1][lc] );
	  atmos->ksca_spectral[isp][iv][lc] /= ( atmos->Z[lc+1] - atmos->Z[lc] );
	  atmos->kabs_spectral[isp][iv][lc] /= ( atmos->Z[lc+1] - atmos->Z[lc] );

	  atmos->ksca_spectral[MCCAOTH_TOT][iv][lc] += atmos->ksca_spectral[isp][iv][lc];
	  atmos->kabs_spectral[MCCAOTH_TOT][iv][lc] += atmos->kabs_spectral[isp][iv][lc];
       
        /*   FIXCE debugging output CHECK!!! */ 
        /*         if (lc==2) */
        /*           if (iv==0) */
        /*   fprintf(stderr, "lc %d lambda %g molabs %g aerabs %g aersca %g \n", lc, atmos->lambda[iv], */
        /*                   alis->tau_molabs[iv][atmos->Nz-1-lc], */
        /*                   alis->tau_aerabs[iv][atmos->Nz-1-lc], */
        /*                   alis->tau_aersca[iv][atmos->Nz-1-lc]); */
          /*   FIXCE debugging output */
          /*           if (lc==2) */
          /*             if (iv==0) */
          /*      fprintf(stderr, "lc %d lambda %g wcabs %g wcsca %g icsca %g \n", lc, atmos->lambda[iv], */
          /*                     alis->tau_wcabs[iv][atmos->Nz-1-lc], */
          /*                     alis->tau_wcsca[iv][atmos->Nz-1-lc], */
          /*                     alis->tau_icsca[iv][atmos->Nz-1-lc]); */
          /*  } */
        }
      }
    }
    
    /* concentration importance sampling (aerosol), to be generalized for caoth*/
    if (sample->concentration_is){
      isp=MCCAOTH_AER;
      for (ic=0; ic<atmos->Nc; ic++){
        
        atmos->kabs_scaled[ic][lc] = alis->aer_scaling_factors[ic] * atmos->kabs ->prof [isp][lc]
          / (atmos->Z[lc+1] - atmos->Z[lc]);
        
        atmos->ksca_scaled[ic][lc] = alis->aer_scaling_factors[ic] *
          atmos->ksca [ids][MCRIS_MODE_NORMAL]->prof [isp][lc]
          / (atmos->Z[lc+1] - atmos->Z[lc]);
      }
    }
  }

  /* now do delta_scaling */
  if ( sample->delta_scaling > -1.0 )
    for (lc=0; lc<atmos->Nz; lc++) {
      ids = MCSC_MODE_DELTA_SCALE;

      for (isp=1; isp<=n_caoth; isp++) /* this necessitates that MCCAOTH_MOL=CAOTH_MOL+1;
					   MCCAOTH_AER+1 */
	atmos->ksca [ids][iris]->prof [isp][lc] = dt_s[isp-1][lc] * om_s[isp-1][lc]
	  * ( 1.0 - ds_s[isp-1][lc] );
    }

  /* now copy into modified (IS) scattering */
  for (lc=0; lc<atmos->Nz; lc++)
    for (ids=0; ids<atmos->nscaDS; ids++)
      for (iris=1; iris<atmos->nscaRIS; iris++)
	for (isp=1; isp<=atmos->ksca [ids][iris]->n_caoth; isp++) /* tot not needed, therefore
								      we start at isp=1 */
	  atmos->ksca [ids][iris]->prof [isp][lc]
	    = atmos->ksca [ids][MCRIS_MODE_NORMAL]->prof [isp][lc];

#if HAVE_LIDAR
  if (sample->LLE_RIS_MAS || sample->RIS_MS) {
    mc_lidar_setup_ris_mas (atmos, sample);
  }
#endif

  if (sample->LLE_jacobian || sample->abs_jacobian)
    sample->Nz_jac = atmos->Nz;

  /* calculate total scattering optical depth */
  for (lc=0; lc<atmos->Nz; lc++) {
    for (ids=0; ids<atmos->nscaDS; ids++)
      for (iris=0; iris<atmos->nscaRIS; iris++) {
	atmos->ksca [ids][iris]->prof [MCCAOTH_TOT][lc] = 0.0;
	for (isp=1; isp<=atmos->ksca [ids][iris]->n_caoth; isp++)
	  atmos->ksca [ids][iris]->prof [MCCAOTH_TOT][lc]
	    += atmos->ksca [ids][iris]->prof [isp][lc];
      }

    atmos->kabs ->prof [MCCAOTH_TOT][lc] = 0.0;
    for (isp=1; isp<=atmos->kabs->n_caoth; isp++)
      atmos->kabs ->prof [MCCAOTH_TOT][lc] += atmos->kabs ->prof [isp][lc];
  }

  /* calculate total optical depths */

  ids = MCSC_MODE_NORMAL;
  iris = MCRIS_MODE_NORMAL;

  if (sample->LLE_D_DIS)
    atmos->ddis_eps = calloc( (size_t) atmos->Nz, sizeof(double));
  else
    atmos->ddis_eps=NULL;

  for (lc=atmos->Nz-1; lc>=0; lc--) {
    for (isp=0; isp<=atmos->ksca [ids][iris]->n_caoth; isp++) {
      atmos->tabstot[isp] += atmos->kabs->prof [isp][lc]; 
      atmos->tscatot[isp] += atmos->ksca [ids][iris]->prof [isp][lc];
    }

    if (sample->LLE_D_DIS)
      atmos->ddis_eps[lc]=1.0 - atmos->tscatot[MCCAOTH_TOT]/10.;
  }
  
  /* normalize with 1/dz */
  for (lc=0; lc<atmos->Nz; lc++) {
    dz_inv = 1.0 / (atmos->Z[lc+1] - atmos->Z[lc]);

    /* scattering */
    for (ids=0; ids<atmos->nscaDS; ids++)
      for (iris=0; iris<atmos->nscaRIS; iris++)
	for (isp=0; isp<=atmos->ksca [ids][iris]->n_caoth; isp++)
	  atmos->ksca [ids][iris]->prof [isp][lc] *= dz_inv;

    /* absorption */
    for (isp=0; isp<=atmos->kabs->n_caoth; isp++)
      atmos->kabs->prof [isp][lc] *= dz_inv;
  }

  /* setup extinction coefficient profile; Virtual IS may not be set up before here! */
  for (ids=0; ids<atmos->nscaDS; ids++)
    for (iris=0; iris<atmos->nscaRIS; iris++)
      for (ivs=0; ivs<atmos->nscaVIS; ivs++)
	for (lc=0; lc<atmos->Nz; lc++)
	  atmos->kext [ids][iris][ivs]->prof [MCCAOTH_TOT][lc]
	    = atmos->ksca [ids][iris]->prof [MCCAOTH_TOT][lc];

  /* convert asymmetry factor to double */
  for (lc=0; lc<atmos->Nz; lc++) {
    for (isp=MCCAOTH_AER; isp<=atmos->g1->n_caoth; isp++) {
      atmos->g1  ->prof [isp][lc] = g1_s [isp-1][lc];
      atmos->g2  ->prof [isp][lc] = g2_s [isp-1][lc];
      atmos->ff  ->prof [isp][lc] = f_s  [isp-1][lc];
      atmos->reff->prof [isp][lc] = re_s [isp-1][lc];
    }
  }

  return 0;
}


/***********************************************************************************/
/* Function: calloc_1D_atmos                                              @62_30i@ */
/* Description:                                                                    */
/*  Allocate memory for the 1D part of atmosphere_struct *atmos.                   */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int calloc_1D_atmos ( atmosphere_struct *atmos,
			     int                nlyr,
			     int                n_caoth,
			     int                nlambda_abs,
			     int                nlambda_ref )
{
  int ids=0, iris=0, ivs=0;
  int iv=0;

  atmos->Nz = nlyr;
  atmos->n_caoth = n_caoth;  

  /* allocate memory for profiles */
  atmos->Z = calloc ((size_t) (atmos->Nz+1), sizeof(double));
  if (atmos->Z==NULL)
    return mem_err_out ( "atmos->Z", ERROR_POSITION );

  atmos->kabs   = calloc_profile (n_caoth, atmos->Nz);
  if (atmos->kabs==NULL)
    return mem_err_out ( "atmos->kabs", ERROR_POSITION );
  
  if (nlambda_abs>0){
    atmos->nlambda_abs = nlambda_abs;
    /* FIXCE CHECK!!! */
    /* may be profile structure should be extended for spectral absorption ???*/
    atmos->lambda = calloc ( nlambda_abs, sizeof(float)); 
    
    atmos->kabs_spectral = calloc (n_caoth+1, sizeof(double**)); 
    atmos->ksca_spectral = calloc (n_caoth+1, sizeof(double**));
    for (ids=0; ids<=n_caoth; ids++){
      atmos->kabs_spectral[ids] = calloc ( nlambda_abs, sizeof(double*)); 
      atmos->ksca_spectral[ids] = calloc ( nlambda_abs, sizeof(double*)); 
      for (iv=0; iv<nlambda_abs; iv++){
        atmos->kabs_spectral[ids][iv] =  calloc (nlyr, sizeof(double));
        atmos->ksca_spectral[ids][iv] =  calloc (nlyr, sizeof(double));
      }
    }
  }

  /* FIXCE: Quick fix CHECK!!! */
  else{
    atmos->nlambda_abs=1;
    atmos->lambda = calloc ( atmos->nlambda_abs, sizeof(float)); 
    atmos->lambda[0]=0;
  }
    
  if (nlambda_ref>0){
    atmos->nlambda_ref = nlambda_ref;
    atmos->ilambda_ref = calloc ( nlambda_ref, sizeof(int) );
  }

  atmos->ksca   = calloc (atmos->nscaDS, sizeof(profile *));
  for (ids=0; ids<atmos->nscaDS; ids++)
    atmos->ksca[ids] = calloc (atmos->nscaRIS, sizeof(profile));
  for (ids=0; ids<atmos->nscaDS; ids++)
    for (iris=0; iris<atmos->nscaRIS; iris++) {
      atmos->ksca[ids][iris]  = calloc_profile (n_caoth, atmos->Nz);
      if (atmos->ksca[ids][iris]==NULL)
	return mem_err_out ( "atmos->ksca[][]", ERROR_POSITION );
    }

  atmos->kext   = calloc (atmos->nscaDS, sizeof(profile **));
  for (ids=0; ids<atmos->nscaDS; ids++)
    atmos->kext[ids] = calloc (atmos->nscaRIS, sizeof(profile *));
  for (ids=0; ids<atmos->nscaDS; ids++)
    for (iris=0; iris<atmos->nscaRIS; iris++)
      atmos->kext[ids][iris] = calloc (atmos->nscaVIS, sizeof(profile));
  for (ids=0; ids<atmos->nscaDS; ids++)
    for (iris=0; iris<atmos->nscaRIS; iris++)
      for (ivs=0; ivs<atmos->nscaVIS; ivs++) {
	atmos->kext[ids][iris][ivs] = calloc_profile (0, atmos->Nz);
	if (atmos->kext[ids][iris][ivs]==NULL)
	  return mem_err_out ( "atmos->kext[][][]", ERROR_POSITION );
    }

  atmos->g1     = calloc_profile (n_caoth, atmos->Nz);
  if (atmos->g1==NULL)
    return mem_err_out ( "atmos->g1", ERROR_POSITION );

  atmos->g2     = calloc_profile (n_caoth, atmos->Nz);
  if (atmos->g2==NULL)
    return mem_err_out ( "atmos->g2", ERROR_POSITION );

  atmos->ff     = calloc_profile (n_caoth, atmos->Nz);
  if (atmos->ff==NULL)
    return mem_err_out ( "atmos->ff", ERROR_POSITION );

  /* ??? set ff to 1 ??? CHECK!!! */

  atmos->reff   = calloc_profile (n_caoth, atmos->Nz);
  if (atmos->reff==NULL)
    return mem_err_out ( "atmos->reff", ERROR_POSITION );

  atmos->tabstot = calloc((size_t) n_caoth+1, sizeof(double));
  if (atmos->tabstot==NULL)
    return mem_err_out ( "atmos->tabstot", ERROR_POSITION );

  atmos->tscatot = calloc((size_t) n_caoth+1, sizeof(double));
  if (atmos->tscatot==NULL)
    return mem_err_out ( "atmos->tscatot", ERROR_POSITION );

  return 0;
}


/***********************************************************************************/
/* Function: calloc_profile                                               @62_30i@ */
/* Description:                                                                    */
/*  Allocate memory for struct profile and initialize structure.                   */ 
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static profile *calloc_profile (int n_caoth, int n)
{
  int isp=0;

  profile *p = calloc(1, sizeof(profile));
  if (p==NULL)
    return NULL;

  p->n_caoth = n_caoth;
  p->n = n;

  p->prof =  calloc((size_t) n_caoth+1, sizeof(double *));
  if (p->prof==NULL)
    return NULL;

  for (isp=0; isp<n_caoth+1; isp++) {
    p->prof[isp] = calloc((size_t) n, sizeof(double));

    if (p->prof[isp]==NULL)
      return NULL;
  }
  
  return p;
}


/***********************************************************************************/
/* Function: calloc_profile3D                                             @62_30i@ */
/* Description:                                                                    */
/*  Allocate memory for struct profile3D and initialize structure.                 */ 
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static profile3D *calloc_profile3D ( int   n_caoth,
				     int  *tocalloc,
				     int   Nz,
				     int   Nx,
				     int   Ny,
				     int **threed )
{
  int kc=0, ic=0, isp=0;

  profile3D *p = calloc(1, sizeof(profile3D));
  if (p==NULL)
    return NULL;

  p->n_caoth = n_caoth;

  /* this vector tells us which caoth need 3d allocation (incl. tot) */
  p->tocalloc = calloc((size_t) n_caoth+1, sizeof(int));
  if (p->tocalloc==NULL)
    return NULL;

  for (isp=0; isp<=n_caoth; isp++)
    p->tocalloc[isp] = tocalloc[isp];

  p->Nz = Nz;
  p->Nx = Nx;
  p->Ny = Ny;

  p->nthreed = calloc ((size_t) n_caoth+1, sizeof(int));
  if (p->nthreed==NULL)
    return NULL;
  p->threed = calloc ((size_t) n_caoth+1, sizeof(int *));
  if (p->threed==NULL)
    return NULL;
  for (isp=0; isp<=n_caoth; isp++) {
    p->threed[isp] = calloc ((size_t) Nz, sizeof(int));
    if (p->threed[isp]==NULL)
      return NULL;
  }

  for (isp=0; isp<=n_caoth; isp++)
    for (kc=0; kc<Nz; kc++) {
      p->threed[isp][kc] = threed[isp][kc];
      if (p->threed[isp][kc]>=1)
	p->nthreed[isp]++;
    }

  p->prof =  calloc((size_t) n_caoth+1, sizeof(float ***));
  if (p->prof==NULL)
    return NULL;

  for (isp=0; isp<=n_caoth; isp++) {
    if (p->tocalloc[isp]) {
      p->prof[isp] = calloc((size_t) Nz, sizeof(float **));

      if (p->prof[isp]==NULL)
	return NULL;
    }
    else
      p->prof[isp]=NULL;
  }

  for (isp=0; isp<=n_caoth; isp++)
    if (p->tocalloc[isp]) {
      for (kc=0; kc<Nz; kc++)
	if (p->threed[isp][kc]>=1) {
	  p->prof[isp][kc] = calloc((size_t) Nx, sizeof(float *));
	  for (ic=0; ic<Nx; ic++)
	    p->prof[isp][kc][ic] = calloc((size_t) Ny, sizeof(float));
	}
    }

  return p; 
}


/***********************************************************************************/
/* Functions: v_mult, v_sub, vn_mult                                      @62_30i@ */
/* Description:                                                                    */
/*   short cuts for vector*vector -> number, vector - vector -> vector,            */
/*                  vector*number -> vector                                        */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int v_cross_product_norm (double *a, double *b, double *c)
{
  double n=0.0;
  int i=0;

  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];

  /* Normalization*/
  for(i=0; i<3; i++)
    n += c[i] * c[i];

  /* a and b are parallel */
  if (n==0.0)
    return 1;

  n=sqrt(n);

  for(i=0; i<3; i++)
    c[i]/=n; 

  return 0;
}

void mat_v_mult (double **a, double *b, double *c)
{
  int i=0, j=0;

  for (i=0; i<4; i++) {
    c[i]=0.0;
    for (j=0; j<4; j++)
      c[i]+=a[i][j]*b[j];
  }
}

void v_mult (double *a, double *b, double *c)
{
  *c = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void v_sub (double *a, double *b, double *c)
{
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2];
}

void v_add (double *a, double *b, double *c)
{
  c[0] = a[0] + b[0];
  c[1] = a[1] + b[1];
  c[2] = a[2] + b[2];
}

void vn_mult (double *a, double *b, double *c)
{
  c[0] = a[0] * *b;
  c[1] = a[1] * *b;
  c[2] = a[2] * *b;
}

void v_cross_product (double *a, double *b, double *c)
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

void v_mult_mu (double *a, double *b, double *c)
{
  *c = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

  if (*c > 1.0) {
    if (*c > 1.0 + MC_EPSILON)
      *c = NOT_A_NUMBER;
    *c = 1.0;
  }

  if (*c < -1.0) {
    if (*c < -1.0 - MC_EPSILON)
      *c = NOT_A_NUMBER;
    *c = -1.0;
  }

  if (*c == NOT_A_NUMBER)
    fprintf(stderr,"Error! v_mult_mu() has unphysical result %e \n",*c);

}

void v_fabs (double *a, double *c)
{
  *c = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

void v_neg (double *a)
{
  a[0] = -a[0];
  a[1] = -a[1];
  a[2] = -a[2];
}

double cosd ( double angle )
{
  return cos ( angle * PI / 180.0 );
}

double sind ( double angle )
{
  return sin ( angle * PI / 180.0 );
}

double tand ( double angle )
{
  return tan ( angle * PI / 180.0 );
}

double acosd ( double rad )
{
  return acos ( rad ) * 180.0 / PI;
}

double asind ( double rad )
{
  return asin ( rad ) * 180.0 / PI;
}

double atand ( double rad )
{
  return atan ( rad ) * 180.0 / PI;
}


/***********************************************************************************/
/* Functions: get_ksca_tot, get_kabs_tot, get_ksca, etc                   @62_30i@ */
/* Description:                                                                    */
/*   return the value of an atmospheric quantity at the location of the photon     */
/*   depends of whether atmos->threed, and whether to deltascale                   */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

/* total absorption optical depth */
static inline double get_kabs_tot (atmosphere_struct *atmos, photon_struct *p)
{
  if (atmos->threed[MCCAOTH_TOT][p->kc]>=1)
    return (atmos->kabs3D->prof [MCCAOTH_TOT])[p->kc][p->ic][p->jc];
  else
    if (p->iv_alis!=-1)
      return atmos->kabs_spectral[MCCAOTH_TOT][p->iv_alis][p->kc];
    else
      return (atmos->kabs->prof [MCCAOTH_TOT])[p->kc];
}

static inline double get_g1 (atmosphere_struct *atmos, photon_struct *p, int isp)
{
  if (atmos->threed[isp][p->kc]>=1)
    return (atmos->g1_3D->prof [isp]) [p->kc] [p->ic] [p->jc];
  else
    return (atmos->g1->prof [isp])[p->kc];
}

static inline double get_g2 (atmosphere_struct *atmos, photon_struct *p, int isp)
{
  if (atmos->threed[isp][p->kc]>=1)
    return (atmos->g2_3D->prof [isp]) [p->kc] [p->ic] [p->jc];
  else
    return (atmos->g2->prof [isp])[p->kc];
}

static inline double get_ff (atmosphere_struct *atmos, photon_struct *p, int isp)
{
  if (atmos->threed[isp][p->kc]>=1)
    return (atmos->ff_3D->prof [isp]) [p->kc] [p->ic] [p->jc];
  else
    return (atmos->ff->prof [isp])[p->kc];
}

static inline double get_reff (atmosphere_struct *atmos, photon_struct *p, int isp)
{
  if (atmos->threed[isp][p->kc]>=1)
    return (atmos->reff_3D->prof [isp]) [p->kc] [p->ic] [p->jc];
  else
    return (atmos->reff->prof [isp])[p->kc];
}

#ifdef CLDPRP
static double inline get_dxlwc (atmosphere_struct *atmos, photon_struct *p, int isp)
{
  if (atmos->threed[isp][p->kc]>=1)
    return (atmos->dxlwc_3D->prof [isp]) [p->kc] [p->ic] [p->jc];
  else
    return 0.0;
}

static double inline get_dylwc (atmosphere_struct *atmos, photon_struct *p, int isp)
{
  if (atmos->threed[isp][p->kc]>=1)
    return (atmos->dylwc_3D->prof [isp]) [p->kc] [p->ic] [p->jc];
  else
    return 0.0;
}

static double inline get_dzlwc (atmosphere_struct *atmos, photon_struct *p, int isp)
{
  if (atmos->threed[isp][p->kc]>=1)
    return (atmos->dzlwc_3D->prof [isp]) [p->kc] [p->ic] [p->jc];
  else
    return 0.0;
}
#endif

/* scattering optical depth of choice */
double get_ksca (atmosphere_struct *atmos, photon_struct *p, int isp)
{
  if ( isp < 0 || isp > atmos->n_caoth ) {
    fprintf(stderr,"Error in get_ksca, no such type %d of scattering!!!\n", isp);
    return NOT_A_NUMBER;
  }
  if ( atmos->threed[isp][p->kc]>=1 )
    return atmos->ksca3D [p->SC_mode][MCRIS_MODE_NORMAL]->prof [isp] [p->kc] [p->ic] [p->jc];
  else
    if (p->iv_alis!=-1)
      return atmos->ksca_spectral[isp][p->iv_alis][p->kc];
    else
      return atmos->ksca   [p->SC_mode][MCRIS_MODE_NORMAL]->prof [isp] [p->kc];
}

/* scattering optical depth of choice */
static inline double get_ksca_spectral (atmosphere_struct *atmos, photon_struct *p, int isp, int iv)
{
  if ( isp < 0 || isp > atmos->n_caoth ) {
    fprintf(stderr,"Error in get_ksca, no such type %d of scattering!!!\n", isp);
    return NOT_A_NUMBER;
  }
  if ( atmos->threed[isp][p->kc]>=1 )
    return atmos->ksca3D [p->SC_mode][MCRIS_MODE_NORMAL]->prof [isp] [p->kc] [p->ic] [p->jc];
  else
    return atmos->ksca_spectral[isp][iv][p->kc];
}


#ifdef CLDPRP
/* absorption optical depth of choice */
double get_kabs (atmosphere_struct *atmos, photon_struct *p, int isp)
{
  if ( isp < 0 || isp > atmos->n_caoth ) {
    fprintf(stderr,"Error in get_kabs, no such type %d of absorption!!!\n", isp);
    return 0.0/0.0;
  }
  if ( atmos->threed[isp][p->kc]>=1 )
    return atmos->kabs3D->prof [isp] [p->kc] [p->ic] [p->jc];
  else
    if (p->iv_alis!=-1)
      return atmos->kabs_spectral[isp][p->iv_alis][p->kc];
    else
      return atmos->kabs->prof [isp] [p->kc];
}
#endif

#ifdef NEWRISQIDD
/* modified scattering optical depth of choice */
double get_kscaIS (atmosphere_struct *atmos, photon_struct *p, int isp)
{
  double ksca=0.0, ksca2=0.0;

  if ( isp < 0 || isp > atmos->n_caoth ) {
    fprintf(stderr,"Error in get_kscaIS, no such type %d of scattering!!!\n", isp);
    return NOT_A_NUMBER;
  }
  if ( atmos->threed[isp][p->kc]>=1 )
    ksca = atmos->ksca3D [p->SC_mode][p->RIS_mode]->prof [isp] [p->kc] [p->ic] [p->jc];
  else
    if (p->iv_alis!=-1)
      ksca = atmos->ksca_spectral[isp][p->iv_alis][p->kc];
    else
      ksca = atmos->ksca   [p->SC_mode][p->RIS_mode]->prof [isp] [p->kc];

  if (p->risqidd_beta != 0.0) { /* here probably error  */
    if ( atmos->threed[isp][p->kc]>=1 )
      ksca2 = atmos->ksca3D [p->SC_mode][p->RIS_mode]->prof [MCCAOTH_TOT] [p->kc] [p->ic] [p->jc];
    else
      ksca2 = atmos->ksca   [p->SC_mode][p->RIS_mode]->prof [MCCAOTH_TOT] [p->kc];

    if (p->risqidd_beta > ksca2) /* this should always be fulfilled */
      ksca *= p->risqidd_beta / ksca2;
  }

  return ksca*atmos->ris_factor;
}
#else
/* modified scattering optical depth of choice */
double get_kscaIS (atmosphere_struct *atmos, photon_struct *p, int isp)
{
  if ( isp < 0 || isp > atmos->n_caoth ) {
    fprintf(stderr,"Error in get_kscaIS, no such type %d of scattering!!!\n", isp);
    return NOT_A_NUMBER;
  }
  if ( atmos->threed[isp][p->kc]>=1 )
    return atmos->ris_factor*atmos->ksca3D [p->SC_mode][p->RIS_mode]->prof [isp] [p->kc] [p->ic] [p->jc];
  else
    if (p->iv_alis!=-1)
      return atmos->ris_factor*atmos->ksca_spectral[isp][p->iv_alis][p->kc];
    else
      return atmos->ris_factor*atmos->ksca   [p->SC_mode][p->RIS_mode]->prof [isp] [p->kc];
}
#endif

/* modified extinction scattering optical depth */
/* this includes virtual scattering rates       */
static inline double get_kext (atmosphere_struct *atmos, photon_struct *p)
{
  double ksca=0.0;
  if (atmos->threed[MCCAOTH_TOT][p->kc]>=1)
    ksca = atmos->kext3D [p->SC_mode][p->RIS_mode][p->VIS_mode]->prof [MCCAOTH_TOT] [p->kc] [p->ic] [p->jc];
  else
    if (p->iv_alis!=-1)
      ksca = atmos->ksca_spectral[MCCAOTH_TOT][p->iv_alis][p->kc];
    else
      ksca = atmos->kext   [p->SC_mode][p->RIS_mode][p->VIS_mode]->prof [MCCAOTH_TOT] [p->kc];

  if (p->vis_beta != 0.0 && ksca > 0. && ksca < p->vis_beta)
      ksca = p->vis_beta;

  /* special case for VIS-FOD */
  if (!p->lest.behind_detector && ksca > 0. && ksca < p->lest.vis_fod_kext)
      ksca = p->lest.vis_fod_kext;

  return ksca*atmos->ris_factor;
}


/* Assign an altitude to a photon and determine the respective layer number; */
/* required by generate_photon. This function is used by generate_photon     */
/* and MUST BE CALLED AFTER THE PHOTON DIRECTION HAS BEEN ASSIGNED BECAUSE   */
/* THE LAYER NUMBER MIGHT DEPEND ON THE PHOTON DIRECTION!                    */

int set_photon_z (float z, atmosphere_struct *atmos, photon_struct *p)
{
  int hit=0;
  static int jlo=0;
  float zact=0;

  float *zz = atmos->Z;
  int n = atmos->Nz+1;

  zact = z;

  if (zz[0]!=0) {
    fprintf (stderr, "Error, assuming that the profile starts at z=0\n");
    return -1;
  }

  if (n<2) {
    fprintf (stderr, "Error, assuming at least two vertical levels\n");
    return -1;
  }

  /* restrict start position to available range */
  if (zact < 0.0)
    zact = 0.0;
  else {
    if (zact > zz[n-1])
      zact = zz[n-1];
  }

  /* an upward photon at top of atmosphere is moved downward a bit */
  /* in order to avoid numerical problems; this is required e.g.   */
  /* to calculate downward irradiance at TOA in backward mode      */
  if (zact == zz[n-1] && p->dir.dx[2]>0)
    zact *= (1.0-MC_EPSILON);

  /* a downward photon at the lowest level is moved upward a bit */
  /* in order to avoid numerical problems; this is required e.g. */
  /* to calculate surface upward irradiance in backward mode     */
  if (zact == 0.0 && p->dir.dx[2]<0)
    zact += MC_EPSILON;

  hunt_modified (zz, n, zact, &jlo, &hit);

  p->kc = jlo;
  p->x[2] = zact;   /* p->x[2] is a double while all other z's are float - does that make sense? */

  /* downward moving photons starting at a level are assigned the index of the lower box */
  if (hit && p->dir.dx[2]<0)
    p->kc--;

  return 0;
}


/*****************************************************************/
/* slightly modified version of hunt() from numerical recipes;   */
/* the layer index kc is assigned so that zz[kc] <= z < zz[kc+1] */
/*****************************************************************/

static void hunt_modified (float *xx, int n, float x, int *jlo, int *hit)
/* Given an array xx[1..n], and given a value x, returns a value jlo such that
   x is between xx[jlo] and xx[jlo+1]. xx[1..n] must be monotonic, either
   increasing or decreasing. jlo=0 or jlo=n is returned to indicate that x is out
   of range. jlo on input is taken as the initial guess for jlo on output. */
{
  int jm=0, jhi=0, inc=0;
  int ascnd=0;

  *hit=0;
  xx-=1;

  ascnd=(xx[n] >= xx[1]); 
  if (*jlo <= 0 || *jlo > n) { 
    *jlo=0;
    jhi=n+1;
  } else {
    inc=1;
    /* hint: setting brackets around first logical, i.e. ((x >= xx[*jlo] ) == ascnd), will fix the compiler warnings */ 
    if ( (x >= xx[*jlo]) == ascnd ) {
      if (*jlo == n) return;
      jhi=(*jlo)+1;
      while ( (x >= xx[jhi]) == ascnd ) {
	*jlo=jhi;
	inc += inc;
	jhi=(*jlo)+inc;
	if (jhi > n) {
	  jhi=n+1;
	  break;
	}
      }
    } else {
      if (*jlo == 1) {
	*jlo=0;
	return;
      }
      jhi=(*jlo)--;
      while ( (x < xx[*jlo]) == ascnd ) {
	jhi=(*jlo);
	inc <<= 1;
	if (inc >= jhi) {
	  *jlo=0;
	  break;
	}
	else *jlo=jhi-inc;
      }
    }
  }
  while (jhi-(*jlo) != 1) {
    jm=(jhi+(*jlo)) >> 1;
    if ( (x >= xx[jm]) == ascnd )
      *jlo=jm;
    else
      jhi=jm;
  }
  if (x == xx[n]) 
      /* changed that, BM */
      /* *jlo=n-1;        */
    *jlo=n;

    
  if (x == xx[1]) 
    *jlo=1;

  *jlo-=1;

  xx+=1;

  if (xx[*jlo] == x) 
    *hit = 1;

  return;
}


/***********************************************************************************/
/* Function: stokes_vector_sca                                            @62_30i@ */
/* Description:                                                                    */
/*  Calculates the scattered Stokes (weight) vector.                               */
/*  The Stokes (weight) vector is used to compute polarized radiative transfer.    */
/*  (this function works only for forward tracing and might be removed later.)     */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Claudia Emde                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/
static void stokes_vector_sca(double *stokes_vector,
			      double *P,
			      double *dx_sca,
			      double *dx_inc,
                              double  phi_sun,
                              double  phi_det)
{
  
  int i=0, j=0, status=0;
  double stokes_vector_inc[4];
  double N=0.0; 
  double n_z[3] = { 0, 0, 1 };
  double n_ort_sca[3], n_ort_ref[3] = { 0, 0, 0 }, n_ort_sin[3];
  double Z[4][4];
  double cos2alpha1=1.0, sin2alpha1=0.0, cosalpha=0.0;
  double cos2alpha2=1.0, sin2alpha2=0.0, sinalpha=0.0;
  double n_ort_sca_z=0.0; 
  
  /* Copy Stokes vector */
  for (i=0; i<4; i++)
    stokes_vector_inc[i]=stokes_vector[i];

  /* Calculate vector orthogonal to scattering plane */
  status = v_cross_product_norm ( dx_inc, dx_sca, n_ort_sca );

  n_ort_sca_z=n_ort_sca[2];
  
  /*incoming and scattered direction are not parallel */
  if(!status){
    /* orthogonal vector to p_inc-z-plane*/
    status = v_cross_product_norm ( dx_inc, n_z, n_ort_ref );
    
    if (status) {
      /* orthogonal vector to (z, phi_esc)-plane*/
      n_ort_ref[0] = -cosd(phi_sun);
      n_ort_ref[1] =  sind(phi_sun);
      n_ort_ref[2] =  0.0;
    
      n_ort_sca_z=-sind(phi_sun)*dx_sca[1]+cosd(phi_sun)*dx_sca[0];
    }

    /* Calculate rotation angle 1 */
    v_mult_mu ( n_ort_sca, n_ort_ref, &cosalpha );
    
    /* Calculate sinalpha using vector product */
    v_cross_product ( n_ort_sca, n_ort_ref, n_ort_sin );
    v_mult ( n_ort_sin, n_ort_sin, &N );
    sinalpha = - sqrt(N);
    
    cos2alpha1=cosalpha*cosalpha-sinalpha*sinalpha;
    sin2alpha1=2.0*sinalpha*cosalpha;
    
    /* orthogonal to p_sca-z-plane*/
    status = v_cross_product_norm ( dx_sca, n_z, n_ort_ref );

    if (status) {
      /* orthogonal vector to (z, phi_esc)-plane*/
      n_ort_ref[0] = -cosd(phi_det);
      n_ort_ref[1] =  sind(phi_det);
      n_ort_ref[2] =  0.0;
    
      n_ort_sca_z = -dx_inc[0]*cosd(phi_det)+dx_inc[1]*sind(phi_det);
    }
    
    /* Calculate rotation angle 2 */
    v_mult_mu ( n_ort_sca, n_ort_ref, &cosalpha );

    /* Calculate sinalpha using vector product */
    v_cross_product ( n_ort_sca, n_ort_ref, n_ort_sin );
    v_mult ( n_ort_sin, n_ort_sin, &N );
    sinalpha = + sqrt(N);
    
    cos2alpha2=cosalpha*cosalpha-sinalpha*sinalpha;
    sin2alpha2=2.0*sinalpha*cosalpha;
  }
  
  
  /* Exact forward or backward scattering, no rotation required, as initialized  */
  /*      cos2alpha1=1.0; */
  /*      cos2alpha2=1.0;  */
  /*      sin2alpha1=0.0;  */
  /*      sin2alpha2=0.0;  */

  /* The sign of the following elements depends on orientation */
  /* The z-component of the cross product of dx_inc and dx_sca is the same as the z-component */
  /* of the cross product of dx_inc_horz and dx_sca_horz (vectors projected in x-y-plane.) */
  /* The orientation of the normal vector to dx_inc_horz and dx_sca_horz is equivalent to  */
  /* whether the azimutal difference between the two vectors in the x-y plane is smaller  */
  /* or larger 180 degrees. */
  if(n_ort_sca_z<0.) {
    sin2alpha1=-sin2alpha1;
    sin2alpha2=-sin2alpha2;
  }

  Z[0][0] =  P[0];
  Z[0][1] =  P[1]*cos2alpha1; 
  Z[0][2] = -P[1]*sin2alpha1; 
  Z[0][3] =  0.0;
  Z[1][0] =  P[1]*cos2alpha2; 
  Z[1][1] =  P[4]*cos2alpha1*cos2alpha2 - P[2]*sin2alpha1*sin2alpha2;
  Z[1][2] = -P[4]*sin2alpha1*cos2alpha2 - P[2]*cos2alpha1*sin2alpha2;
  Z[1][3] = -P[3]*sin2alpha2;  
  Z[2][0] =  P[1]*sin2alpha2;
  Z[2][1] =  P[4]*cos2alpha1*sin2alpha2 + P[2]*sin2alpha1*cos2alpha2;
  Z[2][2] = -P[4]*sin2alpha1*sin2alpha2 + P[2]*cos2alpha1*cos2alpha2;
  Z[2][3] =  P[3]*cos2alpha2; 
  Z[3][0] =  0.0;
  Z[3][1] = -P[3]*sin2alpha1;
  Z[3][2] = -P[3]*cos2alpha1;
  Z[3][3] =  P[5];
  
  /* scattering matrix multiplication */
  for (i=0; i<4; i++){
    stokes_vector[i]=0.0;
    for (j=0; j<4; j++){
      stokes_vector[i]+=Z[i][j]*stokes_vector_inc[j];
    }
    stokes_vector[i]/=P[0];
  }
}

/****************************************************************************/
/* Function: phase_matrix_mult                                     @62_30i@ */
/* Description:                                                             */
/*  Multiplies the phase matrix of a scattering event with the photon       */
/*  phase matrix (product of all phase matrices along the photon path).     */
/* Parameters:                                                              */
/* Return value:                                                            */
/* Example:                                                                 */
/* Files:                                                                   */
/* Known bugs:                                                              */
/* Author: Claudia Emde                                                     */
/*                                                                 @i62_30@ */
/****************************************************************************/
void phase_matrix_mult ( double          **phase_matrix,
			 double           *P,
			 double           *dx_sca,
			 double           *dx_inc,
			 double            phi_sun,
			 double            phi_det,
			 int               escape,
			 int               scattercounter,
			 int               backward )
{
  
  int i=0, j=0, k=0, status=0;
  double phase_matrix_inc[4][4];
  double N=0.0;
  double n_ort_sca[3], n_ort_ref[3], n_ort_sin[3];
  double n_z[3] = { 0, 0, 1 };
  double Z[4][4];
  double cos2alpha1=1.0, sin2alpha1=0.0, cosalpha=0.0;
  double cos2alpha2=1.0, sin2alpha2=0.0, sinalpha=0.0;
  double *dx1=NULL, *dx2=NULL;
  double n_ort_sca_z=0.0;

  /* Copy photon phase matrix */
  for (i=0; i<4; i++)
    for(j=0; j<4; j++)
      phase_matrix_inc[i][j]=phase_matrix[i][j];

  if(backward){
    /*exchange incoming and scattered direction for backward */
    dx2=dx_inc;
    dx1=dx_sca;
    n_z[2]=-1;
  }
  else {
    dx1=dx_inc;
    dx2=dx_sca;
  }

  /* Calculate vector orthogonal to scattering plane */
  status = v_cross_product_norm ( dx1, dx2, n_ort_sca );

  n_ort_sca_z=n_ort_sca[2];

  /* incoming and scattered direction are not parallel */
  if (!status) {
    /* orthogonal vector to dx1-z-plane*/
    status = v_cross_product_norm ( dx1, n_z, n_ort_ref );

    if (status) {
      /* Incoming direction || z-axis happens if SZA=0 */
      /* Polarization plane (z, phi0)                  */
      if(!(backward && escape) && !(!backward && scattercounter==0) ){
        fprintf(stderr, "scattercounter %d \n", scattercounter); 
	fprintf(stderr, "\n");
        fprintf(stderr, "*** Warning, polarization plane of photon lost, because \n");
        fprintf(stderr, "*** photon is scattered exactly in z-direction.\n"); 
        fprintf(stderr, "*** If this happens too often please contact Claudia.\n"); 
	fprintf(stderr, "\n");
      }
      
      /* orthogonal vector to (z, phi_esc)-plane*/
      n_ort_ref[0] = -cosd(phi_sun);
      n_ort_ref[1] =  sind(phi_sun);
      n_ort_ref[2] =  0.0;

      n_ort_sca_z=-sind(phi_sun)*dx2[1]+cosd(phi_sun)*dx2[0];
      if(backward)
         n_ort_sca_z*=-1;
    }

    /* Calculate rotation angle 1 */
    v_mult_mu ( n_ort_sca, n_ort_ref, &cosalpha );

    /* Calculate sinalpha using vector product */
    v_cross_product ( n_ort_sca, n_ort_ref, n_ort_sin );
    v_mult ( n_ort_sin, n_ort_sin, &N );
    sinalpha = - sqrt(N);

    cos2alpha1=cosalpha*cosalpha-sinalpha*sinalpha;
    sin2alpha1=2.0*sinalpha*cosalpha;

    /* orthogonal to dx2z-plane*/
    status = v_cross_product_norm ( dx2, n_z, n_ort_ref );

    if (status) {
      /* In case of umu=+/- 1 we need to consider this special case: */ 
      /* In forward mode and umu +/-1 the last rotation needs to be in the plane */
      /* defined by (z, phi_sensor). */
      /* In backward mode the first rotation.*/
      if(!(!backward && escape) && !(backward && scattercounter==0) ){
        fprintf(stderr, "scattercounter %d \n", scattercounter); 
	fprintf (stderr, "\n");
        fprintf(stderr, "*** Warning, polarization plane of photon lost, because \n");
        fprintf(stderr, "*** photon is scattered exactly in z-direction.\n"); 
        fprintf(stderr, "*** If this happens too often please contact Claudia.\n"); 
	fprintf (stderr, "\n");
      }

      /* orthogonal vector to (z, phi_esc)-plane*/
      n_ort_ref[0] = -cosd(phi_det);
      n_ort_ref[1] =  sind(phi_det);
      n_ort_ref[2] =  0.0;

      n_ort_sca_z = -dx1[0]*cosd(phi_det)+dx1[1]*sind(phi_det);
      
      if(backward)
        n_ort_sca_z*=-1;
    }

    /* Calculate rotation angle 2 */
    v_mult_mu ( n_ort_sca, n_ort_ref, &cosalpha );

    /* Calculate sinalpha using vector product */
    v_cross_product ( n_ort_sca, n_ort_ref, n_ort_sin );
    v_mult ( n_ort_sin, n_ort_sin, &N );
    sinalpha = + sqrt(N);

    cos2alpha2=cosalpha*cosalpha-sinalpha*sinalpha;
    sin2alpha2=2.0*sinalpha*cosalpha;
  }
  
  /* Exact forward or backward scattering, no rotation required, as initialized  */
  /*      cos2alpha1=1.0; */
  /*      cos2alpha2=1.0;  */
  /*      sin2alpha1=0.0;  */
  /*      sin2alpha2=0.0;  */

  /* The sign of the following elements depends on orientation */
  /* The z-component of the cross product of dx_inc and dx_sca is the same as the z-component */
  /* of the cross product of dx_inc_horz and dx_sca_horz (vectors projected in x-y-plane.) */
  /* The orientation of the normal vector to dx_inc_horz and dx_sca_horz is equivalent to  */
  /* whether the azimutal difference between the two vectors in the x-y plane is smaller  */
  /* or larger 180 degrees. */
  if(n_ort_sca_z<0.) {
    sin2alpha1=-sin2alpha1;
    sin2alpha2=-sin2alpha2;
  }

  Z[0][0] =  P[0];
  Z[0][1] =  P[1]*cos2alpha1; 
  Z[0][2] = -P[1]*sin2alpha1; 
  Z[0][3] =  0.0;
  Z[1][0] =  P[1]*cos2alpha2; 
  Z[1][1] =  P[4]*cos2alpha1*cos2alpha2 - P[2]*sin2alpha1*sin2alpha2;
  Z[1][2] = -P[4]*sin2alpha1*cos2alpha2 - P[2]*cos2alpha1*sin2alpha2;
  Z[1][3] = -P[3]*sin2alpha2;  
  Z[2][0] =  P[1]*sin2alpha2;
  Z[2][1] =  P[4]*cos2alpha1*sin2alpha2 + P[2]*sin2alpha1*cos2alpha2;
  Z[2][2] = -P[4]*sin2alpha1*sin2alpha2 + P[2]*cos2alpha1*cos2alpha2;
  Z[2][3] =  P[3]*cos2alpha2; 
  Z[3][0] =  0.0;
  Z[3][1] = -P[3]*sin2alpha1;
  Z[3][2] = -P[3]*cos2alpha1;
  Z[3][3] =  P[5];
 
  /* scattering matrix multiplication */
  for (i=0; i<4; i++){
    for (j=0; j<4; j++){
      phase_matrix[i][j]=0.0;
      for (k=0; k<4; k++){
        if (backward)
          phase_matrix[i][j]+=phase_matrix_inc[i][k]*Z[k][j];
        else
          phase_matrix[i][j]+=Z[i][k]*phase_matrix_inc[k][j];
      }
      phase_matrix[i][j]/=P[0];
    }
  }

}

/****************************************************************************/
/* Function: reflect_matrix_rotate                                 @62_30i@ */
/* Description:                                                             */
/*  Rotates the phase matrix of a reflection event and multiplies it with   */
/*  the photon phase matrix (product of all p. m.s along the photon path).  */
/* Parameters:                                                              */
/* Return value:                                                            */
/* Example:                                                                 */
/* Files:                                                                   */
/* Known bugs:                                                              */
/* Author: Robert Buras                                                     */
/*                                                                 @i62_30@ */
/****************************************************************************/

void reflect_matrix_rotate ( double          **phase_matrix,
			     double          **R,
			     double           *dx_sca,
			     double           *dx_inc,
			     double           *n_hor,
			     double            phi_source,
			     double            phi_target,
			     int               escape,
			     int               scattercounter,
			     int               backward,
			     double           *rotmat1,//XXX
			     double           *rotmat2 )//XXX
{
  
  int i=0, j=0, k=0, status=0;
  double phase_matrix_inc[4][4];
  double N=0.0;
  double n_ort_in[3], n_ort_out[3], n_ort_ref[3], n_ort_sin[3];
  double n_z[3] = { 0, 0, 1 };
  /* trick for making it possible to send matrices Z1,2 to subroutine via pointer */
  double Z1m[4][4], Z2m[4][4];
  double *Z1v[4] = { Z1m[0], Z1m[1], Z1m[2], Z1m[3] }, *Z2v[4]  = { Z2m[0], Z2m[1], Z2m[2], Z2m[3] };
  double **Z1 = Z1v, **Z2 = Z2v;
  double cos2alpha1=1.0, sin2alpha1=0.0, cosalpha=0.0;
  double cos2alpha2=1.0, sin2alpha2=0.0, sinalpha=0.0;
  double scr=0.0;
  double *dx1=NULL, *dx2=NULL;

  /* Copy photon phase matrix */
  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
      phase_matrix_inc[i][j]=phase_matrix[i][j];

  if(backward){
    /*exchange incoming and scattered direction for backward */
    dx2=dx_inc;
    dx1=dx_sca;
    scr=phi_target;
    phi_target=phi_source;
    phi_source=scr;
    n_z[2]=-1;
  }
  else {
    dx1=dx_inc;
    dx2=dx_sca;
  }

  if (n_hor[2]!=1.0) { /* do not rotate if n is already 0,0,1 */

    /* Calculate vector orthogonal to incoming reflection plane */
    status = v_cross_product_norm ( dx1, n_hor, n_ort_in );

    /* incoming and reflect direction are not parallel */
    if (!status) {
      /* incoming */
      status = v_cross_product_norm ( dx1, n_z, n_ort_ref );

      if (status) {
	/* Incoming direction || z-axis happens if SZA=0 */
	/* Polarization plane (z, phi0)                  */
	if(!(backward && escape) && !(!backward && scattercounter==0) ){
	  fprintf(stderr, "scattercounter %d \n", scattercounter); 
	  fprintf (stderr, "\n");
	  fprintf(stderr, "*** Warning, polarization plane of photon lost, because \n");
	  fprintf(stderr, "*** photon is scattered exactly in z-direction.\n"); 
	  fprintf(stderr, "*** If this happens too often please contact Robert.\n"); 
	  fprintf (stderr, "\n");
	}

	/* orthogonal vector to (z, phi_esc)-plane*/
	n_ort_ref[0] = -cosd(phi_source);
	n_ort_ref[1] =  sind(phi_source);
	n_ort_ref[2] =  0.0;

      }

      /* Calculate rotation angle 1 */
      v_mult_mu ( n_ort_in, n_ort_ref, &cosalpha );

      /* Calculate sinalpha using vector product */
      v_cross_product ( n_ort_in, n_ort_ref, n_ort_sin );
      v_mult ( n_ort_sin, n_ort_sin, &N );
      sinalpha = - sqrt(N);

      cos2alpha1=cosalpha*cosalpha-sinalpha*sinalpha;
      sin2alpha1=2.0*sinalpha*cosalpha;
      /* The sign of the following elements depends on phi:*/
      if(n_ort_in[2]<0.)
	sin2alpha1=-sin2alpha1;

      phase_matr_rot ( R, cos2alpha1, sin2alpha1, Z1 );
    }
    else
      /* no rotation required, cos2alpha=1 by default */
      Z1 = R;

    /* Calculate vector orthogonal to outgoing reflection plane */
    status = v_cross_product_norm ( dx2, n_hor, n_ort_out );

    /* incoming and reflect direction are not parallel */
    if (!status) {
      /* outgoing */
      status = v_cross_product_norm ( dx2, n_z, n_ort_ref );

      if (status) {
	/* In case of umu=+/- 1 we need to consider this special case: */ 
	/* In forward mode and umu +/-1 the last rotation needs to be in the plane */
	/* defined by (z, phi_sensor). */
	/* In backward mode the first rotation.*/
	if(!(!backward && escape) && !(backward && scattercounter==0) ){
	  fprintf(stderr, "scattercounter %d \n", scattercounter); 
	  fprintf (stderr, "\n");
	  fprintf(stderr, "*** Warning, polarization plane of photon lost, because \n");
	  fprintf(stderr, "*** photon is scattered exactly in z-direction.\n"); 
	  fprintf(stderr, "*** If this happens too often please contact Claudia.\n"); 
	  fprintf (stderr, "\n");
	}

	/* orthogonal vector to (z, phi_esc)-plane*/
	n_ort_ref[0] = -cosd(phi_target);
	n_ort_ref[1] =  sind(phi_target);
	n_ort_ref[2] =  0.0;
      }

      /* Calculate rotation angle 2 */
      v_mult_mu ( n_ort_out, n_ort_ref, &cosalpha );

      /* Calculate sinalpha using vector product */
      v_cross_product ( n_ort_out, n_ort_ref, n_ort_sin );
      v_mult ( n_ort_sin, n_ort_sin, &N );
      sinalpha = + sqrt(N);

      cos2alpha2=cosalpha*cosalpha-sinalpha*sinalpha;
      sin2alpha2=2.0*sinalpha*cosalpha;

      /* The sign of the following elements depends on phi:*/
      if(n_ort_out[2]<0.)
	sin2alpha2=-sin2alpha2;
    
      /* fprintf(stderr, "no %d escape %d angle %g\n", scattercounter, escape, acos(cosalpha)*180/PI); */

      rot_phase_matr ( Z1, cos2alpha2, sin2alpha2, Z2 );
    }
    else
    /* no rotation required, cos2alpha=1 by default */
      Z2 = Z1;
  }
  else
    /* no rotation required, cos2alpha=1 by default */
    Z2 = R;

  /* scattering matrix multiplication */
  /* fprintf(stderr,"matrix sigmas %e %e %e %e \n",cos2alpha1,sin2alpha1,cos2alpha2,sin2alpha2); */
  for (i=0; i<4; i++){
    for (j=0; j<4; j++){
      /* fprintf(stderr," %e(%e)",R[i][j],Z2[i][j]); */
      phase_matrix[i][j]=0.0;
      for (k=0; k<4; k++){
        if (backward)
          phase_matrix[i][j]+=phase_matrix_inc[i][k]*Z2[k][j];
        else
          phase_matrix[i][j]+=Z2[i][k]*phase_matrix_inc[k][j];
      }
      /* XXX RPB: The following line is commented out. I think it is correct
	 this way. However, you could choose to do this in analogy to
	 scattering, then you would have to multiply R[0][0] to the photon
	 weight. Reason for doing so is that after many reflections, the
	 phase_matrix could become extremely small or large, leading to
	 numerical problems. It could maybe also affect VROOM, such that it
	 becomes extremely slow or very imprecise (the Stokes vector is only
	 partly taken into account for splitting). */
      /* phase_matrix[i][j]/=R[0][0]; CHECK!!! */ 
    }
    /* fprintf(stderr,"\n"); */
  }

}


void phase_matr_rot ( double **Z_in, double cos2alpha, double sin2alpha, double **Z_out )
{
  int i=0;
  for (i=0;i<4;i++) {
    Z_out[i][0] =   Z_in[i][0];
    Z_out[i][1] =   Z_in[i][1] * cos2alpha
                  + Z_in[i][2] * sin2alpha;
    Z_out[i][2] = - Z_in[i][1] * sin2alpha
                  + Z_in[i][2] * cos2alpha;
    Z_out[i][3] =   Z_in[i][3];
  }
}



void rot_phase_matr ( double **Z_in, double cos2alpha, double sin2alpha, double **Z_out )
{
  int j=0;
  for (j=0;j<4;j++) {
    Z_out[0][j] =   Z_in[0][j];
    Z_out[1][j] =   Z_in[1][j] * cos2alpha
                  - Z_in[2][j] * sin2alpha;
    Z_out[2][j] = + Z_in[1][j] * sin2alpha
                  + Z_in[2][j] * cos2alpha;
    Z_out[3][j] =   Z_in[3][j];
  }
}


/***********************************************************************************/
/* Function: get_phase_matrix_caoth                                       @62_30i@ */
/* Description:                                                                    */
/*  Obtain the phase matrix for caoth isp.                                         */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Claudia Emde & Robert Buras                                             */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int get_phase_matrix_caoth ( atmosphere_struct *atmos,
				    photon_struct     *p,
				    int                isp,
				    double             mu,
				    int                np,
				    /* Output */
				    double            *phase_matrix )
{
  int status=0;

  switch(atmos->scatter_type[isp]) {
  case MCSCAT_MOL:    /* Rayleigh scattering */
    get_phase_matrix_rayleigh ( mu, atmos->rayleigh_depol, np, phase_matrix );
    break;
  case MCSCAT_AER:   /* Aerosol scattering   */
    status = get_phase_matrix_pft ( &(atmos->phase_aer[p->kc]), mu, 0, np, phase_matrix );
    if (status)
      return fct_err_out (status, "get_phase_matrix_pft", ERROR_POSITION);
    break;
  case MCSCAT_HG1:
    if (np>1)
      /* Henyey-Greenstein phase function can not be applied for polarization*/
      return err_out ("Error %d returned\n The Henyey Greenstein phase function can not be used to calculate polarized radiances with clouds.\n Please provide scattering phase matrices.", -1);
    phase_matrix[0] = HG ( get_g1(atmos, p, isp), mu );
    break;
  case MCSCAT_HG2:
    if (np>1)
      /* Henyey-Greenstein phase function can not be applied for polarization*/
      return err_out ("Error %d returned\n The Henyey Greenstein phase function can not be used to calculate polarized radiances with clouds.\n Please provide scattering phase matrices.", -1);
    phase_matrix[0] = HG2 ( get_g1(atmos, p, isp),
			    get_g2(atmos, p, isp),
			    get_ff(atmos, p, isp),
			    mu );
    break;
  case MCSCAT_PFT:
    status = get_phase_matrix_pft_interpol_reff ( mu,
						  get_reff(atmos, p, isp),
						  atmos->phase[isp],
						  p->SC_mode,
						  np,
						  phase_matrix );
    if (status)
      return fct_err_out (status, "get_phase_matrix_pft_interpol_reff", ERROR_POSITION);
    break;
  default:    
    fprintf(stderr,"Error, no such type %d of scattering!!!\n", atmos->scatter_type[isp]);
    return -1;
  }
  return 0;
}


/***********************************************************************************/
/* Function: calc_stokes_vector                                           @62_30i@ */
/* Description:                                                                    */
/*  Pre-Calculate weight vector for the sensor direction                           */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Claudia Emde (put into subroutine by Robert Buras)                      */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static inline int calc_stokes_vector_for_escape ( sample_struct     *sample,
						  atmosphere_struct *atmos,
						  double            *phase_matrix,
						  double             mu2,
						  direction          dir_inc,
						  double             phi_esc,
						  direction          dir_esc,
						  int                polmat,
						  /* In-/Output */
						  photon_struct     *p )
{
  int ip=0, j=0;
  double phase_matrix_trivial[6]={ 1, 0, 1, 0, 1, 1 };
  double tmp_dir[3];
  
  /* calculate phase matrix for whole atmosphere and multiply it
     with I_0. This works in forward and backward mode but might
     take more CPU time */
  if (polmat){
    phase_matrix_mult( p->phamat, phase_matrix,
		       dir_esc.dx, dir_inc.dx,
                       p->fw_phi0, p->fw_phi, 1,
		       p->scattercounter, sample->backward );
    
    /* turn polarisation in front of MFOV lidar detector */
    if (sample->LidarLocEst) {

      tmp_dir[0]= - sample->lidar[sample->ili].dir.dx[0];
      tmp_dir[1]= - sample->lidar[sample->ili].dir.dx[1];
      tmp_dir[2]= - sample->lidar[sample->ili].dir.dx[2];

      phase_matrix_mult ( p->phamat, phase_matrix_trivial,
			  tmp_dir, dir_esc.dx,
                          p->fw_phi0, p->fw_phi, 1,
			  p->scattercounter, sample->backward );
    }

    /*Z_tot*I_0 (resulting phase matrix times I0, corresponds
      to first column of Z for unpolarized incoming radiation
      (1,0,0,0))*/
    for(ip=0; ip<4; ip++) {
      p->stokes[ip] = 0.0;
      for(j=0; j<4; j++)
	p->stokes[ip] += p->phamat[ip][j]*p->stokes0[j];
    }
  }
  else{
    stokes_vector_sca( p->stokes, phase_matrix, dir_esc.dx, dir_inc.dx, 
                       p->fw_phi0, p->fw_phi);
  }

  return 0;
}



/***********************************************************************************/
/* Function: derive_deltaphi_horz                                         @62_30i@ */
/* Description:                                                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Claudia Emde (put into subroutine by Robert Buras)                      */
/*                                                                        @i62_30@ */
/***********************************************************************************/

double derive_deltaphi_horz ( double *dx,
			      double   phi_esc,
			      int      scattercounter,
			      int      LidarLocEst,
			      double   phi_init )
{
  double p_phi=0.0, phi=0.0;

  /* Calculate azimuth of photon incoming direction */
  if ( scattercounter !=0 || LidarLocEst )
    p_phi = calc_phi_horz(dx, NULL);
  else
    p_phi = phi_init;

  /* Azimuth difference - very important to compute phase matrix */
  phi = p_phi - phi_esc;

  while (phi < 0.)
    phi += 360.;

  return phi;
}


/***********************************************************************************/
/* Function: get_reflection_probability_matrix                                     */
/* Description:                                                                    */
/*   Calculate Stokes weight vector after a surface reflection taking              */
/*   into account the reflection matrix. This weight vector is then used in        */
/*   in the "local estimate" (escape_probability() ).                              */
/*   So far only the BPDF by Tsang/Mishchenko is available, but other reflection   */
/*   matrices may be included in this function.                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Claudia Emde                                                            */
/*                                                                                 */
/***********************************************************************************/

int get_reflection_probability_matrix ( albedo_struct *albedo,
					photon_struct *p,
					int            nstokes,
					int            backward,
					double         mu_inc,
					double         phi_inc,
					double         mu_esc,
					double         phi_esc, 
					int            ia,
					int            ja,
					int            il,
					float          wvnmlo,
					float          wvnmhi,
					int            iv,
					/* Output */
					double      ***refl_mat )
{
  int status=0;
  float wavelength=0.0;

  switch (albedo->method) {

  case MCALB_HAPKE:
  case MCALB_ROSSLI:
  case MCALB_ROSSLI2D:
  case MCALB_COXANDMUNK:
    fprintf (stderr, "Error, albedo->method %d not yet implemented with polarisation!\n",
             albedo->method);
    return -1;
    
  case MCALB_RPV:
  case MCALB_RPV2D_SPECTRAL:
   /* CE: Why is Tsang here included for RPV2D?? */
    ASCII_calloc_double(refl_mat, 4, 4);

    if (il==0 && albedo->method == MCALB_RPV2D_SPECTRAL) {
      /* Wavelength in micrometer */
      wavelength=0.5*1e4*(1.0/wvnmlo + 1.0/wvnmhi);
    
      if (backward)
	status=bpdf_tsang(albedo->u10, wavelength, mu_esc, phi_esc, mu_inc, phi_inc, refl_mat);
      else
	status=bpdf_tsang(albedo->u10, wavelength, mu_inc, phi_inc, mu_esc, phi_esc, refl_mat);
   
      if (status)
	return fct_err_out (status, "bpdf_tsang", ERROR_POSITION);
    }
    else
      (*refl_mat)[0][0] = rpv_brdf (albedo->rpv[il].rho0, albedo->rpv[il].k,
				    albedo->rpv[il].theta,
				    albedo->rpv[il].scale, albedo->rpv[il].sigma,
				    albedo->rpv[il].t1, albedo->rpv[il].t2,
				    mu_inc, mu_esc,
				    fabs((phi_esc-phi_inc)*180.0/PI));
    break;

  case MCALB_LAM:
    ASCII_calloc_double(refl_mat, 4, 4);
    if (iv==-1)
      (*refl_mat)[0][0]=albedo->albedo;
    else
      (*refl_mat)[0][0]=albedo->spectral_albedo[iv];
    break;
    
  case MCALB_LAM2D:
    ASCII_calloc_double(refl_mat, 4, 4);
    (*refl_mat)[0][0]=albedo->albedo2D[ia][ja];
    break; 

  case MCALB_LAM2D_SPECTRAL:
    ASCII_calloc_double(refl_mat, 4, 4);
    (*refl_mat)[0][0]=albedo->alb_type[(int) albedo->rpv_index[ia][ja]];
    break; 
    
  case MCALB_TSANG:
    ASCII_calloc_double(refl_mat, 4, 4);
    
    /* Wavelength in micrometer */
    wavelength=0.5*1e4*(1.0/wvnmlo + 1.0/wvnmhi);
    
    if (backward)
      status=bpdf_tsang(albedo->u10, wavelength, mu_esc, phi_esc, mu_inc, phi_inc, refl_mat);
    else
      status=bpdf_tsang(albedo->u10, wavelength, mu_inc, phi_inc, mu_esc, phi_esc, refl_mat);
   
    if (status)
      return fct_err_out (status, "bpdf_tsang", ERROR_POSITION);
    
    break;
    
  default:
    fprintf (stderr, "Error, albedo->method %d not yet implemented!\n", albedo->method);
    status = -1;
    return NOT_A_NUMBER;
  }

  return 0;
}
  
/***********************************************************************************/
/* Function: reflection_polarized                                                  */
/* Description:                                                                    */
/*   Calculate Stokes weight vector after a surface reflection taking              */
/*   into account the reflection matrix. This weight vector is then used in        */
/*   in the "local estimate" (escape_probability() ).                              */
/*   So far only the BPDF by Tsang/Mishchenko is available, but other reflection   */
/*   matrices may be included in this function.                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Claudia Emde                                                            */
/*                                                                                 */
/***********************************************************************************/

static double reflection_polarized ( albedo_struct *albedo,
				     photon_struct *p,
				     int            nstokes,
				     int            backward,
				     double        *dx_inc,
				     double        *dx_out,
				     double        *n_hor,
				     double         mu_inc,
				     double         phi_inc,
				     double         mu_esc,
				     double         phi_esc, 
				     double         phi_target,
				     int            ia,
				     int            ja,
				     int            il,
				     float          wvnmlo,
				     float          wvnmhi,
				     int            escape,
				     int            spherical3D, 
				     int            spectral_is,
				     int            nlambda,
				     float         *lambda,
				     int           *status )
{
  int iv=0;
  double **refl_mat=NULL;
  
#ifdef MUCHOUT
  int ip;
#endif

  *status=0;

  if (mu_esc <= 0.0)
    return 0.0;
  
  if (!(albedo->reflectalways)){
    fprintf (stderr, "Error, polarisation requires option mc_surface_reflectalways!\n");
    *status=-1; 
    return NOT_A_NUMBER;
  }

  *status = get_reflection_probability_matrix ( albedo, p, nstokes, backward,
						mu_inc, phi_inc, mu_esc, phi_esc,
						ia, ja, il, wvnmlo, wvnmhi, -1,
						&refl_mat );
  if (*status)
    return fct_err_out (NOT_A_NUMBER, "get_reflection_probability_matrix", ERROR_POSITION);
  
  
  /* Only Lambertian spectral albedo is considered for ALIS so far, for */
  /* bpdf_tsang wavelength dependance is only due to changes in the */
  /* refractive index of water which is very small. For the whole */
  /* spectral range the refractive index of the calculation wavelength */
  /* is assumed. Calling get_reflection_probability_matrix in the loop over*/
  /* wavelengths is very expensive.*/
 
  if (spectral_is){
    switch (albedo->method) {
      
      /* For the first cases we do not expect that ALIS uses spectral albedo,  */
      /* 				  so we do not throw a warning.  */
    case MCALB_HAPKE:
    case MCALB_ROSSLI:
    case MCALB_ROSSLI2D:
    case MCALB_COXANDMUNK:
    case MCALB_RPV:
    case MCALB_LAM2D:
    case MCALB_RPV2D_SPECTRAL:
    case MCALB_TSANG:
      break;
    case MCALB_LAM:
      if (albedo->spectral_albedo!=NULL && refl_mat[0][0]!=0.0){ 
	for (iv=0; iv<nlambda; iv++){
	  p->q_albedo_spectral[iv] *=  albedo->spectral_albedo[iv]/refl_mat[0][0];
	}
      }
      break; 
    case MCALB_LAM2D_SPECTRAL:
      if (refl_mat[0][0]!=0.0){ 
	for (iv=0; iv<nlambda; iv++){
	  if(spherical3D){
#ifdef HAVE_SPHER
	    coord_spherical3D (p,
			       albedo->Nx, albedo->X,
			       albedo->Ny, albedo->Y,
			       0, NULL,
			       0, NULL,
			       0,
			       &ia, &ja, NULL, NULL,
			       0);
#else
	    fprintf(stderr,"Error! you are not allowed to use spherical 3D!\n");
	    *status=-1; 
	    return NOT_A_NUMBER;
#endif
	  }
	  else
	    albedo_coord (p, albedo, &ia, &ja);	
	  
	  p->q_albedo_spectral[iv] *=  
	    albedo->spectral_alb_type[(int) albedo->rpv_index[ia][ja]][iv] / refl_mat[0][0];
	}
      }
      break; 
    default:
      fprintf (stderr, "Error, albedo->method %d not yet implemented!\n", albedo->method);
      *status = -1;
      return NOT_A_NUMBER;
    }
  }
  
  reflect_matrix_rotate ( p->phamat, refl_mat,
			  dx_out, dx_inc, n_hor,
			  p->phi0, phi_target,
			  escape, p->scattercounter, backward,
			  NULL, NULL );

  /* scattering matrix multiplication */
  mat_v_mult ( p->phamat, p->stokes0, p->stokes );

#ifdef MUCHOUT
  if(p->muchoutcounter==MUCHOUTPHOTON) {
    for (ip=0; ip<nstokes; ip++)
      fprintf(stderr,"refl_mat %e %e %e %e\n",refl_mat[ip][0],refl_mat[ip][1],refl_mat[ip][2],refl_mat[ip][3]);
    for (ip=0; ip<nstokes; ip++)
      fprintf(stderr,"phamat %e %e %e %e\n",p->phamat[ip][0],p->phamat[ip][1],p->phamat[ip][2],p->phamat[ip][3]);
    fprintf(stderr,"stokes %e %e %e %e\n",p->stokes[0],p->stokes[1],p->stokes[2],p->stokes[3]);
  }
#endif
  
  ASCII_free_double(refl_mat, 4);
  
  /* Isotropic scattering, BPDF or Lambertian albedo included in Stokes weight vector */
  return 1.0;
  }

/***********************************************************************************/
/* Function: correct_escape_dir_refract                                   @62_30i@ */
/* Description:                                                                    */
/*  Correct the photon direction for the calculation of the escape probability     */
/*  the point of scattering towards the sun.                                       */ 
/*  Note: This is the first not yet fully verified version, which might be changed */
/*        later on.                                                                */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Claudia Emde                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int correct_escape_dir_refract(double *dir_sensor, 
				      photon_struct *photon, 
                                      atmosphere_struct *atmos, 
                                      sample_struct *sample, 
                                      float *refind)
{
  
  double step=0.0;
  int id=0, i=0, status=0; 
  int counter=0; 

  photon_struct *p = calloc_photon ( sample,
				     0,
				     MCSRC_NONE,
				     atmos->nlambda_abs,
                                     atmos->Nc,
				     atmos->n_caoth );
  double diff[3], N=0.0;
  double diff_abs=999., diff_abs_old=999., a_old=0.0, a_new=0.0, a=0.0; 
  /* double dir_scatpoint[3];  */

  /* double dir_sza;  */
  double e_r[3];
  double epsilon=1e-6;
  
  /* work on a copy of the photon, not the original one */
  cp_photon_struct (p, photon, sample, atmos->n_caoth);
  
  /* the direction to be traced is the radiance direction */
  cp_direction (&(p->dir), &(sample->rad[id].dir));
    
  /* initialize status; if unchanged at end of travel step, the photon crosses a grid boundary */
  p->photon_status = MCSTATUS_TRAVEL;

  /* fprintf(stderr, "Correct direction for local estimate ... \n"); */
  /* dir_start=p->dir.dx[2]; */
  /* dir_sza=p->dir.dx[2]; */

  e_r[0]=p->x[0]-atmos->xmax/2.;
  e_r[1]=p->x[1]-atmos->ymax/2.;
  e_r[2]=p->x[2]+atmos->r_earth;
  
  N=0; 
  for (i=0; i<3; i++){
    N+=e_r[i]*e_r[i];
  }
  for (i=0; i<3; i++)
    e_r[i]/=sqrt(N);

  /* fprintf(stderr, "e_r %g %g %g \n", e_r[0],e_r[1],e_r[2]); */
  /* how accurate do we need the escape radiance angle? Is 1e-5 required?*/
  
  while ( diff_abs> 1e-4) {
    
    for(i=0;i<3;i++)
      p->dir.dx[i]+=a*e_r[i];
    
    /* Normalize direction vector */
    N=sqrt(p->dir.dx[0]*p->dir.dx[0] + p->dir.dx[1]*p->dir.dx[1] + p->dir.dx[2]*p->dir.dx[2]);
    for (i=0; i<3; i++)
      p->dir.dx[i]/=N; 
    
    /* fprintf(stderr, "at scattering point p->dir.dx %g  %g  %g \n", p->dir.dx[0], p->dir.dx[1], p->dir.dx[2]); */
    /* for (i=0; i<3; i++) */
    /*   dir_scatpoint[i]= p->dir.dx[i]; */
    
    while (1==1){
            
      /*******************************************************************/
      /* 1. Calculate step length to next boundary                       */
      /*******************************************************************/
      status = intersection1D_spherical (p, atmos, sample->refraction, 
                                         refind, &step);
      /* Photon is absorbed by the surface.*/
      if (status<0)
        break;
      
      if (step<0)
        step=0;
      
      /*******************************************************************/
      /* 2. apply step length, move photon, and exit if photon path ends */
      /*******************************************************************/
      
      /* sum up absorption and scattering optical depths */
      /* step is the escape distance within each box between two box boundaries (or between initital photon position and box boundary) depending on the escape direction */
      p->tauabs.tot += step * 
        ( get_ksca (atmos, p, MCCAOTH_TOT) + get_kabs_tot (atmos, p) );
      
      p->pathlength += step;
      
      status = step1D (p, atmos, step, sample->bcond, 0, 0);
      
      /* end photon travel unless it simply crosses the boundary */
      if (p->photon_status != MCSTATUS_TRAVEL )
        break;
      
      /**********************************/
      /* 3. check if photon leaves grid */
      /**********************************/
      
      /* finally check if photon arrived at the lower surface */
      /* or at the top of the atmosphere                      */
      
      if (p->kc<0) {
        /* fprintf(stderr, "at BOA p->dir.dx[2] %g \n", p->dir.dx[2]); */
        p->photon_status = MCSTATUS_BOUNDARY_LOWER;
        break;
      }
      
      if (p->kc>atmos->Nz-1) {
        /* fprintf(stderr, "at TOA p->dir.dx[2] %g \n", p->dir.dx[2]); */
        p->photon_status = MCSTATUS_BOUNDARY_UPPER;
        break;
      }
      
    }
    diff_abs_old=diff_abs;
    
    N=0;
    for (i=0; i<3; i++){
      diff[i]=p->dir.dx[i]-sample->rad[id].dir.dx[i]; 
      N+=diff[i]*diff[i]; 
    }
    diff_abs=sqrt(N);
     
    /* Initialize new photon */
    cp_photon_struct (p, photon, sample, atmos->n_caoth);
    cp_direction (&(p->dir), &(sample->rad[id].dir));
    p->photon_status = MCSTATUS_TRAVEL;
    
    counter+=1; 
    
    /*   fprintf(stderr, "iteration %d diff %g %g %g diffabs %g \n", counter, diff[0], diff[1], diff[2], diff_abs); */
    /*   fprintf(stderr, "diff_abs_old %g a_old %g diff_abs %g, a_new %g \n",                                       */
    /*        diff_abs_old, a_old, diff_abs, a);                                                                    */
    
    /* Factor a to be found in iteration */
    /* First guess */
    if (counter==1){
      a=0.01;
      a_new=0.01; 
      a_old=0.0;
    }
    else{
      a_new=(diff_abs*a_old+diff_abs_old*a)/(diff_abs_old + diff_abs);
      a_old=a;
      a=a_new;
      /* a+=0.0001; */
    }
    
    if (fabs(diff_abs-diff_abs_old)<epsilon){
      /* fprintf(stderr, "Warning, refraction accurracy %g %g \n", diff_abs_old, diff_abs); */
      break;
    }

    /*    fprintf(stderr, "a_new %g \n", a);   */
    
    
  }
  
  for (i=0; i<3; i++)
    dir_sensor[i]=p->dir.dx[i];
  
  /* free memory */
  destroy_photon (p, atmos->n_caoth);
  
  return 0;
}


/***********************************************************************************/
/* Function: set_panorama_dir                                             @62_30i@ */
/* Description:                                                                    */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Bernhard Reinhard, Robert Buras                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int set_panorama_dir ( sample_struct *sample, float *sza, float *phi0  )
{
  double dtheta=0.0, theta_min=0.0, theta_max=0.0;
  double phi_min=0.0, phi_max=0.0, dphi=0.0;
  double upper_theta=0.0, lower_theta=0.0;
  double addhalf=0.0;
  double sinalpha=0.0, mu2=0.0;
  double xc[3]={0.,0.,0.};

  if (sample->pan_distr_photons_over_pixel && sample->pan_no_pixel){
    fprintf(stderr,"ERROR: pan_no_pixel and pan_distr_photons_over_pixel may not be used together!\n");
    exit(1);
  }    

  /* assuming you want to simulate a CCD-camera                            */
  /* Either by the simple method of sending all photons from the center of */
  /* the individual CCD pixels, or by the more realistic option of         */
  /* distributing the photons over the finite solid angle covered by the   */
  /* individual pixels                                                     */

  /* if pan_no pixel set, distribute angles in such manner that       */
  /* theta_min/theta_max; phi_min/phi_max as provided by user input file   */
  /* are directly faced */

  theta_min = sample->pan_theta_min;
  theta_max = sample->pan_theta_max;
  phi_min   = sample->pan_phi_min;
  phi_max   = sample->pan_phi_max;

  if (sample->pan_umu_min < -2.0) {
    /* "cartesian" coordinates (fish-eye) */
    /* names need cleaning up !!! BCA */
    sinalpha = sind(phi_max);
    xc[0] = sinalpha * ( 2.0 * ( (double) sample->backward_is + 0.5 )
			 / (double) sample->Nx - 1.0 );
    xc[1] = sinalpha * ( 2.0 * ( (double) sample->backward_js + 0.5 )
			 / (double) sample->Ny - 1.0 );
    xc[2] = 0.0;

    mu2 = xc[0] * xc[0] + xc[1] * xc[1];
    if (mu2>sinalpha*sinalpha)
      /* this angle does not exist, or should not be calculated */
      return 1;
    *sza = asind ( sqrt ( mu2 ) );
    *phi0 = calc_phi_horz ( xc, NULL ); /* check if +-90/180; or toggle 0][1 */

    /* error messages for modes not working BCA */

    return 0;
  }

  /* define dtheta */
  if (sample->Ny - sample->pan_no_pixel > 0) {
    dtheta = ( theta_max - theta_min )
      / (double) ( sample->Ny - sample->pan_no_pixel );
  }
  else {
    dtheta = 0.0;
    if (theta_max != theta_min) {
      fprintf(stderr, "%s %s %s %s",
	      "Panorama_no_pixel: Theta_max != Theta_min but setup only",
	      "calculation at one specific theta. Will do the calculation",
	      "at theta_min. Maybe you want to increase the numbers at",
	      "mc_sample_grid?!\n");
      return -1;
    }
  }

  /* define dphi */
  if (sample->Nx - sample->pan_no_pixel > 0)
    dphi = ( sample->pan_phi_max - sample->pan_phi_min )
      / (double) ( sample->Nx - sample->pan_no_pixel );
  else {
    dphi = 0.;
    if (phi_max != phi_min){
      fprintf(stderr, "%s %s %s %s",
	      "Panorama_no_pixel: Phi_max != Phi_min but setup only",
	      "calculation at one specific phi. Will do the calculation",
	      "at phi_min. Maybe you want to increase the numbers at",
	      "mc_sample_grid?!\n");
      return -1;
    }
  }

  if (sample->pan_distr_photons_over_pixel){
    lower_theta = ( sample->pan_theta_min
		    + dtheta * (double) sample->backward_js );
    upper_theta = lower_theta + dtheta;
    sample->pan_FOV = - dphi / 180.0 * PI
      * ( cosd ( upper_theta ) - cosd ( lower_theta ) );
    addhalf=0.0;
  }
  else
    addhalf=0.5;

  /* calculate sza = minimum + delta times number of pixel*/
  *sza = theta_min + dtheta * ( (double) sample->backward_js + addhalf );
  /* calculate phi0 */
  *phi0 = phi_min + dphi * ( (double) sample->backward_is + addhalf );

  sample->pan_dphi = dphi;
  sample->pan_dtheta = dtheta;

  return 0;
}  


/***********************************************************************************/
/* Function: std_noisy                                                    @62_30i@ */
/* Description:                                                                    */
/*   Calculate standard deviation, taking into account that the input might suffer */
/*   from computational noise and set to 0 instead of sqrt(-tiny number) = NaN     */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static inline double std_noisy (double var, double avg, double dmcphotons)
{
  if (fabs(var/avg/avg-1.0) >= MC_EPSILON*MC_EPSILON)
    return sqrt ((var-avg*avg) / dmcphotons);
  else 
    return 0.0;
}

/***********************************************************************************/
/* Function: Fresnel                                                      @62_30i@ */
/* Description:                                                                    */
/*   Stolen from Cox & Munk                                                        */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

double Fresnel ( double nr,
		 double ni,
		 double coschi )
{
  /*
    C to compute the Fresnel's coefficient of reflection (see for
    C example M. Born and E. Wolf, Principles of Optics, Pergamon Press, 
    C fifth edition, 1975, pp 628; seventh edition, 1999, pp 753
    C input parameters: nr=index of refraction of the sea water
    C                   ni=extinction coefficient of the sea water
    C                   coschi & sinchi=cosine and sine of the incident radiation 
    C                                   with respect of the wave facet normal.
    C output parameter: R1=Fresnel's coefficient for reflection
  */

  double nr2=0.0, ni2=0.0, a1=0.0, a2=0.0, u=0.0, v=0.0;
  double Rr2=0.0, b1=0.0, b2=0.0, Rl2=0.0;
  double sinchi = sqrt ( 1.0 - coschi * coschi );

  /* absolute value for a1 to get v=0 when ni=0 */
  nr2 = nr * nr;
  ni2 = ni * ni;

  a1= nr2 - ni2 - sinchi * sinchi;
  a1 = a1 > 0.0 ? a1 : -a1;
  a2 = sqrt ( a1 * a1 + 4.0 * nr2 * ni2 );

  /* bm:
     added the abs because a1 quite often equals a2,
     except for very small rounding errors that cause
     nonsense because the sqrt of a negative number is NAN */
  u = a1 + a2;
  u = u > 0.0 ? u : -u;
  u = sqrt ( 0.5 * u );
  v = -a1 + a2;
  v = v > 0.0 ? v : -v;
  v = sqrt ( 0.5 * v );
  Rr2 = ( ( coschi - u ) * ( coschi - u ) + v * v )
      / ( ( coschi + u ) * ( coschi + u ) + v * v );

  b1 = ( nr2 - ni2 ) * coschi;
  b2 = 2 * nr * ni * coschi;
  /* bm:
     changed the following line which was buggy in the original 
     6S code */
  Rl2 = ( ( b1 - u ) * ( b1 - u ) + ( b2 - v ) * ( b2 - v ) )
      / ( ( b1 + u ) * ( b1 + u ) + ( b2 + v ) * ( b2 + v ) );

  return  ( Rr2 + Rl2 ) / 2.;
}

void index_water ( double wl,
		   double xsal,
		   double *nr,
		   double *ni )
{
  /*
    C input parameters:  wl=wavelength (in micrometers)
    C                    xsal=salinity (in ppt), if xsal<0 then 34.3ppt by default
    C output parameters: nr=index of refraction of sea water
    C                    ni=extinction coefficient of sea water
  */

  int i=0;
  double xwl=0.0, yr=0.0, yi=0.0;
  double nrc=0.0, nic=0.0;

  /*
    C Indices of refraction for pure water from Hale and Querry, 
    C Applied Optique, March 1973, Vol. 12,  No. 3, pp. 555-563
  */
  double twl[] = {
    0.250,0.275,0.300,0.325,0.345,0.375,0.400,0.425,0.445,0.475,
    0.500,0.525,0.550,0.575,0.600,0.625,0.650,0.675,0.700,0.725,
    0.750,0.775,0.800,0.825,0.850,0.875,0.900,0.925,0.950,0.975,
    1.000,1.200,1.400,1.600,1.800,2.000,2.200,2.400,2.600,2.650,
    2.700,2.750,2.800,2.850,2.900,2.950,3.000,3.050,3.100,3.150,
    3.200,3.250,3.300,3.350,3.400,3.450,3.500,3.600,3.700,3.800,
    3.900,4.000 };
  /* double tnr[] = { */
  /*   1.362,1.354,1.349,1.346,1.343,1.341,1.339,1.338,1.337,1.336, */
  /*   1.000,1.100,1.200,1.300,1.400,1.500,1.600,1.700,1.800,1.900, */
  /*   1.330,1.330,1.329,1.329,1.329,1.328,1.328,1.328,1.327,1.327, */
  /*   1.327,1.324,1.321,1.317,1.312,1.306,1.296,1.279,1.242,1.219, */
  /*   1.188,1.157,1.142,1.149,1.201,1.292,1.371,1.426,1.467,1.483, */
  /*   1.478,1.467,1.450,1.432,1.420,1.410,1.400,1.385,1.374,1.364, */
  /*   1.357,1.351 }; */
  /* double tni[] = { */
  /*   3.35E-08,2.35E-08,1.60E-08,1.08E-08,6.50E-09, */
  /*   3.50E-09,1.86E-09,1.30E-09,1.02E-09,9.35E-10, */
  /*   1.00E-09,1.32E-09,1.96E-09,3.60E-09,1.09E-08, */
  /*   1.39E-08,1.64E-08,2.23E-08,3.35E-08,9.15E-08, */
  /*   1.56E-05,1.48E-03,1.25E-01,1.82E-07,2.93E-07, */
  /*   3.91E-07,4.86E-07,1.06E-06,2.93E-06,3.48E-06, */
  /*   2.89E-06,9.89E-06,1.38E-04,8.55E-05,1.15E-04, */
  /*   1.10E-03,2.89E-04,9.56E-04,3.17E-03,6.70E-03, */
  /*   1.90E-02,5.90E-02,1.15E-01,1.85E-01,2.68E-01, */
  /*   2.98E-01,2.72E-01,2.40E-01,1.92E-01,1.35E-01, */
  /*   9.24E-02,6.10E-02,3.68E-02,2.61E-02,1.95E-02, */
  /*   1.32E-02,9.40E-03,5.15E-03,3.60E-03,3.40E-03, */
  /*   3.80E-03,4.60E-03 }; */
  /* double tnr[] = { */
  /*   1.362,1.354,1.349,1.346,1.343,1.341,1.339,1.338,1.337,1.336, */
  /*   1.335,1.334,1.333,1.333,1.332,1.332,1.331,1.331,1.331,1.330, */
  /*   1.330,1.330,1.329,1.329,1.329,1.328,1.328,1.328,1.327,1.327, */
  /*   1.327,1.324,1.321,1.317,1.312,1.306,1.296,1.279,1.242,1.219, */
  /*   1.188,1.157,1.142,1.149,1.201,1.292,1.371,1.426,1.467,1.483, */
  /*   1.478,1.467,1.450,1.432,1.420,1.410,1.400,1.385,1.374,1.364, */
  /*   1.357,1.351 }; */
  /* double tni[] = { */
  /*   3.35E-08,2.35E-08,1.60E-08,1.08E-08,6.50E-09, */
  /*   3.50E-09,1.86E-09,1.30E-09,1.02E-09,9.35E-10, */
  /*   1.00E-15,1.00E-14,1.00E-13,1.00E-12,1.00E-11, */
  /*   1.00E-10,1.00E-09,1.00E-08,1.00E-07,1.00E-06, */
  /*   1.56E-07,1.48E-07,1.25E-07,1.82E-07,2.93E-07, */
  /*   3.91E-07,4.86E-07,1.06E-06,2.93E-06,3.48E-06, */
  /*   2.89E-06,9.89E-06,1.38E-04,8.55E-05,1.15E-04, */
  /*   1.10E-03,2.89E-04,9.56E-04,3.17E-03,6.70E-03, */
  /*   1.90E-02,5.90E-02,1.15E-01,1.85E-01,2.68E-01, */
  /*   2.98E-01,2.72E-01,2.40E-01,1.92E-01,1.35E-01, */
  /*   9.24E-02,6.10E-02,3.68E-02,2.61E-02,1.95E-02, */
  /*   1.32E-02,9.40E-03,5.15E-03,3.60E-03,3.40E-03, */
  /*   3.80E-03,4.60E-03 }; */
  double tnr[] = {
    1.362,1.354,1.349,1.346,1.343,1.341,1.339,1.338,1.337,1.336,
    1.335,1.334,1.333,1.333,1.332,1.332,1.331,1.331,1.331,1.330,
    1.330,1.330,1.329,1.329,1.329,1.328,1.328,1.328,1.327,1.327,
    1.327,1.324,1.321,1.317,1.312,1.306,1.296,1.279,1.242,1.219,
    1.188,1.157,1.142,1.149,1.201,1.292,1.371,1.426,1.467,1.483,
    1.478,1.467,1.450,1.432,1.420,1.410,1.400,1.385,1.374,1.364,
    1.357,1.351 };
  double tni[] = {
    3.35E-08,2.35E-08,1.60E-08,1.08E-08,6.50E-09,
    3.50E-09,1.86E-09,1.30E-09,1.02E-09,9.35E-10,
    1.00E-09,1.32E-09,1.96E-09,3.60E-09,1.09E-08,
    1.39E-08,1.64E-08,2.23E-08,3.35E-08,9.15E-08,
    1.56E-07,1.48E-07,1.25E-07,1.82E-07,2.93E-07,
    3.91E-07,4.86E-07,1.06E-06,2.93E-06,3.48E-06,
    2.89E-06,9.89E-06,1.38E-04,8.55E-05,1.15E-04,
    1.10E-03,2.89E-04,9.56E-04,3.17E-03,6.70E-03,
    1.90E-02,5.90E-02,1.15E-01,1.85E-01,2.68E-01,
    2.98E-01,2.72E-01,2.40E-01,1.92E-01,1.35E-01,
    9.24E-02,6.10E-02,3.68E-02,2.61E-02,1.95E-02,
    1.32E-02,9.40E-03,5.15E-03,3.60E-03,3.40E-03,
    3.80E-03,4.60E-03 };

  /* locate */
  i=1;
  while ( wl >= twl[i] && i<61 )
    i++;

  xwl = twl[i] - twl[i-1];
  yr  = tnr[i] - tnr[i-1];
  yi  = tni[i] - tni[i-1];

  *nr = tnr[i-1] + ( wl - twl[i-1] ) * yr / xwl;
  *ni = tni[i-1] + ( wl - twl[i-1] ) * yi / xwl;

  /*
    c Correction to be applied to the index of refraction and to the extinction 
    c coefficients of the pure water to obtain the ocean water one (see for 
    c example Friedman). By default, a typical sea water is assumed 
    c (Salinity=34.3ppt, Chlorinity=19ppt) as reported by Sverdrup. 
    c In that case there is no correction for the extinction coefficient between 
    c 0.25 and 4 microns. For the index of refraction, a correction of +0.006 
    c has to be applied (McLellan). For a chlorinity of 19.0ppt the correction 
    c is a linear function of the salt concentration. Then, in 6S users are able 
    c to enter the salt concentration (in ppt).
    c REFERENCES:
    c Friedman D., Applied Optics, 1969, Vol.8, No.10, pp.2073-2078.
    c McLellan H.J., Elements of physical Oceanography, Pergamon Press, Inc.,
    c        New-York, 1965, p 129.
    c Sverdrup H.V. et al., The Oceans (Prentice-Hall, Inc., Englewood Cliffs,
    c        N.J., 1942, p 173.
  */

  nrc = 0.006;
  nic = 0.000;
  *nr = *nr + nrc * ( xsal / 34.3 );
  *ni = *ni + nic * ( xsal / 34.3 );

}


/***********************************************************************************/
/* Function: find_tau_ext                                                 @62_30i@ */
/* Description:                                                                    */
/* Find extinction between two given points                                        */
/* Code mostly taken from escape_probability()                                     */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Christian Pause                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int find_tau_ext ( photon_struct     *photon, 
		   atmosphere_struct *atmos,
		   sample_struct     *sample,
		   elevation_struct  *elev,
		   int                ic,
		   int                jc,
		   int                kc,
		   float             *refind,
		   double            *dx,
		   double            *start,
		   double             path,
		   double            *tau_save )
{
  double step=0.0, gamma=0.0;
  int status=0;
  int is_threed=0;
  int j=0;
  /* double pt[3]; */

#ifdef HAVE_SPHER
  int first=1;
#endif

  /* prepare things */
  photon_struct *p = calloc_photon ( sample,
				     atmos->Nz,
				     MCSRC_NONE,
				     atmos->nlambda_abs,
                                     atmos->Nc,
				     atmos->n_caoth );

  
  /* work on a copy of the photon, not the original one */
  cp_photon_struct (p, photon, sample, atmos->n_caoth);

  /* set direction and starting point */   
  for (j=0;j<3;j++) {
    p->dir.dx[j] = dx[j];
    p->x[j] = start[j];
  }

  p->ic = ic;
  p->jc = jc;
  p->kc = kc;
  p->pathlength = 0.;
  p->tauext_tot = 0.;
  p->tauabs.tot = 0.;

  /* adjust p->dir.hit if neccessary */
  for (j=0;j<3;j++) { 
    if (p->dir.dx[j] > 0.0)
      p->dir.hit[j] = 1;
    else if (p->dir.dx[j] < 0.0)
      p->dir.hit[j] = 0;
    else
      p->dir.hit[j] = -1;
  }    

  /* initialize status; if unchanged at end of travel step, the photon crosses a grid boundary */
  p->photon_status = MCSTATUS_TRAVEL;
  p->RIS_mode = MCRIS_MODE_NORMAL;
  p->VIS_mode = MCVIS_MODE_NORMAL;
  
  while (1==1) {  /* until photon reaches the surface or TOA */
    
    is_threed = (atmos->threed[MCCAOTH_TOT][p->kc]>=1); /* 3D layer */

    /**************************************************************/
    /* 1. find out what happens next with photon, and step length */
    /**************************************************************/

    /* calculate next intersection point with the box boundaries   */
    if (sample->spherical3D) {
#ifdef HAVE_SPHER
      status = intersection3D_spherical (p, atmos,
					 sample->spherical3D_scene,
					 &first, &step);
      if (status)
	return err_out("Error %d returned by intersection3D_spherical\n", status);
#else
      fprintf(stderr,"Error! you are not allowed to use spherical 3D!\n");
      return -1;
#endif
    }
    else {
      if (is_threed==1) { /* 3D layer */
#if HAVE_MYSTIC3D
        status = intersection3D (p, atmos, 10000, 0, &step);
	if (status)
	  return err_out("Error %d returned by intersection3D\n", status);
#else
	fprintf(stderr,"Error! you are not allowed to use mystic 3D!\n");
	return -1;
#endif
      }
      else {
        if (sample->spherical){
          status = intersection1D_spherical (p, atmos,
      				       sample->refraction, refind,  &step);
	  /* Photon is absorbed by the surface.*/
	  if (status<0)
	    break;
	}
	else {
	  intersection1D (p, atmos, 10000, 0, &step, 0);
	}
      }
    }

    if (step<0)
      step=0;
    
    /* check if photon crossed the 2D surface */
    if (elev->elev2D)  {
      
      status = cross_surface (p, step,
                              elev,
                              sample->bcond,
                              &gamma);
      
      
      /* if cross_surface ran into some inconsistency */
      if (status<0) {
        fprintf (stderr, "Something terrible has happened! cross_surface() returned\n");
        fprintf (stderr, "status %d\n", status);
        return -1;
      }
      
      if (status>0) {
        step = gamma;
        p->photon_status = MCSTATUS_SURFACE;
      }
    }

    if (p->pathlength + step >= path) {
      step = path - p->pathlength;
      p->photon_status = MCSTATUS_DEFAULT;
    }

    /*******************************************************************/
    /* 2. apply step length, move photon, and exit if photon path ends */
    /*******************************************************************/

    /* sum up absorption optical depth */
    /*tau_ext += step * (get_kabs_tot (atmos, p) + get_ksca (atmos, p, MCCAOTH_TOT));*/
    p->pathlength += step;
    mc_add_optical_depth(atmos, p, step, sample->coherent_backscatter);
 
    /* for (j=0; j<3; j++) */
    /*   pt[j] = p->x[j] + step * p->dir.dx[0]; */
 
    /* move photon to box boundary, to next scattering point,  */
    /* or to the 2D surface                                    */
      
    if (is_threed==1 || sample->spherical3D){ /* 3D layer */
#if HAVE_MYSTIC3D
      status = step3D (p, atmos, step, sample->bcond, 0, sample->spherical3D, 0);
#else
      fprintf(stderr,"Error! you are not allowed to use mystic 3D!\n");
      return -1;
#endif
    }
    else
      status = step1D (p, atmos, step, sample->bcond, 0, 0);
  
    if (status<0) {
      fprintf (stderr, "Error %d returned by step%dD()\n", status,
               (is_threed==1?3:1));
      return status;
    }
      
    /* end photon travel unless it simply crosses the boundary */
    if (p->photon_status != MCSTATUS_TRAVEL )
      break;

    /**********************************/
    /* 3. check if photon leaves grid */
    /**********************************/

    /* finally check if photon arrived at the lower surface */
    /* or at the top of the atmosphere                      */
      
    if (p->kc<0) {
      p->photon_status = MCSTATUS_BOUNDARY_LOWER;
        
      if (elev->elev2D) {
        fprintf (stderr, "FATAL error! Escape photon reached the 1D boundary in a\n");
        fprintf (stderr, "             2D elevation case (%dD layer)!\n\n",
                 (is_threed==1?3:1));
        return -1;
      }
        
      break;
    }
      
    if (p->kc>atmos->Nz-1) {
      p->photon_status = MCSTATUS_BOUNDARY_UPPER;
      break;
    }
  } /* end step loop */
   
  tau_save[0] = p->tauext_tot;
  tau_save[1] = p->tauabs.tot;
 
  /* free memory */
  destroy_photon (p, atmos->n_caoth);

  return 1;
}

/***********************************************************************************/
/* Function: find_ris_factor                                              @62_30i@ */
/* Description:                                                                    */
/* Find out the ris-factor needed for a scattering-optical depth of 1 (or choice)  */
/* Code partly taken from escape_probability()                                     */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Christian Pause                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int find_ris_factor ( photon_struct *photon, 
		      atmosphere_struct *atmos,
		      sample_struct *sample,
		      elevation_struct *elev,
		      float *refind )
{
  double step=0.0, gamma=0.0;
  int status=0;
  int is_threed=0;
  double tau_sca=0.0;
  double opan=0.0;
  double radi=0.0;
  double tsca[6];
  double tautemp[5];
  double taumin=0.0;
  double taumax=0.0;
  double theta=0.0;
  double phi=0.0;
  int i=0;
  int j=0;
  int cloudy=0;
  double psqrt=0;
  int numberOfShots = 0;
  double aadd[6][3] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
  double arad[6][3] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}};

#ifdef HAVE_SPHER
  int first=1;
#endif

  /* prepare things */
  photon_struct *p = calloc_photon ( sample,
				     atmos->Nz,
				     MCSRC_NONE,
				     atmos->nlambda_abs,
                                     atmos->Nc,
				     atmos->n_caoth );

  if (sample->LidarLocEst) {
  /* send out six photons */
  /* first one straight with only molecular scattering to determine magnitude of mol sca */
  /* second one straight with full scattering */
  /* third to sixth from four extremal points on the detector along the maximum opening angles */
  /* if any of these five photons sees a cloud, a ris-factor is calculated */
    numberOfShots = 6;
    theta = sample->laser[sample->ili].theta;
    phi = sample->laser[sample->ili].phi;
    opan = sin(sample->laser[sample->ili].alpha);
    radi = sample->laser[sample->ili].radius;
  /* copy an array blockwise */

/* The following constructs do not seem to be compatible with every compiler:
    memcpy(aadd,(double[6][3]){{0,0,0},
                               {0,0,0},
                               {opan*cosd(phi),opan*sind(phi),0},
                               {-opan*cosd(phi),-opan*sind(phi),0},
                               {opan*sind(phi)*cosd(theta),opan*cosd(phi)*cosd(theta),opan*sind(theta)},
                               {-opan*sind(phi)*cosd(theta),-opan*cosd(phi)*cosd(theta),-opan*sind(theta)}}, sizeof aadd);
    
    memcpy(arad,(double[6][3]){{0,0,0},
                               {0,0,0},
                               {radi*cosd(phi),radi*sind(phi),0},
                               {-radi*cosd(phi),-radi*sind(phi),0},
                               {radi*sind(phi)*cosd(theta),radi*cosd(phi)*cosd(theta),radi*sind(theta)},
                               {-radi*sind(phi)*cosd(theta),-radi*cosd(phi)*cosd(theta),-radi*sind(theta)}}, sizeof arad);
    Doing it by hand instead: */

    aadd[2][0] = opan*cosd(phi);
    aadd[2][1] = opan*sind(phi);
    aadd[3][0] = -opan*cosd(phi);
    aadd[3][1] = -opan*sind(phi);
    aadd[4][0] = opan*sind(phi)*cosd(theta);
    aadd[4][1] = opan*cosd(phi)*cosd(theta);
    aadd[4][2] = opan*sind(theta);
    aadd[5][0] = -opan*sind(phi)*cosd(theta);
    aadd[5][1] = -opan*cosd(phi)*cosd(theta);
    aadd[5][2] = -opan*sind(theta);

    arad[2][0] = radi*cosd(phi);
    arad[2][1] = radi*sind(phi);
    arad[3][0] = -radi*cosd(phi);
    arad[3][1] = -radi*sind(phi);
    arad[4][0] = radi*sind(phi)*cosd(theta);
    arad[4][1] = radi*cosd(phi)*cosd(theta);
    arad[4][2] = radi*sind(theta);
    arad[5][0] = -radi*sind(phi)*cosd(theta);
    arad[5][1] = -radi*cosd(phi)*cosd(theta);
    arad[5][2] = -radi*sind(theta);
  }
  else {
  /* Standard MYSTIC, only one path is needed */
  /* we have to loop two times because the first loop only considers */
  /* molecular scattering, its result is disregarded */
    numberOfShots = 2;
  }

  for (i=0; i<numberOfShots; i++) {

    /* reset tau_sca for each run */
    tau_sca = 0.0;

    /* work on a copy of the photon, not the original one */
    cp_photon_struct (p, photon, sample, atmos->n_caoth);

    /* set direction and starting point on the detector */   
    if (sample->LidarLocEst) {
      cp_direction (&(p->dir), &(sample->laser[sample->ili].dir));
      for (j=0;j<3;j++) { 
        p->dir.dx[j] = p->dir.dx[j]+aadd[i][j];
        p->x[j]=sample->laser[sample->ili].x[j]+arad[i][j];
      }
     
  
      /* normalize new direction and adjust p->dir.hit if neccessary */
      psqrt = sqrt(p->dir.dx[0] * p->dir.dx[0] + p->dir.dx[1] * p->dir.dx[1] + p->dir.dx[2] * p->dir.dx[2]);
      for (j=0;j<3;j++) { 
        p->dir.dx[j] = p->dir.dx[j]/psqrt;
        if (p->dir.dx[j] > 0.0)
          p->dir.hit[j] = 1;
        else if (p->dir.dx[j] < 0.0)
          p->dir.hit[j] = 0;
        else
          p->dir.hit[j] = -1;
      }    
  
      /* check if photon has entered another box after offset on the detector area and adjust if neccessary */
      for (j=0;j<atmos->Nx;j++) {
        if (atmos->X[j] >= p->x[0]) {
          p->ic = j-1;
          break;
        }
      }
      for (j=0;j<atmos->Ny;j++) {
        if (atmos->Y[j] >= p->x[1]) {
          p->jc = j-1;
          break;
        }
      }
      for (j=0;j<atmos->Nz;j++) {
        if (atmos->Z[j] >= p->x[2]) {
          p->kc = j-1;
          break;
        }
      }
    }

    /* initialize status; if unchanged at end of travel step, the photon crosses a grid boundary */
    p->photon_status = MCSTATUS_TRAVEL;
    p->RIS_mode = MCRIS_MODE_NORMAL;
    p->VIS_mode = MCVIS_MODE_NORMAL;
    
    while (1==1) {  /* until photon reaches the surface or TOA */
      
      is_threed = (atmos->threed[MCCAOTH_TOT][p->kc]>=1); /* 3D layer */
  
      /**************************************************************/
      /* 1. find out what happens next with photon, and step length */
      /**************************************************************/
  
      /* calculate next intersection point with the box boundaries   */
      if (sample->spherical3D) {
#ifdef HAVE_SPHER
        status = intersection3D_spherical (p, atmos,
  					 sample->spherical3D_scene,
  					 &first, &step);
        if (status)
  	return err_out("Error %d returned by intersection3D_spherical\n", status);
#else
        fprintf(stderr,"Error! you are not allowed to use spherical 3D!\n");
        return -1;
#endif
      }
      else {
        if (is_threed==1) { /* 3D layer */
#if HAVE_MYSTIC3D
	  status = intersection3D (p, atmos, 10000, 0, &step);
	  if (status)
	    return err_out("Error %d returned by intersection3D\n", status);
#else
	  fprintf(stderr,"Error! you are not allowed to use mystic 3D!\n");
	  return -1;
#endif
        }
        else {
	  if (sample->spherical){
	    status = intersection1D_spherical (p, atmos,
					       sample->refraction, refind,  &step);
  	  /* Photon is absorbed by the surface.*/
  	  if (status<0)
  	    break;
  	}
  	else {
  	  intersection1D (p, atmos, 10000, 0, &step, 0);
  	}
        }
      }
  
      if (step<0)
        step=0;
      
      /* check if photon crossed the 2D surface */
      if (elev->elev2D)  {
        
        status = cross_surface (p, step,
                                elev,
                                sample->bcond,
                                &gamma);
        
        
        /* if cross_surface ran into some inconsistency */
        if (status<0) {
          fprintf (stderr, "Something terrible has happened! cross_surface() returned\n");
          fprintf (stderr, "status %d\n", status);
          return -1;
        }
        
        if (status>0) {
          step = gamma;
          p->photon_status = MCSTATUS_SURFACE;
        }
      }
  
      /*******************************************************************/
      /* 2. apply step length, move photon, and exit if photon path ends */
      /*******************************************************************/
  
      /* sum up scattering optical depth */
      if (i == 0)
        tau_sca += step * get_ksca (atmos, p, MCCAOTH_MOL);
      else
        tau_sca += step * get_ksca (atmos, p, MCCAOTH_TOT);
      p->pathlength += step;
     
      /* move photon to box boundary, to next scattering point,  */
      /* or to the 2D surface                                    */
        
      if (is_threed==1 || sample->spherical3D) { /* 3D layer */
#if HAVE_MYSTIC3D
        status = step3D (p, atmos, step, sample->bcond, 0, sample->spherical3D, 0);
#else
	fprintf(stderr,"Error! you are not allowed to use mystic 3D!\n");
	return -1;
#endif
      }
      else
        status = step1D (p, atmos, step, sample->bcond, 0, 0);
    
      if (status<0) {
        fprintf (stderr, "Error %d returned by step%dD()\n", status,
                 (is_threed==1?3:1));
        return status;
      }
        
      /* end photon travel unless it simply crosses the boundary */
      if (p->photon_status != MCSTATUS_TRAVEL )
        break;
  
      /**********************************/
      /* 3. check if photon leaves grid */
      /**********************************/
  
      /* finally check if photon arrived at the lower surface */
      /* or at the top of the atmosphere                      */
        
      if (p->kc<0) {
        p->photon_status = MCSTATUS_BOUNDARY_LOWER;
          
        if (elev->elev2D) {
          fprintf (stderr, "FATAL error! Escape photon reached the 1D boundary in a\n");
          fprintf (stderr, "             2D elevation case (%dD layer)!\n\n",
                   (is_threed==1?3:1));
          return -1;
        }
          
        break;
      }
        
      if (p->kc>atmos->Nz-1) {
        p->photon_status = MCSTATUS_BOUNDARY_UPPER;
        break;
      }
    } /* end step loop */
    tsca[i]=tau_sca;
  } /* end loop over all directions */
  
    /*********************/
    /* 4. set ris-factor */
    /*********************/

/* check if there is one or more cloudy paths */
/* XXX factor 100 is arbitrary, maybe increase it */
  if (sample->LidarLocEst) {
    for (j=1;j<numberOfShots;j++) {
      if (tsca[j] >= 100*tsca[0]) {
        tautemp[j-1]=tsca[j];
        cloudy=1;
      }
      else
        tautemp[j-1]=-1;
    }
/* if there is one or more cloudy paths, find smallest and largest tau_sca among them */
    if (cloudy) {
      for (j=0;j<5;j++) {
        if (tautemp[j] > 0) {
          taumin=tautemp[j];
          taumax=tautemp[j];
          i=j;
          break;
        }
      }
      for (j=i;j<5;j++) {
        if (tautemp[j] > 0) {
  	if (tautemp[j] < taumin) 
  	  taumin = tautemp[j];
  	else if (tautemp[j] > taumax) 
  	  taumax = tautemp[j];
        }
      }
      tau_sca=0.5*(taumin+taumax);
    }
    else 
      atmos->ris_factor = 1.0;
  }
  /* set tau_sca as average between highest and lowest optical depth */
  /* XXX not sure if this is the best possible way */
  if (tau_sca > sample->ris_optical_depth) {
    atmos->ris_factor = 1.0;
  }
  else {
    atmos->ris_factor = sample->ris_optical_depth/tau_sca;
  }

  sample->pathlength = p->pathlength;

  if (atmos->ris_factor > 1.e10)
    atmos->ris_factor = 1.e10;

  /* free memory */
  destroy_photon (p, atmos->n_caoth);

  return 0;
}



/***********************************************************************************/
/* Function: summarize_result_mishchenko_cb                               @62_30i@ */
/* Description: outputs stuff related to CB into stdout                            */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Christian Pause                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int summarize_result_mishchenko_cb (mishchenko_cb_struct* mish, char* basename )
{
  double I=0., Id=0., Q=0., Qd=0., U=0., Ud=0., V=0., Vd=0.;
  int i,j;

  I  = fiveKahanSum(mish->M_Elements_S0[0][0] , mish->M_Elements_S0[0][1] ,     mish->M_Elements_S0[0][2]
      , mish->M_Elements_S0[1][3] , 2 * mish->M_Elements_S0[1][4]);
  Id = fourKahanSum(mish->M_Elements_S0[0][0] , mish->M_Elements_S0[0][1]
      , mish->M_Elements_S0[1][3] ,     mish->M_Elements_S0[1][4]);

  Q  = fiveKahanSum(mish->M_Elements_S0[1][0] , mish->M_Elements_S0[1][1] ,     mish->M_Elements_S0[1][2]
      , mish->M_Elements_S0[0][3] , 2 * mish->M_Elements_S0[0][4]);
  Qd = fourKahanSum(mish->M_Elements_S0[1][0] , mish->M_Elements_S0[1][1]
      , mish->M_Elements_S0[0][3] ,     mish->M_Elements_S0[0][4]);

  U  = fiveKahanSum(mish->M_Elements_S0[2][0] , mish->M_Elements_S0[2][1] ,     mish->M_Elements_S0[2][2]
      , mish->M_Elements_S0[3][3] , 2 * mish->M_Elements_S0[3][4]);
  Ud = fourKahanSum(mish->M_Elements_S0[2][0] , mish->M_Elements_S0[2][1]
      , mish->M_Elements_S0[3][3] ,     mish->M_Elements_S0[3][4]);

  V  = fiveKahanSum(mish->M_Elements_S0[3][0] , mish->M_Elements_S0[3][1] ,     mish->M_Elements_S0[3][2]
      ,- mish->M_Elements_S0[2][3] ,- 2 * mish->M_Elements_S0[2][4]);
  Vd = fourKahanSum(mish->M_Elements_S0[3][0] , mish->M_Elements_S0[3][1]
      ,- mish->M_Elements_S0[2][3] ,-     mish->M_Elements_S0[2][4]);

  char mishchenkoCBFilename [FILENAME_MAX] = "";
  FILE *mishchenkoCBFile = NULL;
  strcpy (mishchenkoCBFilename, basename);
  strcat (mishchenkoCBFilename, ".mish.cb");
  if ((mishchenkoCBFile = fopen (mishchenkoCBFilename, "w")) == NULL) {
        fprintf (stderr, "Error opening %s for writing\n",mishchenkoCBFilename);
        return -1;
  }

  fprintf(mishchenkoCBFile, "Backscatter Enhancements computed by Mishchenko's formulas:\n");
  fprintf(mishchenkoCBFile, "backscatter enhancement  I: %e\n", I/Id );
  fprintf(mishchenkoCBFile, "backscatter enhancement  Q: %e\n", Q/Qd );
  fprintf(mishchenkoCBFile, "backscatter enhancement  U: %e\n", U/Ud );
  fprintf(mishchenkoCBFile, "backscatter enhancement  V: %e\n", V/Vd );
  fprintf(mishchenkoCBFile, "backscatter enhancement QP: %e\n", (I+Q)/(Id+Qd) );
  fprintf(mishchenkoCBFile, "backscatter enhancement QX: %e\n", (I-Q)/(Id-Qd) );
  fprintf(mishchenkoCBFile, "backscatter enhancement UP: %e\n", (I+U)/(Id+Ud) );
  fprintf(mishchenkoCBFile, "backscatter enhancement UX: %e\n", (I-U)/(Id-Ud) );
  fprintf(mishchenkoCBFile, "backscatter enhancement VP: %e\n", (I+V)/(Id+Vd) );
  fprintf(mishchenkoCBFile, "backscatter enhancement VX: %e\n", (I-V)/(Id-Vd) );

  for(i=0;i<4;i++) {
    for(j=0;j<5;j++) {
      mish->M_Elements_S0[i][j] = 0.;
    }
  }

  (void) fclose (mishchenkoCBFile);

  return 0;
}

/***********************************************************************************/
/* Function: save_phase_matrix_elements_temp                              @62_30i@ */
/* Description: save phase temporary matrix elements for Mishchenko matrices       */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Christian Pause                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void save_phase_matrix_elements_temp ( mishchenko_cb_struct *mish, double **phamat,
                                       int scattercounter, double* phaseMatrix )
{ 
  if (scattercounter == 0) {
    mish->R1t[0] = phamat[0][0];
    mish->R1t[1] = phamat[0][1];
    mish->R1t[2] = phamat[2][2];
    mish->R1t[3] = phamat[2][3];
    mish->R1t[4] = phamat[1][1];
    mish->R1t[5] = phamat[3][3];
  }
  else {
    mish->RMt[0] = phamat[0][0];
    mish->RMt[1] = phamat[0][1];
    mish->RMt[2] = phamat[2][2];
    mish->RMt[3] = phamat[2][3];
    mish->RMt[4] = phamat[1][1];
    mish->RMt[5] = phamat[3][3];
  }
}


/***********************************************************************************/
/* Function: save_phase_matrix_elements                                   @62_30i@ */
/* Description: save phase matrix elements for Mishchenko matrices                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Christian Pause                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void save_phase_matrix_elements ( mishchenko_cb_struct *mish, double totweight,
                                  int scattercounter, double *stokes0 )
{ 
  double R1[6], RM[6], RC[4] = {0.,0.,0.,0.};
  int i;
  
  if (scattercounter == 0) {
    for(i=0;i<6;i++) {
      R1[i] = mish->R1t[i] * totweight;
      mish->R1t[i] = 0.;
    }
    mish->R1C += 1;
    
    mish->M_Elements_S0[0][0] += R1[0] * stokes0[0];
    mish->M_Elements_S0[0][3] += R1[1] * stokes0[0];
    mish->M_Elements_S0[1][0] += R1[4] * stokes0[1];
    mish->M_Elements_S0[1][3] += R1[1] * stokes0[1];
    mish->M_Elements_S0[2][0] += R1[2] * stokes0[2];
    mish->M_Elements_S0[2][3] += R1[3] * stokes0[2];
    mish->M_Elements_S0[3][0] += R1[5] * stokes0[3];
    mish->M_Elements_S0[3][3] += R1[3] * stokes0[3];
  
  }
  else {
    for(i=0;i<6;i++) {
      RM[i] = mish->RMt[i] * totweight;
      mish->RMt[i] = 0.;
    }
    mish->RMC += 1;
    
    RC[0] += 0.5 * fourKahanSum( RM[0] , RM[4] ,- RM[2] , RM[5] );
    RC[1] += 0.5 * fourKahanSum( RM[0] , RM[4] , RM[2] ,- RM[5] );
    RC[2] += 0.5 * fourKahanSum(-RM[0] , RM[4] , RM[2] , RM[5] );
    RC[3] += 0.5 * fourKahanSum( RM[0] ,- RM[4] , RM[2] , RM[5] );
    
    mish->M_Elements_S0[0][1] += RM[0] * stokes0[0];
    mish->M_Elements_S0[0][2] += RC[0] * stokes0[0];
    mish->M_Elements_S0[0][4] += RM[1] * stokes0[0];
    mish->M_Elements_S0[1][1] += RM[4] * stokes0[1];
    mish->M_Elements_S0[1][2] += RC[1] * stokes0[1];
    mish->M_Elements_S0[1][4] += RM[1] * stokes0[1];
    mish->M_Elements_S0[2][1] += RM[2] * stokes0[2];
    mish->M_Elements_S0[2][2] += RC[2] * stokes0[2];
    mish->M_Elements_S0[2][4] += RM[3] * stokes0[2];
    mish->M_Elements_S0[3][1] += RM[5] * stokes0[3];
    mish->M_Elements_S0[3][2] += RC[3] * stokes0[3];
    mish->M_Elements_S0[3][4] += RM[3] * stokes0[3];
  }
}


/***********************************************************************************/
/* Function: kahan summation                                              @62_30i@ */
/* Description: Implementation of the Kahan summation algorithm, taken from        */
/*              http://en.wikipedia.org/wiki/Kahan_summation_algorithm             */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Christian Pause                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/
inline double kahanSum(double *f,int N)
{
  double sum = f[0];
  double c = 0.0,y,t;
  int i;
  for (i=1;i<N;i++) {
    y = f[i] - c;
    t = sum + y;
    c = (t - sum) - y;
    sum = t;
  }
  return sum;
}

double fiveKahanSum(double a, double b, double c, double d, double e) {
  double* f = NULL;
  double s;
  f = calloc(5,sizeof(double));
  f[0] = a;
  f[1] = b;
  f[2] = c;
  f[3] = d;
  f[4] = e;
  s = kahanSum(f,5);
  free(f);
  return s;
}

double fourKahanSum(double a, double b, double c, double d) {
  double* f = NULL;
  double s;
  f = calloc(4,sizeof(double));
  f[0] = a;
  f[1] = b;
  f[2] = c;
  f[3] = d;
  s = kahanSum(f,4);
  free(f);
  return s;
}

double threeKahanSum(double a, double b, double c) {
  double* f = NULL;
  double s;
  f = calloc(3,sizeof(double));
  f[0] = a;
  f[1] = b;
  f[2] = c;
  s = kahanSum(f,3);
  free(f);
  return s;
}
         


/***********************************************************************************/
/* Function: find_max                                                     @62_30i@ */
/* Description:                                                                    */
/* Find maximum value of 3 given variables                                         */
/* Parameters:                                                                     */
/* Return value: max                                                               */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Carolin Klinger, 2013.09.27                                             */
/*                                                                        @i62_30@ */
/***********************************************************************************/
static double find_max(double x, double y, double z)
{
  double max = x;

  if (max<=z)
    max = z;

  if (max<=y)
    max = y;

  return max;
}


/***********************************************************************************/
/* Function: find_min                                                     @62_30i@ */
/* Description:                                                                    */
/* Find minimum value of 3 given variables                                         */
/* Parameters:                                                                     */
/* Return value: min                                                               */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Carolin Klinger, 2013.09.27                                             */
/*                                                                        @i62_30@ */
/***********************************************************************************/
static double find_min(double x, double y, double z)
{
  double min = x;

  if (min>=z)
    min = z;

  if (min>=y)
    min = y;

  return min;
}

/*****************************************************************************************/
/* Function: pst                                                                @62_30i@ */
/* Description:                                                                          */
/* Estimate escape probability (backward heating rates, thermal; Klinger and Mayer,2013) */
/* Parameters:                                                                           */
/* Return value: escape probability                                                      */
/* Example:                                                                              */
/* Files:                                                                                */
/* Known bugs:                                                                           */
/* Author: Carolin Klinger, 2013.09.27                                                   */
/*                                                                              @i62_30@ */
/*****************************************************************************************/

static double pst(double mu, int i, double dtau)
{
  if (dtau==0 || i==0)
    return 0.5*mu*mu;
  
  if (mu==0)
    return 0; 

  /* if the optical thickness is too large we approximate the result with 0 to avoid numerical under/overflow */
  if ((double) i * dtau > 200)
    return 0;

#if HAVE_LIBGSL 
  // we need the exponential integral Ei() from the gsl 
  return 0.5*(mu*exp(-(double) i * dtau/mu)*(mu-(double) i*dtau)-((double) i*dtau)*((double) i*dtau)*gsl_sf_expint_Ei(-(double)i*dtau/mu));
#else
  fprintf (stderr, "Error, for heating rate calculations with EMABS_OPT we need the gsl.\n");
  fprintf (stderr, "Please install and recompile!\n");
  return 0.0/0.0;
#endif
}


/***********************************************************************************/
/* Function: emabsopt_location                                            @62_30i@ */
/* Description:                                                                    */
/* Calculates additional weightings and starting locations for an optimized photon */
/* distribution for the thermal backward heating rate EMABSOPT method.             */
/* Function is called in generate_photon.                                          */
/* An exact description of the method can be found in Klinger and Mayer (2013).    */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Carolin Klinger, 2013.09.27                                             */
/*                                                                        @i62_30@ */
/***********************************************************************************/
static int emabsopt_location (photon_struct *p,
			      atmosphere_struct *atmos,
			      sample_struct     *sample,
			      long int           nphotons)
{

  /* Definition of variables */

  /* General variables */
  static double kabs_save = -1.0;
  static double ***w_xy = NULL, ***w_xz = NULL, ***w_yz = NULL;
  static double *w_z = NULL;
  static double *p_sum2  = NULL;
  static double delta_tau_x = 0.0, delta_tau_y = 0.0, delta_tau_z = 0.0;
  static double tau_x = 0.0, tau_y = 0.0, tau_z  = 0.0;
  static double phot_frac = 0.0;
  static double w_in = 0.0;

  static int N_sub = 0;
  static int N_sub_cutoff = 0;
  static int nx = 0, ny = 0, nz = 0;
  static int n_xy = 0, n_yz = 0, n_xz = 0;
  static int n_max = 0;	

  int ix = 0, iy = 0, iz = 0;

  double delta_z = 0.0;
  double position = 0.0;
  double r_location = 0.0, r_frac = 0.0;
  double delta_tau_orig = 0.0;
  double tau_cutoff = 0.0, tau_min = 0.0;
  double N0 = 0.0;
  double mug = 0.0;
  double p_sum_tot = 0.0;
  double p_avg = 0.0, p_avg_out = 0.0;
  double p_in = 0.0, p_in_avg = 0.0;
  double *p_z = NULL, *p_an_z = NULL;
  double phot_dist_in = 0.0;
      
  /* ** 3d variables** */
  double N0_xy = 0.0, N0_xz = 0.0, N0_yz = 0.0;
  double p_sum_xy = 0.0, p_sum_xz = 0.0, p_sum_yz = 0.0;
  double *p_x = NULL, *p_y = NULL;
  double *p_an_x = NULL, *p_an_y = NULL;
  double ***prob_xy = NULL, ***prob_xz = NULL, ***prob_yz = NULL;
  double *p_sum1_xy = NULL, *p_sum1_xz = NULL, *p_sum1_yz = NULL;
  double ***phot_dist_xy = NULL, ***phot_dist_xz = NULL, ***phot_dist_yz = NULL;
  
  /* ** 1d variables ** */
  double N0_z = 0.0;
  double p_sum_z = 0.0;
  double *prob_z = NULL;
  double *p_sum1_z = NULL;
  double *phot_dist_z = NULL;
	 
  
  /* ************************************************************************* */
  /* In the following section, all variables are only calculated for the first */
  /* photon of aspecific optical thickness. The calculation is perfomed again, */
  /* if the optical thickness of the grid box changes, e.g. for a new          */
  /* correlated-k subband or different atmospheric/cloud properties.           */
  /* The cutoff optical thickness of a grid box is set to 10, the optical      */
  /* thickness of a sub grid box is set to 1, mug is set to 0.                 */
  /* If the cutoff optical thickness is smaller then the half of the total     */
  /* optical thickness, a inner cuboid is left - variance reduction, avoiding  */
  /* spikes!                                                                   */
  /* ************************************************************************* */
  /* The optical thickness in all 3 spatial direction is calculated and the    */
  /* minimum optical thickness of the three is found.                          */
  /* Depending on this minimum optical thickness, the number of sub grid boxes */
  /* is calculated and the escape probabilty is estimated, as well as          */
  /* additional weights.                                                       */
  /* ************************************************************************* */

  if (atmos->kabs3D->threed[MCCAOTH_TOT][p->kc]>=1){   /* if 3D */
    if (atmos->kabs3D->prof [MCCAOTH_TOT][p->kc][p->ic][p->jc]!= kabs_save){

      mug = 0;
      delta_tau_orig = 1;  /* ** Optical thicknes of a subbox */
      tau_cutoff = 10;        /* Cutoff optical thickness - from the edge to the middle of the gridbox ** */

      /* if (n_max>0){ */
	/* free memory of static variable if already allocated w_xy w_yz w_xz */
	if (w_xy != NULL) {
	  for (ix=0; ix<nx; ix++) {
	    for (iy=0; iy<ny; iy++)
	      free (w_xy[ix][iy]);
	    free (w_xy[ix]);
	  }
	  free (w_xy);
	}
	      
	if (w_xz != NULL) {
	  for (ix=0; ix<nx; ix++) {
	    for (iy=0; iy<n_max; iy++)
	      free (w_xz[ix][iy]);
	    free (w_xz[ix]);
	  }
	  free (w_xz);
	}
	      
	if (w_yz != NULL) {
	  for (ix=0; ix<n_max; ix++) {
	    for (iy=0; iy<ny; iy++)
	      free (w_yz[ix][iy]);
	    free (w_yz[ix]);
	  }
	  free (w_yz);
	}
      /* } */

  
 
      /* calculate tau of x,y and z-direction and find mimimum */
      tau_x = atmos->kabs3D->prof [MCCAOTH_TOT][p->kc][p->ic][p->jc] * ( atmos->X[p->ic+1] - atmos->X[p->ic] );
      tau_y = atmos->kabs3D->prof [MCCAOTH_TOT][p->kc][p->ic][p->jc] * ( atmos->Y[p->jc+1] - atmos->Y[p->jc] );
      tau_z = atmos->kabs3D->prof [MCCAOTH_TOT][p->kc][p->ic][p->jc] * ( atmos->Z[p->kc+1] - atmos->Z[p->kc] );
      tau_min=find_min(tau_x,tau_y,tau_z);

      /* maximum number of sub grid boxes within tau_cutoff */
      n_max = (int)(tau_cutoff/delta_tau_orig);

      /* If the minimum optical thickness is smaller then the standard "tau_cutoff=10", find a new tau_cutoff  */
      if (2*tau_cutoff>=tau_min)
	n_max = (int)(tau_min/(2*delta_tau_orig));

      if(n_max<=0){
	fprintf (stderr, "Ooops! n_max ( %d ) should never be zero or smaller than zero! Some error occured! \n", n_max);
	return -1;
      }

      
      /* calculate number of sub cubes / number of sub-distances in one direction */
      nx = (int)(tau_x/delta_tau_orig)+1;
      ny = (int)(tau_y/delta_tau_orig)+1;
      nz = (int)(tau_z/delta_tau_orig)+1;

      N_sub = nx*ny*nz;  /* total number of sub grid boxes */
      N_sub_cutoff = 2*(nx*ny*n_max+nx*nz*n_max+ny*nz*n_max);  /* total number of sub grid boxes within tau_cutoff*/
      delta_tau_x = tau_x/nx; /* optical thickness in x,y,z for a sub grid box*/
      delta_tau_y = tau_y/ny;   
      delta_tau_z = tau_z/nz;	      

      /* number of sub grid boxes on one grid box face within tau_cutoff*/
      n_xy = nx*ny*n_max;
      n_xz = nx*(nz-2*n_max)*n_max;
      n_yz = (ny-2*n_max)*(nz-2*n_max)*n_max;

      /* calculate original photon distribution per sub-cube */
      N0 = (double)nphotons/(double)N_sub;    /* !! Attention! nphotons must be double here!! */
      N0_xy =2*N0*n_xy;  /* photons started in xy plane. Factor 2, because each face exists twice! */
      N0_xz =2*N0*n_xz;  /* photons started in xz plane */
      N0_yz =2*N0*n_yz;  /* photons started in yz plane */


      /* If there is a inner cuboid/cube, perform following calculation*/
      /* if (n_max>0) { */
	/* allocate memory */
	p_x = calloc(nx, sizeof(double));
	p_y = calloc(ny, sizeof(double));
	p_z = calloc(nz, sizeof(double));
	p_an_x = calloc(nx/2, sizeof(double));
	p_an_y = calloc(ny/2, sizeof(double));
	p_an_z = calloc(nz/2, sizeof(double));
	p_sum1_xy = calloc(nx*ny*n_max, sizeof(double));
	p_sum1_xz = calloc(nx*nz*n_max, sizeof(double));
	p_sum1_yz = calloc(ny*nz*n_max, sizeof(double));

	/* free memory of static variable p_sum2 if already allocated*/
	if (p_sum2 != NULL)
	  free(p_sum2);

	p_sum2 = calloc(N_sub_cutoff/2+1, sizeof(double));
	 
	prob_xy = calloc(nx, sizeof(double**));
	for (ix=0; ix<nx; ix++) {
	  prob_xy[ix] = calloc(ny, sizeof(double*));
	  for (iy=0; iy<ny; iy++)
	    prob_xy[ix][iy] = calloc(n_max, sizeof(double)); 
	}

	prob_xz = calloc(nx, sizeof(double**));
	for (ix=0; ix<nx; ix++) {
	  prob_xz[ix] = calloc(n_max, sizeof(double*));
	  for (iy=0; iy<n_max; iy++)
	    prob_xz[ix][iy] = calloc(nz, sizeof(double));
	}  

	prob_yz = calloc(n_max, sizeof(double**));
	for (ix=0; ix<n_max; ix++) {
	  prob_yz[ix] = calloc(ny, sizeof(double*));
	  for (iy=0; iy<ny; iy++)
	    prob_yz[ix][iy] = calloc(nz, sizeof(double));
	}  

	 
	w_xy = calloc(nx, sizeof(double**));
	for (ix=0; ix<nx; ix++) {
	  w_xy[ix] = calloc(ny, sizeof(double*));
	  for (iy=0; iy<ny; iy++)
	    w_xy[ix][iy] = calloc(n_max, sizeof(double));
	}

	w_xz = calloc(nx, sizeof(double**));
	for (ix=0; ix<nx; ix++) {
	  w_xz[ix] = calloc(n_max, sizeof(double*));
	  for (iy=0; iy<n_max; iy++)
	    w_xz[ix][iy] = calloc(nz, sizeof(double));
	}  
	if (w_xz==NULL)
	  fprintf (stderr, "Fatal Error! Out of memory \n");

	w_yz = calloc(n_max, sizeof(double**));
	for (ix=0; ix<n_max; ix++) {
	  w_yz[ix] = calloc(ny, sizeof(double*));
	  for (iy=0; iy<ny; iy++)
	    w_yz[ix][iy] = calloc(nz, sizeof(double));
	} 

	phot_dist_xy = calloc(nx, sizeof(double**));
	for (ix=0; ix<nx; ix++) {
	  phot_dist_xy[ix] = calloc(ny, sizeof(double*));
	  for (iy=0; iy<ny; iy++)
	    phot_dist_xy[ix][iy] = calloc(n_max, sizeof(double));
	}

	phot_dist_xz = calloc(nx, sizeof(double**));
	for (ix=0; ix<nx; ix++) {
	  phot_dist_xz[ix] = calloc(n_max, sizeof(double*));
	  for (iy=0; iy<n_max; iy++)
	    phot_dist_xz[ix][iy] = calloc(nz, sizeof(double));
	}  

	phot_dist_yz = calloc(n_max, sizeof(double**));
	for (ix=0; ix<n_max; ix++) {
	  phot_dist_yz[ix] = calloc(ny, sizeof(double*));
	  for (iy=0; iy<ny; iy++) 
	    phot_dist_yz[ix][iy] = calloc(nz, sizeof(double));
	} 
	/* save optical thickness of current grid box. Necessary to check if optical */
	/* thickness has changed for next photon and if this calculation has to be   */
	/* performed again                                                           */
	kabs_save = atmos->kabs3D->prof [MCCAOTH_TOT][p->kc][p->ic][p->jc];


#if !HAVE_LIBGSL 
	fprintf (stderr, "Error, for heating rate calculations with EMABS_OPT we need the gsl.\n");
	fprintf (stderr, "Please install and recompile!\n");
	return -1;
#endif


	/* Calculate escape probabilty for 3 spatial direction within tau_cutoff */
	for (ix=0;ix<(nx/2);ix++)
	  p_an_x[ix] = ((pst(1,ix,delta_tau_x) - pst(1,(ix+1),delta_tau_x)) - (pst(mug,ix,delta_tau_x) - pst(mug,(ix+1),delta_tau_x))) / delta_tau_x /(1-mug);
	for (iy=0;iy<(ny/2);iy++)
	  p_an_y[iy] = ((pst(1,iy,delta_tau_y) - pst(1,(iy+1),delta_tau_y)) - (pst(mug,iy,delta_tau_y) - pst(mug,(iy+1),delta_tau_y))) / delta_tau_y /(1-mug);
	for (iz=0;iz<(nz/2);iz++)
	  p_an_z[iz] = ((pst(1,iz,delta_tau_z) - pst(1,(iz+1),delta_tau_z)) - (pst(mug,iz,delta_tau_z) - pst(mug,(iz+1),delta_tau_z))) / delta_tau_z /(1-mug);
	

	/* Escape probabilty is the same for the opposit spatial direction, therefore use above calculated escepa probability */
	for (ix=0;ix<(nx/2);ix++) {
	  p_x[ix] = p_an_x[ix];
	  p_x[nx-ix-1] = p_an_x[ix];
	}
	for (iy=0;iy<(ny/2);iy++) {
	  p_y[iy] = p_an_y[iy];
	  p_y[ny-iy-1] = p_an_y[iy];
	}
	for (iz=0;iz<(nz/2);iz++) {
	  p_z[iz] = p_an_z[iz];
	  p_z[nz-iz-1] = p_an_z[iz];
	}
		
		
	/* averaged probability */		
	for (ix=0;ix<nx;ix++)
	  for (iy=0;iy<ny;iy++)
	    for (iz=0;iz<n_max;iz++) {
	      prob_xy[ix][iy][iz] =  find_max(p_x[ix],p_y[iy],p_z[iz]);
	      p_avg_out += (2*N0*sqrt(prob_xy[ix][iy][iz]*(1-prob_xy[ix][iy][iz])));
	    }
				
	for (ix=0;ix<nx;ix++)
	  for (iy=0;iy<n_max;iy++)
	    for (iz=n_max;iz<(nz-n_max);iz++) {
	      prob_xz[ix][iy][iz] =  find_max(p_x[ix],p_y[iy],p_z[iz]);
	      p_avg_out += (2*N0*sqrt(prob_xz[ix][iy][iz]*(1-prob_xz[ix][iy][iz])));
	    }
	  	
	for (ix=0;ix<n_max;ix++)
	  for (iy=n_max;iy<(ny-n_max);iy++)
	    for (iz=n_max;iz<(nz-n_max);iz++){
	      prob_yz[ix][iy][iz] =  find_max(p_x[ix],p_y[iy],p_z[iz]);
	      p_avg_out += (2*N0*sqrt(prob_yz[ix][iy][iz]*(1-prob_yz[ix][iy][iz])));
	    }
	
	p_in = p_x[n_max+1];
	p_in_avg = (nphotons-(N0_yz+N0_xz+N0_xy))*sqrt(p_in*(1-p_in));
	p_avg = (p_avg_out+p_in_avg)/nphotons; 

	for (ix=0;ix<nx;ix++)
	  for (iy=0;iy<ny;iy++)
	    for (iz=0;iz<n_max;iz++) {
	      w_xy[ix][iy][iz] = (p_avg/sqrt(prob_xy[ix][iy][iz]*(1-prob_xy[ix][iy][iz])));
	      phot_dist_xy[ix][iy][iz] = 2*N0/w_xy[ix][iy][iz];
	      phot_frac+=phot_dist_xy[ix][iy][iz];
	      p_sum1_xy[ix*ny*n_max+iy*n_max+iz] = phot_dist_xy[ix][iy][iz]/(N0_xy+N0_xz+N0_yz);
	      p_sum_xy += phot_dist_xy[ix][iy][iz]/(N0_xy+N0_xz+N0_yz);
	    }

	for (ix=0;ix<nx;ix++)
	  for (iy=0;iy<n_max;iy++)
	    for (iz=n_max;iz<(nz-n_max);iz++) {
	      w_xz[ix][iy][iz] = (p_avg/sqrt(prob_xz[ix][iy][iz]*(1-prob_xz[ix][iy][iz])));
	      phot_dist_xz[ix][iy][iz] = 2*N0/w_xz[ix][iy][iz];
	      phot_frac+=phot_dist_xz[ix][iy][iz];
	      p_sum1_xz[ix*nz*n_max+iy*nz+iz] = phot_dist_xz[ix][iy][iz]/(N0_xy+N0_xz+N0_yz);
	      p_sum_xz += phot_dist_xz[ix][iy][iz]/(N0_xy+N0_xz+N0_yz);
	    }

	for (ix=0;ix<n_max;ix++)
	  for (iy=n_max;iy<(ny-n_max);iy++)
	    for (iz=n_max;iz<(nz-n_max);iz++) {
	      w_yz[ix][iy][iz] = (p_avg/sqrt(prob_yz[ix][iy][iz]*(1-prob_yz[ix][iy][iz])));
	      phot_dist_yz[ix][iy][iz] = 2*N0/w_yz[ix][iy][iz];
	      phot_frac+=phot_dist_yz[ix][iy][iz];
	      p_sum1_yz[ix*ny*nz+iy*nz+iz] = phot_dist_yz[ix][iy][iz]/(N0_xy+N0_xz+N0_yz);
	      p_sum_yz += phot_dist_yz[ix][iy][iz]/(N0_xy+N0_xz+N0_yz);
	    }
	 
 
	w_in = (p_avg/sqrt(p_in*(1-p_in)));


	phot_dist_in = (nphotons-(N0_yz+N0_xz+N0_xy))/w_in;
	p_sum_tot=p_sum_xy+p_sum_xz+p_sum_yz;
	p_sum2[0] = 0.0;
	phot_frac=phot_frac/(phot_frac+phot_dist_in);

	for(ix=1;ix<=(nx*ny*n_max);ix++)
	  p_sum2[ix] = p_sum2[ix-1]+p_sum1_xy[ix-1]/p_sum_tot;

	for(ix=(nx*ny*n_max+1);ix<=(nx*ny*n_max+nx*nz*n_max);ix++)
	  p_sum2[ix] = p_sum2[ix-1]+p_sum1_xz[ix-1-nx*ny*n_max]/p_sum_tot;

	for(ix=(nx*ny*n_max+nx*nz*n_max+1);ix<=(nx*ny*n_max+nx*nz*n_max+n_max*ny*nz);ix++)
	  p_sum2[ix] = p_sum2[ix-1]+p_sum1_yz[ix-1-(nx*ny*n_max+nx*nz*n_max)]/p_sum_tot;


	/* **** free memory **** */
	free(p_x);
	free(p_y);
	free(p_z);
	free(p_an_x);
	free(p_an_y);
	free(p_an_z);
	free(p_sum1_xy);
	free(p_sum1_xz);
	free(p_sum1_yz);

	/* free memory of static variable prob */
	for (ix=0; ix<(nx); ix++) {
	  for (iy=0; iy<(ny); iy++)
	    free (prob_xy[ix][iy]);
	  free (prob_xy[ix]);
	}
	free (prob_xy);

	for (ix=0; ix<(nx); ix++) {
	  for (iy=0; iy<(n_max); iy++)
	    free (prob_xz[ix][iy]);
	  free (prob_xz[ix]);
	}
	free (prob_xz);

	for (ix=0; ix<(n_max); ix++) {
	  for (iy=0; iy<ny; iy++)
	    free (prob_yz[ix][iy]);
	  free (prob_yz[ix]);
	}
	free (prob_yz);

	/* free memory of static variable phot_dist */
	for (ix=0; ix<(nx); ix++) {
	  for (iy=0; iy<(ny); iy++)
	    free (phot_dist_xy[ix][iy]);
	  free (phot_dist_xy[ix]);
	}
	free (phot_dist_xy);
	    
	for (ix=0; ix<(nx); ix++) {
	  for (iy=0; iy<(n_max); iy++)
	    free (phot_dist_xz[ix][iy]);
	  free (phot_dist_xz[ix]);
	}
	free (phot_dist_xz);
	    
	for (ix=0; ix<(n_max); ix++) {
	  for (iy=0; iy<ny; iy++)
	    free (phot_dist_yz[ix][iy]);
	  free (phot_dist_yz[ix]);
	}
	free (phot_dist_yz);

    } /* close kabs_save */    
  }  /* end 3D */
	  
	  
  /* If there is a inner cuboid/cube, perform following calculation*/	  
  else {   /* if 1D */

    /* The whole calculation from the 3D case is performed only in z-direction */

    if (atmos->kabs->prof [MCCAOTH_TOT][p->kc]!= kabs_save){

      mug = 0;
      delta_tau_orig = 1;
      tau_cutoff = 10;

      /* calculate tau of x,y and z-direction */
      tau_z = atmos->kabs->prof [MCCAOTH_TOT][p->kc] * ( atmos->Z[p->kc+1] - atmos->Z[p->kc] );

      /* maximum number of sub grid boxes within tau_cutoff */
      n_max = (int)(tau_cutoff/delta_tau_orig);

     /* If the minimum optical thickness is smaller then the standard "tau_cutoff=10", find a new tau_cutoff  */	      
      if (2*tau_cutoff >= tau_z)
	n_max = (int)(tau_z/(2*delta_tau_orig));
	      
      if(n_max<=0){
	fprintf (stderr, "Ooops! n_max ( %d ) should never be zero or smaller than zero! Some error occured! \n", n_max);
	return -1;
      }

      /* calculate number of sub cubes / number of sub-distances in one direction */
      nz = (int)(tau_z/delta_tau_orig)+1;
      N_sub = nz;
      N_sub_cutoff = 2*n_max;
      delta_tau_z = tau_z/nz;	      

      /* calculate original photon distribution per sub-cube */
      N0 = nphotons/N_sub;
      N0_z =N0*N_sub_cutoff;  

      /* If there is a inner cuboid/cube, perform following calculation*/
      /* if (n_max>0) { */

	/* allocate memory */
	p_z = calloc(nz, sizeof(double));
	p_an_z = calloc(nz/2, sizeof(double));
	p_sum1_z = calloc(n_max, sizeof(double));
	prob_z = calloc(nz, sizeof(double));
	w_z = calloc(nz, sizeof(double));
	phot_dist_z = calloc(nz, sizeof(double));


	/* free memory of static variable p_sum2 */
	if (p_sum2 != NULL)
	  free(p_sum2);

	p_sum2 = calloc(N_sub_cutoff/2+1, sizeof(double));

	/* save optical thickness of current grid box. Necessary to check if optical */
	/* thickness has changed for next photon and if this calculation has to be   */
	/* performed again */
	kabs_save = atmos->kabs->prof [MCCAOTH_TOT][p->kc];

	/* estimate escape probability */
	for (iz=0;iz<(nz/2);iz++)
	  p_an_z[iz] = ((pst(1,iz,delta_tau_z) - pst(1,(iz+1),delta_tau_z)) - (pst(mug,iz,delta_tau_z) - pst(mug,(iz+1),delta_tau_z))) / delta_tau_z /(1-mug);


	/* calculate total p for each sub-part in z */
	for (iz=0;iz<(nz/2);iz++) {
	  p_z[iz] = p_an_z[iz];
	  p_z[nz-iz-1] = p_an_z[iz];
	} 

	/* averaged probability */	      
	for (iz=0;iz<n_max;iz++) {
	  prob_z[iz] = p_z[iz];
	  p_avg_out += (2*N0*sqrt(prob_z[iz]*(1-prob_z[iz])));
	}
	  
	p_in = p_z[n_max+1];
	p_in_avg = (nphotons-(N0_z))*sqrt(p_in*(1-p_in));
	p_avg = (p_avg_out+p_in_avg)/nphotons; 

	/* calculate weight */
	for (iz=0;iz<n_max;iz++) {
	  w_z[iz] = (p_avg/sqrt(prob_z[iz]*(1-prob_z[iz])));
	  phot_dist_z[iz] = 2*N0/w_z[iz];
	  phot_frac+=phot_dist_z[iz];
	  p_sum1_z[iz] = phot_dist_z[iz]/N0_z;
	  p_sum_z += phot_dist_z[iz]/N0_z;
	}

	w_in = (p_avg/sqrt(p_in*(1-p_in)));

	phot_dist_in = (nphotons-(N0_z))/w_in;
	p_sum_tot=p_sum_z;
	p_sum2[0] = 0.0;
	phot_frac=phot_frac/(phot_frac+phot_dist_in);

	for (iz=1;iz<=(n_max);iz++)
	  p_sum2[iz] = p_sum2[iz-1]+p_sum1_z[iz-1]/p_sum_tot;
	


	/* **** free memory **** */
	free(p_z);
	free(p_an_z);
	free(p_sum1_z);
	free(prob_z);
	free(phot_dist_z);

    } /* end kabs_save */
  } /* end 1D */

  /* ************************************************************************* */
  /* End of calcualtion for a changing optical thickness                       */
  /* ************************************************************************* */


  /* ************************************************************************* */
  /* In the following - starting locations of the photons within the sub grid  */
  /* boxes are choosen randomly and the corrsponding weight is applied         */
  /* ************************************************************************* */

  if (atmos->kabs3D->threed[MCCAOTH_TOT][p->kc]>=1) { /* if 3D */	  

    /* if(n_max>0){ */
    r_frac=uvspec_random();

    if (r_frac<phot_frac) {  /* if photon is within cutoff optical thickness, if not, see calulation for inner cubiod/cube */
	    
      /* getting position/intervall of cumulative probability */
      r_location=uvspec_random();
      position = locate(p_sum2, N_sub_cutoff/2, r_location);

      /* deriving x,y and z coordinate from culumaltive probability */
      /* Find exact starting position and apply weight */
      if (position<nx*ny*n_max) {
	ix = position / (ny*n_max);
	iy = (position - ix*ny*n_max)/n_max;
	iz = position - ix*ny*n_max - iy*n_max;
	p->weight *= w_xy[ix][iy][iz];
	p->weight_emis *= w_xy[ix][iy][iz]; 
	r_location=uvspec_random();
	if (r_location<0.5) {
	  p->x[0] = ((double) sample->backward_is + ((double) ix + uvspec_random()) / (double)nx)*sample->delX; 
	  p->x[1] = ((double) sample->backward_js + ((double) iy + uvspec_random()) / (double)ny)*sample->delY; 
	  delta_z = ( atmos->Z[p->kc+1] - atmos->Z[p->kc] )/nz * (iz + uvspec_random());
	}
	else {
	  p->x[0] = ((double) sample->backward_is + ((double) ix + uvspec_random()) / (double)nx)*sample->delX; 
	  p->x[1] = ((double) sample->backward_js + ((double) iy + uvspec_random()) / (double)ny)*sample->delY; 
	  delta_z = ( atmos->Z[p->kc+1] - atmos->Z[p->kc] )/nz * ((nz-iz-1) + uvspec_random());
	}
      }	    

      else if ((position<(nx*ny*n_max+nx*nz*n_max))&&(position>=(nx*ny*n_max))) {
	position=position-nx*ny*n_max;
	ix = position / (nz*n_max);
	iy = (position - ix*nz*n_max)/nz;
	iz = position - ix*nz*n_max - iy*nz;
	p->weight *= w_xz[ix][iy][iz];
	p->weight_emis *= w_xz[ix][iy][iz];
	r_location=uvspec_random();
	if (r_location<0.5) {
	  p->x[0] = ((double) sample->backward_is + ((double) ix + uvspec_random()) / (double)nx)*sample->delX; 
	  p->x[1] = ((double) sample->backward_js + ((double) iy + uvspec_random()) / (double)ny)*sample->delY; 
	  delta_z = ( atmos->Z[p->kc+1] - atmos->Z[p->kc] )/nz * (iz + uvspec_random());
	}
	else {
	  p->x[0] = ((double) sample->backward_is + ((double) ix + uvspec_random()) / (double)nx)*sample->delX; 
	  p->x[1] = ((double) sample->backward_js + ((double) (ny-iy-1) + uvspec_random()) / (double)ny)*sample->delY; 
	  delta_z = ( atmos->Z[p->kc+1] - atmos->Z[p->kc] )/nz * (iz + uvspec_random());
	} 
      }	

      else {
	position=position-nx*ny*n_max-nx*nz*n_max;
	ix = position / (ny*nz);
	iy = (position - ix*ny*nz)/nz;
	iz = position - ix*ny*nz - iy*nz;
	p->weight *= w_yz[ix][iy][iz];
	p->weight_emis *= w_yz[ix][iy][iz];
	r_location=uvspec_random();
	if (r_location<0.5) {
	  p->x[0] = ((double) sample->backward_is + ((double) ix + uvspec_random()) / (double)nx)*sample->delX; 
	  p->x[1] = ((double) sample->backward_js + ((double) iy + uvspec_random()) / (double)ny)*sample->delY; 
	  delta_z = ( atmos->Z[p->kc+1] - atmos->Z[p->kc] )/nz * (iz + uvspec_random());
	}
	else {
	  p->x[0] = ((double) sample->backward_is + ((double) (nx-ix-1) + uvspec_random()) / (double)nx)*sample->delX; 
	  p->x[1] = ((double) sample->backward_js + ((double) iy + uvspec_random()) / (double)ny)*sample->delY; 
	  delta_z = ( atmos->Z[p->kc+1] - atmos->Z[p->kc] )/nz * (iz + uvspec_random());
	}
      }

      if (!sample->spherical3D)
	atmos_coord (p, atmos, &(p->ic), &(p->jc));
	    
      set_photon_z (p->x[2] + delta_z, atmos, p);  

    } /* end phot_frac */


    else { /* if photon is started in inner cuboid/cube */
      p->x[0] = ((double) sample->backward_is + (((double) n_max + (double)(nx-2*n_max) * uvspec_random())) / (double)nx)*sample->delX; 
      p->x[1] = ((double) sample->backward_js + (((double) n_max + (double)(ny-2*n_max) * uvspec_random())) / (double)ny)*sample->delY; 
      delta_z = (atmos->Z[p->kc+1]-atmos->Z[p->kc])/nz * (n_max + ((nz-2*n_max) * uvspec_random()));
      p->weight *= w_in;
      p->weight_emis *= w_in;
	    
      if (!sample->spherical3D)
	atmos_coord (p, atmos, &(p->ic), &(p->jc));
	    
      set_photon_z (p->x[2] + delta_z, atmos, p);
    }
  }   /* end 3D */


  else {     /* if 1d */
    /* if(n_max>0){ */
    r_frac=uvspec_random();

    if (r_frac<phot_frac) { /* if photon is within cutoff optical thickness, if not, see calulation for inner cubiod/cube */
	    
      /* getting position/intervall of cumulative probability */
      r_location=uvspec_random();
      position = locate(p_sum2, N_sub_cutoff/2, r_location);
	      
      /* deriving x,y and z coordinate from culumaltive probability */
      /* Find exact starting position and apply weight */
      iz = position;
      p->weight *= w_z[iz];
      p->weight_emis *= w_z[iz]; 
	      
      r_location=uvspec_random();
	      
      if (r_location<0.5)
	delta_z = ( atmos->Z[p->kc+1] - atmos->Z[p->kc] )/nz * (iz + uvspec_random());
      else
	delta_z = ( atmos->Z[p->kc+1] - atmos->Z[p->kc] )/nz * ((nz-iz-1) + uvspec_random());
	      
      if (!sample->spherical3D)
	atmos_coord (p, atmos, &(p->ic), &(p->jc));
	      
      set_photon_z (p->x[2] + delta_z, atmos, p);
    }

    else { /* if photon is started in inner cuboid/cube */
      delta_z = (atmos->Z[p->kc+1]-atmos->Z[p->kc])/nz * (n_max + ((nz-2*n_max) * uvspec_random()));
      p->weight *= w_in;
      p->weight_emis *= w_in;
	      
      if (!sample->spherical3D)
	atmos_coord (p, atmos, &(p->ic), &(p->jc));
	      
      set_photon_z (p->x[2] + delta_z, atmos, p);
    }
  }       /* end if 1d */

  return 0;
} /* End function emabsopt_location */


/***********************************************************************************/
/* Function: denet_weight                                                 @62_30i@ */
/* Description:                                                                    */
/* Calculates additional weightings according to starting direction                */
/* for dEnet method. Counts photons (counter1)                                     */
/* An exact description of the method can be found in Klinger and Mayer (2013).    */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Carolin Klinger, 2013.09.27                                             */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int denet_weight (photon_struct *p,
			 atmosphere_struct *atmos,
			 sample_struct *sample,
			 long int nphotons,
			 int *counter1)
{
   
  if (atmos->kabs3D->threed[MCCAOTH_TOT][p->kc]>=1) {    /* 3D case */
    
    p->weight = sample->weight_heat;

    /* upper and lower area  - set z new, x and y are constant*/
    if ((*counter1 < (sample->n_xy/4)) || (*counter1 < (3*sample->n_xy/4) && *counter1 >= (sample->n_xy/2))) {
      /* set location to the edge of the layer (lower edge)*/
      set_photon_z (atmos->Z[p->kc], atmos, p);
      p->n_dEnet = sample->n_xy/4;     

      if (p->dir.dx[2]<0)
	p->dEnet_component = 0;
      if (p->dir.dx[2]>0){
	p->weight *= -1.0;
	p->dEnet_component = 1;   
      }
    }

    else if ((*counter1 < (sample->n_xy/2) && *counter1 >= (sample->n_xy/4)) || (*counter1 < (sample->n_xy) && *counter1 >= (3*sample->n_xy/4))) {
      /* set location to the edge of the layer (upper edge)*/
      set_photon_z ((atmos->Z[p->kc+1]), atmos, p);
      p->n_dEnet = sample->n_xy/4;        

      if (p->dir.dx[2]>0)
	p->dEnet_component = 2;
      if (p->dir.dx[2]<0){
	p->weight *= -1.0;
	p->dEnet_component = 3;
      }  
    }

    /* left and right area - set x new, y and z are constant*/
    else if ((*counter1 < (sample->n_xy+(sample->n_yz/4)) && *counter1 >= (sample->n_xy)) || (*counter1 < (sample->n_xy+(3*sample->n_yz/4)) && *counter1 >= (sample->n_xy+(sample->n_yz/2)))) {
      /* set location to the edge of the layer (left)*/
      p->x[0] = (double) sample->backward_is*sample->delX;
      p->n_dEnet = sample->n_yz/4;     

      if (p->dir.dx[0]<0)
	p->dEnet_component = 4;
      if (p->dir.dx[0]>0){
	p->weight *= -1.0;
	p->dEnet_component = 5;
      }
    }

    else if ((*counter1 < (sample->n_xy+(sample->n_yz/2)) && *counter1 >= (sample->n_xy+(sample->n_yz/4))) || (*counter1 < (sample->n_xy+(sample->n_yz)) && *counter1 >= (sample->n_xy+(3*sample->n_yz/4)))) {
      /* set location to the edge of the layer (right)*/;
      p->x[0] = (double) (sample->backward_is+1)*sample->delX;
      p->n_dEnet = sample->n_yz/4;     

      if (p->dir.dx[0]>0)
	p->dEnet_component = 6;
      if (p->dir.dx[0]<0){
	p->weight *= -1.0;
	p->dEnet_component = 7;
      }
    }

    /* front and back area - set y new, x and z are constant*/
    else if ((*counter1 < (sample->n_xy+sample->n_yz+(sample->n_xz/4)) && *counter1 >= (sample->n_xy+(sample->n_yz))) || (*counter1 < (sample->n_xy+sample->n_yz+(3*sample->n_xz/4)) && *counter1 >= (sample->n_xy+sample->n_yz+(sample->n_xz/2)))) {
      /* set location to the edge of the layer (forward)*/
      p->x[1] = (double) sample->backward_js*sample->delY;
      p->n_dEnet = sample->n_xz/4;     

      if (p->dir.dx[1]<0)
	p->dEnet_component = 8;
      if (p->dir.dx[1]>0){
	p->weight *= -1.0;
	p->dEnet_component = 9;
      }
    }

    else { 
      /* set location to the edge of the layer (backward)*/
      p->x[1] = (double) (sample->backward_js + 1)*sample->delY;
      p->n_dEnet = sample->n_xz/4; 
    
      if (p->dir.dx[1]>0)
	p->dEnet_component = 10;
      if (p->dir.dx[1]<0){
	p->weight *= -1.0;
	p->dEnet_component = 11;  
      }
    }
  }

  else { /* 1d case */
    /* upper and lower area  - set z new, x and y are constant*/
    //   if ((*counter1 < (nphotons/4)) || (*counter1 < (3*nphotons/4) && *counter1 >= (2*nphotons/4))) {
    if ((*counter1 < (sample->n_xy/4)) || (*counter1 < (3*sample->n_xy/4) && *counter1 >= (sample->n_xy/2))) {
      /* set location to the edge of the layer (lower edge)*/
      set_photon_z (atmos->Z[p->kc], atmos, p);
      p->n_dEnet = sample->n_xy/4;     

      if (p->dir.dx[2]<0)
	p->dEnet_component = 0;
      if (p->dir.dx[2]>0) {
	p->weight *= -1.0;
	p->dEnet_component = 1;
      }
    }
    //    else if ((*counter1 < (2*nphotons/4) && *counter1 >= (nphotons/4)) || (*counter1 < (4*nphotons/4) && *counter1 >= (3*nphotons/4))) {
    else { //((*counter1 < (2*sample->n_xy/4) && *counter1 >= (1*sample->n_xy/4)) || (*counter1 < (4*sample->n_xy/4) && *counter1 >= (3*sample->n_xy/4))) {
      /* set location to the edge of the layer (upper edge)*/
      set_photon_z ((atmos->Z[p->kc+1]), atmos, p);
      p->n_dEnet = sample->n_xy/4;   

      if (p->dir.dx[2]>0)
	p->dEnet_component = 2;
      if (p->dir.dx[2]<0) {
	p->weight *= -1.0;
	p->dEnet_component = 3;
      }
    }
  }	    

  atmos_coord (p, atmos, &(p->ic), &(p->jc));
	  
  (*counter1)++;

  if (*counter1==nphotons)
    (*counter1)=0;	    
	    

  return 0;
}


/***********************************************************************************/
/* Function: denet_set_direction                                          @62_30i@ */
/* Description:                                                                    */
/* Sets diretion for dEnet method (photon direction at each face of grid box.      */
/* An exact description of the method can be found in Klinger and Mayer (2013).    */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Carolin Klinger, 2013.09.27                                             */
/*                                                                        @i62_30@ */
/***********************************************************************************/
static int denet_set_direction (photon_struct *p,
				atmosphere_struct *atmos,
				sample_struct *sample,
				long int nphotons,
				int *counter1,
				double *norm)
{
  int i=0;

  if (atmos->kabs3D->threed[MCCAOTH_TOT][p->kc]>=1) {   /* 3D */

    if (*counter1 < (sample->n_xy/2))
      /* downward */
      random_Lambertian_normal (&(p->dir), norm);
    else if (*counter1 < (2*sample->n_xy/2)){
      /* upward */
      for (i=0; i<3; i++)
	norm[i] = - norm[i];
      random_Lambertian_normal (&(p->dir), norm);
    }

    else if (*counter1 >= (2*sample->n_xy/2) && *counter1 < (sample->n_xy+(sample->n_yz/2))) {
      /* right */
      norm[0] = 1.0;
      norm[1] = 0.0;
      norm[2] = 0.0;
      random_Lambertian_normal (&(p->dir), norm);
    }

    else if (*counter1 >= (sample->n_xy+(sample->n_yz/2)) && *counter1 < (sample->n_xy+(2*sample->n_yz/2))) {
      /* left */
      norm[0] = -1.0;
      norm[1] = 0.0;
      norm[2] = 0.0;
      random_Lambertian_normal (&(p->dir), norm);
    }

    else if (*counter1 >= (sample->n_xy+(2*sample->n_yz/2)) && *counter1 < (sample->n_xy+sample->n_yz+(sample->n_xz/2))) {
      /* forward */
      norm[0] = 0.0;
      norm[1] = -1.0;
      norm[2] = 0.0;
      random_Lambertian_normal (&(p->dir), norm);
    }

    else {
      /* backward */
      norm[0] = 0.0;
      norm[1] = 1.0;
      norm[2] = 0.0;
      random_Lambertian_normal (&(p->dir), norm);
    }
  }

  else { /* 1D */
    //  if (*counter1 < (nphotons/2)) {
    if (*counter1 < (sample->n_xy/2)){
      /* downward */
      random_Lambertian_normal (&(p->dir), norm);
    }
    else {
      /* upward */
      for (i=0; i<3; i++)
	norm[i] = - norm[i];
      random_Lambertian_normal (&(p->dir), norm);
    }
  }

  return 0;
}

/***********************************************************************************/
/* Function: dist_photon_direction_by_phasemax                                     */
/* Description:                                                                    */
/*  Distribute initial backward photon directions according to phase_max.          */
/*  This used for the option mc_panorama circumsolar_var_red                       */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs: Extent of sundisk or extraterrestrial sunshape is not considered.   */
/*             Implementation should be easy. Just fold phase_max with sunshape.   */
/* Author: Bernhard Reinhardt, 2014.03.27                                          */
/*                                                                        @i62_30@ */
/***********************************************************************************/
static int dist_photon_direction_by_phasemax (photon_struct *p, 
					sample_struct *sample,
					double *sza,
					double *cossza,
					double *sinsza)
{
  double F_max=0.0, F_min=0.0;
  double mu_min=0.0, mu_max=0.0;
  double phase_max_value = 0.0;

  /* Distribute photon directions according to
     phase_max. Store upper and lower boundaries of
     pixel*/
  mu_min    = *cossza;
  mu_max    = cosd ( *sza + sample->pan_dtheta );
  /* Get cummulative values F of phase_max for upper and
     lower boundary of pixel */

  F_min = get_F( sample->phase_max[0], -mu_min, 0 );
  F_max = get_F( sample->phase_max[0], -mu_max, 0 );

  /* Choose a random starting mu according to phase
     function phase_max for the photon by abusing the
     scattering routine sc_mu */
  *cossza = - sc_mu ( sample->phase_max[0], 1, 0, F_max, F_min );
  *sza = acosd ( *cossza );
  *sinsza = sqrt ( 1.0 - (*cossza) * (*cossza) );

  /* Now adjust photon weight for the "unphysical"
     choice of direction: Divide photon-weight by
     probability from phase_max for sending it into this
     direction */
  phase_max_value = get_phase_max ( sample->phase_max, sample->n_phase_max, -1.0* *cossza, p->DDIS_SC_mode );

  p->weight /= 2.0 / fabs ( F_max - F_min ) * phase_max_value;
  /* Now multiply photon weight by the photons
     probability it would have been sent to this
     direction by geometrical considerations */
  p->weight *= 2.0 / fabs( - mu_max + mu_min );

  /* Handling of input file errors */
  if( sample->pan_alignment != MCPAN_ALIGNMENT_SUN ){
    fprintf(stderr,"%s %s",
	    "\n\n mc_panorama circumsolar_var_red can only be used together with",
	    "mc_panorama_alignment sun!\n\n");
    return -1;
  }  
  return 0;
}
