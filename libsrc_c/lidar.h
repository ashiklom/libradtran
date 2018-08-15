/************************************************************************
 * $Id: lidar.h 3107 2015-05-12 12:17:19Z Claudia.Emde $
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

#ifndef __lidar_h
#define __lidar_h

#undef MOON

/* local estimator type */
#define MCLIDAR_NONE              0
#define MCLIDAR_SPACE             1
#define MCLIDAR_FALCON            2
#define MCLIDAR_ASCOPE            3
#define MCLIDAR_POLARIZE          4
#define MCLIDAR_SIMPLE            5
#define MCLIDAR_TEST              6
#define MCLIDAR_NODDIS            7
#define MCLIDAR_PILOT             8
#define MCLIDAR_MOON              9
#define MCRADAR                  10

/* additional channels for Lidar */
#define LIDAR_CHANNEL_NORMAL      0
#define LIDAR_CHANNEL_RAMAN       1
#define LIDAR_CHANNEL_HSRL        2

typedef struct qidd_struct {
  int vis;
  int ris;
  int willvis;
  int willris;
  int recalc;

  double vis_rfac;
  double vis_r;
  double vis_step;
  double vis_step2;
  double vis_dist;
  double vis_cosalpha;

  double ris_r;
  double ris_step;
} qidd_struct;

/* Initialization functions */

int set_lidar_settings (int locest, sample_struct *sample);

int mc_lidar_fix_taumax (sample_struct *sample, atmosphere_struct *atmos,
			 result_struct *result, elevation_struct *elev,
			 photon_struct *p, float *refind,
			 int quiet);

void mc_lidar_setup_ris_mas (atmosphere_struct *atmos, sample_struct *sample);

lidar_result_field *calloc_lidar_result_field ( int Nli,
						int Nt,
						int n_caoth,
						int ntypes,
						int norders,
						int Np,
						int No,
						int Na,
						int std,
						int jacobian,
						int Nz );

jacobian_result_field *calloc_jacobian_result_field (int std, int Nx, int Ny, int Nz, int n_caoth);

void mc_lidar_setup_jac_exp (sample_struct *sample, atmosphere_struct *atmos);

void calloc_photon_lidar_part ( photon_struct *p,
				sample_struct *sample,
				int            Nz,
				int            source,
				int            n_caoth );

/* finishing functions */

void destroy_photon_lidar_part (photon_struct *p, int n_caoth);

void free_lidar_result_field (lidar_result_field *res, int Nli, int Nt, int std);
void free_jacobian_result_field (jacobian_result_field *res);

/* functions outside photon journey */

int generate_photon_lidar (photon_struct *p, sample_struct *sample, atmosphere_struct *atmos);
int initialize_vis_fod (sample_struct *sample, photon_struct *p);

void lidar_collect_jacobian_stddev ( lidar_result_field *result_lidar,
				     int                 Nt,
				     int                 ili,
				     int                 n_caoth );

void abs_collect_jacobian_stddev (jacobian_result_field *result_jacobian, int is, int js);

/* travel_tau functions */

void mc_lidar_prepare_travel_tau (sample_struct *sample, photon_struct *p,
				  int *moving_towards_detector);


int mc_lidar_prepare_travel_tau_step (sample_struct *sample, photon_struct *p,
				      qidd_struct *qidd, atmosphere_struct *atmos,
				      double dtau,
				      int moving_towards_detector);

int mc_lidar_travel_tau_action (photon_struct *p, double *step, qidd_struct *qidd);

double ksca_for_tausca (photon_struct *p, double *step, double *ksca,
			qidd_struct *qidd);

void mc_lidar_travel_tau_addition (atmosphere_struct *atmos, sample_struct *sample,
				   photon_struct *p,
				   double *step, qidd_struct *qidd, double *ksca);

void mc_lidar_travel_special_cases (sample_struct *sample, photon_struct *p,
				    qidd_struct *qidd, double *ksca, double dtau);

/* local estimator functions */

int mc_lidar_locest ( result_struct     *result,
		      sample_struct     *sample,
		      atmosphere_struct *atmos,
		      photon_struct     *p,
		      int                itl,
		      int                LoWeRR_crit );

int lidar_detector_kernel ( lidar_detector  lidar,
			    double         *x_p,
			    /* Output */
			    locest_struct  *lest );

/* helper functions */

int calc_locest_connection ( double        *x_p,
			     double        *x_det,
			     double        *x_end,
			     /* Output */
			     locest_struct *lest_out );

int derive_cphi (double *dx1, double *dx2, double stheta, double *cphi);

/* scattering functions */

int calc_distance_to_cone (double c_th, double s_th, double cphi2,
			   double y_sc_cc, double z_sc_cc, double l_sc_cc_inv,
			   double t_det,
			   locest_struct *lest);

int check_lidar_in_atmos (sample_struct     *sample,
			  atmosphere_struct *atmos );

/* non-lidary functions */

int intersection3D_spherical (photon_struct     *p,
			      atmosphere_struct *atmos,
			      int                spherical3D_scene,
			      int               *first,
			      double            *step);

int coord_spherical3D (photon_struct *p,
		       int            Nx,
		       double        *X,
		       int            Ny,
		       double        *Y,
		       int            Nz,
		       float         *Z,
		       int            Nyg,
		       double        *Yg,
		       double         r_earth,
		       int           *ic,
		       int           *jc,
		       int           *kc,
		       int           *jcg,
		       int            warn);

int setup_spherical3D (char    *type,
		       int      Nx,
		       int      Ny,
		       int      spherical3D_scene,
		       double   spherical3D_scene_lon_min,
		       double   spherical3D_scene_lon_max,
		       double   spherical3D_scene_lat_min,
		       double   spherical3D_scene_lat_max,
		       double  *delX,
		       double  *delY,
		       double  *xmax,
		       double  *ymax,
		       double **X,
		       double **Y,
		       int     *Nyg,
		       double **Yg,
		       int    **realY,
		       double **sin2_lat,
		       double **tan2_lat,
		       int     *jc_eq,
		       int      quiet);

void cart_to_spher (double *x,
		    double *dx,
		    double *r,
		    double *costheta,
		    double *sintheta,
		    double *e_r);

int move_satellite_photon_to_toa ( photon_struct     *p,
				   atmosphere_struct *atmos );

/* Functions for coherent backscattering CB and other stuff */
void add_coord_to_cohebasca ( photon_struct *p );
void clear_cohebasca (cohebasca_struct *cb);
void clear_pss (p_save_struct *pss);
void copy_pss (p_save_struct *pss, p_save_struct *pss_copy);
void copy_cohebasca (cohebasca_struct *cb, cohebasca_struct *cb_copy);
void calc_stokes_cb (cohebasca_struct *cb, double *stokes );
void add_tauext_to_cohebasca ( photon_struct     *p,
			       sample_struct     *sample,
			       lidar_cb_struct   *lidcb);
int calc_stokes_alt ( sample_struct     *sample,
		      atmosphere_struct *atmos,
		      photon_struct     *p,
                      elevation_struct  *elev,
		      float             *refind );
void calc_connection (double *S, double *M, double *n, double *MK, double *SKn);
int check_opan_cb ( photon_struct *p, photon_struct *p_escape, sample_struct *sample, atmosphere_struct *atmos );
double plane_wave_locest (photon_struct *p, sample_struct *sample);
int summarize_result_lidar_cb (mishchenko_cb_struct* mish, lidar_cb_struct* lidcb, sample_struct *sample, double ris_factor, 
                               long int nphotons, char* basename );
double float_zeros ( float f, int acc );
void write_sr_lidar_output_part(FILE* file, sample_struct* sample, double factor, int ili, int io, int ia, int nn);
int summarize_result_lidar (sample_struct* sample,
                            atmosphere_struct* atmos,
                            result_struct* result,
                            char* basename,
                            int write_files,
                            int write_output_as_netcdf,
                            double factor );
int write_mmclx_locest_file (char *mmclxlocfilename, sample_struct *sample,
                             result_struct *result, atmosphere_struct *atmos);

int write_netcdf_locest_file (char *netcdflocfilename, sample_struct *sample,
                              result_struct *result, atmosphere_struct *atmos,
                              double nphotons);

#endif
