/*--------------------------------------------------------------------
 * $Id: uvspec.c 3279 2017-07-07 20:35:28Z Claudia.Emde $
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

#include "uvspec.h"
#include "numeric.h"


/* prototypes of internal functions */
void pmesg (char *message, int verbose);



/**************************************************************/
/* The main uvspec function.                                  */
/**************************************************************/

int uvspec (input_struct input, output_struct *output)
{
  int status = 0;
  char function_name[]="uvspec";
  char file_name[]="uvspec.c";


  /**** Check the model input data ****/
  pmesg (" ... calling uvspec_check(), checking model input data\n", input.verbose);
  status = uvspec_check(input); /* in this file: uvspec.c */
  if (status!=0) {
    fprintf (stderr, "%d error(s) in uvspec input-file, aborting\n", abs(status));
    return status;
  }


  /**** Set wavelength grid for the transmittance calculation ****/
  pmesg (" ... calling setup_wlgrid(), generating transmittance wavelength grid\n", input.verbose);
  status = setup_wlgrid (input, output);  /* in ancillary.c */
  if (status!=0) {
    fprintf (stderr, "Error %d setting up transmittance wavelength grid in %s (%s)\n", status, function_name, file_name);
    return status;
  } 


  /**** Set wavelength grid for the radiative transfer calculation ****/
  pmesg (" ... calling setup_rte_wlgrid(), generating radiative transfer wavelength grid\n", input.verbose);
  status = setup_rte_wlgrid (input, output);  /* in ancillary.c */
  if (status!=0) {
    fprintf (stderr, "Error %d setting up transmittance wavelength grid in %s (%s)\n", status, function_name, file_name);
    return status;
  } 


  /**** Setup sza and phi0 ****/
  pmesg (" ... calling setup_sza(), generating model solar zenith and azimuth\n", input.verbose);
  status = setup_sza (input, output); /* in sza.c */
  if (status!=0) {
    fprintf (stderr, "Error %d setting up sza and phi0 in %s (%s)\n", status, function_name, file_name);
    return status;
  }


  /**** Independent pixel stuff ****/
  pmesg (" ... calling setup_ipa(), generating independent pixels\n", input.verbose);
  status = setup_ipa (input, output); /* in ipa.c */
  if (status!=0) {
    fprintf (stderr, "Error %d setting up IPA properties in %s (%s)\n", status, function_name, file_name);
    return status;
  }

  /**** 3D molecular atmosphere ****/
  pmesg(" ... setup 3D molecular atmosphere\n", input.verbose);
  setup_molecular3d(input,output); /* in molecular3d.c */
  if (status!=0) {
    fprintf (stderr, "Error %d setting up 3D molecular atmosphere in %s (%s)\n", status, function_name, file_name);
    return status;
  }
  
  
  /**** get surface altitude ****/
  pmesg (" ... calling setup_altitude(), determine surface elevation \n", input.verbose);
  status = setup_altitude (input, output); /* in atmosphere.c */
  if (status!=0) {
    fprintf (stderr, "Error %d determining surface altitude in %s (%s) \n", status, function_name, file_name);
    return status;
  }


  /**** Caoth stuff (1D) ****/
  pmesg (" ... calling setup_all_caoth(), generating all clouds etc.\n", input.verbose);
  status = setup_all_caoth (input, output); /* in cloud.c */ 
  if (status!=0) {
    fprintf (stderr, "Error %d setting up all cloud properties in %s (%s)\n", status, function_name, file_name);
    return status;
  }


  /**** Cloudless sky atmosphere ****/
  pmesg (" ... calling setup_atmosphere(), generating model atmosphere\n", input.verbose);
  status = setup_atmosphere (input, output); /* in atmosphere.c */
  if (status!=0) {
    fprintf (stderr, "Error %d setting up atmosphere in %s (%s)\n", status, function_name, file_name);
    return status;
  }  


  /**** Setup temperature at computing levels ****/
  pmesg (" ... calling setup_temperature(), generating model temperature profile\n", input.verbose);
  status = setup_temperature (input, output); /* in atmosphere.c */
  if (status!=0) {
    fprintf (stderr, "Error %d setting up temperature in %s (%s)\n", status, function_name, file_name);
    return status;
  }


  /**** Read cross section files, like O3, BrO, OclO, NO2, etc. ****/
  pmesg (" ... calling setup_crs(), reading cross section files\n", input.verbose);
  status = setup_crs (input, output); /* in molecular.c */
  if (status!=0) {
    fprintf (stderr, "Error %d setting up absorption cross sections in %s (%s)\n", status, function_name, file_name);
    return status;
  }


  /**** Calculate Rayleigh scattering cross section ****/
  pmesg (" ... calling setup_rayleigh(), calculating Rayleigh scattering\n", input.verbose);
  status = setup_rayleigh (input, output); /* in molecular.c */ 
  if (status!=0) {
    fprintf (stderr, "Error %d setting up Rayleigh cross sections in %s (%s)\n", status, function_name, file_name);
    return status;
  }


  /**** Calculate Raman scattering cross section ****/
  pmesg (" ... calling setup_raman(), calculating Raman scattering\n", input.verbose);
  status = setup_raman (input, output); /* in molecular.c */
  if (status!=0) {
    fprintf (stderr, "Error %d setting up Raman cross sections in %s (%s)\n", status, function_name, file_name);
    return status;
  }


  /****  Optical properties of trace gases ****/
  pmesg (" ... calling setup_gases(), generating optical properties of trace gases\n", input.verbose);
  status = setup_gases (input, output);  /* in atmosphere.c */
  if (status!=0) {
    fprintf (stderr, "Error %d setting up optical properties of trace gases in %s (%s)\n", status, function_name, file_name);
    return status;
  }


  /****   Aerosol stuff   ****/
  pmesg (" ... calling setup_aerosols(), generating aerosols\n", input.verbose);
  status = setup_aerosol (input, output); /* in aerosol.c */
  if (status!=0) {
    fprintf (stderr, "Error %d setting up aerosol properties in %s (%s)\n", status, function_name, file_name);
    return status;
  }

  /**** read 3D caoth for MYSTIC ****/
  pmesg (" ... calling setup_caoth3D(), generating 3D clouds etc.\n", input.verbose);
  status = setup_caoth3D (input, output); /* in cloud3d.c */
  if (status!=0) {
    fprintf (stderr, "Error %d setting up 3D clouds etc. in %s (%s)\n", status, function_name, file_name);
    return status;
  }

  /**** Redistribute profiles to required vertical resolution ****/
  pmesg (" ... calling setup_redistribute(), redistributing vertical profiles\n", input.verbose);
  status = setup_redistribute (input, output); /* in redistribute.c */
  if (status!=0) {
    fprintf (stderr, "Error %d doing vertical redistribution of optical properties in %s (%s)\n", status, function_name, file_name);
    return status;
  }
  
  /**** 3D atmosphere for MYSTIC - to be called after redistribution! ****/
  pmesg (" ... calling setup_atmosphere3D(), generating 3D atmosphere\n", input.verbose);
  status = setup_atmosphere3D (input, output); /* in cloud3d.c */
  if (status!=0) {
    fprintf (stderr, "Error %d setting up 3D atmosphere in %s (%s)\n", status, function_name, file_name);
    return status;
  }
  
  /**** Calculate 3D molecular optical properties ****/
  pmesg (" ... setup_optprop_molecular3D(), generating 3D gas absorption and Rayleigh scattering coefficients, include in caoth3d\n", input.verbose);
  status = setup_optprop_molecular3d (input, output); /* in molecular3d.c */
  if (status!=0) {
    fprintf (stderr, "Error %d setting up 3D optical properties for molecules %s (%s)\n", status, function_name, file_name);
    return status;
  }
  

  /****  Read extraterrestrial spectrum and correct for Earth-Sun distance ****/
  pmesg (" ... calling setup_extraterrestrial(), reading extraterrestrial spectrum\n", input.verbose);
  status = setup_extraterrestrial (input, output); /* in extraterrestrial.c */
  if (status!=0) {
    fprintf (stderr, "Error %d setting up extraterrestrial irradiance in %s (%s)\n", status, function_name, file_name);
    return status;
  }


  /**** Setup wavelength dependent albedo ****/
  pmesg (" ... calling setup_albedo(), generating surface albedo\n", input.verbose);
  status = setup_albedo (input, output); /* in albedo.c */
  if (status!=0) {
    fprintf (stderr, "Error %d setting up wavelength dependent albedo in %s (%s)\n", status, function_name, file_name);
    return status;
  }


  /**** Setup wavelength dependent fluorescence ****/
  pmesg (" ... calling setup_fluorescence(), generating surface fluorescence\n", input.verbose);
  status = setup_fluorescence (input, output); /* in fluorescence.c */
  if (status!=0) {
    fprintf (stderr, "Error %d setting up wavelength dependent fluorescence in %s (%s)\n", status, function_name, file_name);
    return status;
  }  

  /****     Solve RTE     ****/
  pmesg (" ... calling solve_rte(), solving RTE\n", input.verbose);
  status = solve_rte (input, output); /* in solve_rte.c */ 
  if (status!=0) {
    fprintf (stderr, "Error %d solving RTE in %s (%s)\n", status, function_name, file_name);
    return status;
  }


  /**** Distribute RTE results from internal wavelength grid to transmittance wavelength grid ****/
  pmesg (" ... calling internal_to_transmittance_grid(), \n", input.verbose);
  status = internal_to_transmittance_grid (input, output); /* in ancillary.c */
  if (status!=0) {
    fprintf (stderr, "Error %d distributing internal wavelength grid to transmittance wavelength grid in %s (%s)\n", status, function_name, file_name);
    return status;
  }


  /**** Interpolate transmittance to high resolution wavelength grid ****/
  pmesg (" ... calling interpolate_transmittance(), interpolating transmittance\n", input.verbose);
  status = interpolate_transmittance (input, output); /* in ancillary.c */
  if (status!=0) {
    fprintf (stderr, "Error %d interpolating model output to hires wavelength grid in %s (%s)\n", status, function_name, file_name);
    return status;
  }


  /**** Multiply with extraterrestrial irradiance ****/
  pmesg (" ... calling multiply_extraterrestrial()\n", input.verbose);
  status = multiply_extraterrestrial (input, output); /* in ancillary.c */
  if (status!=0) {
    fprintf (stderr, "Error %d multiplying with extraterrestrial irradiance in %s (%s)\n", status, function_name, file_name);
    return status;
  }


  /**** Convolve with slit function ****/
  if (input.convolve) {
    pmesg (" ... calling convolve(), convoluting with slit function\n", input.verbose);
    status = convolve (input, output); /* in ancillary.c */
    if (status!=0) {
      fprintf (stderr, "Error %d convolving hires spectrum with slit function in %s (%s)\n", status, function_name, file_name);
      return status;
    }
  }


  /**** Finally interpolate to output grid ****/
  if (input.spline) {
    pmesg (" ... calling spline_interpolate(), interpolating to output grid\n", input.verbose);
    status = spline_interpolate (input, output); /* in ancillary.c */
    if (status!=0) {
      fprintf (stderr, "Error %d interpolating spectrum to output grid in %s (%s)\n", status, function_name, file_name);
      return status;
    }
  }

  
  /**** Process 3D output (wavelength integration, ...) ****/
  if (output->mc.sample.passback3D) {
    pmesg (" ... calling processing3D(), post processing of 3D fields\n", input.verbose);
    status = processing3D (input, output); /* in ancillary.c */
    if (status!=0) {
      fprintf (stderr, "Error %d integrating 3D output over wavelength in %s (%s)\n", status, function_name, file_name);
      return status;
    }
  }
 

  /**** Process 1D output (wavelength integration, ...) ****/
  pmesg (" ... calling processing1D(), post processing\n", input.verbose);
  status = processing1D (input, output);  /* in ancillary.c */
  if (status!=0) {
    fprintf (stderr, "Error %d processing 1D output in %s (%s)\n", status, function_name, file_name);
    return status;
  }

  return 0;
}



/**************************************************************/
/* Check the model input for consistency.                     */
/**************************************************************/

int uvspec_check (input_struct Input)
{

#define EPSILON 1.E-6        /*relative rounding error for floats */

  int i=0, ichk=0, status=0;
  int s1=0, s2=0;
  FILE *f=NULL;
  int lev=0;
  float dz=0;

  char function_name[]="uvspec_check";
  char file_name[]="uvspec.c";

  /* ***FIXME*** Still need to check (also in uvspec function):
     aerosol_haze
     aerosol_scale_ssa
     aerosol_set_gg
     aerosol_set_ssa
     aerosol_season
     aerosol_visibility
     aerosol_vulcan
     angstrom
     atmosphere_file
     data_files_path
     nstr
     rte_solver
     solar_file
     sza
     sza_file
     umu
     wc_file
     wc_set_gg
     wc_set_ssa
     wc_layer
     wvn
   */

  status-= caoth_check(Input);

  /* check latitude and longitude */
  /* position inside normal longitude latitude boundaries ? */
  if ( fabs (Input.latitude - NOT_DEFINED_FLOAT)  > EPSILON && Input.latitude  < - 90.0) {
    fprintf (stderr, "Error in %s (%s), latitude is %7.2f < -90 S\n", function_name, file_name, Input.latitude);
    return -1;
  }
  if ( fabs (Input.latitude - NOT_DEFINED_FLOAT)  > EPSILON && Input.latitude  >   90.0) {
    fprintf (stderr, "Error in %s (%s), latitude is %7.2f > +90 N\n", function_name, file_name, Input.latitude);
    return -1;
  }
  if ( fabs (Input.longitude - NOT_DEFINED_FLOAT) > EPSILON && Input.longitude < -360.0) {
    fprintf (stderr, "Error in %s (%s), longitude is %7.2f < -360 W\n", function_name, file_name, Input.longitude);
    return -1;
  }
  if ( fabs (Input.longitude - NOT_DEFINED_FLOAT) > EPSILON && Input.longitude >  360.0) {
    fprintf (stderr, "Error in %s (%s), longitude is %7.2f > +360 E\n", function_name, file_name, Input.longitude);
    return -1;
  }

  /* time check */
  if ( Input.UTC.tm_year != NOT_DEFINED_INTEGER && (Input.UTC.tm_year < 1902.-1900. || 2038.-1900. <= Input.UTC.tm_year) ) {
    fprintf (stderr, "Error, sorry the mktime-function of C does only work for the time range\n");
    fprintf (stderr, "       1902 ... 2038. If anybody has a function with broader range, \n");
    fprintf (stderr, "       please email it to the libRadtran programmers. Thanks. \n");
    return -1;
  } 

  /* Albedo checks */
  /* Albedo value range for direct input */
  if (Input.alb.albedo > 1.0) {
    fprintf (stderr, "Error in %s (%s), surface albedo %f > 1 \n", function_name, file_name, Input.alb.albedo);
    status--;
  }
  else if (Input.alb.albedo < 0.0 && Input.alb.source == ALBEDO_CONSTANT) {
    fprintf (stderr, "Error in %s (%s), surface albedo %f < 0\n", function_name, file_name, Input.alb.albedo);
    status--;
  }

  if ( Input.alb.surf_type_map == TRUE )
    if ( Input.alb.source != ALBEDO_USER_LIBRARY && Input.alb.source != ALBEDO_IGBP_LIBRARY) {
      fprintf (stderr, "Error, surface_type_map was specified, but also another albedo information in the input file.\n");
      switch(Input.alb.source) {
      case ALBEDO_CONSTANT:
        fprintf (stderr, "       constant albedo\n");
        break;
      case ALBEDO_FROM_ALBEDO_MAP:
        fprintf (stderr, "       albedo_map \n");
        break;
      case ALBEDO_FROM_EMISSIVITY_MAP:
        fprintf (stderr, "       emissivity_map \n");
        break;
      case ALBEDO_FROM_ALBEDO_FILE:
        fprintf (stderr, "       albedo_file \n");
        break;
      case ALBEDO_USER_LIBRARY:
      case ALBEDO_IGBP_LIBRARY:
        /* does never happen */
        break;
      default:
        fprintf (stderr, "Error: unsupported albedo source in %s (%s)\n", function_name, file_name);
        return -1;
      }
      fprintf (stderr, "       please specify only one albedo information in the input file! \n");
      status--;
    }

  if (Input.alb.source == ALBEDO_FROM_EMISSIVITY_MAP && Input.source != SRC_THERMAL) {
    fprintf (stderr, "Error, emissivity_map only possible for thermal simulations in %s (%s)\n", function_name, file_name);
    fprintf (stderr, "       please use albedo_map instead! \n");
    status--;
  }

  /* Day of year checks */

  if (Input.UTC.tm_yday > 366) {
    fprintf (stderr, "Error in %s (%s), day of year %d > 366\n", function_name, file_name, Input.UTC.tm_yday);
    status--;
  }
  else if (Input.UTC.tm_yday < 0 && Input.UTC.tm_yday != -999) {
    fprintf (stderr, "Error in %s (%s), day of year %d < 0 \n", function_name, file_name, Input.UTC.tm_yday);
    status--;
  }

  /* ECMWF data + time checks */
  if (Input.UTC.tm_mday <= 0) {
    ichk=0;
    if (Input.atm.rs_source == RS_FROM_ECMWF)
      ichk=1;
    if (Input.i_wc!=-1)
      if ( Input.caoth[Input.i_wc].source == CAOTH_FROM_ECMWF )
	ichk=1;
    if (Input.i_ic!=-1)
      if ( Input.caoth[Input.i_ic].source == CAOTH_FROM_ECMWF )
	ichk=1;
    if (ichk==1) {
      fprintf (stderr, "Error, if you use ECMWF data, it is necessary to specify also the time;\n");
      fprintf (stderr, "       please use the 'time' option in the input file\n");
      status--;
    }
  }

  /* mol_modify vs. mixing_ratio checks */

  if (Input.atm.column[MOL_O2] >= 0.0 && Input.mixing_ratio[MX_O2] >= 0.0) {
    fprintf (stderr, "Error, either 'mol_modify O2' or 'mixing_ratio O2' can be specified\n");
    fprintf (stderr, "in the input file but not both. Please remove one of these lines!\n");
    return -1;
  }

  if (Input.atm.column[MOL_H2O] >= 0.0 && Input.mixing_ratio[MX_H2O] >= 0.0) {
    fprintf (stderr, "Error, either 'mol_modify H2O' or 'mixing_ratio H2O' can be specified\n");
    fprintf (stderr, "in the input file but not both. Please remove one of these lines!\n");
    return -1;
  }

  if (Input.atm.column[MOL_CO2] >= 0.0 && Input.mixing_ratio[MX_CO2] >= 0.0) {
    fprintf (stderr, "Error, either 'mol_modify CO2' or 'mixing_ratio CO2' can be specified\n");
    fprintf (stderr, "in the input file but not both. Please remove one of these lines!\n");
    return -1;
  }

  if (Input.atm.column[MOL_NO2] >= 0.0 && Input.mixing_ratio[MX_NO2] >= 0.0) {
    fprintf (stderr, "Error, either 'mol_modify NO2' or 'mixing_ratio NO2' can be specified\n");
    fprintf (stderr, "in the input file but not both. Please remove one of these lines!\n");
    return -1;
  }

  if (Input.atm.column[MOL_CH4] >= 0.0 && Input.mixing_ratio[MX_CH4] >= 0.0) {
    fprintf (stderr, "Error, either 'mol_modify CH4' or 'mixing_ratio CH4' can be specified\n");
    fprintf (stderr, "in the input file but not both. Please remove one of these lines!\n");
    return -1;
  }

  if (Input.atm.column[MOL_N2O] >= 0.0 && Input.mixing_ratio[MX_N2O] >= 0.0) {
    fprintf (stderr, "Error, either 'mol_modify N2O' or 'mixing_ratio N2O' can be specified\n");
    fprintf (stderr, "in the input file but not both. Please remove one of these lines!\n");
    return -1;
  }


  /* Aerosol checks */

  /* if any of the following four variables is set, all four must be set; */
  /* otherwise the information is incomplete!                             */
  if (Input.aer.seasn>0 || Input.aer.vulcan>0 || Input.aer.haze>0 || Input.aer.visibility>0) {

    if (Input.aer.seasn<=0) {
      fprintf (stderr, "Error in %s (%s), aerosol information incomplete; specify aerosol_season or aerosol_default!\n", function_name, file_name);
      status--;
    }

    if (Input.aer.vulcan<=0) {
      fprintf (stderr, "Error in %s (%s), aerosol information incomplete; specify aerosol_vulcan or aerosol_default!\n", function_name, file_name);
      status--;
    }

    if (Input.aer.haze<=0) {
      fprintf (stderr, "Error in %s (%s), aerosol information incomplete; specify aerosol_haze or aerosol_default!\n", function_name, file_name);
      status--;
    }

    if (Input.aer.visibility<=0) {
      fprintf (stderr, "Error in %s (%s), aerosol information incomplete; specify aerosol_visibility or aerosol_default!\n", function_name, file_name);
      status--;
    }
  }
  
  /* aerosol_default does not work with polarization */
  if (Input.aer.seasn>0 && Input.aer.n_species == NOT_DEFINED_INTEGER  && 
      ((Input.rte.solver == SOLVER_POLRADTRAN &&
        Input.rte.polradtran[POLRADTRAN_NSTOKES]>1) || 
       Input.rte.mc.polarisation == 1 )) {
    fprintf (stderr, "Error in %s (%s), default aerosol is not implemented for polarized radiance calculations!\nPlease use option aerosol_species_file (e.g., OPAC aerosol).\n", function_name, file_name);
    status--;
    }


  /* if no standard aerosol is defined, we don't want any of the aerosol varibles set */
  if (!Input.aer.standard && Input.aer.spec) {  
    fprintf (stderr, "Error in %s (%s), aerosol information incomplete; specify aerosol_default or each\n", function_name, file_name);
    fprintf (stderr, "of aerosol_vulcan, aerosol_haze, and aerosol_visibility!\n"); 
    status--;
  }


  if ((Input.aer.seasn < 0 && Input.aer.seasn!=-999) || Input.aer.seasn >= AEROSOL_SEASON_NN) {
    fprintf (stderr, "Error in %s (%s), aerosol season out of range.\n", function_name, file_name);
    status--;
  }

  if ((Input.aer.haze < 0 && Input.aer.haze!=-999) || Input.aer.haze >= AEROSOL_HAZE_NN) {
    fprintf (stderr, "Error in %s (%s), aerosol haze out of range.\n", function_name, file_name);
    status--;
  }

  if ((Input.aer.vulcan < 0 && Input.aer.vulcan!=-999) || Input.aer.vulcan >= AEROSOL_VULCAN_NN) {
    fprintf (stderr, "Error in %s (%s), aerosol vulcan out of range.\n", function_name, file_name);
    status--;
  }
    

  s1=strlen(Input.aer.filename[FN_AER_SIZ]);
  s2=strlen(Input.aer.filename[FN_AER_REF]);

  if (s1>0 && (s2==0 && Input.aer.re==0)) {
    fprintf (stderr, "Error in %s (%s), aerosol information incomplete; specify aerosol_refrac_index \n", function_name, file_name);
    fprintf (stderr, "or aerosol_refrac_file or omit the size distribution\n");
    status--;
  }

  if (s1==0 && (s2>0 || Input.aer.re>0)) {
    fprintf (stderr, "Error in %s (%s), aerosol information incomplete; specify aerosol_sizedist_file\n", function_name, file_name);
    fprintf (stderr, "or omit the refractive index!\n");
    status--;
  }


  /* RTE checks */

  switch (Input.rte.solver) {
  case SOLVER_SPSDISORT:
    if (Input.source==SRC_THERMAL) {
      fprintf (stderr, "Error, spsdisort does not include thermal emission\n");
      status--;
    }

    break;

  case SOLVER_SDISORT:

    if (Input.source==SRC_THERMAL) {
      fprintf (stderr, "Error, sdisort does not include thermal emission\n");
      status--;
    }

    if (Input.rte.sdisort[SDISORT_NSCAT] < 0 || Input.rte.sdisort[SDISORT_NSCAT] > 2) {
      fprintf (stderr, "Error, sdisort nscat %d is out of bounds.\n", Input.rte.sdisort[SDISORT_NSCAT]);
      fprintf (stderr, "Only 1 or 2 is allowed in combination with the sdisort solver.\n");
      status--;
    }

    if (Input.rte.sdisort[SDISORT_NREFRAC] < 0 || Input.rte.sdisort[SDISORT_NREFRAC] > 2) {
      fprintf (stderr, "Error, sdisort nrefrac = %d is out of bounds\n", Input.rte.sdisort[SDISORT_NREFRAC]);
      status--;
    }
    
    break;

  case SOLVER_TWOSTR:
  case SOLVER_FTWOSTR:

    if (Input.rte.nprndis>0)
      for (i=0;i<Input.rte.nprndis;i++)
	if (!(Input.rte.prndis[i] == 1 || Input.rte.prndis[i] == 2)) { 
	  fprintf (stderr, "Error, twostr can only handle prndis 1 and/or 2.\n");
	  status--;
	}

    break;

  case SOLVER_RODENTS:
  case SOLVER_TWOSTREBE:
  case SOLVER_TWOMAXRND:
  case SOLVER_SSLIDAR:
  case SOLVER_FDISORT1:
  case SOLVER_SOS:
  case SOLVER_MONTECARLO:
  case SOLVER_POLRADTRAN:
  case SOLVER_FDISORT2:
  case SOLVER_TZS:
  case SOLVER_SSS:
  case SOLVER_SSSI:
  case SOLVER_NULL:
  case SOLVER_DISORT:
    break;

  default:
    fprintf (stderr, "Error, unknown rte_solver %d in uvspec_check\n", Input.rte.solver);
    return -1;
  }


  if (Input.rte.mc.absorption==MCFORWARD_ABS_EMISSION && Input.source != SRC_THERMAL) {
    fprintf (stderr, "Error, mc_emission does only make sense with \"source thermal\"\n");
    status--;
  }

  /* either specify a sensor direction or calculate surface-parallel */
  if (Input.rte.mc.surfaceparallel && Input.rte.mc.sensordirection) {
    fprintf (stderr, "Error, mc_surfaceparallel together with  mc_sensordirection\n");
    fprintf (stderr, "does not make sense. Please remove one of them!\n");
    status--;
  }

  /* sensor direction does only work in backward mode */
  if (Input.rte.mc.sensordirection && !Input.rte.mc.backward.yes) {
    fprintf (stderr, "Error, mc_sensordirection only possible in backward mode!\n");
    status--;
  }

  /* sensor position does only work in backward mode */
  if (Input.rte.mc.sensorposition && !Input.rte.mc.backward.yes) {
    fprintf (stderr, "Error, mc_sensorposition only possible in backward mode!\n");
    status--;
  }

  /* check if absorption/emission was consistently defined in forward/backward mode; */
  /* if both are defined, unneccesarily large amounts of memory might be allocated.  */
  if (Input.rte.mc.absorption && Input.rte.mc.backward.absorption) {
    fprintf (stderr, "Error, specified both forward and backward calculation\n");
    fprintf (stderr, "of absorption/emission/actinic flux. Please choose either\n");
    fprintf (stderr, "forward (mc_absorption, ...) or backward (mc_backward_output abs).\n");
    fprintf (stderr, "The programmers are sorry for any inconvenience this might have caused ...\n");
    status--;
  }
    
  /* importance sampling checks */
  if (Input.rte.mc.spectral_is || Input.rte.mc.concentration_is){
    if(Input.source == SRC_THERMAL && !Input.rte.mc.backward.yes){
      fprintf(stderr, "Error, thermal calculations with importance sampling ('mc_spectral_is' and 'mc_aeris') work only \n");
      fprintf(stderr, "       in backward mode. Please specify option 'mc_backward'. \n");
    }
  }

  /* check if z_atm_forced are sorted in ascending (descending) order */
  /* the order of the levels is reversed just direct after the input  */
  if (Input.atm.nz_atm_forced_sea != NOT_DEFINED_INTEGER) {
    for (i=1; i<Input.atm.nz_atm_forced_sea; i++) {
      if (Input.atm.z_atm_forced_sea[i] - Input.atm.z_atm_forced_sea[i-1] >= 0 ) {
        fprintf (stderr, "Error, 'atmosphere_zgrid' levels not sorted in ascending order\n");
        fprintf (stderr, "       checking atmosphere_zgrid-levels in %s (%s)\n", function_name, file_name);
        for (lev=0; lev<Input.atm.nz_atm_forced_sea; lev++)
          fprintf (stderr, "    zgrid[%3d] = %8.3f\n",lev,Input.atm.z_atm_forced_sea[Input.atm.nz_atm_forced_sea-1-lev]);
        status--;
      }
    }
  }

  /* forbidden combination z_atm_forced and altitude with 2 arguments */
  if (Input.atm.nz_atm_forced_sea != NOT_DEFINED_INTEGER && Input.alt.altitude_dz > 0.0 ) {
    fprintf (stderr, "Error, combination 'atmosphere_zgrid' and 'altitude with second argument' is not allowed\n");
    status--;
  }
  /* forbidden combination z_atm_forced and molecular_tau_file */
  if (Input.atm.nz_atm_forced_sea != NOT_DEFINED_INTEGER && strlen(Input.filename[FN_MOL_TAU_ABS]) > 0) {
    fprintf (stderr, "Error, combination 'atmosphere_zgrid' and 'molecular_tau_file' is not allowed\n");
    status--;
  }
  /* forbidden combination z_atm_forced and rayleigh_tau_file */
  if (Input.atm.nz_atm_forced_sea != NOT_DEFINED_INTEGER && strlen(Input.filename[FN_MOL_TAU_SCA]) > 0) {
    fprintf (stderr, "Error, combination 'atmosphere_zgrid' and 'rayleigh_tau_file' is not allowed\n");
    status--;
  }

  /* check if z_atm_forced are sorted in ascending (descending) order */
  /* the order of the levels is reversed just direct after the input  */
  if (Input.atm.nz_atm_forced_sea != NOT_DEFINED_INTEGER) {
    for (i=1; i<Input.atm.nz_atm_forced_sea; i++) {
      if (Input.atm.z_atm_forced_sea[i] - Input.atm.z_atm_forced_sea[i-1] >= 0 ) {
        fprintf (stderr, "Error, 'atmosphere_zgrid' levels not sorted in ascending order\n");
        fprintf (stderr, "       checking atmosphere_zgrid-levels in %s (%s)\n", function_name, file_name);
        for (lev=0; lev<Input.atm.nz_atm_forced_sea; lev++)
          fprintf (stderr, "    zgrid[%3d] = %8.3f\n",lev,Input.atm.z_atm_forced_sea[Input.atm.nz_atm_forced_sea-1-lev]);
        status--;
      }
    }
  }

  /* forbidden combination z_atm_forced and altitude with 2 arguments */
  if (Input.atm.nz_atm_forced_sea != NOT_DEFINED_INTEGER && Input.alt.altitude_dz > 0.0 ) {
    fprintf (stderr, "Error, combination 'atmosphere_zgrid' and 'altitude with second argument' is not allowed\n");
    status--;
  }
  /* forbidden combination z_atm_forced and molecular_tau_file */
  if (Input.atm.nz_atm_forced_sea != NOT_DEFINED_INTEGER && strlen(Input.filename[FN_MOL_TAU_ABS]) > 0) {
    fprintf (stderr, "Error, combination 'atmosphere_zgrid' and 'molecular_tau_file' is not allowed\n");
    status--;
  }
  /* forbidden combination z_atm_forced and rayleigh_tau_file */
  if (Input.atm.nz_atm_forced_sea != NOT_DEFINED_INTEGER && strlen(Input.filename[FN_MOL_TAU_SCA]) > 0) {
    fprintf (stderr, "Error, combination 'atmosphere_zgrid' and 'rayleigh_tau_file' is not allowed\n");
    status--;
  }

  /* check if zout is sorted in ascending order */
  switch (Input.atm.zout_source) {
  case OUTLEVEL_ZOUT_ABOVE_SUR:
    for (i=1; i<Input.atm.nzout; i++) {
      if (Input.atm.zout_sur[i] - Input.atm.zout_sur[i-1] <= 0 ) {
        fprintf (stderr, "Error, zout (above surface) not sorted in ascending order\n");
        fprintf (stderr, "       checking zout_sur-level in %s (%s)\n", function_name, file_name);
        for (lev=0; lev<Input.atm.nzout; lev++)
          fprintf (stderr, "    zout[%3d] = %8.3f\n",lev,Input.atm.zout_sur[lev]);
        status--;
      }
    }
    break;
  case OUTLEVEL_ZOUT_ABOVE_SEA:
    for (i=1; i<Input.atm.nzout; i++) {
      if (Input.atm.zout_sea[i] - Input.atm.zout_sea[i-1] <= 0 ) {
        fprintf (stderr, "Error, zout (above sea level) not sorted in ascending order\n");
        fprintf (stderr, "       checking zout_sea-level in %s (%s)\n", function_name, file_name);
        for (lev=0; lev<Input.atm.nzout; lev++)
          fprintf (stderr, "    zout[%3d] = %8.3f\n",lev,Input.atm.zout_sea[lev]);
        status--;
      }
    }
    break;
  case OUTLEVEL_PRESS:
    for (i=1; i<Input.atm.nzout; i++) {
      if ( Input.atm.press_zout[i] - Input.atm.press_zout[i-1] >= 0
           && Input.atm.press_zout[i-1] != ZOUT_SURFACE
           && Input.atm.press_zout[i]   != ZOUT_TOA
           && !( Input.atm.press_zout[i] == ZOUT_CPT && i==Input.atm.nzout-1 ) ) {
        if ( (Input.atm.press_zout[i] == ZOUT_CPT || Input.atm.press_zout[i-1] == ZOUT_CPT) && i!=Input.atm.nzout-1) {
          fprintf (stderr, "Error, outputlevel 'CPT' must be the last level.\n");
        }
        else {
          fprintf (stderr, "Error, outputlevel 'pressure_out' not sorted in descending order %f > %f \n",
                                                    Input.atm.press_zout[i],Input.atm.press_zout[i-1] );
        }
        fprintf (stderr, "       checking pressure_out-level in %s (%s)\n", function_name, file_name);
        for (lev=0; lev<Input.atm.nzout; lev++)
          fprintf (stderr, "    pressure_out[%3d] = %8.3f\n",lev,Input.atm.press_zout[lev]);
        status--;
      }
    }
    break;
  case OUTLEVEL_ATM_LEVELS:
  case OUTLEVEL_ALL_LEVELS:
  case OUTLEVEL_MODEL_LEVELS:
  case OUTLEVEL_MODEL_LAYERS:
  case OUTLEVEL_MODEL_LEVELS_AND_LAYERS:
    /* no checks here */
    break;
  default:
    fprintf (stderr, "Error, unknown output level type %d (function %s in %s) \n", Input.atm.zout_source, __func__, __FILE__ );
    fprintf (stderr, "This is a program bug, please report this bug to the programmers and include the input file. Thanx.\n");
    return -1;
  }

  /* transmittance_wl_file and thermal_bands_file are mutually exclusive */
  if ((strlen(Input.filename[FN_WLTRANS])>0) &&
      (strlen(Input.filename[FN_WLBANDS])>0)) {
    fprintf (stderr, "Error, it does not make sense to specify both transmittance_wl_file and thermal_bands_file!\n");
    status--;
  }

  /* solar_file and thermal_bands_file are mutually exclusive */
  if ((strlen(Input.filename[FN_EXTRATERRESTRIAL])>0) &&
      (strlen(Input.filename[FN_WLBANDS])>0)) {
    fprintf (stderr, "Error, it does not make sense to specify both\n");
    fprintf (stderr, "       solar_file and thermal_bands_file!\n");
    status--;
  }

  /* Need to define a wavelength range if no correlated_k and no transmittance_wl_file */
  if ((Input.wl.start<0 || Input.wl.end<0) && 
      (Input.ck_scheme == CK_CRS || Input.ck_scheme == CK_RAMAN || Input.ck_scheme == CK_LOWTRAN) &&
      (strlen(Input.filename[FN_WLTRANS])==0) &&
      (strlen(Input.filename[FN_WLBANDS])==0) &&
      (strlen(Input.filename[FN_MOL_TAU_ABS])==0)) {
    fprintf (stderr, "Error, need to define a wavelength range with 'wavelength lower upper'!\n");
    status--;
  }

  /* wvn_index only makes sense with explicitely defined wavelength indices */
  if ((Input.wl.start_index>0 || Input.wl.end_index>0) && 
      ((Input.ck_scheme == CK_CRS || Input.ck_scheme == CK_RAMAN || Input.ck_scheme == CK_LOWTRAN) &&
       (strlen(Input.filename[FN_WLTRANS])==0) &&
       (strlen(Input.filename[FN_WLBANDS])==0) &&
       (strlen(Input.filename[FN_MOL_TAU_ABS])==0))) {
    fprintf (stderr, "Error, wavelength_index can only be used together with\n");
    fprintf (stderr, "explicitely defined wavelength grids!\n");
    status--;
  }

  /* wvn and wvn_index are mutually exclusive */
  if ((Input.wl.start_index>0 || Input.wl.end_index>0) && 
      (Input.wl.start>0 || Input.wl.end>0)) {
    fprintf (stderr, "Error, it does not make sense to define both 'wavelength' and 'wavelength_index'!\n");
    status--;
  }

  if (Input.spline && !(strlen(Input.filename[FN_SPLINE])>0)) {
    if (Input.wl.start>0 && Input.wl.start > Input.spline_lambda_0)  {
      fprintf (stderr, "Error, spline wavelength %f is smaller than wavelength %f\n",
	       Input.spline_lambda_0, Input.wl.start);
      status--;
    }
    if (Input.wl.end>0 && Input.wl.end < Input.spline_lambda_1)  {
      fprintf (stderr, "Error, spline wavelength %f is greater than wavelength %f\n",
	       Input.spline_lambda_1, Input.wl.end);
      status--;
    }
  }

  /* Check if slit_function file exists */
  if (Input.convolve) {
    if ((f = fopen(Input.filename[FN_SLITFUNCTION], "r")) == NULL) {
      fprintf (stderr, "Error, cannot read slit function file %s\n", Input.filename[FN_SLITFUNCTION]);
      status--;
    }
    else 
      fclose(f);
  }

  
  /* ck and libnetcdf */ 
#if HAVE_LIBNETCDF
#else
  if (Input.ck_scheme == CK_KATO || Input.ck_scheme == CK_KATO2 || Input.ck_scheme == CK_KATO2_96 || 
      Input.ck_scheme == CK_REPTRAN) {
      /*Input.ck_scheme == CK_GENERIC || Input.ck_scheme == CK_REPTRAN) {*/
    fprintf (stderr, "***********************************************************************\n");
    fprintf (stderr, "* Error, you have built uvspec without libnetcdf and hence may not    *\n");
    fprintf (stderr, "* use the 'mol_abs_param kato' or 'mol_abs_param reptran' options.    *\n");
    fprintf (stderr, "* Please get netcdf and rebuild.                                      *\n");
    fprintf (stderr, "***********************************************************************\n");
    status--;
  }
#endif


  /* The MYSTIC solver */ 
#if HAVE_MYSTIC
#else
  if (Input.rte.solver == SOLVER_MONTECARLO) {
    fprintf (stderr, "*********************************************************************\n");
    fprintf (stderr, "* Error, you have built uvspec without the three-dimensional MYSTIC *\n");
    fprintf (stderr, "* module. MYSTIC is currently not released to the public. If you    *\n");
    fprintf (stderr, "* have an interest in three-dimensional calculations and if you     *\n");
    fprintf (stderr, "* have a reasonable project in mind, please contact                 *\n");
    fprintf (stderr, "* bernhard.mayer@dlr.de.                                            *\n");
    fprintf (stderr, "*********************************************************************\n");
    status--;
  }
#endif


  /* Fu and Liou */ 
#if HAVE_FULIOU
#else
  if (Input.ck_scheme == CK_FU) {
    fprintf (stderr, "**********************************************************************\n");
    fprintf (stderr, "* Error, you have built uvspec without support for the Fu and Liou   *\n");
    fprintf (stderr, "* parameterization. If you are interested, please contact the        *\n");
    fprintf (stderr, "* libRadtran authors.                                                *\n");
    fprintf (stderr, "**********************************************************************\n");
    status--;
  }
#endif
  

  /* LOWTRAN/SBDART */ 
#if HAVE_LOWTRAN
  if (Input.ck_scheme == CK_LOWTRAN && !Input.quiet) {
    fprintf (stderr, "*************************************************************************\n");
    fprintf (stderr, "* When using the LOWTRAN/SBDART gas absorption for any result that goes *\n");
    fprintf (stderr, "* into a publication, please don't forget to add a reference like:      *\n");
    fprintf (stderr, "* \"Molecular absorption was parameterized with the LOWTRAN band        *\n");
    fprintf (stderr, "*   model (Pierluissi and Peng, 1985), as adopted from the SBDART code  *\n");
    fprintf (stderr, "*   (Ricchiazzi et al., 1998).                                          *\n");
    fprintf (stderr, "*************************************************************************\n");
  }
#else
  if (Input.ck_scheme == CK_LOWTRAN) {
    fprintf (stderr, "***********************************************************************\n");
    fprintf (stderr, "* Error, you have built uvspec without support for the LOWTRAN/SBDART *\n");
    fprintf (stderr, "* gas parameterization. If you are interested, please contact the     *\n");
    fprintf (stderr, "* libRadtran authors.                                                 *\n");
    fprintf (stderr, "***********************************************************************\n");
    status--;
  }
#endif

  /* consistent?  output_unit (W/(m2 nm) or W/(m2 cm-1) or W/(m2 band)) and processing ('output sum' or 'output integrate')*/
  switch(Input.output_unit) {
  case UNIT_PER_NM:
    if (Input.processing == PROCESS_SUM && !Input.quiet) {
      fprintf (stderr, "*** Warning, the result of the radiative transfer calculation is given\n");
      fprintf (stderr, "*** per nanometer: W/(m**2 nm) because 'output_process per_nm' was chosen.\n");
      fprintf (stderr, "*** In combination with 'output_process sum', the sum of the spectral results is only equal\n");
      fprintf (stderr, "*** to the integrated result, if the wavelength calculation grid \n");
      fprintf (stderr, "*** (defined by the solar_file) has 1nm steps.\n\n");
    }
    break;
  case UNIT_PER_CM_1:
    if (Input.processing == PROCESS_SUM && !Input.quiet) {
      fprintf (stderr, "*** Warning, the result of the radiative transfer calculation is given\n");
      fprintf (stderr, "*** per wavenumber interval: W/(m**2 cm-1) because 'output_process per_cm-1' was chosen.\n");
      fprintf (stderr, "*** In combination with 'output_process sum', the sum of the spectral results is only equal\n");
      fprintf (stderr, "*** to the integrated result, if the wavelength calculation grid (defined by wavelength_grid_file)\n");
      fprintf (stderr, "*** has 1cm-1 steps and the bandwidth of each calculation is 1cm-1 (defined by thermal_bandwidth)\n");
    }
    break;
  case UNIT_PER_BAND:
    if (Input.processing == PROCESS_INT) {
      fprintf (stderr, "Error, the result of the radiative transfer calculation is given\n");
      fprintf (stderr, "per wavelength band: W/(m2 band) because 'output_process per_band' was chosen.\n");
      fprintf (stderr, "This does not make sense in combination with 'output_process integrate'!\n");
      fprintf (stderr, "Please use 'output_process sum' instead of 'output_process integrate' in the input file.\n");
      return -1;
    }
    break;
  case UNIT_NOT_DEFINED:
    /* not defined, that's OK, the user has to know what he is doing */
    break;
  default:
    fprintf (stderr, "Error, unknown output option %d. (function %s in %s) \n", Input.output_unit, __func__, __FILE__ );
    fprintf (stderr, "This is a program bug, please contact the programmers!\n");
    return -1;
  }

  /* Output checks */
  if (Input.heating != HEAT_NONE) {
    switch (Input.rte.solver) {
    case SOLVER_NULL:
    case SOLVER_FTWOSTR:
    case SOLVER_TWOSTR:
    case SOLVER_RODENTS:
    case SOLVER_TWOSTREBE:
    case SOLVER_TWOMAXRND:
    case SOLVER_FDISORT1:
    case SOLVER_DISORT:
    case SOLVER_SDISORT:
    case SOLVER_FDISORT2:
    case SOLVER_SPSDISORT:
    case SOLVER_MONTECARLO:
      break;
    case SOLVER_POLRADTRAN:
      if (Input.heating == HEAT_LOCAL) {
        fprintf (stderr, "Error, 'heating_rate local' cannot be used with polradtran.\n");
        fprintf (stderr, "Please choose heating_rate layer_cd, layer_fd, or a different rte_solver.\n");
        status--;
      }
      break;
    default:
      fprintf (stderr, "Error, the chosen solver %d does not provide heating rates.\n", Input.rte.solver);
      status--;
    }

    if ( Input.heating == HEAT_LOCAL &&  Input.rte.solver == SOLVER_MONTECARLO ) {
      fprintf (stderr, "Error, combination 'heating_rate local' and 'rte_solver MONTECARLO' \n");
      fprintf (stderr, "       is not allowed, as local heating rates require additional z-levels.\n");
      fprintf (stderr, "       Please use: 'mc_absorption K_per_day' instead \n");
      status--;
    }

    if (Input.heating == HEAT_LAYER_CD || Input.heating == HEAT_LAYER_FD) {
      if ( Input.atm.nzout < 2 && 
           Input.atm.zout_source != OUTLEVEL_ATM_LEVELS && 
           Input.atm.zout_source != OUTLEVEL_ALL_LEVELS &&
           Input.atm.zout_source != OUTLEVEL_MODEL_LEVELS &&
           Input.atm.zout_source != OUTLEVEL_MODEL_LAYERS &&
           Input.atm.zout_source != OUTLEVEL_MODEL_LEVELS_AND_LAYERS ) {
        fprintf (stderr, "Error, to calculate heating rates, at least two output levels are required.\n");
        fprintf (stderr, "Please use the zout command to specify them!\n");
        status--;
      }
    }
    if (Input.heating == HEAT_LAYER_CD) {

      switch (Input.atm.zout_source) {
      case OUTLEVEL_ZOUT_ABOVE_SUR:
        dz = fabs(Input.atm.zout_sur[1]-Input.atm.zout_sur[0]);
        for (lev=1; lev<Input.atm.nzout-1; lev++)      {
	  if (fabs(dz - (Input.atm.zout_sur[lev+1]-Input.atm.zout_sur[lev]))/dz > 0.1 && !Input.quiet) {
	    fprintf (stderr, "*** Warning, zout layers should be equally thick;\n");
            fprintf (stderr, "*** otherwise output might be inaccurate!\n"); 
            break;
          }
	}
	break;
      case OUTLEVEL_ZOUT_ABOVE_SEA:
        dz = fabs(Input.atm.zout_sea[1]-Input.atm.zout_sea[0]);
        for (lev=1; lev<Input.atm.nzout-1; lev++)      {
	  if (fabs(dz - (Input.atm.zout_sea[lev+1]-Input.atm.zout_sea[lev]))/dz > 0.1 && !Input.quiet) {
	    fprintf (stderr, "*** Warning, zout layers should be equally thick;\n");
            fprintf (stderr, "*** otherwise output might be inaccurate!\n"); 
            break;
          }
	}
	break;
      case OUTLEVEL_PRESS:
      case OUTLEVEL_ATM_LEVELS:
      case OUTLEVEL_ALL_LEVELS:
      case OUTLEVEL_MODEL_LEVELS:
      case OUTLEVEL_MODEL_LAYERS:
      case OUTLEVEL_MODEL_LEVELS_AND_LAYERS:
        /* no warning in this case, as it is not possible at this moment jet */
	break;
      default:
        fprintf (stderr, "Error, unknown output level type %d. (function %s in %s) \n", Input.atm.zout_source, __func__, __FILE__);
        fprintf (stderr, "This is a program bug, please report this bug to the programmers and include the input file. Thanx.\n");
        return -1;
      }
    }
  }

  if (Input.rte.mc.backward.yes && Input.rte.solver == SOLVER_MONTECARLO)
    if (Input.rte.mc.absorption==MCFORWARD_ABS_HEATING && Input.source == SRC_SOLAR) {
      fprintf (stderr, "Error, solar heating rates not implemented in backward Monte Carlo mode.\n");
	fprintf (stderr, "This could be done but probably forward is more efficient than backward\n");
	fprintf (stderr, "because photons usually contribute to many layers (and heating rates are\n");
	fprintf (stderr, "usually calculated for many layers rather than only one).\n");
	status--;
    }
  




  /* can only simulate wind output, if there is a wind map defined */
  for ( i=0; i<Input.n_output_user; i++ ) 
    if (Input.output_user[i]==OUTPUT_USER_WIND_U || Input.output_user[i]==OUTPUT_USER_WIND_V || Input.output_user[i]==OUTPUT_USER_WIND_W || 
        Input.output_user[i]==OUTPUT_USER_HEAT_AD_X || Input.output_user[i]==OUTPUT_USER_HEAT_AD_Y || Input.output_user[i]==OUTPUT_USER_HEAT_AD_Z || 
        Input.output_user[i]==OUTPUT_USER_HEAT_AD ) {
      if (strlen(Input.filename[FN_ECMWF_WIND_MAP])==0) {
        fprintf (stderr, "Error, in order to get wind or advection information \n");
        fprintf (stderr, "       it is crutial to define an ECMWF_wind_file. \n");
        return -1;
      }
    }
 
  return status;
}

 /*Quick fix for new option name s no_absorption and no_scattering*/
int get_caothoff_index ( caothoff_inp_struct **caothoff,
		      int               *n_caothoff,
		      char              *name )
{
  int isp=0;
  caothoff_inp_struct *newcaoth;

  for (isp=0; isp<*n_caothoff;isp++)
    if ( strcasecmp( (*caothoff)[isp].name, name ) == 0 )
      break;

  /* new caoth, initialize */
  if (isp==*n_caothoff) {

    /* increase number of caoth */
    (*n_caothoff)++;

    /* allocate new caoth struct array with one element more */
    newcaoth=calloc(*n_caothoff, sizeof(caothoff_inp_struct));
    if (newcaoth==NULL) {
      fprintf(stderr,"Error allocating caothoff structure,  (line %d, function %s in %s)\n",
	      __LINE__, __func__, __FILE__);
      exit(1);
    }

    /* copy caoth from old struct array to new struct array */
    for (isp=0; isp<(*n_caothoff)-1;isp++) {
      newcaoth[isp]=(*caothoff)[isp];
      newcaoth[isp].name          = (*caothoff)[isp].name;
    }

    /* free old caoth struct array */ 
    free(*caothoff);

    /* set pointer to new caoth struct array */
    *caothoff = newcaoth;

    /* initialize new caoth */
    (*caothoff)[isp].no_absorption = 0;
    (*caothoff)[isp].no_scattering = 0;

    /* set caoth name */
    (*caothoff)[isp].name = (char *) calloc (strlen(name)+1, sizeof (char));
    strcpy ((*caothoff)[isp].name, name);

  }

  return isp;
}

int get_caoth_index ( caoth_inp_struct **caoth,
		      int               *n_caoth,
		      char              *name,
		      int                old_input )
{
  int isp=0;
  int i,j;
  caoth_inp_struct *newcaoth;

  for (isp=0; isp<*n_caoth;isp++)
    if ( strncasecmp( (*caoth)[isp].name, name, strlen(name) ) == 0 )
      break;

  /* new caoth, initialize */
  if (isp==*n_caoth) {

    /* increase number of caoth */
    (*n_caoth)++;

    /* allocate new caoth struct array with one element more */
    newcaoth=calloc(*n_caoth, sizeof(caoth_inp_struct));
    if (newcaoth==NULL) {
      fprintf(stderr,"Error allocating caoth structure,  (line %d, function %s in %s)\n",
	      __LINE__, __func__, __FILE__);
      exit(1);
    }

    /* copy caoth from old struct array to new struct array */
    for (isp=0; isp<(*n_caoth)-1;isp++) {
      newcaoth[isp]=(*caoth)[isp];
      newcaoth[isp].name          = (*caoth)[isp].name;
      newcaoth[isp].name1         = (*caoth)[isp].name1;
      newcaoth[isp].name2         = (*caoth)[isp].name2;
      newcaoth[isp].fullname      = (*caoth)[isp].fullname;
      newcaoth[isp].filename      = (*caoth)[isp].filename;
      newcaoth[isp].properties_filename = (*caoth)[isp].properties_filename;
    }

    /* free old caoth struct array */
    free(*caoth);

    /* set pointer to new caoth struct array */
    *caoth = newcaoth;

    /* initialize new caoth */
    (*caoth)[isp].cloudcover     = NOT_DEFINED_FLOAT;
    (*caoth)[isp].ipa            = 0;
    (*caoth)[isp].source         = NO_CAOTH;
    (*caoth)[isp].layer          = TRUE;
    (*caoth)[isp].properties     = PROP_NONE;
    (*caoth)[isp].unscaled       = 1;
    (*caoth)[isp].interpolate    = 0;
    (*caoth)[isp].no_scattering  = FALSE;

    (*caoth)[isp].habit          = IC_HABIT_SOLID_COLUMN;
    (*caoth)[isp].reff_prop      = REFF_OU;         /* up to now only for ECMWF */
    (*caoth)[isp].reff           = 20.0;            /* up to now only for ECMWF */
    (*caoth)[isp].fu2yang        = 1;               /* Yang definition of Reff */

    int ic_fu_nn = 2;
    (*caoth)[isp].ic_fu = (int *) calloc ( ic_fu_nn , sizeof (int));
    for (i = 0; i< ic_fu_nn; i++ ) {
	(*caoth)[isp].ic_fu[i] = 0;
    }

    (*caoth)[isp].modify = (float **) calloc (MODIFY_VAR_NN, sizeof (float *));
    for (i=0; i< MODIFY_VAR_NN; i++) {
	(*caoth)[isp].modify[i] = (float*) calloc (MODIFY_TYPE_NN, sizeof (float));
	for (j=0; j< MODIFY_TYPE_NN; j++) {
		(*caoth)[isp].modify[i][j] = NOT_DEFINED_FLOAT;
	}
    }

    /* set caoth name */
    (*caoth)[isp].name = (char *) calloc (strlen(name)+1, sizeof (char));
    strcpy ((*caoth)[isp].name, name);
    if ( old_input ) {
      (*caoth)[isp].name1 = (char *) calloc (strlen(name)+1, sizeof (char));
      strcpy ((*caoth)[isp].name1, name);
      (*caoth)[isp].name2 = (char *) calloc (1, sizeof (char));
      strcpy ((*caoth)[isp].name2, "");
    }
    else {
      (*caoth)[isp].name1 = (char *) calloc (8, sizeof (char));
      strcpy ((*caoth)[isp].name1, "profile");
      (*caoth)[isp].name2 = (char *) calloc (strlen(name)+2, sizeof (char));
      strcpy ((*caoth)[isp].name2, " ");
      strcat ((*caoth)[isp].name2, name);
    }
    (*caoth)[isp].fullname = (char *)
      calloc ( strlen((*caoth)[isp].name1)
	       + strlen((*caoth)[isp].name2)+1, sizeof (char));
    strcpy ((*caoth)[isp].fullname, (*caoth)[isp].name1);
    strcat ((*caoth)[isp].fullname, (*caoth)[isp].name2);

    /* allocate filename char arrays */
    (*caoth)[isp].filename       = calloc (FILENAME_MAX+1, sizeof(char));
    (*caoth)[isp].properties_filename  = calloc (FILENAME_MAX+1, sizeof(char));
  }

  return isp;
}


void pmesg (char *message, int verbose)
{
  if (verbose) {
    fprintf (stderr, "%s", message);
    fflush (stderr);
  }
}


int caoth_check( input_struct Input )
{
  /* Caoth checks */
  int isp=0, status=0;
  double cloudcover=0.0;

  char function_name[]="caoth_check";
  char file_name[]="uvspec.c";

  cloudcover=0.0;

  for (isp=0; isp<Input.n_caoth; isp++) {
    if ( Input.caoth[isp].modify[MODIFY_VAR_TAU   ][MODIFY_TYPE_SET]>=0.0 &&
	 Input.caoth[isp].modify[MODIFY_VAR_TAU550][MODIFY_TYPE_SET]>=0.0 ) {
	fprintf (stderr, "Error in %s (%s), it doesn't make sense to specify both modify tau and modify tau550 for %s;\n",
		 function_name, file_name,Input.caoth[isp].fullname);
      	fprintf (stderr, "please remove either one and retry!\n");
      	status--;
    }

    /* caoth and polarization works only for mie files */
    if ( ( Input.rte.solver == SOLVER_POLRADTRAN && Input.rte.polradtran[POLRADTRAN_NSTOKES]>1 ) || 
	 Input.rte.mc.polarisation == 1 ) {
      switch (Input.caoth[isp].properties) {
      case PROP_FILE:
      case PROP_HEY:
      case PROP_BAUM_V36:
      case PROP_YANG2013:
	break;
      case PROP_MIE:
	/*why did I put this here? RPB:	if (strncasecmp("wc",Input.caoth[isp].name,2)==0) */
	break;
      default:
	fprintf (stderr, "Error in %s (%s), for polarized radiance calculations and %s \n",
		 function_name, file_name, Input.caoth[isp].fullname);
	fprintf(stderr, "you have to use '%s_properties%s XXX', where XXX must contain\n",
		Input.caoth[isp].name1, Input.caoth[isp].name2);
	fprintf(stderr, "phase function matrices (e.g. mie for wc, hey or baum_v36 for ic, or your own optical properties data file)!\n");
	status--;
      }
    }

    /* need to specify a file for the given caoth properties */
    if ( ( Input.caoth[isp].properties != PROP_NONE &&
	   Input.caoth[isp].properties != PROP_EXPLICIT ) &&
	 strlen(Input.caoth[isp].filename)==0 ) {
      fprintf (stderr, "Error, %s_properties%s only makes sense together with '%s_file%s'!\n",
	       Input.caoth[isp].name1,Input.caoth[isp].name2,
	       Input.caoth[isp].name1,Input.caoth[isp].name2);
      status--;
    }

    if (strlen(Input.caoth[isp].properties_filename)>0 && !Input.caoth[isp].layer) {
      fprintf (stderr, "Error, species_properties %s requires optical properties defined per layer for 1D clouds.\n",
	       Input.caoth[isp].name);
      fprintf (stderr, "Please specify species_layer %s and retry if that is what you mean!\n",
	       Input.caoth[isp].name);
      status--;
    }

    if ( ( Input.caoth[isp].properties==PROP_YANG ||
	   Input.caoth[isp].properties==PROP_KEY )
	&& Input.caoth[isp].habit < 0 ) {
      fprintf (stderr, "Error, unknown ice crystal habit for %s. Check habit!\n",
	       Input.caoth[isp].fullname);
      status--;
    }

    /* need to specify a Raytracing file if Raytracing properties were specified */
    if ((Input.caoth[isp].properties == PROP_RAYTRACING) && strlen(Input.filename[FN_RAYTRACING]) == 0) {
	fprintf (stderr, "Error, 'ic_properties RAYTRACING' requires a filename definition with 'ic_raytracing_file filename'!\n");
	status--;
      }
	
    /*check caoth cloudcover*/
    if (Input.caoth[isp].cloudcover != NOT_DEFINED_FLOAT) {
      	if (Input.caoth[isp].cloudcover < 0.0 || Input.caoth[isp].cloudcover > 1.0) {
		fprintf (stderr, "Error, cloudcover %f for %s not between 0 and 1\n",
		Input.caoth[isp].cloudcover, Input.caoth[isp].fullname);
		status--;
	}

      	if ( cloudcover == 0.0 )
		cloudcover = Input.caoth[isp].cloudcover;

      	if ( cloudcover != Input.caoth[isp].cloudcover ) {
		fprintf (stderr, "Error, all cloudcovers must be equal!");
		status--;
      	}
    	if ( strlen(Input.caoth[isp].filename)==0 ||  Input.caoth[isp].source != CAOTH_FROM_1D ) {
		fprintf (stderr, "Error, cloudcover for %s must be accompanied by cloud file\n",
		Input.caoth[isp].fullname);
		status--;
    	}

	if ( Input.caoth[isp].source == CAOTH_FROM_IPA_FILES ) {
		fprintf (stderr, "Error, cloudcover does not work together with ipa_files.\n"); 
		fprintf (stderr, "Please choose either cloudcover OR ipa_files!\n"); 
		status--;
	} 

    }

    /* caoth warnings */
    switch (Input.caoth[isp].properties) {
	    case PROP_MIE:
	    case PROP_FILE:
	    case PROP_EXPLICIT:
	    case PROP_IC_MIE:
	    case PROP_BAUM:
	    case PROP_BAUM_V36:
            case PROP_HEY:
            case PROP_YANG2013:
	      if (Input.rte.solver==SOLVER_FDISORT1) {
	        fprintf (stderr, "\n");
	        fprintf (stderr, "*** Warning, you specified a tabulated phase function for %s.\n",
			 Input.caoth[isp].fullname);
	        fprintf (stderr, "*** It is strongly recommended to use cdisort in this case because\n");
	        fprintf (stderr, "*** the older disort 1.3 version uses only the first NSTR moments of the\n");
	        fprintf (stderr, "*** phase function which might result in really bad distortions.\n");
	        fprintf (stderr, "*** cdisort consideres all moments correctly.\n");
	        fprintf (stderr, "\n\n");
	      }
	      break;
	    case PROP_YANG:
	#if HAVE_YANG
	#else
	      fprintf (stderr, "***********************************************************************\n");
	      fprintf (stderr, "* Error, you have built uvspec without Yang/Key/Mayer support and     *\n");
	      fprintf (stderr, "* hence cannot use 'ic_properties yang'.                              *\n");
	      fprintf (stderr, "***********************************************************************\n");
	      return -1;
	#endif
	    break;
	    default:
	      break;
    }

    if (strlen(Input.caoth[isp].filename)>0 && Input.caoth[isp].properties==PROP_FU) {
	if ( Input.caoth[isp].modify[MODIFY_VAR_GG    ][MODIFY_TYPE_SCALE] >= -1 ||
	     Input.caoth[isp].modify[MODIFY_VAR_GG    ][MODIFY_TYPE_SET  ] >= -1 ||
	     Input.caoth[isp].modify[MODIFY_VAR_SSA   ][MODIFY_TYPE_SCALE] >= 0  ||
	     Input.caoth[isp].modify[MODIFY_VAR_SSA   ][MODIFY_TYPE_SET  ] >= 0  ||
	     Input.caoth[isp].modify[MODIFY_VAR_TAU   ][MODIFY_TYPE_SET  ] >  0  ||
	     Input.caoth[isp].modify[MODIFY_VAR_TAU550][MODIFY_TYPE_SET  ] > 0 ) {
	  fprintf (stderr, "*** Warning, you define an ice cloud using the parameterization by Fu (1996/98) *\n");
	  fprintf (stderr, "*** and explicitely modify one of the optical properties, optical thickness,    *\n");
	  fprintf (stderr, "*** asymmetry parameter, or single scattering albedo. Please make sure that you *\n");
	  fprintf (stderr, "*** understand the discussion of delta-scaling in context with the Fu (1996/98) *\n");
	  fprintf (stderr, "*** parameterization in the libRadtran manual; see the description of           *\n");
	  fprintf (stderr, "*** ic_properties, ic_fu_tau, and ic_fu_reff!                                   *\n");
	}
    }

  } /* end loop isp */
  
  return status;

}

int caoth_set( input_struct *Input, output_struct *Output )
{
  /* Caoth settings */
  int isp=0, status=0;

  /*big caoth loop*/
  for (isp=0; isp<(*Input).n_caoth; isp++) {
    if ( strncasecmp((*Input).caoth[isp].name,"ic",2)==0 ) {
      (*Input).i_ic=isp; 

      /* something special: the new Mayer/Key/Yang parameterization is not yet */ 
      /* included in the distribution; if it is not found, Key et al is used   */ 
      /* instead.                                                              */ 
#if HAVE_YANG 
#else 
      if ((*Input).caoth[isp].properties == PROP_YANG) 
	(*Input).caoth[isp].properties = PROP_KEY; 
#endif 
		
      if ((*Input).caoth[isp].ic_fu[IC_FU_REFF_DEF] == SWITCH_ON )
	(*Input).caoth[isp].fu2yang  = 0;
      if ((*Input).caoth[isp].ic_fu[IC_FU_DELTASCALING] == SWITCH_ON )
	(*Input).caoth[isp].unscaled = 0;
    }

    if ( strncasecmp((*Input).caoth[isp].name,"wc",2)==0 )
      (*Input).i_wc=isp; 

    switch((*Input).caoth[isp].source) {
    case CAOTH_FROM_MOMENTS:
      (*Input).caoth[isp].properties = PROP_EXPLICIT;
      break;
    case CAOTH_FROM_IPA_FILES:
      (*Input).caoth[isp].ipa = 1;
      (*Input).ipa = 1;
      break;
    default:
      break;	
    }

    /* caoth default values */
    if ( (*Input).caoth[isp].properties==PROP_NONE ) {
      if ( strncasecmp((*Input).caoth[isp].name,"wc",2)==0 ) {
	(*Input).caoth[isp].properties = PROP_HU;
	if (!(*Input).quiet)
	  fprintf (stderr, " ... using default water cloud properties, Hu and Stamnes (1993)\n");
      }
      else if ( strncasecmp((*Input).caoth[isp].name,"ic",2)==0 ) {
	(*Input).caoth[isp].properties = PROP_FU;
		
	if (!(*Input).quiet)
	  fprintf (stderr, " ... using default ice cloud properties, Fu et al. (1996/98)\n");
      }
      else {
	fprintf(stderr, "Error! You need to specify properties for %s!\n",
		(*Input).caoth[isp].fullname );
	status--;
      }
    }

    switch((*Input).caoth[isp].properties)	{
    case PROP_ECHAM4:
      if ( strncasecmp((*Input).caoth[isp].name,"wc",2) !=0 &&  strncasecmp((*Input).caoth[isp].name,"ic",2) !=0 ) {
	fprintf (stderr, "Error, profile properties for ECHAM4 does not work with profile %s! Need to use 'wc' or 'ic'!\n",(*Input).caoth[isp].name);
	status--;
      }
      break;
    case PROP_MIE:
    case PROP_HU:
      if ( strncasecmp((*Input).caoth[isp].name,"ic",2) == 0 ) {
	fprintf (stderr, "Error, profile_properties in input file is a water cloud property and is not allowed for profile 'ic'! Define other profile name (e.g. 'wc')!\n");
	status--;
      }
      break;
    case PROP_FU:
    case PROP_IC_MIE:
    case PROP_YANG:
    case PROP_KEY:
    case PROP_BAUM_V36:
    case PROP_BAUM:
    case PROP_HEY:
    case PROP_YANG2013:
      if (strncasecmp( (*Input).caoth[isp].name, "wc", 2)==0) {
	fprintf (stderr, "Error, profile_properties in input file is an ice cloud property and is not allowed for profile 'wc'! Define other profile name (e.g. 'ic')!\n");
	status--;
      }
      break;
    default:
      break;
    }

    if ( (*Input).caoth[isp].source == CAOTH_FROM_ECMWF && (*Input).caoth[isp].layer != 1 ) {
      (*Input).caoth[isp].layer=1;
      if ((*Input).verbose)	fprintf (stderr, " ... switching on layer for %s, as ECMWF cloud data is representative for layers \n", (*Input).caoth[isp].fullname);
    }
    /* quick fix for cloud_overlap introduced at SBCA */
    if ((*Input).cloud_overlap != CLOUD_OVERLAP_OFF && (*Input).cloud_overlap != CLOUD_OVERLAP_RAND) {
      (*Input).caoth[isp].ipa=1;
      (*Input).ipa=1; 
      (*Output).nipa=2;
    }


  } /*end caoth loop*/

  return status;
}

