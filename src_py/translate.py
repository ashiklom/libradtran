
def translate_input(inputfile, outputname='',force=False):
	import os

	lines=inputfile.readlines()
	inputfile.close()
	
	#preparser
	preparse_dict = preparse(lines)
	zip_source		= search_for_zip_source(lines)

	if outputname:
		if os.path.exists(outputname):
			if not force:
				import sys
				print """Error occured!\nFilename {0} for new input file already exists.\nPlease call python with another output filename:\n
				python {1} --new_filename new_file.inp\n
				or overwrite file with:\n
				python {1} --new_filename new_file.inp --force""".format(outputname,inputfile.name)
				sys.exit()
	else: outputname = False

	translate(lines, preparse_dict=preparse_dict, zip_source=zip_source, outputfile=outputname)

def preparse(lines):
	import re

	preparse_dict = {
		'wc_properties_interpolate' : False,
		'ic_properties_interpolate' : False,
		'profile_properties_interpolate' : False,
		'raman_original' : False,
		'filter_function_normalize' : False,
		'mc_spectral_is_wvl' : False,
		'mc_hiddenline' : False
	}

	for line in lines:
		match_str=((line.split('#')[0]).strip()).split()
		if not match_str: continue
		for key in preparse_dict.keys():
			if key == 'mc_spectral_is_wvl':
				if match_str[0] == key: preparse_dict[key] = match_str[1]
			else:
				if match_str[0] == key: preparse_dict[key] = True

	return preparse_dict

def search_for_zip_source(lines):
	import re

	source = ''
	solarfile = ''
	unit=''
        for line in lines:
		match_str=((line.split('#')[0]).strip()).split()
                if match_str and match_str[0]=='solar_file': 
			if len(match_str) > 1: solarfile= match_str[1]
			if len(match_str) > 2: unit	= match_str[2]

        for line in lines:
		match_str=((line.split('#')[0]).strip()).split()
                if match_str and match_str[0]=='source': 
			if len(match_str) > 1: source = match_str[1]

	if source and solarfile:
		return 'source {0} {1} {2}'.format(source,solarfile,unit)
	else: return False

def translate(lines,preparse_dict={},zip_source=False,outputfile=''):

	if outputfile:
		try:
			output = file(outputfile,'w')
		except Exception,e:
			import sys
			print "Exception occured when opening file: ",e
			sys.exit()


	for line in lines:
		line_translate = translate_line(line,preparse_dict=preparse_dict,zip_source=zip_source)
		
		if outputfile:
			output.write(line_translate)
		else:
			if line_translate: print line_translate.strip()

	if zip_source:
                if outputfile:
                        output.write(zip_source)
                else:
                        print zip_source


	if outputfile: output.close()

 
def translate_line(line,preparse_dict={},zip_source=False):
	import re
	match_str=((line.split('#')[0]).strip()).split()

	if not match_str:			return line

	#spectral_options
	if match_str[0]=='transmittance_wl_file':	return re.sub('transmittance_wl_file',      'wavelength_grid_file',     line,count=1)	
	if match_str[0]=='wvn':		return re.sub('wvn',		'wavelength',		line,count=1)
	if match_str[0]=='wvn_index':	return re.sub('wvn_index',	'wavelength_index',	line,count=1)

	#source_options
	if match_str[0]=='solar_file':
		if zip_source:	return ''
		else:		return re.sub('solar_file',           'source solar',        line,count=1)
	if match_str[0]=='source':
		if zip_source:	return ''
		else:		return line
	if match_str[0]=='mc_sun_radius':	return re.sub('mc_sun_radius',	'mc_sun_angular_size',  line,count=1)
	
	#generell_atmosphere_options
	if match_str[0]=='wc_level':		return re.sub('wc_level',	'interpret_as_level wc',line,count=1)
	if match_str[0]=='ic_level':		return re.sub('ic_level',	'interpret_as_level ic',line,count=1)
	if match_str[0]=='profile_level':	return re.sub('profile_level',	'interpret_as_level',	line,count=1)
	if match_str[0]=='ic_no_scattering':	return re.sub('ic_no_scattering','no_scattering ic',	line,count=1)
	if match_str[0]=='wc_no_scattering':	return re.sub('wc_no_scattering','no_scattering wc',	line,count=1)
	if match_str[0]=='aerosol_no_scattering':	return re.sub('aerosol_no_scattering','no_scattering aer',	line,count=1)
	if match_str[0]=='no_rayleigh':		return re.sub('no_rayleigh','no_scattering mol',	line,count=1)
	if match_str[0]=='no_molecular_absorption':	return re.sub('no_molecular_absorption','no_absorption mol',	line,count=1)
	if match_str[0]=='reverse':		return re.sub('reverse',	      'reverse_atmosphere',  line,count=1)	

	#molecular_atmosphere_options
	if match_str[0]=='correlated_k':	return re.sub('correlated_k',	      'mol_abs_param',  line,count=1)	
	if match_str[0]=='fu_h2o_continuum':	
		line = re.sub('fu_h2o_continuum',     'ck_fu_h2o_continuum',  line,count=1)
		if   match_str[1]=='2.1': return re.sub('2.1', 'v2.1', line, count=1)
		elif match_str[1]=='2.4': return re.sub('2.4', 'v2.4', line, count=1)
		else: return line
	if match_str[0]=='absorption':		return re.sub('absorption',	      'ck_lowtran_absorption',  line,count=1)	
	if match_str[0]=='molecular_tau_file':	return re.sub('molecular_tau_file',   'mol_tau_file abs',  line,count=1)	
	if match_str[0]=='rayleigh_tau_file':	return re.sub('rayleigh_tau_file',    'mol_tau_file sca',  line,count=1)	
	if match_str[0]=='dens_file':		
		if len(match_str) == 4:	return 'mol_file {1} {3} {2}\n'.format( *match_str )
		else: 			return re.sub('dens_file', 'mol_file', line, count =1) 
	if match_str[0]=='rh_file':	return 'mol_file H2O {1} rh\n'.format( *match_str )
	if match_str[0]=='rayleigh_crs':	return re.sub('rayleigh_crs',	      'crs_model rayleigh',  line,count=1)	
	if match_str[0]=='o3_crs':		return re.sub('o3_crs',	 	      'crs_model O3', 	     line,count=1)	
	if match_str[0]=='no2_crs':		return re.sub('no2_crs',	      'crs_model NO2', 	     line,count=1)	
	if match_str[0]=='ozone_column':	return 'mol_modify O3 {1} DU\n'.format( *match_str )
	if match_str[0]=='h2o_precip':		return 'mol_modify H2O {1} MM\n'.format( *match_str )
	if match_str[0]=='no2_column_du':	return 'mol_modify NO2 {1} DU\n'.format( *match_str )
	if match_str[0]=='no2_column_moleccm-2':return 'mol_modify NO2 {1} CM_2\n'.format( *match_str )
	if match_str[0]=='dens_column':
		if len(match_str) == 4:		return 'mol_modify {1} {2} {3}\n'.format( *match_str )
		if match_str[1]=='o3' or match_str[1]=='O3':
			return 'mol_modify {1} {2} DU\n'.format( *match_str )
		else:				return 'mol_modify {1} {2} CM_2\n'.format ( *match_str )
	if match_str[0]=='f11_mixing_ratio':	return 'mixing_ratio F11 {1}\n'.format( *match_str )
	if match_str[0]=='f12_mixing_ratio':	return 'mixing_ratio F12 {1}\n'.format( *match_str )
	if match_str[0]=='f22_mixing_ratio':	return 'mixing_ratio F22 {1}\n'.format( *match_str )
	if match_str[0]=='ch4_mixing_ratio':	return 'mixing_ratio CH4 {1}\n'.format( *match_str )
	if match_str[0]=='o2_mixing_ratio':	return 'mixing_ratio O2 {1} n'.format( *match_str )
	if match_str[0]=='n2o_mixing_ratio':	return 'mixing_ratio N2O {1}\n'.format( *match_str )
	if match_str[0]=='h2o_mixing_ratio':	return 'mixing_ratio H2O {1}\n'.format( *match_str )
	if match_str[0]=='co2_mixing_ratio':	return 'mixing_ratio CO2 {1}\n'.format( *match_str )

	#aerosol_options
	if match_str[0]=='aerosol_gg_file':	return re.sub('aerosol_gg_file',      'aerosol_file gg',     line,count=1)	
	if match_str[0]=='aerosol_ssa_file':	return re.sub('aerosol_ssa_file',     'aerosol_file ssa',    line,count=1)	
	if match_str[0]=='aerosol_tau_file':	return re.sub('aerosol_tau_file',     'aerosol_file tau',    line,count=1)	
	if match_str[0]=='aerosol_files':	return re.sub('aerosol_files',	      'aerosol_file explicit',     line,count=1)	
	if match_str[0]=='aerosol_moments_file':return re.sub('aerosol_moments_file', 'aerosol_file moments',     line,count=1)	
	if match_str[0]=='aerosol_set_gg':	return re.sub('aerosol_set_gg',	'aerosol_modify gg set',	line,count=1)
	if match_str[0]=='aerosol_set_ssa':	return re.sub('aerosol_set_ssa',	'aerosol_modify ssa set',	line,count=1)
	if match_str[0]=='aerosol_set_tau':	return re.sub('aerosol_set_tau',	'aerosol_modify tau set',	line,count=1)
	if match_str[0]=='aerosol_set_tau550':	return re.sub('aerosol_set_tau550',	'aerosol_modify tau550 set',	line,count=1)
	if match_str[0]=='aerosol_scale_gg':	return re.sub('aerosol_scale_gg',	'aerosol_modify gg scale',	line,count=1)
	if match_str[0]=='aerosol_scale_ssa':	return re.sub('aerosol_scale_ssa',	'aerosol_modify ssa scale',	line,count=1)
	if match_str[0]=='aerosol_scale_tau':	return re.sub('aerosol_scale_tau',	'aerosol_modify tau scale',	line,count=1)
	if match_str[0]=='angstrom':		return re.sub('angstrom',	      'aerosol_angstrom',  line,count=1)	

	#profile_options
	if match_str[0]=='profile_set_gg':	return re.sub('profile_set_gg\s+'+match_str[1],		'profile_modify %s gg set'% (match_str[1]) ,	line,count=1)
	if match_str[0]=='profile_set_ssa':	return re.sub('profile_set_ssa\s+'+match_str[1],	'profile_modify %s ssa set'% (match_str[1]) ,	line,count=1)
	if match_str[0]=='profile_set_tau':	return re.sub('profile_set_tau\s+'+match_str[1],	'profile_modify %s tau set'% (match_str[1]) ,	line,count=1)
	if match_str[0]=='profile_set_tau550':	return re.sub('profile_set_tau550\s+'+match_str[1],	'profile_modify %s tau550 set'% (match_str[1]),	line,count=1)
	if match_str[0]=='profile_scale_gg':	return re.sub('profile_scale_gg\s+'+match_str[1],	'profile_modify %s gg scale'% (match_str[1]) ,	line,count=1)
	if match_str[0]=='profile_scale_ssa':	return re.sub('profile_scale_ssa\s+'+match_str[1],	'profile_modify %s ssa scale'% (match_str[1]) ,	line,count=1)
	if match_str[0]=='profile_scale_tau':	return re.sub('profile_scale_tau\s+'+match_str[1],	'profile_modify %s tau scale'% (match_str[1]) ,	line,count=1)
	if match_str[0]=='profile_properties':
		if preparse_dict['profile_properties_interpolate']:	return re.sub(match_str[1],	'{} interpolate'.format(match_str[1]), line,count=1)
		else:				return line
	if match_str[0]=='profile_properties_interpolate': return ''

	#cloud_options
	if match_str[0]=='wc_file':		
		if match_str[1]=='1D' or match_str[1]=='3D' or match_str[1]=='moments' or match_str[1]=='ipa' or match_str[1]=='ipa_files':		return line
		else:				return re.sub('wc_file',	'wc_file 1D',		line,count=1)
	if match_str[0]=='mc_wcloud_file':	return re.sub('mc_wcloud_file',	'wc_file 3D',		line,count=1)
	if match_str[0]=='wc_files':		return re.sub('wc_files',	'wc_file moments',	line,count=1)
	if match_str[0]=='wc_ipa_files':	return re.sub('wc_ipa_files',	'wc_file ipa_files',	line,count=1)
	if match_str[0]=='ic_file':		
		if match_str[1]=='1D' or match_str[1]=='3D' or match_str[1]=='moments' or match_str[1]=='ipa' or match_str[1]=='ipa_files':		return line
		else:				return re.sub('ic_file',	'ic_file 1D',		line,count=1)
	if match_str[0]=='mc_icloud_file':	return re.sub('mc_icloud_file',	'ic_file 3D',		line,count=1)
	if match_str[0]=='ic_files':		return re.sub('ic_files',	'ic_file moments',	line,count=1)
	if match_str[0]=='ic_ipa_files':	return re.sub('ic_ipa_files',	'ic_file ipa_files',	line,count=1)
	if match_str[0]=='wc_set_gg':		return re.sub('wc_set_gg',	'wc_modify gg set',	line,count=1)
	if match_str[0]=='wc_set_ssa':		return re.sub('wc_set_ssa',	'wc_modify ssa set',	line,count=1)
	if match_str[0]=='wc_set_tau':		return re.sub('wc_set_tau',	'wc_modify tau set',	line,count=1)
	if match_str[0]=='wc_set_tau550':	return re.sub('wc_set_tau550',	'wc_modify tau550 set',	line,count=1)
	if match_str[0]=='wc_scale_gg':		return re.sub('wc_scale_gg',	'wc_modify gg scale',	line,count=1)
	if match_str[0]=='wc_scale_ssa':	return re.sub('wc_scale_ssa',	'wc_modify ssa scale',	line,count=1)
	if match_str[0]=='wc_scale_tau':	return re.sub('wc_scale_tau',	'wc_modify tau scale',	line,count=1)
	if match_str[0]=='wc_cloudcover':	return re.sub('wc_cloudcover',	'cloudcover wc',	line,count=1)
	if match_str[0]=='ic_set_gg':		return re.sub('ic_set_gg',	'ic_modify gg set',	line,count=1)
	if match_str[0]=='ic_set_ssa':		return re.sub('ic_set_ssa',	'ic_modify ssa set',	line,count=1)
	if match_str[0]=='ic_set_tau':		return re.sub('ic_set_tau',	'ic_modify tau set',	line,count=1)
	if match_str[0]=='ic_set_tau550':	return re.sub('ic_set_tau550',	'ic_modify tau550 set',	line,count=1)
	if match_str[0]=='ic_scale_gg':		return re.sub('ic_scale_gg',	'ic_modify gg scale',	line,count=1)
	if match_str[0]=='ic_scale_ssa':	return re.sub('ic_scale_ssa',	'ic_modify ssa scale',	line,count=1)
	if match_str[0]=='ic_scale_tau':	return re.sub('ic_scale_tau',	'ic_modify tau scale',	line,count=1)
	if match_str[0]=='ic_cloudcover':	return re.sub('ic_cloudcover',	'cloudcover ic',	line,count=1)
	if match_str[0]=='ic_properties':
		if match_str[1] == 'mie': match_str[1] = 'ic-mie'
		if preparse_dict['ic_properties_interpolate']:	return re.sub(match_str[1],	'{} interpolate'.format(match_str[1]), line,count=1)
		else:				return line
	if match_str[0]=='wc_properties':
		if preparse_dict['wc_properties_interpolate']:	return re.sub(match_str[1],	'{} interpolate'.format(match_str[1]), line,count=1)
		else:				return line
	if match_str[0]=='wc_properties_interpolate': return ''
	if match_str[0]=='ic_properties_interpolate': return ''
	if match_str[0]=='ic_fu_reff':		return re.sub('ic_fu_reff',	'ic_fu reff_def on',	line,count=1)
	if match_str[0]=='ic_fu_tau':		
		if match_str[1] == 'scaled':		line = re.sub('scaled',	'on',	line,count=1)
		elif match_str[1] == 'unscaled':	line = re.sub('unscaled','off',	line,count=1)
		return re.sub('ic_fu_tau',     'ic_fu deltascaling',	line,count=1)

	#surface_options
	if match_str[0]=='cox_and_munk_solar_wind': return re.sub('cox_and_munk_solar_wind', 'brdf_cam_solar_wind', line,count=1)
	if match_str[0]=='cox_and_munk_sal':	return re.sub('cox_and_munk_sal',	'brdf_cam sal',		line,count=1)
	if match_str[0]=='cox_and_munk_pcl':	return re.sub('cox_and_munk_pcl',	'brdf_cam pcl',		line,count=1)
	if match_str[0]=='cox_and_munk_u10':	return re.sub('cox_and_munk_u10',	'brdf_cam u10',		line,count=1)
	if match_str[0]=='cox_and_munk_uphi':	return re.sub('cox_and_munk_uphi',	'brdf_cam uphi',	line,count=1)
	if match_str[0]=='rpv_library':		return re.sub('rpv_library',		'brdf_rpv_library',	line,count=1)
	if match_str[0]=='rpv_k':		return re.sub('rpv_k',			'brdf_rpv k',		line,count=1)
	if match_str[0]=='rpv_rho0':		return re.sub('rpv_rho0',		'brdf_rpv rho0',	line,count=1)
	if match_str[0]=='rpv_theta':		return re.sub('rpv_theta',		'brdf_rpv theta',	line,count=1)
	if match_str[0]=='rpv_sigma':		return re.sub('rpv_sigma',		'brdf_rpv sigma',	line,count=1)
	if match_str[0]=='rpv_t1':		return re.sub('rpv_t1',			'brdf_rpv t1',		line,count=1)
	if match_str[0]=='rpv_t2':		return re.sub('rpv_t2',			'brdf_rpv t2',		line,count=1)
	if match_str[0]=='rpv_scale':		return re.sub('rpv_scale',		'brdf_rpv scale',	line,count=1)
	if match_str[0]=='rpv_file':		return re.sub('rpv_file',		'brdf_rpv_file',	line,count=1)
	if match_str[0]=='mc_ambrals_file':	return re.sub('mc_ambrals_file',        'mc_rossli_file',	line,count=1)
	if match_str[0]=='surface_type':	return re.sub('surface_type',		'brdf_rpv_type',	line,count=1)
	if match_str[0]=='surface_temperature':	return re.sub('surface_temperature',	'sur_temperature',	line,count=1)
	if match_str[0]=='mc_temperature_file':	return re.sub('mc_temperature_file',	'sur_temperature_file',	line,count=1)

	#solver_options
	if match_str[0]=='rte_solver':
		if match_str[1] == 'ctwostr':	return re.sub('ctwostr',		'twostr',		line,count=1)
		if match_str[1] == 'twostream':	return re.sub('twostream',		'twostr',		line,count=1)
		if match_str[1] == 'twostreampp': return re.sub('twostreampp',		'twostr',		line,count=1)
		if match_str[1] == 'twostrpp':	return re.sub('twostrpp',		'twostr',		line,count=1)
		if match_str[1] == 'cdisort':	return re.sub('cdisort',		'disort',		line,count=1)
		if match_str[1] == 'disort2':	return re.sub('disort2',		'disort',		line,count=1)
	if match_str[0]=='cdisort_spherical_albedo':	return re.sub('cdisort_spherical_albedo', 'disort_spherical_albedo', line,count=1)
	if match_str[0]=='fisot':		return re.sub('fisot',			'isotropic_source_toa',	line,count=1)
	if match_str[0]=='raman':
		if preparse_dict['raman_original'] and len(match_str) == 1 :
			return re.sub('raman',	'raman original', line,count=1)
		else:	return line
	if match_str[0]=='raman_original':	return ''
	if match_str[0]=='nstr':		return re.sub('nstr',			'number_of_streams',	line,count=1)
	if match_str[0]=='disort_icm':		return re.sub('disort_icm',		'disort_intcor',	line,count=1)
	if match_str[0]=='sslidar_area':	return re.sub('sslidar_area',		'sslidar area',		line,count=1)
	if match_str[0]=='sslidar_E0':		return re.sub('sslidar_E0',		'sslidar E0',		line,count=1)
	if match_str[0]=='sslidar_eff':		return re.sub('sslidar_eff',		'sslidar eff',		line,count=1)
	if match_str[0]=='sslidar_position':	return re.sub('sslidar_position',	'sslidar position',	line,count=1)
	if match_str[0]=='sslidar_range':	return re.sub('sslidar_range',		'sslidar range',	line,count=1)
	if match_str[0]=='cdisort_pseudospherical':		return re.sub('cdisort_pseudospherical',	'pseudospherical',line,count=1)
	if match_str[0]=='ctwostr_pseudospherical':		return re.sub('ctwostr_pseudospherical',	'pseudospherical',line,count=1)
	if match_str[0]=='polradtran_aziorder':	return re.sub('polradtran_aziorder',	'polradtran aziorder',	line,count=1)
	if match_str[0]=='polradtran_nstokes':	return re.sub('polradtran_nstokes',	'polradtran nstokes',	line,count=1)
	if match_str[0]=='polradtran_src_code':	return re.sub('polradtran_src_code',	'polradtran src_code',	line,count=1)
	if match_str[0]=='nscat':		return re.sub('nscat',			'sdisort nscat',	line,count=1)
	if match_str[0]=='nrefrac':		return re.sub('nrefrac',		'sdisort nrefrac',	line,count=1)
	if match_str[0]=='ichap':		return re.sub('ichap',			'sdisort ichapman',	line,count=1)

	#mc_options
	if match_str[0]=='mc_ris_optical_depth':	return re.sub('mc_ris_optical_depth',	'mc_ris optical_depth',    line,count=1)
	if match_str[0]=='mc_ris_factor':	
		if len(match_str) > 1:	return re.sub('mc_ris_factor',	'mc_ris factor',    line,count=1)
		else:			return 'mc_ris optical_depth 1\n'
	if match_str[0]=='mc_spectral_is_wvl' : return ''
	if match_str[0]=='mc_spectral_is':
		if preparse_dict['mc_spectral_is_wvl']:	return re.sub('mc_spectral_is', 'mc_spectral_is {}'.format(preparse_dict['mc_spectral_is_wvl']), line,count=1)
		else:					return line
	if match_str[0]=='mc_hiddenline' : return ''
	if match_str[0]=='mc_visualize' :
		if preparse_dict['mc_hiddenline'] :	return re.sub('mc_visualize', 'mc_visualize hiddenline', line,count=1)
		else:					return line

	#geometry_options
	if match_str[0]=='mc_bcond':		return re.sub('absorp',	'absorb', line,count=1)
	if match_str[0]=='mc_panorama_quicklook':			return re.sub('mc_panorama_quicklook',			'mc_panorama quicklook', 		line,count=1)
	if match_str[0]=='mc_panorama_no_pixel':			return re.sub('mc_panorama_no_pixel',			'mc_panorama no_pixel', 		line,count=1)
	if match_str[0]=='mc_panorama_distr_photons_over_pixel':	return re.sub('mc_panorama_distr_photons_over_pixel',	'mc_panorama distr_photons_over_pixel', line,count=1)
	if match_str[0]=='mc_panorama_with_direct_rad':			return re.sub('mc_panorama_with_direct_rad',		'mc_panorama with_direct_rad', 		line,count=1)
	if match_str[0]=='mc_spherical':
		if len(match_str) == 1:	return re.sub('mc_spherical',	'mc_spherical 1D', line,count=1)
		else:			return line
	if match_str[0]=='mc_spherical3D': return re.sub('mc_spherical3D',   'mc_spherical 3D', line,count=1)

	#output_options
	if match_str[0]=='prndis':		return re.sub('prndis',		'print_disort_info',		line,count=1)
	if match_str[0]=='mc_absorption':	return re.sub('mc_absorption',	'mc_forward_output absorption', line,count=1)
	if match_str[0]=='mc_actinic':		return re.sub('mc_actinic',	'mc_forward_output actinic',    line,count=1)
	if match_str[0]=='mc_emission':		return re.sub('mc_emission',	'mc_forward_output emission',   line,count=1)
	if match_str[0]=='mc_heating':		return re.sub('mc_heating',	'mc_forward_output heating',    line,count=1)
	if match_str[0]=='mc_backward_output' and match_str[1]=='f':		return re.sub('f',	'act',    line,count=1)
	if match_str[0]=='output_process' and match_str[1]=='per_ck_band':	return re.sub('per_ck_band',	'per_band',    line,count=1)

	if match_str[0]=='output':		return re.sub('output', 	'output_process', line, count=1)
	if match_str[0]=='flexstor':		return re.sub('flexstor', 	'output_format flexstor', line, count=1)
	if match_str[0]=='brightness':		return re.sub('brightness', 	'output_quantity brightness', line, count=1)
	if match_str[0]=='reflectivity':	return re.sub('reflectivity', 	'output_quantity reflectivity', line, count=1)
	if match_str[0]=='transmittance':	return re.sub('transmittance', 	'output_quantity transmittance', line, count=1)
	if match_str[0]=='filter_function_normalize': return ''
	if match_str[0]=='filter_function_file':
		if preparse_dict['filter_function_normalize']:	return re.sub(match_str[1],	'{} normalize'.format(match_str[1]), line,count=1)
		else:				return line
	
	return line


import argparse
#read input from cmd-line arguments
parser = argparse.ArgumentParser(description='Python programm to translate old uvspec input file to new input files.')
parser.add_argument('inputfile', type=file, 
			help='Name of uvspec input-file to translate from old to new input style.')
parser.add_argument('--new_filename', default='', type=str,
			help='Name of new style uvspec input-file.')
parser.add_argument('--force', action='store_true', default=False,
			help='Overwrite new_filename if file already exists.' )

try: 
	args = parser.parse_args()
except IOError,e:
	import sys
	print "Error occured while parsing cmd-line-arguments: ",e
	sys.exit(1)
except Exception,e:
	print e
	pass

translate_input(args.inputfile, args.new_filename, args.force)
