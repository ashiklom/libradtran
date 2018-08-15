#
# Copyright (C) 1997 Arve Kylling
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 1, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# To obtain a copy of the GNU General Public License write to the
# Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139,
# USA.

#-------------------------------------------------------------------------

# Naming convention: subs that start with an uppercase letter may be
#                    useful on their own and have their own usage
#                    subs. Subs starting with a lowercase letter are
#                    intended to be called by other subs, hence no 
#                    usage subs for them. Uppercase subs comes first,
#                    lowercase subs after.

package UVTools;
use Carp;

sub load_BAND_response_fnc
# loads response function from 2 column file with wavelength (nm) and normalized response
# lines starting with # are ignored.
# This need to be her up in front in order to be recognised by BEGIN

{
    local($file) = @_;
    open(BAND, $file);
    local($wvn,$resp) = 0.0;
  LINE:
    while($line=<BAND>)  {
	$_ = $line;                   # Make m match $line
	if (m/^[\#]/) {  next LINE; } # Skip comment lines
	$line =~ s/^[ ]*//;
	chop $line;
	($wvn,$resp) = split(/[ \t]+/, $line);
	$wvn = sprintf("%7.3f",$wvn);
	$response_fnc{$wvn} = $resp;
    }
    close(BAND);
    %response_fnc;
}
sub load_RB_response_fnc

# This need to be her up in front in order to be recognised by BEGIN

{
    local($file) = @_;
    open(RB, $file);
    local($wvn,$resp) = 0.0;
  LINE:
    while($line=<RB>)  {
	$_ = $line;                   # Make m match $line
	if (m/^[\#]/) {  next LINE; } # Skip comment lines
	$line =~ s/^[ ]*//;
	chop $line;
	($wvn,$resp) = split(/[ \t]+/, $line);
	$wvn = sprintf("%7.3f",$wvn/10.); # Divide by 10 to get from Angstrom to nm
	$response_fnc{$wvn} = $resp;
    }
    close(RB);
    %response_fnc;
}
    


BEGIN { 
    $time_resolution = 1.0;  # in hours
    $toolspath       = "../bin/";
    $sza2time        = $toolspath."sza2time";
    $sza2time_out    = "./sza2time_out";
    $zenith          = $toolspath."zenith";
    $zenith_out      = "./zenith_out";
    
    $conv         = $toolspath."conv";
    $integrate    = $toolspath."integrate -q";
    $spline       = $toolspath."spline -q";
    $noon       = $toolspath."noon";

    # uvspec related initialization
    $uvspec_inp            = "tmp_UVSPEC.INP"; # Name of uvspec input file
    $uvspec_home           = "../";
    $uvspec                = $toolspath."uvspec";
    $uvspec_mpi            = $toolspath."uvspec_mpi";
    $uvspec_data_file_path = $uvspec_home."data/";
    $uvspec_solar_file     = $uvspec_data_file_path."solar_flux/apm_1nm";
#    $uvspec_solar_file     = $uvspec_data_file_path."solar_flux/atlas_plus_modtran";
    $uvspec_wvn1           = 290.0;
    $uvspec_wvn2           = 700.0;
    $uvspec_zout           = -999;
    $uvspec_atm_type       = "1";
    $atmmod_file           = $uvspec_data_file_path."atmmod/afglus.dat";
    $uvspec_o3_crs         = "Molina";
    $uvspec_solver         = "disort";
    %RB_response_fnc       = load_RB_response_fnc("RB_RESP.DAT");
    %BAND_response_fnc       = load_BAND_response_fnc("GUV_spectral_response313.txt");
#    %BAND_response_fnc       = load_BAND_response_fnc("GUV_spectral_response380.txt");
    $aerosol_ssa           = -1.0;
    $aerosol_gg            = -2.0;
} 

1;

sub DailyUVDose
{
    my(@input) = @_;
    $input = @input;
    if ( $input != 14 ) { usage_DailyUVDose();}
    my($longitude, $longitude_standard, $latitude, $altitude, $mday, $month,
	  $albedo, $ozone_column, $alpha, $beta, $wctau, $wc_file, $wc_cloudcover, $dosetypes) = @input;

    # Calculate sunrise and sunset times in standard time.

    my $sza  = 90.0; # This must be changed for locations polewards
                        # of the arctic circles during the midnigth sun period.
    my(@times) = UpsunAndDownsun($latitude, $longitude,
				    $longitude_standard,
				    $month, $mday, $sza);
    my $begin_day = $times[0];
    my $end_day = $times[1];
    my $day_length = $end_day - $begin_day;
    my $number_of_time_steps = 0;
    if ( $day_length<= 0 ) {
	$number_of_time_steps = 0;
    }
    else {
	$number_of_time_steps = sprintf("%d",$day_length/$time_resolution + 1);
    }

    $dailydoses{"CIE"} = 0;
    $dailydoses{"RB"} = 0;
    $dailydoses{"BAND"} = 0;
    $dailydoses{"UVA"} = 0;
    $dailydoses{"UVB"} = 0;
    $dailydoses{"PAR"} = 0;
    $dailydoses{"CODEGG"} = 0;
    $dailydoses{"CALANUS"} = 0;
    $dailydoses{"305"} = 0;
    $doserates{"CIE"} = 0;
    $doserates{"RB"} = 0;
    $doserates{"BAND"} = 0;
    $doserates{"UVA"} = 0;
    $doserates{"UVB"} = 0;
    $doserates{"PAR"} = 0;
    $doserates{"CODEGG"} = 0;
    $doserates{"CALANUS"} = 0;
    $doserates{"305"} = 0;
#    printf STDERR "begin: $begin_day end: $end_day $number_of_time_steps\n";
    for ($it=0;$it<$number_of_time_steps;$it++) {
#    for ($it=0;$it<2;$it++) {

	# Calculate sza for each timestep between sunrise and sunset.

	local($time) = $begin_day + $it*$time_resolution;
	$sza         = SolarZenithAngle($month, $mday, $time, 
				   $latitude, $longitude, $longitude_standard);

#	printf STDERR "it: $it day: $mday, month: $month, sza: $sza,  time: $time";
#	printf STDERR "\n";

	# Calculate spectrum for $sza.
        # Other input parameters to RT model may be set in spectrum

	local($UVfile) = "UVSPEC.OUT";
	local($day_of_year) = MMDD2dayno($month,$mday);

	local($wvn1) = $uvspec_wvn1;
	local($wvn2) = $uvspec_wvn2;
#	local($atmmod_file) = "../atmmod/afglus.dat";
#	local($wc_file) = "";
#	local($wc_tau) = 0.0;
	Spectrum($albedo, $day_of_year, $ozone_column, $alpha, $beta, $sza, 
		 $wvn1, $wvn2, $atmmod_file, $altitude, $UVfile, $wctau, $wc_file, $wc_cloudcover);

	# Calculate the wanted doserates.

	%doserates = DoseRates($UVfile, $dosetypes, %doserates);
#	print STDERR %doserates, "\n";

    }

    # Integrate to get daily doses.

#    print "doserates ",%doserates,"\n";

    $dosetypes =~ s/^[ ]*//;
    @dosetypes = split(/[ ]+/,$dosetypes);
    $dt = $time_resolution*3600.; # Convert from hours to seconds
                                  # doserates is in something per second

#    print "timeinfo ", $time_resolution, " ", $number_of_time_steps,"\n";

#    for ($d=0;$d<=$#dosetypes;$d++) {
#	printf STDERR "     %-6s ",$dosetypes[$d];
#    }
#    printf STDERR "\n";
    $factor = 1./1000.0; # Convert from milliJoules to Joules.
    for ($d=0;$d<=$#dosetypes;$d++) {
	$dailydoses{$dosetypes[$d]} = $doserates{$dosetypes[$d]}*$dt*$factor;
#	printf STDERR " %10.4f ", $dailydoses{$dosetypes[$d]};
    }
#    printf STDERR "\n";

    %dailydoses;
}

sub usage_DailyUVDose
{
    printf STDERR 
"Usage: DailyUVDose(longitude, longitude_standard, latitude, altitude,
		   daynumber, latitude, daynumber, month, albedo, 
		   ozone_column, alpha, beta, wctau, wc_file, dosetypes)\n";
    exit;
}

sub UVIndex
{
    my(@input) = @_;
    $input = @input;
    if ( $input != 10 ) { usage_UVIndex();}
    my($longitude, $longitude_standard, $latitude, $altitude, $day, $month,
	  $albedo, $ozone_column, $alpha, $beta) = @input;

    # Calculate noon solar zenith angle

    my $year = year();
    my($ans) = `noon $day $month -y $year -p -a $latitude -o $longitude -s $longitude_standard`;
    my ($t1,$t2,$sza)   = split(/ +/,$ans);

    # Calculate spectrum for $sza.
    # Other input parameters to RT model may be set in spectrum

    my $UVfile = "UVSPEC.OUT";
    my $day_of_year = MMDD2dayno($month,$day);
    my $wvn1 = $uvspec_wvn1;
    my $wvn2 = $uvspec_wvn2;
    my $wc_file = "";
    my $wc_tau = 0.0;
    Spectrum($albedo, $day_of_year, $ozone_column, $alpha, $beta, $sza, 
	     $wvn1, $wvn2, $atmmod_file, $altitude, $UVfile, $wctau, $wc_file);

# 	# Calculate the wanted doserates.

    my $dosetypes = "CIE";
    my %doserates = {};
    %doserates = DoseRates($UVfile, $dosetypes, %doserates);

    return ($sza,$doserates{"CIE"}*40./1000.); # Divide by 1000 to get from mW/m2 to W/m2
                                               # Multiply with 40 to get UVindex
}

sub usage_UVIndex
{
    printf STDERR 
"Usage: UVIndex(longitude, longitude_standard, latitude, altitude,
		   daynumber, latitude, daynumber, month, albedo, 
		   ozone_column, alpha, beta)\n";
    exit;
}

sub DoseRates
{
    local(@input) = @_;
    $input = @input;
#    print "inpt @input\n";
    if ( $input < 3 ) { usage_DoseRates();}
    local($UVfile, $dosetypes, %doses) = @_;

    $dosetypes =~ s/^[ ]*//;
    @dosetypes = split(/[ ]+/,$dosetypes);

    # Get global irradiance

    local(%irr) = 0; 
    open(IRR, $UVfile);
    local($wvn1,$dwvn) = 0.0;
    while($line=<IRR>) {
	$line =~ s/^[ ]*//;
	chop $line;
	local($wvn,$rfldir,$rfldn) = split(/[ ]+/, $line);
	$wvn = sprintf("%7.3f",$wvn);
	$dwvn = $wvn - $wvn1;
	$wvn1 = $wvn;
	local($totdn) = $rfldir + $rfldn;
	$irr{$wvn} = $totdn;
#	print "wvn $dwvn $wvn @irr{$wvn}\n";
    }
    close(IRR);

    # Calculate the different doses

    for ($d=0;$d<=$#dosetypes;$d++) {
	$dose = $dosetypes[$d];
#	print $d," $dose\n";
      DOSE: {
	  if ($dose eq "CIE" ) {$doses{"CIE"}  += dose_CIE($dwvn,%irr); }
	  if ($dose eq "RB" )  {$doses{"RB"}   += dose_RB($dwvn,%irr); }
	  if ($dose eq "BAND" )  {$doses{"BAND"}   += dose_BAND($dwvn,%irr); }
#; print "RBDOSES @doses{RB}\n";}
	  if ($dose eq "UVA" ) {$doses{"UVA"}  += dose_UVA($dwvn,%irr);}
	  if ($dose eq "UVB" ) {$doses{"UVB"}  += dose_UVB($dwvn,%irr);}
	  if ($dose eq "PAR" ) {
	      $doses{"PAR"}  += dose_PAR($dwvn,%irr);
	      if ($wvn < 680) {
		  print "warning: PAR computed with irradiances at wavelengths only much less than 700nm?, end wl: $wvn \n";
	      }
	  }
	  if ($dose eq "CODEGG" ) {$doses{"CODEGG"}  += dose_CODEGG($dwvn,%irr);}
	  if ($dose eq "CALANUS" ) {$doses{"CALANUS"}  += dose_CALANUS($dwvn,%irr);}
	  if ($dose eq "305" ) {$doses{"305"}  += dose_305($dwvn,%irr);}
      }
    }
#    print "doses ", %doses, "\n";
    %doses;
}

sub daysinmonth
{
    my(@input) = @_;
    $input = @input;
    if ( $input != 1 ) { usage_daysinmonth();}
    my($month) = @input;
    $month = sprintf("%02d",$month);

    my %days;
    $days{"01"} = 31;
    $days{"02"} = 28; 
    $days{"03"} = 31; 
    $days{"04"} = 30; 
    $days{"05"} = 31; 
    $days{"06"} = 30; 
    $days{"07"} = 31; 
    $days{"08"} = 31; 
    $days{"09"} = 30; 
    $days{"10"} = 31;
    $days{"11"} = 30;
    $days{"12"} = 30;

    $days{$month};
}

sub usage_daysinmonth
{
    printf STDERR 
"
Usage:  daysinmonth(month)

Returns the number of days in month.

Does not take into account leap years.

";
    exit;
}

sub MMDD2dayno
{
    my(@input) = @_;
    $input = @input;
    if ( $input != 2 ) { usage_MMDD2dayno();}
    my($month,$mday) = @input;

    my @days=();
    $days[0] = 0;
    $days[1] = 31 + $days[0];
    $days[2] = 28 + $days[1];
    $days[3] = 31 + $days[2];
    $days[4] = 30 + $days[3];
    $days[5] = 31 + $days[4];
    $days[6] = 30 + $days[5];
    $days[7] = 31 + $days[6];
    $days[8] = 31 + $days[7];
    $days[9] = 30 + $days[8];
    $days[10] = 31 + $days[9];
    $days[11] = 30 + $days[10];

    my $dayno = $days[$month-1]+$mday;
}

sub usage_MMDD2dayno
{
    printf STDERR 
"
Usage:  MMDD2dayno(month, mday)

Calculates the number of day in a year given month and mday.

Does not take into account leap years.

";
    exit;
}

sub dayno2MMDD
{
    local(@input) = @_;
    $input = @input;
    if ( $input != 1 ) { usage_MMDD2dayno();}
    local($dayno) = @input;

    $days[0] = 0;
    $days[1] = 31 + $days[0];
    $days[2] = 28 + $days[1];
    $days[3] = 31 + $days[2];
    $days[4] = 30 + $days[3];
    $days[5] = 31 + $days[4];
    $days[6] = 30 + $days[5];
    $days[7] = 31 + $days[6];
    $days[8] = 31 + $days[7];
    $days[9] = 30 + $days[8];
    $days[10] = 31 + $days[9];
    $days[11] = 30 + $days[10];
    $days[12] = 31 + $days[11];

    my $i;
    for ($i=0;$<=$#days;$i++) {
	if ($dayno < $days[$i]+1) { last;}
    }
    my $dd;
    my $mm;
    if ($i==0) {
	$mm = sprintf("%2d",$i);
	$dd = sprintf("%2d",$dayno);
    }
    else {
	$mm = sprintf("%2d",$i);
	$dd = sprintf("%2d",$dayno-$days[$i-1]);
    }
    @ret = ($mm,$dd);
}

sub usage_dayno2MMDD
{
    printf STDERR 
"
Usage:  dayno2MMDD(dayno)

Calculates month and day given the number of day in a year.

Does not take into account leap years.

";
    exit;
}

sub yyyymmdd2DD
{
    local(@input) = @_;
    $input = @input;
    if ( $input != 1 ) { usage_yyyymmdd2MMDD();}
    local($yyyymmdd) = @input;

    my $dd = substr($yyyymmdd,6,2);
    return $dd;
}

sub yyyymmdd2MM
{
    local(@input) = @_;
    $input = @input;
    if ( $input != 1 ) { usage_yyyymmdd2MMDD();}
    local($yyyymmdd) = @input;

    my $mm = substr($yyyymmdd,4,2);
    return $mm;
}

sub SolarZenithAngle
{
    # Calculate sza for a given time and location

    local(@input) = @_;
    $input = @input;
    if ( $input != 6 ) { usage_SolarZenithAngle();}
    local($month, $mday, $time, $latitude, $longitude, $longitude_standard) = @_;

    $hour = sprintf("%d",$time);
    $min  = sprintf("%d",($time-$hour)*60.);
#    if ($min < 0 ) {$min=0;}
    $sec  = sprintf("%d",(($time-$hour)*60.-$min)*60.);
#    if ($sec < 0 ) {$sec=0;}
    my $year = year();
    while ($hour>=24.){$hour-=24.}  # Fix it so hour is always within [0,24>. 
    while ($hour<0.){$hour+=24.}    # Fix due to Ola Engelsen at NILU.no
    my $ret = system("zenith -p -a $latitude -o $longitude -s $longitude_standard -y $year $mday $month $hour $min $sec > $zenith_out");
    if ($ret != 0) {croak "Unable to gabbab call $zenith";} 

    
    open(SZA, $zenith_out);
    $line = <SZA>;
    close(SZA);
    system("rm -f $zenith_out");
    $line =~ s/^[ ]*//;
    chop $line;
    local($dummy, $sza, $azimuth) = split(/[ ]+/, $line);
    $sza;
}

sub usage_SolarZenithAngle
{
    printf STDERR 
"
Usage:  SolarZenithAngle(month, mday, time, 
                         latitude, longitude, longitude_standard)

Calculates the solar zenith angle for a given time and location.

day is day of month [1,28/30/31]
time is in hours CUT
latitude is in degrees, North positive [-90,90]
longitude is in degrees, West positive [-180,180]
";

    exit;
}

sub Global_Spectrum
{
  my (@input) = @_;
  my $input = @input;
  if ( $input < 12 ) { usage_Global_Spectrum();}
  my ($albedo, $day, $ozone_column, $alpha, $altitude, $beta, 
      $sza, $step, $atmmod_file, $UVfile, $wctau, $wc_file, $umu, $phi) = @input;

  if ( $uvspec_zout < 0.0 ) {
    $uvspec_zout = get_zout($atmmod_file);
  }

  my $global = 1;
  write_uvspec_inputfile($albedo, $day, $ozone_column, $alpha, $beta, 
			 $sza, $wvn1, $wvn2, $atmmod_file, $altitude, $wctau, $wc_file,
			 $umu, $phi,
			 $slitfunction, $wvnc, $wvnc, $step, $global);
  run_uvspec($UVfile);
}

sub Convolved_Spectrum
{
    local(@input) = @_;
    $input = @input;
    if ( $input != 18 ) { usage_Convolved_Spectrum();}
    local($albedo, $day, $ozone_column, $alpha, $altitude, $beta, 
	  $sza, $wvn1, $wvn2, $wvnc, $step, $atmmod_file, 
	  $slitfunction, $UVfile, $wctau, $wc_file, $umu, $phi) = @input;

    if ( $uvspec_zout < 0.0 ) {
	$uvspec_zout = get_zout($atmmod_file);
    }

    write_uvspec_inputfile($albedo, $day, $ozone_column, $alpha, $beta, 
			   $sza, $wvn1, $wvn2, $atmmod_file, $altitude, $wctau, $wc_file,
			   $umu, $phi,
			   $slitfunction, $wvnc, $wvnc, $step);

    run_uvspec($UVfile);
}

sub usage_Convolved_Spectrum
{
    printf STDERR 
"
Usage:  Convolved_Spectrum(albedo, day_of_year, ozone_column, alpha, altitude,
                           beta, sza, wvn1, wvn2, wvnc, step, atmmod_file, 
                           slitfunction, UVfile, wctau, wc_file, umu, phi)

Calculates a surface UV spectrum for a given albedo, altitude, day of year,
ozone column and solar zenith angle between wavelengths wvn1 and wvn2. 
The spectrum is convolved with slitfunction and irradiance returned between 
from wvnc in step resolution. The atmosphere file used is atmmod_file. 
Result is output to the file UVfile.

Aerosols are included if beta is larger than zero.


Clouds may be included by specifying a water cloud input file wc_file. The
cloud optical depth may be set by wctau.

If any other changes want to be made, please modify the sub
write_uvspec_inputfile.

";
    exit;
}
sub Spectrum
{
    local(@input) = @_;
    $input = @input;
    if ( $input != 14 ) { usage_Spectrum();}
    local($albedo, $day, $ozone_column, $alpha, $beta, 
	  $sza, $wvn1, $wvn2, $atmmod_file, $altitude, $UVfile, $wctau, $wc_file, $wc_cloudcover) = @input;

    if ( $uvspec_zout < 0.0 ) {
	$uvspec_zout = get_zout($atmmod_file);
    }

    write_uvspec_inputfile($albedo, $day, $ozone_column, $alpha, $beta, 
			   $sza, $wvn1, $wvn2, $atmmod_file, $altitude, $wctau, $wc_file, $wc_cloudcover);
    run_uvspec($UVfile);
}

sub usage_Spectrum
{
    printf STDERR 
"
Usage:  Spectrum(albedo, day_of_year, ozone_column, alpha, 
                 beta, sza, wvn1, wvn2, atmmod_file, altitude, UVfile, wctau, wc_file, $wc_cloudcover)

Calculates a surface UV spectrum for a given albedo, day of year,
ozone column and solar zenith angle between wavelengths wvn1 and wvn2. 
The atmosphere file used is atmmod_file. Result is output to the 
file UVfile.

Aerosols are included if beta is larger than zero.

Clouds may be included by specifying a water cloud input file wc_file. The
cloud optical depth may be set by wctau.

If any other changes want to be made, please modify the sub
write_uvspec_inputfile.

";
    exit;
}

sub UpsunAndDownsun
{

    # Calculate upgoing and downgoing times in standard time
    # for a given sza.

    my ($latitude, $longitude, $longitude_standard, $month, $mday, $sza) = @_;
    my $zen = 180;
    my ($begin_day,$end_day) =0;
    if (abs($latitude) >= 66.5 ) {
	# Check to see if sun below horizon at noon, if so we have 24 hrs darkness.
        my $year = year();
	my($ans) = `$noon $mday $month -y $year -p -a $latitude -o $longitude -s $longitude_standard`;
	($time,$noont,$azi) = split(/ +/,$ans);
	($hh,$mm,$ss) = split(/:/,$time);
	($ans) = `$zenith $mday $month  $hh $mm $ss -y $year -p -a $latitude -o $longitude -s $longitude_standard`;
	($time,$zen,$azi) = split(/ +/,$ans);
	if ($zen > 90 ) {
	    return ($begin_day,$end_day);
	}
	# Check to see if sun above horizon at midnight, if so we have midnight sun.
	$hh = $hh+12;
	while ($hh>=24.){$hh-=24.}  # Fix it so hour is always within [0,24>. 
	while ($hh<0.){$hh+=24.}    # Fix due to Ola Engelsen at NILU.no
	my($ans) = `$zenith $mday $month  $hh $mm $ss -y $year -p -a $latitude -o $longitude -s $longitude_standard`;
	($time,$zen,$azi) = split(/ +/,$ans);
    }
    if (abs($latitude) >= 66.5 && $zen < 90.0   ) {
	$begin_day = 0.0;
	$end_day   = 23.99;
    }
    else {
    	my $ret = system("$sza2time -p -a $latitude -o $longitude -s $longitude_standard $mday $month $sza > $sza2time_out");
	if ($ret != 0) {croak "Unable to call $sza2time ";} 
	open(TIME, $sza2time_out);
	my $line = <TIME>;
	close(TIME);
	system("rm -f $sza2time_out");

	$line =~ s/^[ ]*//;
	chop $line;
	$line =~ s/^[ ]+//;
	($sza, $begin_day,$end_day) = split(/[ ]+/, $line);
	#	if ( $begin_day<0) {$begin_day=0;} # Commented out due to a suggestion by Ola Engelsen
    }
    @times = ($begin_day, $end_day);
}

sub usage_UpsunAndDownsun
{
    printf STDERR 
"
Usage: UpsunAndDownsun(latitude, longitude, longitude_standard,
                         month, mday, sza )

Calculate upgoing and downgoing times in standard time for a given sza.

mday is day of month [1,28/30/31]
latitude is in degrees, North positive [-90,90]
longitude is in degrees, West positive [-180,180]
";

    exit;
}

sub dose_CIE
{

    # Calculate CIE weighted dose

    my ($dwvn,%irr) = @_;

    my $dose = 0.0;
    my @keys = sort(keys(%irr));
    for (my $i=0;$i<=$#keys;$i++) {
        $wvn = $keys[$i];
        if ($wvn < 280 ) { next; }
        if ($wvn > 400 ) { last; }
        $dose += $irr{$wvn}*erythema_weight($wvn);
#       print "CIE $wvn $dwvn $dose @irr{$wvn}\n";
    }
    $dose*$dwvn;

}

sub erythema_weight{

# Returns the erythema action spectrum value for wavelength $wvn.       
# From Table 1 Dahlback et al. 1989, Photochemistry and Photobiology, 49, 621-625.

    my $wvn = shift;
    my $weight;
    if ($wvn <= 298.0 ) {$weight = 1.0;}
    if ($wvn >= 299.0 && $wvn <= 328.0 ) {$weight = exp(0.2164*(298.0-$wvn));}
    if ($wvn >= 329.0 && $wvn <= 400.0 ) {$weight = exp(0.0345*(139.0-$wvn));}
    $weight;
}

sub dose_RB
{

    # Calculate RB weighted dose


    my ($dwvn,%irr) = @_;
    my $dose = 0.0;
    my @rb_keys = sort(keys(%RB_response_fnc));
    for (my $i=0;$i<=$#rb_keys;$i++) {
        my $wvn = $rb_keys[$i];
        $dose += $irr{$wvn}*$RB_response_fnc{$wvn};
#       print "RB $wvn $dwvn $dose @irr{$wvn} @RB_response_fnc{$wvn}\n";
    }
    $dose*$dwvn;
}


sub dose_UVB
{

    # Calculate UVB dose

    my ($dwvn,%irr) = @_;

    my $dose = 0.0;
    my @keys = sort(keys(%irr));
    for (my $i=0;$i<=$#keys;$i++) {
        my $wvn = $keys[$i];
        if ($wvn < 280 ) { next; }
        if ($wvn > 315 ) { last; }
        $dose += $irr{$wvn};
#       print "UVB $wvn $dwvn $dose @irr{$wvn}\n";
    }
    $dose*$dwvn;
}

sub dose_UVA
{

    # Calculate UVA dose

    local($dwvn,%irr) = @_;

    my $dose = 0.0;
    my @keys = sort(keys(%irr));
    for (my $i=0;$i<=$#keys;$i++) {
        $wvn = $keys[$i];
        if ($wvn < 315 ) { next; }
        if ($wvn > 400 ) { last; }
        $dose += $irr{$wvn};
#       print "$wvn $dose @irr{$wvn}\n";
    }
    $dose*$dwvn;
}

sub dose_PAR
{

    # Calculate PAR dose

    my ($dwvn,%irr) = @_;

    my $dose = 0.0;
    my @keys = sort(keys(%irr));
    for (my $i=0;$i<=$#keys;$i++) {
        $wvn = $keys[$i];
        if ($wvn < 400 ) { next; }
        if ($wvn > 700 ) { last; }
        $dose += $irr{$wvn};
#       print "$wvn $dose @irr{$wvn}\n";
    }
    $dose*$dwvn;
}


sub dose_BAND
{

    # Calculate RB weighted dose

    local($dwvn,%irr) = @_;
    local($dose)=0.0;
    local($sum)=0.0;
    local(@BAND_keys) = sort(keys(%BAND_response_fnc));
    for ($i=0;$i<=$#BAND_keys;$i++) {
	$wvn = $BAND_keys[$i];
	$dose += $irr{$wvn}*$BAND_response_fnc{$wvn};
	$sum+=$BAND_response_fnc{$wvn}
    }
    $dose/=$sum;
}

sub dose_CODEGG
{
# Calculate cod egg mortality weighted dose
#(Kouwenberg et al., 1999, Marine Biology).

    my ($dwvn,%irr) = @_;

    my $dose = 0.0;
    my @keys = sort(keys(%irr));
    for (my $i=0;$i<=$#keys;$i++) {
        $wvn = $keys[$i];
        if ($wvn < 280 ) { next; }
        if ($wvn > 400 ) { last; }
        $dose += $irr{$wvn}*codegg_weight($wvn);
    }
    $dose*$dwvn;
}

sub codegg_weight{

# Returns the cod egg mortality action spectrum value for wavelength $wvn.       
#(Kouwenberg et al., 1999, Marine Biology).
    my $wvn = shift;
    my $weight;
    $weight  = exp(-(10.56+0.121*($wvn-290.)));
    $weight;
}

sub dose_CALANUS
{
# Calculate Calanus finmarchicus (zooplankton) mortality weighted dose
#(Kouwenberg et al., 1999, Marine Biology).

    my ($dwvn,%irr) = @_;

    my $dose = 0.0;
    my @keys = sort(keys(%irr));
    for (my $i=0;$i<=$#keys;$i++) {
        $wvn = $keys[$i];
        if ($wvn < 280 ) { next; }
        if ($wvn > 400 ) { last; }
        $dose += $irr{$wvn}*calanus_weight($wvn);
    }
    $dose*$dwvn;

}

sub calanus_weight{

# Returns the Calanus finmarchicus mortality (zooplankton) action spectrum value for wavelength $wvn.       
#(Kouwenberg et al., 1999, Marine Biology).
    my $wvn = shift;
    my $weight;
    $weight  = exp(-(7.083+0.134*($wvn-290.)));
    $weight;
}

sub dose_305
{

    # Calculate 305nm dose

    my ($dwvn,%irr) = @_;

    my $dose = 0.0;
    my @keys = sort(keys(%irr));
    for (my $i=0;$i<=$#keys;$i++) {
        my $wvn = $keys[$i];
        if ($wvn < 305 ) { next; }
        if ($wvn > 305 ) { last; }
        $dose += $irr{$wvn};
#       print "305nm $wvn $dwvn $dose @irr{$wvn}\n";
    }
    $dose*$dwvn;
}

sub get_zout {
    local($file) = @_;
    open(FILE,$file) || die "read_bandpass: Could not open file: $file\n";
    local(@lines) = <FILE>;
    close(FILE);
    local($line) = $lines[$#lines];
    $line =~ s/^[ ]*//;
    local($zout, $rest) =  split(/[ \t]+/,$line);
    $zout;
}

sub run_uvspec
{
    local($out_file)  = @_;
    $| = 1;
#    print "run uvspec > $out_file\n";
    $uvspec                = $toolspath."uvspec";
    my $ret = system("$uvspec < $uvspec_inp > $out_file");
    if ($ret != 0) {croak "Unable to call $uvspec";} 
}

sub write_uvspec_inputfile
{

  local($albedo, $day, $ozone_column, $alpha, $beta, 
	$sza, $wvn1, $wvn2, $atmmod_file, $altitude, $wctau, $wc_file, $wc_cloudcover) = @_;
  
  open(OUT, ">$uvspec_inp")  || die "\nwrite_uvspec_inputfile: couldn't open $uvspec_inp\n";
  
  print OUT "solar_file $uvspec_solar_file\n";    
  printf OUT "sza %f\n", $sza;    
  printf OUT "albedo %f\n", $albedo;    
  printf OUT "data_files_path %s\n", $uvspec_data_file_path;
  printf OUT "atmosphere_file %s\n", $atmmod_file;
  printf OUT "altitude %f\n", $altitude;
  printf OUT "ozone_column %f\n", $ozone_column;    
  printf OUT "day_of_year %d\n", $day;
  printf OUT "rte_solver $uvspec_solver\n";
  printf OUT "o3_crs $uvspec_o3_crs\n";
  printf OUT "wvn %f %f\n", $wvn1,  $wvn2;
  printf OUT "zout %f\n", $uvspec_zout;
  #    printf OUT "zout %f\n", $altitude; # always compute surface irradiances (check!!)
  
  if ( $beta > 0  ) { 
    print OUT "aerosol_vulcan 1\n";
    print OUT "aerosol_haze 4\n";
    print OUT "aerosol_season 1\n";
    print OUT "aerosol_visibility 100.0\n";
    print OUT "angstrom $alpha $beta\n";
    if ( $aerosol_ssa > -1.0 ) {
      print OUT "aerosol_set_ssa $aerosol_ssa\n";
    }
    if ( $aerosol_gg > -2.0 ) {
      print OUT "aerosol_set_gg $aerosol_gg\n";
    }
  }
  if ($wc_file) {
    print OUT "wc_file $wc_file\n";
    if ($wctau >= 0.0) {
      printf OUT "wc_set_tau %f\n", $wctau;
    }
    if ($wc_cloudcover > 0.0) {
      printf OUT "wc_cloudcover %f\n", $wc_cloudcover;
    }
  }
  
  if ( $wctau > 0.0 ) {
    print OUT "wc_file $wc_file\n";
    printf OUT "wc_set_tau %f\n", $wctau;
  }
  
  if ( length $slitfunction != '' ) {
    printf OUT "slit_function_file %s\n", $slitfunction;
  }
  
  if ( $wvnlower > 0.0 && length $slitfunction != '' ) {
    printf OUT "spline %f %f %f\n", $wvnlower, $wvnupper, $wvnstep;
  }
  
  # radiance; increase NSTR in order to get reasonable results
  if ( $umu != 0.0 ) {
    printf OUT "nstr 16\n";
    printf OUT "umu %f\n", $umu;
    printf OUT "phi %f\n", $phi;
  }
  
  close(OUT);
}


sub load_response_fnc

{
    local($file) = @_;
    open(RB, $file) || die "Could not open file: $file\n";
    local($wvn,$resp) = 0.0;
    local(%response_fnc);
    undef %response_fnc;
  LINE:
#    print "here $file\n";
    while($line=<RB>)  {
	$_ = $line;                   # Make m match $line
	if (m/^[\#]/) {  next LINE; } # Skip comment lines
	$line =~ s/^[ ]*//;
	chop $line;
	if ( split(/[ \t]+/, $line)==2) { 
	    ($wvn,$resp) = split(/[ \t]+/, $line);
	    $wvn = sprintf("%7.3f",$wvn); 
	    $response_fnc{$wvn} = $resp;
	}
    }
    close(RB);
    %response_fnc;
}

sub year {
  my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);

  $year;
}

    
# Local Variables:
# mode: Perl
# End:
