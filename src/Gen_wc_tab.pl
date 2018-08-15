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

#####
# The documentation going into the libRadtran documentation is found in 
# the file read_o3_tab.c. This so because sdoc do not read perl yet.
#####

use UVTools;
use Getopt::Long;
use strict;

# The water cloud table is a function of solar zenith angles and  cloudy/clear
# global irradiances ratios. The solar zenith angles and ratios in the table are
# set below. The numbers in the table are cloud optical depths for specific
# combinations of solar zenith angle and ratios.

# Solar zenith angles that will be used to calculate the wc table
my @szas = (0.1, 5,10,15,
	 20,22,24,26,28,
	 30,32,34,36,38,
	 40,42,44,46,48,
	 50,51,52,53,54,55,56,57,58,59,
	 60.01,61,62,63,64,65,66,67,68,69,
	 70,71,72,73,74,75,76,77,78,79,80); 

# The cloudy irradiance values (ary named ratio because this was adopted from 
# Gen_o3_tab, but it is not a ratio, rather an irradiance value, so since I did
# not change the ary name I got to write this silly comment instead:-)

my $rat = 0.0;
my $i   = 0;
my $dr  = 1.0;
my @ratio;
while ($rat<10.0) {
    @ratio[$i++] = $rat;
    $rat += $dr;
}
$dr  = 2.5;
while ($rat<100.0) {
    @ratio[$i++] = $rat;
    $rat += $dr;
}
$dr  = 5.;
while ($rat<1000.0) {
    @ratio[$i++] = $rat;
    $rat += $dr;
}
$dr  = 10.;
while ($rat<=6000.0) {
    @ratio[$i++] = $rat;
    $rat += $dr;
}

# $wc_lower_limit = 0.0;
# $wc_upper_limit = 500.0;
# $wc_step        = 10.0;
# $n_wc           = ($wc_upper_limit - $wc_lower_limit) / $wc_step;


my $day          = 1.;    # Day number must be corrected for latter.
my $ozone        = 300.0; # Ozone does not matter as we do this for wavelengths where
                          # ozone absorption does not occurr.

my $UVfile = "tmp_Gen_wc_tab_conv";

# Default command line options 
my $absolute     = 0;
my $actinic      = 0;
my $albedo       = 0.0;
my $alpha        = 0.0;
my $altitude     = 0.0;
my $atmmod_file  = "../data/atmmod/afglus.dat";
my $wc_file      = "";
my $bandpassfile = "";
my $beta         = 0.0;
my $file         = "wc_table";
my $slitfunction = '';
my $lambda       = 380.0;
my $o3_crs       = "Molina";
my $global       = 0;

my $umu          = 0;
my $phi          = 0;

my $input_line = "$0 @ARGV";

my $ret = GetOptions(
		  "--absolute",\$absolute,
		  "--actinic",\$actinic,
		  "--albedo=f",\$albedo,
		  "--alpha=f",\$alpha,
		  "--altitude=f",\$altitude,
		  "--atmmod=s",\$atmmod_file,
		  "--wc_file=s",\$wc_file,
		  "--bandpass=s",\$bandpassfile,
		  "--beta=f",\$beta,
		  "--file=s",\$file,
		  "--global",\$global,
		  "--help",\&usage,
		  "--slitfunction=s",\$slitfunction,
		  "--o3_crs=s",\$o3_crs,
		  "--lambda=f",\$lambda,
		  "<>",\&usage);
if ($ret==0) {
    exit;
}

my $file_clear = $file.".CLEAR";

# set wc optical depth 

my @wcs = ();

if ( $actinic ) {
  @wcs = (3.0, 4.0, 5.0, 7.5, 10.0,
	  15.0, 20.0, 25.0, 30.0, 40., 50., 60., 80., 100.,
	  120., 140., 160., 180., 200., 225., 250., 275.,
	  300., 350., 400., 450., 500., 550., 600., 650.,
	  700., 750., 800., 850., 900., 950., 1000., 1050.,
	  1100., 1200., 1300., 1400., 1500.
	 );
}
else {
  @wcs = (0.0, 0.05, 1.0, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0,
	  15.0, 20.0, 25.0, 30.0, 40., 50., 60., 80., 100.,
	  120., 140., 160., 180., 200., 225., 250., 275.,
	  300., 350., 400., 450., 500., 550., 600., 650.,
	  700., 750., 800., 850., 900., 950., 1000., 1050.,
	  1100., 1200., 1300., 1400., 1500
	 );
}

my $n_wc = $#wcs;

# A little checking of input data will not hurt

$UVTools::uvspec_o3_crs = $o3_crs;
$UVTools::uvspec_solver = "disort";
$UVTools::nstr          = 16;

printf STDERR "Water cloud table is written to the file $file\n";
if ( $global ) {
  printf STDERR "Global calculation\n";
  if ( length($bandpassfile)> 0 ) { 
    printf STDERR "Do not specify both --bandpass and --global.\n"; die;
  }
  if ( $actinic ) { 
    printf STDERR "The --actinic option may not be combined with --global.\n"; die;
  }
}
else {
  printf STDERR "Wavelength: %f\n", $lambda;
}

# Now do the real things.

my $bandpass = 0;
my ($wvn1, $wvn2, $wvncenter, $step);
if ( length($bandpassfile)> 0 ) {
    $bandpass = 1;

    my ($wvn1band,$wvn2band,$wvnstep) = read_bandpass($bandpassfile);
    
    $wvn1 = $wvn1band;
    $wvn2 = $wvn2band;
    if (!$absolute) {
	$wvn1 = $wvn1 + $lambda;
	$wvn2 = $wvn2 + $lambda;
    } 
    $wvncenter = $wvn1band;
    $step  = $wvnstep;
}
elsif ( !$global) {
    my ($wvn_left, $wvn_right) = (0,0);
    if ( length $slitfunction != '' ) {
      ($wvn_left, $wvn_right) = read_slitfunction($slitfunction);
      if ( abs($wvn_left ) < 1 ) { $wvn_left  = -1.0;}
      if ( abs($wvn_right) < 1 ) { $wvn_right = 1.0;}
    }
    $wvn1      = $lambda + $wvn_left;  # Plus sign because $wvn_left is negative
    $wvn2      = $lambda + $wvn_right;
    $wvncenter = ($wvn2+$wvn1)/2.0;
    $step  = 15.;  # Some large number to make spline only produce one number
}

my $rat_file    = "tmp_Gen_wc_tab_rat";
my $rat_wc_file = "tmp_Gen_wc_tab_ratwc";
my $tmptable    = "tmp_Gen_wc_tab_table";

open(RAT,">".$rat_file);
for (my $i=0;$i<=$#ratio;$i++) {
    printf RAT "%f\n", $ratio[$i];
}
close(RAT);

# Do all the water cloud spectra

open(IRR_CLEAR,">".$file_clear);
my @table;
# Loop over all zenith angles
for (my $i=0;$i<=$#szas;$i++) {
  open(RAT_WC,">".$rat_wc_file);
  # Loop over all wc optical depths
  for (my $j=$n_wc;$j>=0;$j--) {
    
    printf STDERR "Now doing sza: %5.2f and wc: %6.3f\n", $szas[$i], $wcs[$j];
    
    my $irr;
    if ( $global ) {
      UVTools::Global_Spectrum($albedo, $day, $ozone, 
			       $alpha, $altitude, $beta, $szas[$i], 
			       $step, $atmmod_file, 
			       $UVfile, $wcs[$j], $wc_file);
      $irr    = read_global_irr($UVfile);
    }
    else {
      UVTools::Convolved_Spectrum($albedo, $day, $ozone, 
				  $alpha, $altitude, $beta, $szas[$i], 
				  $wvn1, $wvn2, $wvncenter, 
				  $step, $atmmod_file, 
				  $slitfunction, $UVfile, $wcs[$j], $wc_file,
				  $umu, $phi);
      $irr = 0;
      if ($bandpass) {
	$irr = multiply_with_bandpass_and_integrate( $UVfile, $bandpassfile, $actinic);
      }
      else {
	$irr    = read_irr($UVfile, $actinic);
      }
    }
    my $ratio_cloudy_clear = $irr;

    printf  "%f %f\n", $ratio_cloudy_clear, $wcs[$j];    
    printf RAT_WC "%f %f\n", $ratio_cloudy_clear, $wcs[$j];    

    if ($wcs[$j] == 0.0 ) {
      printf IRR_CLEAR "%f %f\n", $szas[$i], $irr;
    }
  }

  close(RAT_WC);

  # Interpolate calculated ratios to ratio
  my $ret = system("$UVTools::spline -x $rat_file $rat_wc_file > $tmptable");
  $ret = $ret/256;
  if ( $ret!=0 ) {
    printf STDERR "Spline command failed: $UVTools::spline -x $rat_file $rat_wc_file > $tmptable\n"; 
    die;
  }
  #  die;
  open(TMPTABLE,$tmptable);
  my @tmp = <TMPTABLE>;
  close(TMPTABLE);
  for (my $k=0;$k<=$#tmp;$k++) {
    my $line = @tmp[$k];
    $line =~ s/^[ ]*//;
    my ($rat,$wc) = split(/[ \t]+/,$line);
    $_ = $wc;
    if ( m/[0-9\.\-\e]+/ ) {
      $table[$k][$i] = $wc;
    }
    else {
      $table[$k][$i] = -1;
    }
  }
}
close(IRR_CLEAR);

# Write the table.  
open(TABLE,">".$file);
printf TABLE "\# Generated by the following command:\n";
printf TABLE "\# $input_line\n";
printf TABLE "%11.1f ", -1.0;
for (my $i=0; $i<=$#szas; $i++) {
    printf TABLE "%11.5f ", $szas[$i];
}
printf TABLE "\n";
  
for (my $k=0;$k<=$#ratio;$k++) {
    printf TABLE "%11.5f ", $ratio[$k];
    for (my $i=0; $i<=$#szas; $i++) {
	printf TABLE "%11.5f ",  $table[$k][$i];
    }
    printf TABLE "\n";
}
close(TABLE);

sub read_global_irr 
{
  my $file = shift @_;
  my $totrad = 0;
  open(FILE,$file);
  while (my $l = <FILE> ) {
    $l =~ s/^[ ]*//;
    my ($wvn, $dir, $diff, $up, $actso, $actdn, $actup) = split(/[ \t]+/,$l);
    $totrad += $dir+$diff;
  }
  close(FILE);
  $totrad;
}

sub read_irr 
{
  my $file = shift @_;
  my $actinic = shift @_;
  open(FILE,$file);
  my $line = <FILE>;
  close(FILE);
  $line =~ s/^[ ]*//;
  my ($wvn, $dir, $diff, $up, $actso, $actdn, $actup) = split(/[ \t]+/,$line);
  my $totrad = 0;
  if ( $actinic ) {
    $totrad = $actso + $actdn;
  }
  else {
    $totrad = $dir + $diff;
  }
  $totrad;
}

sub read_bandpass
{
    my $file = shift @_;
    open(FILE,$file) || die "read_bandpass: Could not open file: $file\n";
    my @lines = <FILE>;
    @lines = grep !/^#/, @lines;
    close(FILE);
    my $line = @lines[0];
    $line =~ s/^[ ]*//;
    my ($wvn1band,$rest) = split(/[ \t]+/,$line);
    $line = @lines[$#lines];
    $line =~ s/^[ ]*//;
    my ($wvn2band,$rest) = split(/[ \t]+/,$line);
    my $wvnstep = ($wvn2band-$wvn1band)/$#lines;
    ($wvn1band,$wvn2band, $wvnstep);
}

sub read_slitfunction
{
  # Returns the upper (right) and lower (left) limits of the slitfunction
    my $file = shift @_;
    open(FILE,$file) || die "Could not open file: $file\n";
    my @lines = <FILE>;
    close(FILE);
    @lines = grep (!/\#/, @lines);  # Rm comment lines
    my ($line) = @lines[0];
    $line =~ s/^[ ]*//;
    my ($wvn_left, $rest) = split(/[ \t]+/,$line);
    $line = @lines[$#lines];
    $line =~ s/^[ ]*//;
    my ($wvn_right, $rest) = split(/[ \t]+/,$line);
    ($wvn_left, $wvn_right);
}

sub multiply_with_bandpass_and_integrate
{
  my $irrfile = shift @_;
  my $bpfile  = shift @_;
  my $actinic = shift @_;

  open(FILE,$irrfile) || die "multiply_with_bandpass_and_integrate: Could not open file: $irrfile\n";
  my @irrlines = <FILE>;
  close(FILE);
  open(FILE,$bpfile) || die "multiply_with_bandpass_and_integrate: Could not open file: $bpfile\n";
  my @bplines = <FILE>;
  @bplines = grep !/^\#/, @bplines;
  close(FILE);
  
  # Multiply irradiance with bandpass
  if ($#irrlines != $#bplines)  { 
    die "multiply_with_bandpass_and_integrate: Bandpass and uvfile have different length\n";
  }
  my @irr;
  my @wvn;
  for (my $i=0;$i<=$#irrlines;$i++) {
    my $line = @irrlines[$i];
    chop $line;
    $line =~ s/^[ \t]+//;
    my ($wvn, $dir, $diff, $up, $actso, $actdn, $actup) = split(/[ \t]+/,$line);
    if ( $actinic ) {
      $irr[$i] = $actso + $actdn;
    }
    else {
      $irr[$i] = $dir + $diff;
    }
    $wvn[$i] = $wvn;
    $line = $bplines[$i];
    chop $line;
    $line =~ s/^[ \t]+//;
    my ($wvnbp, $bp, $gabba) = split(/[ \t]+/,$line);
#	print "wvn  $wvn wvnbp: $wvnbp bp: $bp dir: $dir\n";
#	if ($wvn != $wvnbp)  { 
#	    die "multiply_with_bandpass_and_integrate: Bandpass and uvfile wvn differ $wvn $wvnbp\n";
#	}
#    printf "%d %f %f %f %f\n", $i, $wvnbp, $bp, $irr[$i], $irr[$i]*$bp;
    $irr[$i] = $irr[$i]*$bp;
  }

  # Integrate up the irradiance
  
  my $tmpbp = "tmp_Gen_wc_tab_bp";
  my $tmpirr = "tmp_Gen_wc_tab_irr";
  open(BP,">".$tmpbp);
  for (my $i=0;$i<=$#irrlines;$i++) {
    printf BP "%f %f\n", $wvn[$i], $irr[$i];
  }
  close(BP);
  
  system("$UVTools::integrate -p  $tmpbp > $tmpirr");
  
  open(IRR,$tmpirr);
  my $irr = <IRR>;
  chop $irr;
  close(IRR);
  
  #    system("$UVTools::conv $irrfile $bpfile > tmpconv");
  
  #    die "multi $irr die\n";
  $irr;

}

sub usage {
    printf STDERR "\n";
    printf STDERR "Gen_wc_tab generates a table of cloudy/clear sky irradiance\n";
    printf STDERR "ratios versus solar zenith angle for different water cloud\n";
    printf STDERR "optical depth for a given wavelength. The table is\n";
    printf STDERR "read and water cloud optical depth derived by the read_o3_tab program.\n";
    printf STDERR "\n";
    printf STDERR "The irradiance ratios may be calculated in the following ways:\n";
    printf STDERR "(See the libRadtran documentation for examples.)\n";
    printf STDERR "\n";
    printf STDERR "   * Cloudy/clear ratio calculated from single wavelength\n";
    printf STDERR "     measurements.\n";
    printf STDERR "   * Cloudy/clear ratio calculated from measurements multiplied by\n";
    printf STDERR "     a bandpass function and integrated over the bandpass.\n";
    printf STDERR "\n";
    printf STDERR "Gen_wc_tab understands the following options:\n";
    printf STDERR "\n";
    printf STDERR "--absolute               : The wavelengths of the bandpass file are in\n";
    printf STDERR "                           absolute units. Default is relative units.\n";
    printf STDERR "--albedo <value>         : Lambertian surface albedo. Default is 0.0.\n";
    printf STDERR "--alpha <value>          : Angstrom alpha coefficient. Default is 0.0.\n";
    printf STDERR "--altitude <value>       : Altitude of observation site (km). Default is 0.0.\n";
    printf STDERR "--atmmod <name>          : Name of atmosphere file. Default atmmod/afglus.dat.\n";
    printf STDERR "--wc_file <name>         : Name of water cloud file. Default none.\n";
    printf STDERR "--beta <value>           : Angstrom beta coefficient. Default is 0.0.\n";
    printf STDERR "--help                   : Prints this message.\n";
    printf STDERR "--o3_crs <name>          : Name of o3 cross section to use. Default is\n";
    printf STDERR "                           Molina. See libRadtran documentation for other options.\n";
    printf STDERR "--slitfunction <name>    : Name of slitfunction file.\n";
    printf STDERR "--bandpass <name>        : Name of file holding bandpass for wavelength.\n";
    printf STDERR "--file <name>            : Name of file where the table will be stored.\n";
    printf STDERR "--lambda <value>         : Value for wavelength, in nm.\n";
    printf STDERR "\n";
    exit;
}

# Local Variables:
# mode: Perl
# End:
