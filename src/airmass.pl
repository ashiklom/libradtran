#!/usr/bin/perl

# Copyright (C) 2002 Arve Kylling
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

# To see what this thing does try: perl airmass.pl --help

# TODO:
#
#  1) Make it able to handle several wavelength, that is for one sza
#     get airmass for many wavelengths.
#  2) Make it able to handle several umu and phi in one go.             DONE
#  3) Make it able to handle several altitudes.
#  4) Include BrO and OClO                                              DONE
#      - BrO airmass test
#      - Z,p,T file, o3 file, no2 file, bro file.                       DONE
#  5) Make uvspec output integrated optical depth at user levels 



# New files data/atmmod/{afglus_zpt.dat,afglus_o3.dat,afglus_no2.dat,bro.dat,oclo.dat}
#           data/crs/crs_bro_iasb.dat



#  HOWTO add another trace gas into uvspec.
#  
#  1) uvspec.h:  * Increase MOL_NN by one
#                * Add appropriate molecule identifer for trace
#                  gas X:  MOL_X
#                * add X_column and X_scale_factor to output structure
#                * add spectrum identifier (e.g. BRO_IASB and no2_spectrum)
#
#  2) uvspec_lex.l: * Add required coding for input parameters "dens_file" 
#                  "dens_column" and "crs_file".
#                * Add documentation for "dens_file", "dens_column", and "crs_file".
#                * Set default unit "Input.atm.unit[MOL_X] = CM_2;"  
#                * Set default spectrum, e.g. Input.bro_spectrum      = BRO_IASB;
#
#  3) atmosphere.c * Add scaling of profile to function setup_atmosphere
#                * Allocate memory for species in function read_atmosphere
#                * Do stuff in function shift_atmosphere
#                * Do stuff in function interpolate_atmosphere
#  4) ancillary.c* Add MOL_X to calculation of babso in function
#                  optical_properties
#  5) molecular.c* Add function to read cross section for MOL_X
#                * Add call to above function in function setup_crs
#  6)            * Put crs file in data/crs
#                * Put example profile files in data/atmmod
#  7) TEST EVERYTHING MANY TIMES.

use Getopt::Long;
use strict;
use Carp;
use Math::Complex;

my $version = 1.0;    # 16.05.2002, Arve.Kylling@nilu.no

my $sname   = $0;
$sname      =~ s/\.pl//; # rm .pl
$sname      =~ s&.+/&&;  # rm possible path

my $help               = 0;
my $amfX               = "";
my $sza_file           = "";
my $sza0               =-99;
my $sza1               =-99;
my $dsza               =-99;
my $wvl                =-99;
my $pi       = atan2(1,1)*4;

my $ret = GetOptions(
		     '--amfX=s'              => \$amfX,
		     '--dsza=f'              => \$dsza,
		     '--sza0=f'              => \$sza0,
		     '--sza1=f'              => \$sza1,
		     '--sza_file=s'          => \$sza_file,
		     '--wvl=f'               => \$wvl,
                     '--help'                => \$help,
);
if ($ret==0) { exit;}
if ($help)   { usage($version); }

my @szas;
($wvl, @szas) = check_line_options ($sza0, $sza1, $dsza, $sza_file, $wvl, $amfX);

my $header             = 1;
my  @pls = <DATA>;    # Fixed stuff for this run from end of this file
foreach my $sza (@szas) {
  $sza =~ s/\n//; # Get rid of possible new lines from --sza_file option

  # First with absorber
  my $abscol     = 1;
  my $uvinp_file = write_uvspec_input_file($sza, $amfX, $abscol, $wvl);
  my $uvinp_file1 = $uvinp_file."1";
  `mv $uvinp_file $uvinp_file1`;
  my $uvout_file_1 = "tmp_".$sname."_UVSPEC.OUT_1";
  `./uvspec < $uvinp_file1 > $uvout_file_1 2>&1`;
  my ($umu,$uu1,$nphi,$phi,$numu,$tau1,$vcol1,$crs1) = get_rad($uvout_file_1, $amfX);
  my @umu = @$umu;
  my @uu1 = @$uu1;
  my @phi = @$phi;

  # Next without absorber
  $abscol       = 0;
  my $uvinp_file   = write_uvspec_input_file($sza, $amfX, $abscol, $wvl);
  my $uvinp_file0 = $uvinp_file."0";
  `mv $uvinp_file $uvinp_file0`;
  my $uvout_file_0 = "tmp_".$sname."_UVSPEC.OUT_0";
  `./uvspec < $uvinp_file0 > $uvout_file_0 2>&1`;
  my ($umu,$uu0,$nphi,$phi,$numu,$tau0,$vcol0,$crs0) = get_rad($uvout_file_0, $amfX);
  my @umu = @$umu;
  my @uu0 = @$uu0;

  # Calculate airmass
  my $tau           = $tau1-$tau0;

  if ($header) {
    printf "%57s", " ";
    for (my $i=0;$i<=$#umu;$i++) {
      printf "| %22s %12.9f (%6.2f) %22s   ", " ", $umu[$i], 180*acos($umu[$i])/$pi-90, " " ;
    }
    printf "\n";
    $header = 0;
  }

  for (my $k=0;$k<=$#phi;$k++) {
    printf "%5.2f %5.2f %6.2f %12.6e %12.6e %12.6e ", $sza, $phi[$k], $wvl, $tau, $vcol1, $crs1;
    for (my $i=0;$i<=$#umu;$i++) {
      my $amf           = log($uu0[$i][$k]/$uu1[$i][$k])/$tau;
      my $slant_column0 = $amf*$vcol1;
      my $slant_column1 = log($uu0[$i][$k]/$uu1[$i][$k])/$crs1;
      printf " %15.9e %15.9e %12.9f %12.6e %12.6e ", 
	$uu0[$i][$k], $uu1[$i][$k], $amf, $slant_column0, $slant_column1;
    }
    printf "\n";
  }

#  printf "%5.2f %6.2f %12.6e %12.6e %12.9f %12.6e %12.6e %12.9f %12.6e %12.6e\n",
#    $sza, $wvl, $uu0, $uu1, $tau, $vcol1, $crs1, $amf, $slant_column0, $slant_column1;
#  die;

}


sub check_line_options {
  ## Check line options and set various stuff
  my $sza0     = shift @_;
  my $sza1     = shift @_;
  my $dsza     = shift @_;
  my $sza_file = shift @_;
  my $wvl      = shift @_;
  my $amfX     = shift @_;

  my @szas;
  if ($sza0!=-99 && $sza1!=-99 && $dsza!=-99 ) {
    #if ($sza0!=-99 && $sza1!=-99 && $dsza!=-99 && $wvl!=-99 ) {
    my $sza = $sza0;
    my $i   = 0;
    while ($sza <= $sza1) {
      $szas[$i]  = $sza;
      $sza       = $sza0 + $dsza*($i+1);
      $i++;
    }
  }
  elsif ( $sza_file ne "") {
    open(FP,$sza_file) || croak "Could not open file $sza_file";
    @szas = <FP>;
    close(FP);
  }
  else {
    printf stderr "\n";
    printf stderr "The combinations of the options --sza0, --sza1 and --dsza\n";
    printf stderr "or --sza_file are not not correctly specified\n";
  }

  $_ = $amfX;
 SWITCH: {
    /^bro/ && do {
      if ($wvl < 0 ) { $wvl = 355.0; }
      last SWITCH;};
    /^no2/ && do {
      if ($wvl < 0 ) { $wvl = 440.0; }
      last SWITCH;};
    /^o3/ && do {
      if ($wvl < 0 ) { $wvl = 510.0; }
      last SWITCH;};
    /^oclo/ && do {
      if ($wvl < 0 ) { $wvl = 360.0; }
      last SWITCH;};
    /^hcho/ && do {
      if ($wvl < 0 ) { $wvl = 356.0; }
      last SWITCH;};
    {
      croak "$0 can not calculate airmass factor for specie $amfX\n";
    };
  }  
  ($wvl, @szas);
}

sub get_rad {
  my $of   = shift @_;
  my $amfX = shift @_;
  open (FP,$of);
  my @ls   =  <FP>; 
  close(FP);
  my @opt;
  my $crs;
  my $vcol;
  my $crs;
  my $tau;
  if (grep (/ airmass calculations /, @ls)) {
    @opt  = grep (/ airmass calculations /, @ls);
    my @opts = split(/\s+/, @opt[0]);
    $crs   = $opts[$#opts];
    $vcol  = $opts[$#opts-1];
    $tau   = $opts[$#opts-2];
  }
  else {
    @opt  = grep (/AMF/, @ls);
    my @opts = grep (/AMFsum/, @opt);
    my @opts = split(/\|/, @opts[0]);
    my @optc = grep (/crs /, @ls);
    my @optsc = split(/\|/, @optc[0]);
    my @optv = grep (/ column /, @opt);
    my @optsv = split(/\|/, @optv[0]);
    $_ = $amfX;
  SWITCH: {
      /^bro/ && do {
	$tau   = $opts[$#opts-2];
	$crs   = $optsc[$#optsc-2];
	$vcol  = $optsv[$#optsv-2];
	last SWITCH;};
      /^no2/ && do {
	$tau   = $opts[$#opts-3];
	$crs   = $optsc[$#optsc-3];
	$vcol  = $optsv[$#optsv-3];
	last SWITCH;};
      /^o3/ && do {
	$tau   = $opts[$#opts-6];
	$crs   = $optsc[$#optsc-6];
	$vcol  = $optsv[$#optsv-6];
        last SWITCH;};
      /^oclo/ && do {
	$tau   = $opts[$#opts-1];
	$crs   = $optsc[$#optsc-1];
	$vcol  = $optsv[$#optsv-1];
	last SWITCH;};
      /^hcho/ && do {
	$tau   = $opts[$#opts];
	$crs   = $optsc[$#optsc];
	$vcol  = $optsv[$#optsv];
	last SWITCH;};
    }
  }
  
  # Get nphi and numu from input file:
  my @numu = grep(/umu /, @pls);
  $_ = $numu[0]; s/#+.*\n//;  # rm comments 
  $numu[0] = $_;
  my @nphi = grep(/phi /, @pls);
  $_ = $nphi[0]; s/#+.*\n//;
  $nphi[0] = $_;
  my @numu = split(/\s+/,@numu[$#numu]);
  my @nphi = split(/\s+/,@nphi[$#nphi]);
  my $numu = $#numu;
  my $nphi = $#nphi;
  my @phi  = @nphi[1..$#nphi];

  my @umu; my @uu;
  my $ans  =  $ls[$#ls];
  $ans     =~ s/^[ ]+//;
  my (@elem) = split(/\s+/,$ans);
  my $i=$#ls;
  my $j=0;
  while ($#elem == $nphi+1 ) {
    $umu[$j] = $elem[0];
    for (my $k=2;$k<=$#elem;$k++) { $uu[$j][$k-2]  = $elem[$k];}
    my $ans  =  $ls[--$i];
    $ans     =~ s/^[ ]+//;
    (@elem) = split(/\s+/,$ans);
    $numu++;
    $j++;
  }
  $nphi=$#elem-2;

  # my ($umu,$u0u,$uu) = split(/\s+/,$ans);

  (\@umu,\@uu,$nphi,\@phi,$numu,$tau,$vcol,$crs);
}

sub INIT {

  sub write_uvspec_input_file {
    my $sza    = shift @_;
    my $amfX   = shift @_;
    my $abscol = shift @_;
    my $wvl    = shift @_;
    
    my $uvinp_file = "tmp_".$sname."_UVSPEC.INP";
    open(FP,">$uvinp_file")  || croak("Could not open $uvinp_file");
    
    my @plst=@pls;
    if ( $abscol ==0 ) {
      $_ = $amfX;
    SWITCH: {
#	/^bro/ && do {
#	  @plst = grep (!/dens_file[\s]+BRO/,@plst);
#	  last SWITCH;};
#	/^no2/ && do {
#	  @plst = grep (!/dens_file[\s]+NO2/,@plst);
#	  last SWITCH;};
#	/^o3/ && do {
#	  @plst = grep (!/dens_file[\s]+O3/,@plst);
#	  last SWITCH;};
#	/^oclo/ && do {
#	  @plst = grep (!/dens_file[\s]+OCLO/,@plst);
#	  last SWITCH;};
#	/^hcho/ && do {
#	  @plst = grep (!/dens_file[\s]+HCHO/,@plst);
#	  last SWITCH;};
      }
    }

    print  FP @plst;
    
    printf FP "###########################################################\n";
    printf FP "### Everything following these comment lines is changed ###\n";
    printf FP "### by $0 for each call to uvspec. Note that    ###\n";
    printf FP "### for options specified both above and below these    ###\n";
    printf FP "### lines, it is the one below that uvspec finally uses.###\n";
    printf FP "###########################################################\n";
    printf FP "\n";
    printf FP "sza $sza\n";
    printf FP "wvn $wvl $wvl\n";
    
    if ( $abscol ==0 ) {
      $_ = $amfX;
    SWITCH: {
	/^bro/ && do {
	  printf FP "dens_column BRO 0.0\n";
	  last SWITCH;};
	/^no2/ && do {
	  printf FP "dens_column NO2 0.0\n";
	  last SWITCH;};
	/^o3/ && do {
	  printf FP "dens_column O3 0.0\n";
	  last SWITCH;};
	/^oclo/ && do {
	  printf FP "dens_column OCLO 0.0\n";
	  last SWITCH;};
	/^hcho/ && do {
	  printf FP "dens_column HCHO 0.0\n";
	  last SWITCH;};
      }
    }
    close(FP);
    
    $uvinp_file;
  }
}

sub usage {
  my $version = shift @_;
  printf STDERR "\n";
  printf STDERR "This is $0, version $version.\n";
  printf STDERR "$0 does the following:\n";
  printf STDERR " \n";
  printf STDERR "Usage:                                                             \n";
  printf STDERR "        perl $0 [options]                                          \n";
  printf STDERR " \n";
  printf STDERR "$0 requires uvspec libradtran package in the path.                 \n";
  printf STDERR " \n";
  printf STDERR "Output is to stdout and comes in two varieties:                    \n";
  printf STDERR "\n";
  printf STDERR "$0 understands the following options:                              \n";
  printf STDERR "\n";
  printf STDERR "--amfX         : Specify for which absorber to calculate AMF       \n";
  printf STDERR "                 The following is currently implemented:           \n";
  printf STDERR "                 --amfX=o3      (ozone)                            \n";
  printf STDERR "                 --amfX=no2     (nitrogen dioxide)                 \n";
  printf STDERR "                 --amfX=bro     (bromine oxide)                    \n";
  printf STDERR "                 --amfX=oclo                                       \n";
  printf STDERR "--help         : Prints this message.                              \n";
  printf STDERR "\n";
  exit;
}

# Local Variables:
# mode: Perl
# End:


__END__

##############################################################
# This stuff is set at the bottom of  $0

data_files_path ../data/ # Location of internal uvspec data files
                         # Location of atmospheric profile file. 
atmosphere_file ../data/atmmod/afglus_zpt.dat
                         # Location of the extraterrestrial spectrum
solar_file ../data/solar_flux/apm_1nm

albedo 0.2              # Surface albedo

# Radiative transfer solver stuff
nscat 2
rte_solver sdisort       # Radiative transfer equation solver
deltam  on               # delta-M scaling on
nstr  16                 # Number of streams

phi 0                    # Azimuth angle of line of sight
umu -1                   # Cosine of polar angle of line of sight

#dens_file   O3 ../data/atmmod/afglus_o3.dat
#dens_column O3  300         DU   # O3 column in Dobson units

#dens_file   NO2 ../data/atmmod/afglus_no2.dat
#dens_column NO2 5.57104e+15 CM_2 # NO2 column in molecules per cm**2
#dens_column NO2 0.0 CM_2 # NO2 column in molecules per cm**2

#dens_file   BRO ../data/atmmod/bro.dat
#dens_column BRO 0.0 CM_2 # BrO column in molecules per cm**2
#dens_column BRO 1.61722e+13 CM_2 # BrO column in molecules per cm**2

#dens_file   OCLO ../data/atmmod/oclo.dat
#dens_column OCLO 1.8e+13 CM_2 # OClO column in molecules per cm**2

verbose                  # To get total vertical absorption optical depth of atmosphere
transmittance

dens_file   BRO ../examples/AMF_bro_mat_am_1a_uvspec.dat
dens_file   O3  ../examples/AMF_o3_pro_1a_uvspec.dat
dens_file   NO2 ../examples/AMF_no2_pro_1a_uvspec.dat

crs_file    BRO ../examples/AMF_bro_025.xs
crs_file    O3  ../examples/AMF_o3_025.xs
crs_file    NO2 ../examples/AMF_no2_025.xs

atmosphere_file ../examples/AMF_pres_temp_1a_uvspec.dat
