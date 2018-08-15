#!/usr/bin/perl

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

use Carp;
use strict;
use Getopt::Long;
use lib "../flexstor";
use Flexstor;
use Spectrum;

# Restrictions:
#                May only be used for one altitude

my %phi;        # Quantum yields as functions of wavelength
                # Wavelength dependence differs for each phi
my %crs;        # Cross sections as functions of wavelength
                # Wavelength dependence differs for each crs
my %Jrat;       # Photolysis frequencies

my $temp = 293.17;
$phi{"o3(1d)"} = phi_o3_talukdar($temp);
#$phi{"o3(1d)"}->print();

my $type = "4pi";
# Process input options
my $ret = GetOptions(
		     "--type=s",\$type,
		     "--help",\&usage,
);
if ($ret==0) {
    exit;
}
my $uvspec_file = shift;
#my ($hh,$mm)    = split /:/,$uvspec_file;
#$hh = substr($hh,length($hh)-2,2);
#$mm = substr($mm,0,2);
#my $hhmm_dec = $hh + $mm/60.;

# Read actinic flux from uvspec output file
my $S_act  = get_act($uvspec_file,$type);
$temp   = $S_act->temper();
# Output actinic flux resolution to file
my $tmp_wl = "tmp_J_wl";
$S_act->write_key_to_file($tmp_wl);

# Read cross sections and quantum yields

$crs{"o3"}      = crs_o3_bass_and_paur($temp);
#$phi{"o3(1d)"}  = phi_o3_talukdar($temp);
$phi{"o3(1d)"}  = phi_o3_matsumi($temp);
$crs{"no2"}     = crs_no2_davidson($temp);
#$crs{"no2"}     = crs_no2_schneider();
$phi{"no2"}     = phi_no2_jpl();
$crs{"hno2"}    = crs_hno2_bongartz();
$phi{"hno2"}    = phi_hno2();

### Interpolate cross sections and quantum yields to actinic flux resolution
my $key;
foreach $key (keys %phi) {
    $phi{$key} = $phi{$key}->linear_interpolate($tmp_wl);
}

foreach $key (keys %crs) {
    $crs{$key} = $crs{$key}->linear_interpolate($tmp_wl);
}    
#$phi{"o3(1d)"}->print();
#$crs{"no2"}->print(); die;
#$phi{"hno2"}->print();

# Multiply actinic flux, cross section and quantum yield for each J-value and integrate
my $S_mult      = new Spectrum("mult");

# O3 -> O(1d)
$S_mult         = Spectrum::multiply3spectra($crs{"o3"},$phi{"o3(1d)"},$S_act);
$Jrat{"o3(1d)"} = $S_mult->integrate();

# NO2
$S_mult         = Spectrum::multiply3spectra($crs{"no2"},$phi{"no2"},$S_act);
$Jrat{"no2"}    = $S_mult->integrate();

# HNO2
$S_mult         = Spectrum::multiply3spectra($crs{"hno2"},$phi{"hno2"},$S_act);
$Jrat{"hno2"}   = $S_mult->integrate();

# Output
#rintf "%f %e %e\n", $hhmm_dec, $Jrat{"o3(1d)"}, $Jrat{"no2"};
printf "%e %e %e\n", $Jrat{"o3(1d)"}, $Jrat{"no2"}, $Jrat{"hno2"};

exit;

sub get_act {
    my $file = shift;
    my $type = shift; 
  
    my $f        = new Flexstor($file);
    my @table    = $f->table('data');
    my @columns  = $table[0]->columns();
    my @wl       = $columns[0]->data();
    my @sza      = $columns[1]->data();
    my @uavgso   = $columns[6]->data();
    my @uavgdn   = $columns[7]->data();
    my @uavgup   = $columns[8]->data();
    my $temp     = Flexstor::get_sgl($file,"self_surface_temperature");
    my $i;
    my $pi       = 3.14159;

    my %uavg = ();
    for ($i=0;$i<=$#wl;$i++) {
      # Convert from mw m-2 nm -1 to quanta s-1 cm-2 nm -1
      # h*c = 6.6260E-34*2.9979E+08 = 1.986409E-25 Jm, 1 nm = 1.E-09 m.
      # 1 W = 1000 mW, 1 m^2 = 10000 cm^2
      my $fact = ($wl[$i] / 1.9864E-16) / (1000.*10000.); 
      $uavgso[$i] *= $fact*4.*$pi;
      $uavgdn[$i] *= $fact*4.*$pi;
      $uavgup[$i] *= $fact*4.*$pi;
      my $lambda = sprintf("%7.3f",$wl[$i]);
      if ( $type eq "4pi") {
	$uavg{$lambda} = $uavgso[$i]+$uavgdn[$i]+$uavgup[$i];
      }
      elsif ( $type eq "down") {
	$uavg{$lambda} = $uavgso[$i]+$uavgdn[$i];
      }
    }
    
    my $spectrum  = $type;
    my $S         = new Spectrum($type);
    $S->spectrum(%uavg);
    $S->temper($temp);
    return $S;
}

sub crs_hno2_bongartz {
    my $file = "../data/crs/crs_hno2_007.dat";
    open(DAT,$file) || croak "Could not open file $file";
    my @lines = <DAT>;
    close(DAT);
    @lines = grep (!/\#/, @lines);

    my $iv;
    my %crs;
    for ($iv=0;$iv<$#lines;$iv++) {
	my $line = $lines[$iv];
	$line    =~ s/^[ ]*//;
	$line    =~ s/[\n]*//;
	my ($wvl ,$sigma) = split(/[ \t\n]+/,$line);
	my $w = sprintf("%0007.3f",$wvl);
	$crs{$w} = $sigma;
    }

    my $S_crs = new Spectrum("hno2 crs, Bongartz");
    $S_crs->spectrum(%crs);
    return $S_crs;
}

sub phi_hno2 {
    my $file = "../data/crs/qy_hno2.dat";
    open(DAT,$file) || croak "Could not open file $file";
    my @lines = <DAT>;
    close(DAT);
    @lines = grep (!/\#/, @lines);

    my $iv;
    my %phi;
    for ($iv=0;$iv<$#lines;$iv++) {
	my $line = $lines[$iv];
	$line    =~ s/^[ ]*//;
	$line    =~ s/[\n]*//;
	my ($wvl ,$qy) = split(/[ \t\n]+/,$line);
	my $w = sprintf("%0007.3f",$wvl);
	$phi{$w} = $qy;
    }

    my $S_crs = new Spectrum("hno2 qy");
    $S_crs->spectrum(%phi);
    return $S_crs;
}

sub crs_no2_davidson {
    my $temp = shift;
    my $file = "../data/crs/crs_no2_davidson.dat";
    open(DAT,$file) || croak "Could not open file $file";
    my @lines = <DAT>;
    close(DAT);
    @lines = grep (!/\#/, @lines);

    my $tempC = $temp - 273.15; # Convert from Kelvin to Celsius
    my $iv;
    my %crs;
    for ($iv=0;$iv<$#lines;$iv++) {
	my $line = $lines[$iv];
	$line    =~ s/^[ ]*//;
	$line    =~ s/[\n]*//;
	my ($wvl1, $wvl2, $sigma, $a) = split(/[ \t\n]+/,$line);
	my $w = sprintf("%0007.3f",($wvl1+$wvl2)*0.5);
	$crs{$w} = $sigma*1.e-19 + $a*$tempC*1.e-22;
    }

    my $S_crs = new Spectrum("no2 crs, Davidson");
    $S_crs->spectrum(%crs);
    return $S_crs;
}

sub crs_no2_schneider {
    my $file = "../data/crs/crs_no2_012.dat";
    open(DAT,$file) || croak "Could not open file $file";
    my @lines = <DAT>;
    close(DAT);
    @lines = grep (!/\#/, @lines);

    my $iv;
    my %crs;
    for ($iv=0;$iv<$#lines;$iv++) {
	my $line = $lines[$iv];
	$line    =~ s/^[ ]*//;
	$line    =~ s/[\n]*//;
	my ($wvl ,$sigma) = split(/[ \t\n]+/,$line);
	my $w = sprintf("%0007.3f",$wvl);
	$crs{$w} = $sigma;
    }

    my $S_crs = new Spectrum("no2 crs, Schneider");
    $S_crs->spectrum(%crs);
    return $S_crs;
}

sub phi_no2_jpl {
    my $file = "../data/crs/qy_no2_004.dat";
    open(DAT,$file) || croak "Could not open file $file";
    my @lines = <DAT>;
    close(DAT);
    @lines = grep (!/\#/, @lines);

    my $iv;
    my %phi;
    for ($iv=0;$iv<$#lines;$iv++) {
	my $line = $lines[$iv];
	$line    =~ s/^[ ]*//;
	$line    =~ s/[\n]*//;
	my ($wvl ,$qy) = split(/[ \t\n]+/,$line);
	my $w = sprintf("%0007.3f",$wvl);
	$phi{$w} = $qy;
    }

    my $S_crs = new Spectrum("no2 qy, jpl");
    $S_crs->spectrum(%phi);
    return $S_crs;
}

sub crs_o3_bass_and_paur {
#    print "# Using Bass and Paur cross sections\n";
    my $temp = shift;
    my $file = "../data/crs/crs_o3_pab_cf.dat";
    open(DAT,$file) || croak "Could not open file $file";
    my @lines = <DAT>;
    close(DAT);
    @lines = grep (!/\#/, @lines);

    my $iv;
    my @wvl;
    my @C0;
    my @C1;
    my @C2;
    my $T0 = 273.13;
    my %crs;
    my $tdiff = $temp-$T0;
    for ($iv=0;$iv<$#lines;$iv++) {
	my $line = $lines[$iv];
	$line    =~ s/^[ ]*//;
	$line    =~ s/[\n]*//;
	($wvl[$iv],$C0[$iv],$C1[$iv],$C2[$iv]) = split(/[ \t\n]+/,$line);
	my $w = sprintf("%0007.3f",$wvl[$iv]);
	$crs{$w} = ($C0[$iv]+$C1[$iv]*$tdiff+$C2[$iv]*$tdiff*$tdiff)*1.e-20;
    }

    my $S_crs = new Spectrum("o3 crs, Bass and Paur");
    $S_crs->spectrum(%crs);
    return $S_crs;
}

sub phi_o3_matsumi {
#    print "# Using Matusumi o3 quantum yield\n";
# Matusumi et al., JGR, 107, 10.1029/2001JD000510, 2002
  my $temp = shift;
  
  #280-700

  my %phi;

  for (my $iv=280;$iv<=700;$iv++) {
    my $wl = sprintf("%0007.3f",$iv);
    if ( $iv <= 305 ) {
      $phi{$wl} = 0.9;
    }
    elsif ( $iv > 305  && $iv <=328 ) {
      my $X1 = 304.225;
      my $X2 = 314.957;
      my $X3 = 310.737;
      my $w1 = 5.576;
      my $w2 = 6.601;
      my $w3 = 2.187;
      my $A1 = 0.8036;
      my $A2 = 8.9061;
#      my $A2 = .89061;
      my $A3 = 0.1192;
      my $v2 = 825.518;
      my $R  = 0.695;
      my $c  = 0.0765;

      my $q1   = 1;
      my $q2   = exp(-($v2/($R*$temp)));
      my $qrat1 = $q1/($q1+$q2);
      my $qrat2 = $q2/($q1+$q2);

      my $exp1 = (($X1-$wl)/$w1)**4;
      my $exp2 = (($X2-$wl)/$w2)**2;
      my $exp3 = (($X3-$wl)/$w3)**2;

      $phi{$wl} = $qrat1*$A1*exp(-$exp1) + $qrat2*$A2*(($temp/300)**2)*exp(-$exp2)+
	$A3*(($temp/300)**1.5)*exp(-$exp3) + $c;
    }
    elsif ( $iv >  328 && $iv <= 340 ) {
      $phi{$wl} = 0.08;
    }
    elsif ( $iv >  340 ) {
      $phi{$wl} = 0.0;
    }
  }

  my $S_qy = new Spectrum("o3 qy, Matsumi");
  $S_qy->spectrum(%phi);
  return $S_qy;
}

sub phi_o3_talukdar {
#    print "# Using Talukdar o3 quantum yield\n";
    my $temp = shift;
    my $file = "../data/crs/qy_o3_talu.dat";
    open(DAT,$file) || croak "Could not open file $file";
    my @lines = <DAT>;
    close(DAT);
    @lines = grep (!/\#/, @lines);

    my $iv;
    my @wvl;
    my @A;
    my @B;
    my %phi;
    for ($iv=0;$iv<$#lines-1;$iv++) {
	my $line = $lines[$iv];
	$line    =~ s/^[ ]*//;
	$line    =~ s/[\n]*//;
	($wvl[$iv],$A[$iv],$B[$iv]) = split(/[ \t\n]+/,$line);
	my $w = sprintf("%0007.3f",$wvl[$iv]);
	$phi{$w} = 0.06 + $A[$iv] * exp(-$B[$iv]/$temp);
    }
    $iv = $#lines-1;
    my $line = $lines[$iv];
    $line    =~ s/^[ ]+//;
    ($wvl[$iv],$A[$iv],$B[$iv]) = split(/[ ]+/,$line);
    my ($wvl_begin,$wvl_end) = split(/[-]+/,$wvl[$iv]);

    for ($iv=$wvl_begin;$iv<=$wvl_end;$iv++) {
	my $w = sprintf("%0007.3f",$iv);
	$phi{$w} = 0.06;
    }
    $iv = $#lines;
    $line = $lines[$iv];
    $line    =~ s/^[ ]+//;
    ($wvl[$iv],$A[$iv],$B[$iv]) = split(/[ ]+/,$line);
    ($wvl_begin,$wvl_end) = split(/[-]+/,$wvl[$iv]);

    for ($iv=$wvl_begin;$iv<=$wvl_end;$iv++) {
	my $w = sprintf("%0007.3f",$iv);
	$phi{$w} = 0.0;
    }


    my $S_qy = new Spectrum("o3 qy, Talukdar");
    $S_qy->spectrum(%phi);
    return $S_qy;
}

sub usage {
    printf STDERR "\n";
    printf STDERR "Calc_J.pl calculates various J-values (currently only J(O(1D)),    \n";
    printf STDERR "J(NO2) and J(HN02) from a uvspec output file.                      \n";
    printf STDERR "                                                                   \n";
    printf STDERR "Usage:                                                             \n";
    printf STDERR "                                                                   \n";
    printf STDERR "perl Calc_J.pl [options] uvspec_output_file                        \n";
    printf STDERR "                                                                   \n";
    printf STDERR "The following restrictions apply                                   \n";
    printf STDERR "to the uvspec output file:                                         \n";
    printf STDERR "   *) It must be in flexstor format, that is, set the flexstor option.\n";
    printf STDERR "   *) It must only contain results for one altitude.               \n";
    printf STDERR "   *) There must be a line specifying the temperature to be used   \n";
    printf STDERR "      for the calculation. It must be of the form:                 \n";
    printf STDERR "        self_surface_temperature(K) 289.99                         \n";
    printf STDERR "                                                                   \n";
    printf STDERR "Calc_J.pl will output one line, with the following information:    \n";
    printf STDERR "        Field 1: J(O(1D))                                          \n";
    printf STDERR "        Field 2: J(NO2)                                            \n";
    printf STDERR "        Field 3: J(HNO2)                                           \n";
    printf STDERR "                                                                   \n";
    printf STDERR "Calc_J.pl must be placed such that perl will find the Flexstor.pm  \n";
    printf STDERR "and Spectrum.pm files provided in the ../flexstor directory.       \n";
    printf STDERR "                                                                   \n";
    printf STDERR "An example of input and output files to Calc_J is provided in the  \n";
    printf STDERR "examples directory.                                                \n";
    printf STDERR "\n";
    printf STDERR "Calc_J.pl understands the following options:                       \n";
    printf STDERR "\n";
    printf STDERR "--type [type]        : Type may be set to either `4pi` or `down`.  \n";
    printf STDERR "                       `4pi' means that both the upwelling and downwelling\n";
    printf STDERR "                       radiation is used to calculate the J-values. If type\n";
    printf STDERR "                       is set to `down`, only the downwelling radiation\n";
    printf STDERR "                       is used. Default is `4pi`.                  \n";
    printf STDERR "--help               : Prints this message.                        \n";
    printf STDERR "\n";
    exit;
}

# Local Variables:
# mode: Perl
# End:
