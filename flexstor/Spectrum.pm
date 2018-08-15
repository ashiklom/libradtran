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
use lib "../src";
use UVTools;

package Spectrum;
BEGIN {
    $ddd         = "Null";
    $ddd_dhhmm   = "Null";
    $hhmm        = "Null";    
    $mmdd        = "Null";    
    $dd          = "Null";    
    $mm          = "Null";    
    $toolspath    = "../src/";
    $conv         = $toolspath."conv";
    $spline       = $toolspath."spline -q";
    $integrate    = $toolspath."integrate -q";
}
1;

sub new {
    my $that     = shift;
    my $type     = shift;
    my $class    = ref($that) || $that;
    my $self     = { 
	ALBEDO   => $albedo,
	ALPHA    => $alpha,
	BETA     => $beta,
	DD       => $dd,
	DDD      => $ddd,        # Day of year
	DDD_DHHMM=> $ddd_dhhmm,  # Day of year dot decimal hours and minutes
	HHMM     => $hhmm,       # hours and minutes
	MM       => $mm,
	MM_DD    => {($mm,$dd)},
	OZONE    => $ozone,
	SPECTRUM => { %data },
	SZA      => $sza,
	TEMPER   => $temper,
	T_END    => $t_end,
	T_START  => $t_start,
	TYPE     => $type,
	UVA      => $uva,
	UVB      => $uvb,
	YYYYMMDD => $yyyymmdd,
    };
    $self->{TYPE} = $type;
    bless $self, $class;
    return $self;
}

sub print {
    my $self = shift;
    print "print\n";
    foreach $key (sort {$a<=>$b}( keys %{$self->{SPECTRUM}})) {
#	if ($key > 320 && $key < 321) {
#	    printf "%s %8.4f %13.6e\n", $self->{TYPE}, $key, ${ $self->{SPECTRUM}}{$key};
	    print  "$self->{TYPE}, $key, ${ $self->{SPECTRUM}}{$key}\n";
#	}
    }    
}


sub albedo {
    my $self = shift;
    if (@_) {
	my $albedo = shift;
	$self->{ALBEDO} = $albedo;
    }
    else {
	return $self->{ALBEDO};
    }
}

sub alpha {
    my $self = shift;
    if (@_) {
	my $alpha = shift;
	$self->{ALPHA} = $alpha;
    }
    else {
	return $self->{ALPHA};
    }
}

sub beta {
    my $self = shift;
    if (@_) {
	my $beta = shift;
	$self->{BETA} = $beta;
    }
    else {
	return $self->{BETA};
    }
}

sub ddd {
    my $self = shift;
    if (@_) {
	my $ddd = shift;
	$self->{DDD} = $ddd;
    }
    else {
	my $dd = $self->dd();
	my $mm = $self->mm();
	if ( $self->{DDD} eq "Null" ) {
	    $self->{DDD} = UVTools::MMDD2dayno($mm,$dd);
	    return $self->{DDD};
	}
	else {
	    return $self->{DDD};
	}
    }
}

sub ddd_dhhmm {
    my $self = shift;
    if ( $self->{DDD_DHHMM} eq "Null" ) {
	if ( $self->{HHMM} eq "Null" ) {
	    my $ddd = $self->ddd();
	    $self->{DDD_DHHMM} = $ddd.".".sprintf("%004d",10000.*$self->{T_START}/24.);
	}
	else {
	    my $ddd = ddd();
	    my $hhmm = $self->{HHMM};
	    my $hh   = substr($hhmm,0,2);
	    my $mm   = substr($hhmm,2,2);
	    $self->{DDD_DHHMM} = $ddd.sprintf("%5.4f",($hh+($mm/60.))/24.);
	}
	return $self->{DDD_DHHMM};
    }
    else {
	return $self->{DDD_DHHMM};
    }
    return $self->{DDD_DHHMM};
}

sub hhmm {
    my $self = shift;
    if (@_) {
	my $hhmm = shift;
	$self->{HHMM} = $hhmm;
    }
    else {
	return $self->{HHMM};
    }
}

sub dd {
    my $self = shift;
    if (@_) {
	my $dd = shift;
	$self->{DD} = ($mm,$dd);
    }
    else {
	if ( $self->{DD} eq "Null" ) {
	    my $yyyymmdd = $self->{YYYYMMDD};
	    my ($dd) = UVTools::yyyymmdd2DD($yyyymmdd);	    
	    $self->{DD} = $dd;
	}
	return $self->{DD};
    }
}

sub mm {
    my $self = shift;
    if (@_) {
	my $mm = shift;
	$self->{MM} = ($mm,$mm);
    }
    else {
	if ( $self->{MM} eq "Null" ) {
	    my $yyyymmdd = $self->{YYYYMMDD};
	    my ($mm) = UVTools::yyyymmdd2MM($yyyymmdd);	    
	    $self->{MM} = $mm;
	}
	return $self->{MM};
    }
}

sub ozone {
    my $self = shift;
    if (@_) {
	my $ozone = shift;
	$self->{OZONE} = $ozone;
    }
    else {
	return $self->{OZONE};
    }
}

sub sza {
    my $self = shift;
    if (@_) {
	my $sza = shift;
	$self->{SZA} = $sza;
    }
    else {
	return $self->{SZA};
    }
}

sub spectrum{
    my $self = shift;
    if (@_) {
	my %dat = @_;
	$self->{SPECTRUM} = \%dat;
    }
    else {
	return $self;
    }
}

sub temper {
    my $self = shift;
    if (@_) {
	my $temper = shift;
	$self->{TEMPER} = $temper;
    }
    else {
	return $self->{TEMPER};
    }
}

sub data{
    my $self = shift;
    if (@_) {
	my %dat = @_;
	$self->{SPECTRUM} = \%dat;
    }
    else {
	return %{$self->{SPECTRUM}};
    }
}

sub type {
    my $self = shift;
    return $self->{TYPE};
}

sub time {
    my $self = shift;
    if (@_) {
	my $t_start = shift;
	my $t_end   = shift;
	$self->{T_START} = $t_start;
	$self->{T_END}   = $t_end;
    }
    else {
	return ($self->{T_START},$self->{T_END});
    }
}

sub yyyymmdd {
    my $self = shift;
    if (@_) {
	my $yyyymmdd = shift;
	$self->{YYYYMMDD} = $yyyymmdd;
    }
    else {
	return $self->{YYYYMMDD};
    }
}

sub hash {
    my $self = shift;
    return %{$self->{SPECTRUM}};
}

sub set_max {
    my $self = shift;
    my $max = shift;
    my $key;
    foreach $key (sort {$a<=>$b}( keys %{$self->{SPECTRUM}})) {
	if ( ${$self->{SPECTRUM}}{$key} > $max )  { ${$self->{SPECTRUM}}{$key}= $max; }
    }
}

sub set_min {
    my $self = shift;
    my $min = shift;
    my $key;
    foreach $key (sort {$a<=>$b}( keys %{$self->{SPECTRUM}})) {
	if ( ${$self->{SPECTRUM}}{$key} < $min )  { ${$self->{SPECTRUM}}{$key}= $min; }
    }
}

sub x_ary {
    my $self = shift;
    my @x = sort {$a<=>$b}( keys %{$self->{SPECTRUM}});
    return @x;
}

sub y_ary {
    my $self = shift;
    my $i = 0;
    my @y;
    foreach $key (sort {$a<=>$b}( keys %{$self->{SPECTRUM}})) {
 	$y[$i++] = ${$self->{SPECTRUM}}{$key};
    }
    return @y;
}

sub x_max {
    my $self = shift;
    my @x = sort {$a<=>$b}( keys %{$self->{SPECTRUM}});
    return $x[$#x];
}

sub x_min {
    my $self = shift;
    my @x = sort {$a<=>$b}( keys %{$self->{SPECTRUM}});
    return $x[0];
}

sub y_max {
    my $self = shift;
    my @y = sort {$a<=>$b}( values %{$self->{SPECTRUM}});
    return $y[$#y];
}

sub y_max_at_x {
    my $self = shift;
    my @y = sort {$a<=>$b}( values %{$self->{SPECTRUM}});
    my $max = $y[$#y];
    my $xmax=0;
    my $i;
    foreach $key (sort {$a<=>$b}( keys %{$self->{SPECTRUM}})) {
#	printf "key %f %f %f\n", $key, ${$self->{SPECTRUM}}{$key}, $max;
	if ( abs (${$self->{SPECTRUM}}{$key}-$max) < 0.000000001) { $xmax = $key; last; }
    }
    return ($xmax, $max);
}

sub y_at_x {
    my $self = shift;
    my $at_x = shift;
    my $y;
    foreach $key (sort {$a<=>$b}( keys %{$self->{SPECTRUM}})) {
#	printf "key %f %f %f\n", $key, ${$self->{SPECTRUM}}{$key}, $max;
	if ( abs($key-$at_x) > 0.00001 ) { $y = ${$self->{SPECTRUM}}{$key}; last; }
    }
    return $y;
}

sub y_absmax_at_x {
    my $self = shift;
    my @y = sort {$a<=>$b}( values %{$self->{SPECTRUM}});
    my $yy;
    my $max = 0.0;     
    foreach $yy (@y) {
	if (abs($yy) > $max ) { $max = abs($yy); }
    }
    my $xmax=0;
    my $i;
    foreach $key (sort {$a<=>$b}( keys %{$self->{SPECTRUM}})) {
#	printf "key %f %f %f\n", $key, ${$self->{SPECTRUM}}{$key}, $max;
	if ( abs (abs(${$self->{SPECTRUM}}{$key})-$max) < 0.000000001) { $xmax = $key; last; }
    }
    return ($xmax, $max);
}

sub y_min {
    my $self = shift;
    my @y = sort {$a<=>$b}( values %{$self->{SPECTRUM}});
    return $y[0];
}

sub add {
    my $self  = shift;
    my $value   = shift;
    my %h1 = $self->hash();
    foreach $key (sort {$a<=>$b}( keys %h1)) {
	$h1{$key} = $h1{$key} + $value;
    }
    $self->spectrum(%h1);
    return $self;
}

sub add_to_hash {
    my $self    = shift;
    my $key     = shift;
    my $value   = shift;
    my %h1 = $self->hash();
    $h1{$key} = $value;
    $self->spectrum(%h1);
    return $self;
}

sub avg {
    my @spectra = @_;
    my $spectrum;
    my %avg=0;
    my %avg_sq=0;
    my %std=0;
    my $wvn=0;
    foreach $spectrum (@spectra) {
#	printf "type %s\n", $spectrum->type();
	my %h = $spectrum->hash();
	foreach $key (sort {$a<=>$b}( keys %h)) {
	    $avg{$key} += $h{$key};
	    $avg_sq{$key} += $h{$key}*$h{$key};
#	    if (abs($key-320.0)<0.1) { 
#		print "$key $h{$key} $avg{$key} $avg_sq{$key}\n";
#	    }
	}	
    }
    $nspectra = $#spectra+1;
#    print "nspectra $nspectra\n";
    foreach $key (sort {$a<=>$b}( keys %avg)) {
	$avg{$key} = $avg{$key}/$nspectra;
	$avg_sq{$key} = $avg_sq{$key}/$nspectra;
    }	
    
    foreach $key (sort {$a<=>$b}( keys %avg)) {
	$std{$key} = sqrt($avg_sq{$key}-$avg{$key}*$avg{$key});
    }	

    my $S_avg = new Spectrum("average");
    $S_avg->spectrum(%avg);
    my $S_std = new Spectrum("std");
    $S_std->spectrum(%std);
    return ($S_avg,$S_std);
}      

sub fractional_difference {
#    my $self  = shift;
    my $self1 = shift;
    my $self2 = shift;
    my %h1 = $self1->hash();
    my %h2 = $self2->hash();
    my %h3;
    foreach $key (sort {$a<=>$b}( keys %h1)) {
	if (exists $h2{$key}) {
	    if ( $h2{$key} != 0.0 ) { $h3{$key} = ($h1{$key}-$h2{$key})/$h2{$key} }
	    else                  { $h3{$key} = 0.0 }
	}
    }
    my $S_fd = new Spectrum("fractional_difference");
    $S_fd->spectrum(%h3);
    $S_fd->albedo($self2->albedo());
    return $S_fd;
}      

sub linear_estimate_Y {
#    my $self  = shift;
    my @ary = @_;
    my $self1 = shift;
    my $self2 = shift;
    my $alb_file1;
    my $alb_file2;
    my %x1;
    my %x2;
    my %h1 = $self1->hash();
    my %h2 = $self2->hash();
    my %h3;
    my @keys = sort {$a<=>$b}( keys %h1);

#    print " ary: @ary\n";

    if ($#ary > 1) { 
	$alb_file1 = shift;
	$alb_file2 = shift;
	my $S = new Spectrum;
	%x1  = $S->read_spectrum_from_file($alb_file1);
	%x2  = $S->read_spectrum_from_file($alb_file2);

#	print "alb $alb_file1\n";
#	print %x1;
#	exit;
    }
    else {
	my $x_1 = $self1->albedo();
	my $x_2 = $self2->albedo();
	foreach $key (@keys) {
	    $x1{$key} = $x_1;
	    $x2{$key} = $x_2;
	}
    }
    foreach $key (@keys) {
	my ($c1,$c2);
 	if (exists $h2{$key}) {
 	    if ( $h2{$key} != 0.0 ) { 
#		$key, $h2{$key}, $h1{$key}, $x2{$key}, $x1{$key};
 		$c2 = ($h2{$key}-$h1{$key})/($x2{$key}-$x1{$key});
		$c1 = $h1{$key} - $c2*$x1{$key};
		$h3{$key} = - $c1/$c2;
#		$h1{$key}, $h2{$key}, $h3{$key}, $key, $x1{$key}, $x2{$key}, 
#		($h2{$key}-$h1{$key}), ($x2{$key}-$x1{$key}), $c1, $c2;
 	    }
 	    else { 
 		$h3{$key} = ($x1{$key}+$x2{$key})/2.;
 	    }
 	}
    }
    my $S_leY = new Spectrum("linear_estimate_Y");
    $S_leY->spectrum(%h3);
    return $S_leY;
}      

sub mult {
    my $self  = shift;
    my $value   = shift;
    my %h1 = $self->hash();
    foreach $key (sort {$a<=>$b}( keys %h1)) {
	$h1{$key} = $h1{$key} * $value;
    }
    $self->spectrum(%h1);
    return $self;
}

sub mult_cosine_sza {
    my $self  = shift;
    my $sza   = shift;
    my %h1 = $self->hash();
    my $pi = atan2(1,1) * 4;
    my $sza_rad = $sza*$pi/180.;
#    print "$sza $sza_rad\n"; exit;
    foreach $key (sort {$a<=>$b}( keys %h1)) {
#	printf "%f %f %f\n", $h1{$key}, cos($sza_rad), $h1{$key} * cos($sza_rad);
	$h1{$key} = $h1{$key} * cos($sza_rad);
    }
    $self->spectrum(%h1);
    return $self;
}

sub multiply3spectra {
    my $self1 = shift;
    my $self2 = shift;
    my $self3 = shift;
    my %h1 = $self1->hash();
    my %h2 = $self2->hash();
    my %h3 = $self3->hash();
    my %h4;
    foreach $key (sort {$a<=>$b}( keys %h1)) {
	$h4{$key} = $h1{$key}*$h2{$key}*$h3{$key}; 
    }
    my $S = new Spectrum("multiply3");
    $S->spectrum(%h4);
    return $S;
}      

sub ratio {
    my $self1 = shift;
    my $self2 = shift;
#    my $self3 = shift;
    my %h1 = $self1->hash();
    my %h2 = $self2->hash();
    my %h3;
    foreach $key (sort {$a<=>$b}( keys %h1)) {
#	print "ratio $key h1: $h1{$key} h2: $h2{$key} h3: $h3{$key} \n";
	if (exists $h2{$key}) {
	    if ( $h2{$key} != 0.0 ) { $h3{$key} = $h1{$key}/$h2{$key} }
	    else                  { $h3{$key} = 0.0 }
#	    print "ratio $key $h1{$key} $h2{$key} $h3{$key} \n";
	}
    }
    my $S_ratio = new Spectrum("ratio");
    $S_ratio->spectrum(%h3);
    my $albedo = $self2->albedo();
    $S_ratio->albedo($albedo);
    return $S_ratio;
}      

sub read_spectrum_from_file {
    my $self  = shift;
    my $file  = shift;
    my %h;
    open(TMP,"$file");
    while ($line=<TMP>) {
	chop $line;
	$line =~ s/^[ ]+//;
	my($key, $value) = split(/[ ]+/,$line);
	$key = sprintf("%7.3f",$key);
	$h{$key} = $value;
    }
    close(TMP);
    return %h; 
}

sub scale_to_range {
    my $self  = shift;
    my $x1   = shift;
    my $x2   = shift;
    my $wanted_range = $x2-$x1;
    my $max = $self->y_max();
    my $min = $self->y_min();
    my $original_range = $max-$min;
    my $fact = $wanted_range/$original_range;
    my %h1 = $self->hash();
    foreach $key (sort {$a<=>$b}( keys %h1)) {
	$h1{$key} = $h1{$key} - $min;
	$h1{$key} = $h1{$key} * $fact;
	$h1{$key} = $h1{$key} + $x1;
    }
    $self->spectrum(%h1);
    return $self;
}

sub smooth {
    my $self = shift;
    my $filter_file = shift;
    my $tmp_file1 = "./tmp1";
    my $tmp_file2 = "./tmp2";
    $self->write_spectrum_to_file($tmp_file1);
    system("$conv $tmp_file1 $filter_file > $tmp_file2");
    my %h  = $self->read_spectrum_from_file($tmp_file2);
    system("rm $tmp_file1 $tmp_file2");
    my $S_smooth = new Spectrum("smooth");
    $S_smooth->spectrum(%h);
    $S_smooth->albedo($self->albedo());
    return $S_smooth;
}

sub integrate {
    my $self = shift;
    my $filter_file = shift;
    my $tmp_file1 = "./tmp1";
    my $tmp_file2 = "./tmp2";
    $self->write_spectrum_to_file($tmp_file1);
    system("$integrate -p $tmp_file1 > $tmp_file2");
    open(FP,$tmp_file2);
    my $line = <FP>;
    close(FP);
    chop $line;
    my $val = sprintf("%e",$line);
    return $val;
}

sub spline_interpolate {
    my $self = shift;
    my $filter_file = shift;
    my $tmp_file1 = "./tmp1";
    my $tmp_file2 = "./tmp2";
    $self->write_spectrum_to_file($tmp_file1);
    system("$spline -x $filter_file  $tmp_file1 > $tmp_file2");
    my %h  = $self->read_spectrum_from_file($tmp_file2);
#    system("rm $tmp_file1 $tmp_file2");
    my $S_spline = new Spectrum("spline_interpolate");
    $S_spline->spectrum(%h);
    $S_spline->albedo($self->albedo());
    return $S_spline;
}

sub linear_interpolate {
    my $self = shift;
    my $filter_file = shift;
    my $tmp_file1 = "./tmp1";
    my $tmp_file2 = "./tmp2";
    $self->write_spectrum_to_file($tmp_file1);
    system("$spline -l -x $filter_file  $tmp_file1 > $tmp_file2");
    my %h  = $self->read_spectrum_from_file($tmp_file2);
#    system("rm $tmp_file1 $tmp_file2");
    my $S_spline = new Spectrum("linear_interpolate");
    $S_spline->spectrum(%h);
    $S_spline->albedo($self->albedo());
    return $S_spline;
}

sub sub {
    my $self  = shift;
    my $value   = shift;
    my %h1 = $self->hash();
    foreach $key (sort {$a<=>$b}( keys %h1)) {
	$h1{$key} = $h1{$key} - $value;
    }
    $self->spectrum(%h1);
    return $self;
}

sub uva {
    my $self = shift;
    my %h = $self->hash();
    my $uva = 0.0;
    foreach $key (sort {$a<=>$b}( keys %h)) {
	if ( $key >= 320.0 && $key <= 400.0 ) { $uva += $h{$key} }
    }	
#    print "uva $uva\n";
    $self->{UVA} = $uva;
    return $uva;
}      

sub uvb {
    my $self = shift;
    my %h = $self->hash();
    my $uvb = 0.0;
    foreach $key (sort {$a<=>$b}( keys %h)) {
	if ( $key >= 280.0 && $key <= 315.0 ) { $uvb += $h{$key} }
    }	
    $self->{UVB} = $uvb;
    return $uvb;
}      
	      
sub value {
    my $self = shift;
    my $key  = shift;
    if ( !exists $self->{SPECTRUM}{$key} ) {
	fprintf STDERR "Spectrum->value(%f) does not exist\n", $key;
	die;
    }
    my $value = ${ $self->{SPECTRUM} } {$key};
    return $value;
}

sub write_spectrum_to_file {
    my $self  = shift;
    my $file  = shift;
    open(TMP,">$file");
    foreach $key (sort {$a<=>$b}( keys %{$self->{SPECTRUM}})) {
	printf TMP  "$key, ${ $self->{SPECTRUM}}{$key}\n";
    }    
    close(TMP);
}

sub write_key_to_file {
    my $self  = shift;
    my $file  = shift;
    open(TMP,">$file");
    foreach $key (sort {$a<=>$b}( keys %{$self->{SPECTRUM}})) {
	printf TMP  "$key\n";
    }    
    close(TMP);
}

# Local Variables:
# mode: Perl
# End:
