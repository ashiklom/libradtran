#
# Copyright (C) 1998 Arve Kylling
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
# the file snowalbedo.c. This so because sdoc do not read perl yet.
#####

use UVTools;
use Getopt::Long;

# This script generates three tables of the ice qext, ssa and gg 
# as a function of ice sphere radius and wavelength. The tables are
# used by snowalbedo to calculate the albedo of snow.

$smoothfile = "MIE_SMOOTH";

# Ice sphere radii (um)
@radii = (10., 25., 50., 100., 200., 500., 800.); 
for($i=0;$i<=4;$i++) {
    $prog{sprintf("%6.1f",$radii[$i])} = "MIEV0";
}
for($i=5;$i<=$#radii;$i++) {
    $prog{sprintf("%6.1f",$radii[$i])} = "MIEV0";
}
# Default command line options 
$file         = "Stamnes_table";

$ret = GetOptions(
		  "--file=s",\$file,
		  "--help",\&usage,
		  "<>",\&usage);
if ($ret==0) {
    exit;
}

# Now do the real things.

$tmptable    = "tmp_table";

# Loop over all radii
#for ($i=0;$i<=0;$i++) {
for ($i=0;$i<=$#radii;$i++) {

    printf STDERR "Now doing radius: %5.2f\n", $radii[$i];
    
    # Mie calculation for radius
    $radius = sprintf("%6.1f",$radii[$i]);
    $out = "tmp_out";
    $program = $prog{$radius};
#    printf "program $program\n";
    mie_calc($radius, $program, $out);
    # Smooth the Mie results
#    smooth($smoothfile, $out);
    
    # Put in table
    open(TMP,$out);
    @tmp = <TMP>;
    close(TMP);
    for ($k=0;$k<=$#tmp;$k++) {
	$line = @tmp[$k];
	$line =~ s/^[ ]*//;
	local($wlr,$re,$im,$qext,$ssa,$gg, $spike) = split(/[ \t]+/,$line);
	$_ = $qext;
	if ( m/[0-9\.\-\e]+/ ) {
	    $table_qext[$k][$i] = $qext;
	    $table_ssa[$k][$i]  = $ssa;
	    $table_gg[$k][$i]   = $gg;
	    if ($i==0) {
		$wl[$k]             = $wlr;
	    }
	}
	else {
	    $table_qext[$k][$i] = -99;
	    $table_ssa[$k][$i]  = -99;
	    $table_gg[$k][$i]   = -99;
	}
    }
}

# Write the three qext, ssa, and gg tables.  

open(TABLE_QEXT,">".$file.".qext");
printf TABLE_QEXT "%6.1f ", -99.0;
open(TABLE_SSA,">".$file.".ssa");
printf TABLE_SSA "%6.1f ", -99.0;
open(TABLE_GG,">".$file.".gg");
printf TABLE_GG "%6.1f ", -99.0;
for ($i=0; $i<=$#radii; $i++) {
    printf TABLE_QEXT "%12.1f", $radii[$i];
    printf TABLE_SSA "%16.1f", $radii[$i];
    printf TABLE_GG "%12.1f", $radii[$i];
}
printf TABLE_QEXT "\n";
printf TABLE_SSA "\n";
printf TABLE_GG "\n";
  
for ($k=0;$k<=$#wl;$k++) {
    printf TABLE_QEXT "%6.1f ", $wl[$k];
    printf TABLE_SSA "%10.1f ",  $wl[$k];
    printf TABLE_GG "%6.1f ",   $wl[$k];
    for ($i=0; $i<=$#radii; $i++) {
	printf TABLE_QEXT "%12.5e",  $table_qext[$k][$i];
	printf TABLE_SSA "%16.9e",   $table_ssa[$k][$i];
#	printf "%16.9e",   $table_ssa[$k][$i];
	printf TABLE_GG "%12.5e",    $table_gg[$k][$i];
    }
    printf TABLE_QEXT "\n";
    printf TABLE_SSA "\n";
    printf TABLE_GG "\n";
}
close(TABLE_QEXT);
close(TABLE_SSA);
close(TABLE_GG);

sub smooth 
{
    my $smooth_file = shift @_;
    my $out         = shift @_;
    my $tmp = "tmp";
    system("./conv $out $smooth_file > $tmp ");
#    system("diff $tmp $out; mv $tmp $out");
    system("mv $tmp $out");
}

sub mie_calc
{
    my $radius  = shift @_;
    my $program = shift @_;
    my $out     = shift @_;
    my $tmp = "tmp";
    open(TMP,">".$tmp);
    printf TMP "mie_program $program\n";
    printf TMP "mimcut 1.e-09\n";
    printf TMP "refrac ice\n";
    printf TMP "r_mean $radius\n";
    printf TMP "wvn 280.0 700.0";
    printf TMP "wvn_step 5.\n";
    close(TMP);
    system("mie < $tmp > $out");
    system("rm $tmp");
#     open(TMP,$out);
#     @tmp = <TMP>;
#     close(TMP);
#     system("rm $out");
# #    system("cp $out tmpgabba");
#     open(TMP,">".$out);
#     for ($k=0;$k<=$#tmp;$k++) {
# 	$line = @tmp[$k];
# 	chop $line;
# 	$line =~ s/^[ ]*//;
# 	local($wlr,$re,$im,$qext,$ssa,$gg, $spike) = split(/[ \t]+/,$line);
# 	if ($spike > 0.9999 || $program ne "MIEV0" ) {
# 	    printf TMP "$line\n";
# 	}
#     }
#     close(TMP);    
#    exit;
}

sub usage {
    printf STDERR "\n";
    printf STDERR "Gen_snow_tab generates three tables of the \n";
    printf STDERR "ice qext, ssa and gg as a function of ice \n";
    printf STDERR "sphere radius and wavelength. The tables are\n";
    printf STDERR " used by snowalbedo to calculate the albedo of snow.\n";
    printf STDERR "\n";
    printf STDERR "Gen_snow_tab understands the following options:\n";
    printf STDERR "\n";
    printf STDERR "--help                   : Prints this message.\n";
    printf STDERR "--file <name>            : Name of file where the table will be stored.\n";
    printf STDERR "\n";
    exit;
}

# Local Variables:
# mode: Perl
# End:
