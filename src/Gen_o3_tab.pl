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

# The Stamnes table is a function of solar zenith angles and Stamnes
# ratios. The solar zenith angles and Stamnes ratios in the table are
# set below. The numbers in the table is ozone values for specific
# combinations of solar zenith angles and Stamnes ratios.

# Solar zenith angles that will be used to calculate the Stamnes table
@szas = (0.1, 5,10,15,
	 20,22,24,26,28,
	 30,32,34,36,38,
	 40,42,44,46,48,
	 50,51,52,53,54,55,56,57,58,59,
	 60.01,61,62,63,64,65,66,67,68,69,
	 70,71,72,73,74,75,76,77,78,79,80); 

# The Stamnes ratio 
@ratio = (1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,
	  2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,
	  3,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,
	  4,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,
	  5,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,
	  6,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,
	  7,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,
	  8,8.1,8.2,8.3,8.4,8.5,8.6,8.7,8.8,8.9,
	  9,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,
	  10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,
	  15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,
	  20,20.5,21,21.5,22,22.5,23,23.5,24,24.5,
	  25,25.5,26,26.5,27,27.5,28,28.5,29,29.5,
	  30,31,32,33,34,35,36,37,38,39,40,41,42,
	  43,44,45,46,47,48,49,50,55,60,65,70,75,80,
	  85,90,95,100,105,110,115,120,125,130,
	  135,140,145,150,155,160,165,170,175,
	  180,185,190,195,200,210,220,230,240,
	  250,260,270,280,290,300,310,320,330,
	  340,350,360,370,380,390,400,410,420,
	  430,440,450,460,470,480,490,500,550,600);

$ozone_lower_limit = 0.0;
$ozone_upper_limit = 800.0;
$ozone_step        = 50.0;
$n_ozone           = ($ozone_upper_limit - $ozone_lower_limit) / $ozone_step;

# set ozone 
for ($j=0; $j<=$n_ozone; $j++)  {
    $ozones[$j] = $ozone_lower_limit + $j * $ozone_step;
}

$day          = 1.0; # Day number does not matter since we are dealing with ratios.
$UVfile_lower = "tmp_convlower";
$UVfile_upper = "tmp_convupper";

# Default command line options 
$altitude      = 0.0;
$absolute      = 0;
$albedo        = 0.0;
$alpha         = 0.0;
$atmmod_file   = "../data/atmmod/afglus.dat";
$bandpasslower = "";
$bandpassupper = "";
$beta          = 0.0;
$file          = "Stamnes_table";
$slitfunction  = '';
$lower_lambda  = 305.0;
$upper_lambda  = 340.0;
$o3_crs        = "Molina";
$wctau         = 0.0;
$wc_file       = "";
$zenith        = 0;

my $input_line = "$0 @ARGV";

$ret = GetOptions("--absolute",\$absolute,
		  "--albedo=f",\$albedo,
		  "--alpha=f",\$alpha,
		  "--altitude=f",\$altitude,
		  "--atmmod=s",\$atmmod_file,
		  "--bandpasslower=s",\$bandpasslower,
		  "--bandpassupper=s",\$bandpassupper,
		  "--beta=f",\$beta,
		  "--file=s",\$file,
		  "--help",\&usage,
		  "--slitfunction=s",\$slitfunction,
		  "--o3_crs=s",\$o3_crs,
		  "--lower_lambda=f",\$lower_lambda,
		  "--upper_lambda=f",\$upper_lambda,
		  "--zenith",\$zenith,
		  "<>",\&usage);
if ($ret==0)  {
    exit;
}


# A little checking of input data will not hurt

if ($lower_lambda >= $upper_lambda) {
    die "lower_lambda >= upper_lambda ($lower_lambda >= $upper_lambda)\n";
}

$UVTools::uvspec_o3_crs = $o3_crs;

# zenith radiance
if ($zenith>0)  {
    printf STDERR "Calculating zenith radiance!\n";
    $umu = -1.0;
    $phi =  0.0;
}
else  {
    $umu = 0.0;
    $phi = 0.0;
}

printf STDERR "Stamnes table is written to the file $file\n";
printf STDERR "Upper wavelength: %f\n", $upper_lambda;
printf STDERR "Lower wavelength: %f\n", $lower_lambda;

# Now do the real things.

$bandpass = 0;
if (length($bandpasslower)> 0 && length($bandpassupper) > 0) {
    $bandpass = 1;

    ($wvn1bandupper,$wvn2bandupper,$wvnstepupper) = read_bandpass($bandpassupper);
    ($wvn1bandlower,$wvn2bandlower,$wvnsteplower) = read_bandpass($bandpasslower);
    
    $wvn1upper = $wvn1bandupper;
    $wvn2upper = $wvn2bandupper;

    if (!$absolute) {
	$wvn1upper = $wvn1upper + $upper_lambda;
	$wvn2upper = $wvn2upper + $upper_lambda;
    } 

    $wvncenterupper = $wvn1bandupper;
    $stepupper      = $wvnstepupper;
    $wvn1lower      = $wvn1bandlower;
    $wvn2lower      = $wvn2bandlower;

    if (!$absolute) {
	$wvn1lower = $wvn1lower + $lower_lambda;
	$wvn2lower = $wvn2lower + $lower_lambda;
    } 

    $wvncenterlower = $wvn1bandlower;
    $steplower  = $wvnsteplower;
}
else {
    my ($wvn_left, $wvn_right) = (0,0);
    if ( length $slitfunction != '' ) {
      ($wvn_left, $wvn_right) = read_slitfunction($slitfunction);
      if ( abs($wvn_left ) < 1 ) { $wvn_left  = -1.0;}
      if ( abs($wvn_right) < 1 ) { $wvn_right = 1.0;}
    }
    $wvn1upper      = $upper_lambda + $wvn_left;  # Plus sign because $wvn_left is negative
    $wvn2upper      = $upper_lambda + $wvn_right;
    $wvncenterupper = ($wvn2upper+$wvn1upper)/2.0;
    $stepupper      = 15.0;  # Some large number to make spline only produce one number
    $wvn1lower      = $lower_lambda + $wvn_left;  # Plus sign because $wvn_left is negative
    $wvn2lower      = $lower_lambda + $wvn_right;
    $wvncenterlower = ($wvn2lower+$wvn1lower)/2.0;
    $steplower      = 15.0;  # Some large number to make spline only produce one number
}

$rat_file    = "tmp_rat";
$rat_o3_file = "tmp_rato3";
$tmptable    = "tmp_table";

open(RAT,">".$rat_file);
for ($i=0; $i<=$#ratio; $i++) {
    printf RAT "%f\n", $ratio[$i];
}
close(RAT);

# Loop over all zenith angles
for ($i=0;$i<=$#szas;$i++) {
    my $prev_rat = 0;
    open(RAT_O3,">".$rat_o3_file);

    # Loop over all ozone value
    for ($j=0;$j<=$n_ozone;$j++) {
	
	printf STDERR "Now doing sza: %5.2f and o3: %6.3f\n", $szas[$i], $ozones[$j];
	
	# Calculate upper wavelength
      UVTools::Convolved_Spectrum ($albedo, $day, $ozones[$j], 
				   $alpha, $altitude, $beta, $szas[$i], 
				   $wvn1upper, $wvn2upper, $wvncenterupper, 
				   $stepupper, $atmmod_file, 
				   $slitfunction, $UVfile_upper, $wctau, $wc_file,
				   $umu, $phi);
#	die;
	if ($bandpass) {
	    $irr_upper = multiply_with_bandpass_and_integrate ($UVfile_upper, 
							       $bandpassupper, 
							       $zenith);
	}
	
	# Calculate lower wavelength
      UVTools::Convolved_Spectrum($albedo, $day, $ozones[$j], 
				  $alpha, $altitude, $beta, $szas[$i], 
				  $wvn1lower, $wvn2lower, $wvncenterlower, 
				  $steplower, $atmmod_file, 
				  $slitfunction, $UVfile_lower, $wctau, $wc_file,
				  $umu, $phi);


	if ($bandpass) {
	    $irr_lower = multiply_with_bandpass_and_integrate ($UVfile_lower, 
							       $bandpasslower, 
							       $zenith);
	}

	# Calculate ratio
	if (!$bandpass) { 
	    if (!$zenith) {
		$irr_upper    = read_irr($UVfile_upper);
		$irr_lower    = read_irr($UVfile_lower);
	    }
	    else {
		$irr_upper    = read_rad($UVfile_upper);
		$irr_lower    = read_rad($UVfile_lower);
	    }
	}

	$ratio_up_low = $irr_upper/$irr_lower;

	# Do this $prev_rat to avoid problems with posible negative ratios from negative 
	# intensities, sigh...
	if ( $ratio_up_low < $prev_rat) {
	    $ratio_up_low = $prev_rat+1;
	}

	system("rm $UVfile_upper $UVfile_lower");
	printf RAT_O3 "%f %f\n", $ratio_up_low, $ozones[$j];
	$prev_rat = $ratio_up_low;
	
	printf STDERR "This ratio_up_low: %5.2f and prev_rat: %6.3f\n", $ratio_up_low, $prev_rat;
    }
    close(RAT_O3);

    # Interpolate calculated ratios to ratio
    system("$UVTools::spline -x $rat_file $rat_o3_file > $tmptable");
    open(TMPTABLE,$tmptable);
    @tmp = <TMPTABLE>;
    close(TMPTABLE);

    for ($k=0; $k<=$#tmp; $k++) {
	$line = @tmp[$k];
	$line =~ s/^[ ]*//;
	local($rat,$o3) = split(/[ \t]+/,$line);
	$_ = $o3;

	$table[$k][$i] = -1.0;
	if ( m/[0-9\.\-\e]+/ ) {
	  if ( $o3 >= 0 ) {
	    $table[$k][$i] = $o3;
	  }
	}
    }
}

# Write the table.  
open(TABLE,">".$file);
printf TABLE "\# Generated by the following command:\n";
printf TABLE "\# $input_line\n";
printf TABLE "%6.1f ", -1.0;
for ($i=0; $i<=$#szas; $i++) {
    printf TABLE "%5.0f", $szas[$i];
}
printf TABLE "\n";
  
for ($k=0;$k<=$#ratio;$k++) {
    printf TABLE "%6.1f ", $ratio[$k];
    for ($i=0; $i<=$#szas; $i++) {
	printf TABLE "%5d",  $table[$k][$i];
    }
    printf TABLE "\n";
}
close(TABLE);

# read total irradiance from UVSPEC output file
sub read_irr 
{
    local($file) = shift @_;
    open(FILE,$file);
    local($line) = <FILE>;
    close(FILE);
    $line =~ s/^[ ]*//;
    local($wvn,$dir,$diff,$rest) = split(/[ \t]+/,$line);
    local($totdir) = $dir+$diff;
    $totdir;
}

# read zenith radiance from UVSPEC output file
sub read_rad
{
    local($file) = shift @_;
    local($line) = "";
    open(FILE,$file);
    $line = <FILE>;
    $line = <FILE>;
    $line = <FILE>;  # we want the third line
    close(FILE);
    $line =~ s/^[ ]*//;
    local($wvn,$rest,$rad) = split(/[ \t]+/,$line);
    $rad;
}


sub read_bandpass
{
    my $file = shift @_;
    open(FILE,$file) || die "Could not open file: $file\n";
    my @lines = <FILE>;
    close(FILE);
    @lines = grep (!/\#/, @lines);  # Rm comment lines
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
    local($irrfile) = shift @_;
    local($bpfile)  = shift @_;
    local($zenith)  = shift @_;

    local($i, $j, $nlines);

    open(FILE,$irrfile) || die "Could not open file: $irrfile\n";
    local(@irrlines) = <FILE>;
    close(FILE);

    $nlines = $#irrlines;
    if ($zenith)  {
	$nlines = (($nlines + 1) / 3 - 1);
    }

    open(FILE,$bpfile)  || die "Could not open file: $bpfile\n";
    local(@bplines) = <FILE>;
    @bplines = grep (!/\#/, @bplines);  # Rm comment lines
    close(FILE);


    # Multiply irradiance with bandpass
    if ($nlines != $#bplines) {
	print "nlines = ", $nlines, "bplines = ", $#bplines;
	die "multiply_with_bandpass_and_integrate: Bandpass and uvfile have different length\n";
    }

    for ($i=0; $i<=$nlines; $i++) {
	# read radiation data
	if (!$zenith)  {
	    # irradiance
	    $line = @irrlines[$i];
	    chop $line;
	    $line =~ s/^[ \t]+//;
	    local($wvn, $dir, $diff) = split(/[ \t]+/,$line);
	    $irr[$i] = $dir + $diff;
	    $wvn[$i] = $wvn;
	}
	else {
	    # radiance, need to interpret 1st line (wavelength) 
	    # and 3rd line (radiance)
	    $line = @irrlines[3*$i];
	    chop $line;
	    $line =~ s/^[ \t]+//;
	    local($wvn) = split(/[ \t]+/,$line);
	    $wvn[$i] = $wvn;
	    
	    $line = @irrlines[3*$i+2];
	    chop $line;
	    $line =~ s/^[ \t]+//;
	    local($temp1, $temp2, $rad) = split(/[ \t]+/,$line);
	    $irr[$i] = $rad;
	}

	# read bandpass data
	$line = $bplines[$i];
	chop $line;
	$line =~ s/^[ \t]+//;
	local($wvnbp, $bp, $gabba) = split(/[ \t]+/,$line);

#	print "wvn  $wvn wvnbp: $wvnbp bp: $bp dir: $dir\n";
#	if ($wvn != $wvnbp)  { 
#	    die "multiply_with_bandpass_and_integrate: Bandpass and uvfile wvn differ $wvn $wvnbp\n";
#	}

	$irr[$i] = $irr[$i]*$bp;
    }

    # Integrate the bandpass-weighted irradiance
    
    $tmpbp = "tmpbp";
    $tmpirr = "tmpirr";
    open(BP,">".$tmpbp);
    for ($i=0;$i<=$nlines;$i++) {
	printf BP "%f %f\n", $wvn[$i], $irr[$i];
    }
    close(BP);
    
    system("$UVTools::integrate -p  $tmpbp > $tmpirr");
    
    open(IRR,$tmpirr);
    $irr = <IRR>;
    chop $irr;
    close(IRR);

#    system("$UVTools::conv $irrfile $bpfile > tmpconv");

#    die "multi $irr die\n";
    $irr;

}

sub usage {
    printf STDERR "\n";
    printf STDERR "Gen_o3_tab generates a table of some wavelength ratio versus\n";
    printf STDERR "solar zenith angles for different ozone amounts. The table is\n";
    printf STDERR "read and ozone columns derived by the read_o3_tab program.\n";
    printf STDERR "\n";
    printf STDERR "The following wavelength ratios are supported:\n";
    printf STDERR "(See the libRadtran documentation for examples.)\n";
    printf STDERR "\n";
    printf STDERR "   * Ratios between two wavelengths where each is \n";
    printf STDERR "     calculated from single wavelength measurements.\n";
    printf STDERR "   * Ratios between two wavelengths where each is \n";
    printf STDERR "     calculated from wavelength measurements multiplied by\n";
    printf STDERR "     a bandpass function and integrated over the bandpass.\n";
    printf STDERR "\n";
    printf STDERR "Gen_o3_tab understands the following options:\n";
    printf STDERR "\n";
    printf STDERR "--absolute               : The wavelengths of the bandpass files are in\n";
    printf STDERR "                           absolute units. Default is relative units.\n";
    printf STDERR "--albedo <value>         : Lambertian surface albedo. Default is 0.0.\n";
    printf STDERR "--alpha <value>          : Angstrom alpha coefficient. Default is 0.0.\n";
    printf STDERR "--altitude <value>       : Altitude above sea level [km]. Default is 0.0.\n";
    printf STDERR "--atmmod <name>          : Name of atmosphere file. Default atmmod/afglus.dat.\n";
    printf STDERR "--beta <value>           : Angstrom beta coefficient. Default is 0.0.\n";
    printf STDERR "--help                   : Prints this message.\n";
    printf STDERR "--o3_crs <name>          : Name of o3 cross section to use. Default is Molina.\n";
    printf STDERR "                           See libRadtran documentation for other options.\n";
    printf STDERR "--slitfunction <name>    : Name of slitfunction file.\n";
    printf STDERR "--bandpasslower <name>   : Name of file holding bandpass for lower wavelength.\n";
    printf STDERR "--bandpassupper <name>   : Name of file holding bandpass for upper wavelength.\n";
    printf STDERR "--file <name>            : Name of file where the table will be stored.\n";
    printf STDERR "--lower_lambda <value>   : Value for lower wavelength, in nm.\n";
    printf STDERR "--upper_lambda <value>   : Value for upper wavelength, in nm.\n";
    printf STDERR "--zenith                 : Calculate zenith sky radiance table.\n";
    printf STDERR "\n";
    exit;
}

# Local Variables:
# mode: Perl
# End:
