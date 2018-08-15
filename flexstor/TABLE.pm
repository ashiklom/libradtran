package TABLE;
use strict;
use Carp;
use Column;
use SGL;

1;

sub require_version {
}

sub new {
    croak("Usage: new table( lines )")
      unless(@_ >= 2);

    my $that  = shift;
    my $class = ref($that) || $that;
    my @lines  = @_;
    my ($identifier, @columns) = TABLE_get(@lines);
    my $self  = {
        COLUMNS     => \@columns,
	IDENTIFIER  => $identifier,
    };
#    print "identifier $identifier\n";
    bless $self, $class;
    return $self;
}

sub TABLE_get {

    # One table is currently made up of many columns,
    # thus we put table stuff into package Column

    my @lines = @_;
    
    # Get table identifier from first line
    my $l = 0;
    my $line = $lines[$l];
    $line =~ s/^[ ]+//;
    $line =~ s/\n//;
    my ($gabba, $gabba, $identifier) = split(" ",$line);

    # Get data description line
    # Should be first line after "begin table $identifier"
    $l++;
    $line = $lines[$l];
    $line =~ s/^[ ]+//;
    $line =~ s/\n//;
    my(@data_description) = split(/\)/,$line);
    my $n_cols = $#data_description;

    # Get possible ancillary line
    $l++;
    $line = $lines[$l];
    my $anc_action="";
    my @anc_data = ();
    if (grep /^_/,$line) {
	$line =~ s/^[ ]+//;
	$line =~ s/\n//;
	($anc_action, @anc_data) = split(/[ ]+/,$line);
#	print "anc $anc_action \n";
    }
    else { $l--; }
    
    # Get data columns
    my @data = ();
    my $ll=0;
    for (++$l;$l<$#lines;$l++) {       
	$line = $lines[$l];
	$line =~ s/^[ ]+//;
	$line =~ s/\n//;
	my (@cols) = split(/[ ]+/,$line);
	my $ic;
	for ($ic=0;$ic<=$n_cols;$ic++) {	    
	    if ($anc_action eq "_div") {
		$data[$ic][$ll] = $cols[$ic]/$anc_data[$ic];
	    }
	    elsif ($anc_action eq "_mult") {
                $data[$ic][$ll] = $cols[$ic]*$anc_data[$ic];
            }
            elsif ($anc_action eq "_plus") {
                $data[$ic][$ll] = $cols[$ic]+$anc_data[$ic];
            }
            elsif ($anc_action eq "_minus") {
                $data[$ic][$ll] = $cols[$ic]-$anc_data[$ic];
            }
            elsif ($anc_action eq "_miss") {
                $data[$ic][$ll] = -99999.999;
            }
	    else {
		$data[$ic][$ll] = $cols[$ic];
	    }
	}
#	print "\n";
	$ll++;
    }    

#    print "data :", @{$data[0]}, "\n";

    # Finally fill all Columns
    my @columns=();
    my $ic;
    for ($ic=0;$ic<=$n_cols;$ic++) {
	# check if symbol name contains unit description or not
	# and split accordingly
	my ($symbol_name, $unit) = split(/[\(]/,$data_description[$ic]);
	$columns[$ic] = new Column($symbol_name, $unit);
	$columns[$ic]->data(@{$data[$ic]});
    }

    return ($identifier, @columns);
}

sub identifier {
    my $self = shift;
    if (@_) {
        my $identifier = shift;
        $self->{IDENTIFIER} = $identifier;
    }
    else {
        return $self->{IDENTIFIER};
    }
}

sub symbol_name {
    my $self = shift;
    if (@_) {
        my $symbol_name = shift;
        $self->{SYMBOL_NAME} = $symbol_name;
    }
    else {
        return $self->{SYMBOL_NAME};
    }
}

sub columns {
    my $self = shift;
    if (@_) {
        my $columns = shift;
        $self->{COLUMNS} = $columns;
    }
    else {
        return @{$self->{COLUMNS}};
    }
}
