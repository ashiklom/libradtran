package SGL;
use strict;
use Carp;

1;

sub new {
    croak("Usage: new SGL( line )")
      unless(@_ == 2);

    my $that  = shift;
    my $class = ref($that) || $that;
    my $line  = shift;
    my ($symbol_name, $unit, @data ) = SGL_line_get($line);
    my $self  = {
        DATA        => \@data,
	UNIT        => $unit,
	SYMBOL_NAME => $symbol_name,
    };
    bless $self, $class;
    return $self;
}

sub SGL_line_get {
    my $line = shift;
    $line =~ s/^[ ]+//;  # rm whitespace in beginning of line
    my ($symbol_name_and_unit, $data ) = {};
    
    # check if symbol name contains unit description or not
    # and split accordingly
    if ( grep /\)/, $line ) {
	($symbol_name_and_unit, $data ) = split(/\)/,$line,2);
    }
    else {
	($symbol_name_and_unit, $data ) = split(/[ ]+/,$line,2);
    }
    my ($symbol_name, $unit) = split(/[\(]/,$symbol_name_and_unit);

    # data in $data ary may need to be regrouped
    # Check number of fields in unit
    my (@fields) = split (/\s+/,$unit);
    my $unit_fields = split (/\s+/,$unit); 
    $data =~ s/^[ ]+//;  # rm whitespace in beginning of data
    my $no_data_pnts;
    my @data = ();
    if ( $unit_fields > 0 ) {
	$no_data_pnts = split (/\s+/,$data);
	if ( $unit_fields == $no_data_pnts ) {
	    @data = ($data);
	    return ($symbol_name, $unit, @data);
	}
	else {
	    my $i;
	    my @data_vals = (); 
	    @data_vals = split (/[ ]+/,$data);
	    my $n_elem = $no_data_pnts/$unit_fields;
	    for ($i=0;$i<$n_elem;$i++) {
		my $j;
		for ($j=0;$j<$unit_fields;$j++) {
		    $data[$i] = $data[$i]." ".$data_vals[($i)*$unit_fields+$j];
		}	    
	    }
	    return ($symbol_name, $unit, @data);
	}
    }
    else {
	# In this case we have a text string
	$no_data_pnts = split (/ "/,$data);
	my $n_elem = $no_data_pnts;
	my @data_vals = split (/ "/,$data);
        my $i;
	for ($i=0;$i<$n_elem;$i++) {
	    $data[$i] = $data_vals[$i];
	    $data[$i]  =~ s/"//g;
	}
	return ($symbol_name, $unit, @data);
    }
}

sub unit {
    my $self = shift;
    if (@_) {
        my $unit = shift;
        $self->{UNIT} = $unit;
    }
    else {
        return $self->{UNIT};
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

sub data {
    my $self = shift;
    if (@_) {
        my $data = shift;
        $self->{DATA} = $data;
    }
    else {
        return @{$self->{DATA}};
    }

sub value {
    my $self = shift;
    if (@_) {
        my $data = shift;
        $self->{DATA} = $data;
    }
    else {
        return @{$self->{DATA}}[0];
    }
}

}
