package Column;
use strict;
use Carp;

1;

sub new {
    croak("Usage: new Column( line )")
      unless(@_ == 3);

    my $that  = shift;
    my $class = ref($that) || $that;
    my $symbol_name = shift;
    my $unit        = shift;
    my @data = ();
    my $self  = {
        DATA        => \@data,
	UNIT        => $unit,
	SYMBOL_NAME => $symbol_name,
    };
    bless $self, $class;
    return $self;
}

# sub Column_line_get {
#     my $line = shift;
#     $line =~ s/^[ ]+//;  # rm whitespace in beginning of line
#     my ($symbol_name_and_unit, $data ) = {};
    
#     # check if symbol name contains unit description or not
#     # and split accordingly
#     if ( grep /\)/, $line ) {
# 	($symbol_name_and_unit, $data ) = split(/\)/,$line,2);
#     }
#     else {
# 	($symbol_name_and_unit, $data ) = split(/[ ]+/,$line,2);
#     }
#     my ($symbol_name, $unit) = split(/[\(]/,$symbol_name_and_unit);

#     # data in $data ary may need to be regrouped
#     # Check number of fields in unit
#     my $unit_fields = split (/[ ]+/,$unit);
#     $data =~ s/^[ ]+//;  # rm whitespace in beginning of data
#     my $no_data_pnts;
#     my @data = ();
#     if ( $unit_fields > 0 ) {
# 	$no_data_pnts = split (/[ ]+/,$data);
# 	if ( $unit_fields == $no_data_pnts ) {
# 	    @data = ($data);
# 	    return ($symbol_name, $unit, @data);
# 	}
# 	else {
# 	    my $i;
# 	    my @data_vals = (); 
# 	    @data_vals = split (/[ ]+/,$data);
# 	    my $n_elem = $no_data_pnts/$unit_fields;
# 	    for ($i=0;$i<$n_elem;$i++) {
# 		my $j;
# 		for ($j=0;$j<$unit_fields;$j++) {
# 		    $data[$i] = $data[$i]." ".$data_vals[($i)*$unit_fields+$j];
# 		}	    
# 	    }
# 	    return ($symbol_name, $unit, @data);
# 	}
#     }
#     else {
# 	# In this case we have a text string
# 	$no_data_pnts = split (/ "/,$data);
# 	my $n_elem = $no_data_pnts;
# 	my @data_vals = split (/ "/,$data);
#         my $i;
# 	for ($i=0;$i<$n_elem;$i++) {
# 	    $data[$i] = $data_vals[$i];
# 	    $data[$i]  =~ s/"//g;
# 	}
# 	return ($symbol_name, $unit, @data);
#     }
# }

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
        my @data = @_;
        $self->{DATA} = \@data;
    }
    else {
        return @{$self->{DATA}};
    }
}
