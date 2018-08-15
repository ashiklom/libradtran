package Flexstor;
use strict;
use vars qw($VERSION);

use Carp;
use SGL;
use TABLE


$VERSION = '1.0';
BEGIN {
    $Flexstor::file = "Null";
}
1;

sub VERSION {
    # Version of Flexstor
    return $Flexstor::VERSION;
}

sub new {
    croak("Usage: new Flexstor( file )")
      unless(@_ == 2);

    my $that  = shift;
    my $class = ref($that) || $that;
    my $file  = shift;
    if (!-e $file) { 
	warn("File $file not found");
	return -1;
    }
    my @content = read_Flexstor_file($file);
    my $self  = {
        CONTENT => \@content,
	FILE    => $file,
    };
    bless $self, $class;
    return $self;
}

sub read_Flexstor_file {
    my $file = shift;
    my (@items) = split("/",$file);
    my $filepath="";
    my $i;
    for ($i=0;$i<$#items;$i++) {
	$filepath = $filepath.$items[$i]."/";
    }
    open(FP,$file);
    my @content = <FP>;
    close(FP);
    my @lines = grep /include/, @content;
    for ($i=0;$i<=$#lines;$i++) {
	my $line = $lines[$i];
	$line =~ s/^[ ]+//;
	my ($gabba,$file) = split(" ",$line);
	$file =~ s/\"//g;
	$file = $filepath.$file;
	open(FP,$file) || warn("Can not open include file $file");
	my @inc_content = <FP>;
	close(FP);
	@content = (@content,@inc_content);
    }
    return @content;
}

sub sgl
{
    my $self = shift;
    my $sgl_symbol = shift;
    my @lines = ();
    @lines = grep /$sgl_symbol/, @{$self->{CONTENT}};
    if ($#lines<0) {
	warn("Single line symbol $sgl_symbol not found");
    }
    my $i;
    my @sgl = ();
    for ($i=0;$i<=$#lines;$i++) {
	my $line = $lines[$i];
	$line =~ s/\n//;
	$sgl[$i] = new SGL($line);
    }
    return @sgl;
}

sub sgl_value {
    my $self = shift;
    my $sgl_symbol = shift;
    my @lines = ();
    @lines = grep /$sgl_symbol/, @{$self->{CONTENT}};
    if ($#lines<0) {
	warn("Single line symbol $sgl_symbol not found");
    }
    my $i;
    my $sgl;
    my $line = $lines[0];
    $line =~ s/\n//;
    $sgl = new SGL($line);
    return $sgl->value();
}

sub table
{
    my $self = shift;
    my $table_identifier = shift;
    my @content = @{$self->{CONTENT}};
    my $no_tables = grep /$table_identifier/, @content;
    if ($no_tables<0) {
	warn("No tables with identifier $table_identifier found");
    }
    my $i;
    my $j;
    my $i_table;
    my @tables = ();
    $no_tables = $no_tables/2; #There is both begin and end so divide
                               # by 2 to get number of tables.
    for ($i=0;$i<=$#content;$i++) {
	my $line = $content[$i];
	$line =~ s/\n//;
	if ( grep /table $table_identifier/, $line ) {
	    # Found begin table
	    my $i_begin = $i;
	    my $i_end = $i;
	    # Get end table
	    for ($j=$i+1;$j<=$#content;$j++) {
		my $line = $content[$j];
		if ( grep /table $table_identifier/, $line ) { 
		    $i_end = $j;
		    last;
		}
	    }
	    $tables[$i_table++] = new TABLE(@content[$i_begin..$i_end]);
	    $i = $i_end;   # Advance $i beyond this table
	} 
    }
    return @tables;
}

sub table_identifiers
{
    my $self = shift;
    my @content = @{$self->{CONTENT}};
    my $no_tables = grep /begin[ ]+table/, @content;
    if ($no_tables<0) {
	warn("No tables found");
    }
    my @table_identifiers = ();
    my $j = 0;
    for (my $i=0;$i<=$#content;$i++) {
 	my $line = $content[$i];
 	$line =~ s/\n//;
 	if ( grep /begin[ ]+table/, $line ) {
	    $line =~ s/begin//;
	    $line =~ s/table//;
	    $line =~ s/[ ]+//;
	    $table_identifiers[$j++] = $line;
	}
    }
    return @table_identifiers;
}

sub get_sgl {
    my $file      = shift;
    my $sgl_symbol= shift;
    my $f         = new Flexstor($file);
    my @sgl       = $f->sgl($sgl_symbol);
    my @sgl_value = $sgl[0]->data();
    return $sgl_value[0];
}


__END__

=head1 NAME

Flexstor - Query and/or read Flexstor files

=head1 SYNOPSIS

DOCUMENTATION NOT YET FINISHED...................

  	use Flexstor;
        # The rest is too complicated for a synopsis; keep reading.

=head1 DESCRIPTION

my $f = new Flexstor($file);  # Open and reads Flexstor 
                              # file $file.

# Single line information

my $sgl = $f->sgl($single_line_symbol); # $sgl is the 
                           # number of single lines of 
                           # type  $single_line_symbol

$sgl    = $f->sgl_value($single_line_symbol); # $sgl 
                           # contains the value of the 
                           # first line of type 
                           # $single_line_symbol

my @sgl = $f->sgl($single_line_symbol); # @sgl contains 
                           # all lines of type 
                           # $single_line_symbol. Appears 
                           # in the order which they are
                           # encountered in  the file. 

my @table_identifiers = $f->table_identifiers(); # @table_identifiers 
                           # hold all the table identifiers flexstor 
                           # file pointed to by $f.

my @table  = $f->table($table_identifier); # @table 
                           # contains all tables of 
                           # type $table_identifier

my @columns = $table[0]->columns(); @columns holds all
                           # data from the first table 
                           # in the file

my $symbol = @columns[0]->symbol_name(); # The symbol name 
                           #(from the 
                           # data_description_line)
                           # of the first column

my $unit   = @columns[0]->unit(); # The unit of the first 
                           # column (from the 
                           # data_description_line)

my @x   = $columns[0]->data(); # @x is an ary holding data 
                           # from the first column

Note that the possible action on data stated in an ancillary line is performed 
before data is read into table.

=head1 BUGS

Included files are not included in their original position when
the Flexstor file is read.

Only the '_div' action of the ancillary lines is accounted for. Other
options, e.g. '_mult', is not implemented.

=head1 AUTHOR

Flexstor was written by Arve Kylling (arve.kylling@nilu.no).

=cut

