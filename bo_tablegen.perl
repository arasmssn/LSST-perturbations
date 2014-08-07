#!/usr/bin/perl
use warnings;
use strict;

# script to generate a table representation for each of the relevant 
# camera perturbation terms
# This is performed by repeatedly calling extract_from_fits.perl
# to reverse what was done in the packaging script:
# CAMERA_FEA/generate_misalignments_package.perl
# Bo needs the surface distortion data to be provided in grid format.
# this necessarily means that there must be one table per surface per 
# driving condition: 
# L{1,2,3}_S{1,2}_Z{0,90}_CR{0,90}

my $extract_rows = [[ 1..3 ],
		    [ 4..11 ]];

my @output_filesuffix=("gravity","thermalsoak");

foreach my $run ( 0..$#output_filesuffix ) {

    my %outputdata=("L1_S1_distortion" => [],"L1_S2_distortion" => [],
		    "L2_S1_distortion" => [],"L2_S2_distortion" => [],
		    "L3_S1_distortion" => [],"L3_S2_distortion" => []);

    foreach my $row ( @{ $extract_rows->[$run] } ) {

	`./extract_from_fits.perl $row`;

	# read in the 6 files 
	# tmp_L{1,2,3}_S{1,2}_distortion.tnt
	# ignore the solidbody files for now; Bo presumably gets these from 
	# the multi-extension binary FITS tables.

	foreach my $lens (1,2,3) {
	    foreach my $side (1,2) {
		my $key=sprintf("L%1d_S%1d_distortion",$lens,$side);
		my $file=sprintf("tmp_%s.tnt",$key);
		push($outputdata{$key},[file2list($file)]);
	    }
	}
    }
    foreach my $lens (1,2,3) {
	foreach my $side (1,2) {
	    my $key=sprintf("L%1d_S%1d_distortion",$lens,$side);
	    my $file=sprintf("%s_%s.multicol",$key,$output_filesuffix[$run]);
	    # open output file
	    open(F,">$file") || die;
	    # go through available lists, line by line to produce a 
	    # concatenated file like what Bo wants
	    foreach my $row ( 2 .. $#{$outputdata{$key}->[0]} ) {
		foreach my $col ( 0..$#{$outputdata{$key}} ) {
		    if ($col == 0 ) {
			printf F "%s ",$outputdata{$key}->[$col]->[$row];
		    } else {
	printf F "%s ",(split(' ',$outputdata{$key}->[$col]->[$row]))[2];
		    }
		}
		printf F "\n";
	    }
	    close(F);
	    push($outputdata{$key},[file2list($file)]);
	}
    }
}

sub file2list {
    my($filename)=@_;
    open(F,$filename) || die;
    my @list=<F>;
    close(F);
    chomp @list;
    @list;
}

