#!/usr/bin/perl
# script to extract from the packaged up camera perturbation file
# a single instance of the surface distortions. THIS DOES NOT INCLUDE
# THE SOLID BODY TRANSFORMATIONS NOR THE ZERNIKE REPRESENTATIONS FOR THE
# DISTORTIONS, BUT ONLY ACCESSES THE FEA CALCULATION SUPERPOSITIONS DESCRIBED
# IN THE DOCUMENTATION NARRATIVES.
# USAGE: 
# extract_from_fits.perl <row number>
# 
# output: 6 files, named 
# tmp_{L1,L2,L3}_{S1,S2}_distortion.tnt
#
# which may be plotted using hippodraw or any other such ntuple plotter.
#
# arasmus@slac.stanford.edu, 140730.
use strict;
use warnings;
# script to repeatedly call CAMERA_perturbations.per
# and to package up the output into FITS and excel files
use POSIX;
use Astro::FITS::CFITSIO qw ( :longnames :constants );

my $entry=387; # starting value
my $status=0;
my $fptr=Astro::FITS::CFITSIO::open_file("output_fitstable.fits",
					 READONLY(),
					 $status);

my $data=[[],[],[],[],[],[]];

if ($#ARGV>=0) {
    $entry=$ARGV[0];
}
printf STDERR "entry=$entry\n";

foreach my $label ("L1_S1_distortion",		   "L1_S2_distortion",
		   "L2_S1_distortion",		   "L2_S2_distortion",
		   "L3_S1_distortion",		   "L3_S2_distortion") {

    $fptr->movnam_hdu(BINARY_TBL,$label,0,$status);
    check_error($status);
    # get number of columns
    my ($ncols,$colname);
    $fptr->get_num_cols($ncols,$status);
    printf STDERR "ncols $ncols\n";
    # get column names

    my @colnames=();
    my @coln=();
    printf STDERR "label=$label:\n";
    for (my $col=0;$col<$ncols;$col++) {
	$fptr->get_colname(CASEINSEN,"*",$colname,$col+1,$status);
	push(@colnames,$colname);
    }
    printf STDERR "colnames: %s\n",join(':',@colnames);
    $status=0;

    my @types=();
    for (my $col=0;$col<$ncols;$col++) {

	my ($type,$repeat,$width,$offset);
	$fptr->get_coltype($col+1,$type,$repeat,$width,$status);
	printf STDERR "col $col type $type repeat $repeat width $width\n";

	push(@types,$type);

	if ($type < 0) { # variable length column..
	    $fptr->read_descript($col+1,$entry,$repeat,$offset,$status);
	    printf STDERR "col $col repeat $repeat\n";
	}

	$fptr->read_col_flt($col+1,$entry,1,$repeat,[],$data->[$col],[],$status);
	check_error($status);
    }


    # string this together as a .tnt file and output..
    my $datlabel="";
    my @datcols=();
    for (my $col=0;$col<$ncols;$col++) {
	if ($types[$col]>0) {
	    # append the data label
	    $datlabel.=sprintf("%s=%f ",$colnames[$col],$data->[$col]->[0]);
	} else {
	    push(@coln,$colnames[$col]);
	    push(@datcols,$data->[$col]);
	}
    }
    open(GG,">tmp_".$label.".tnt") || die;
    printf GG "%s\n",$datlabel;
    printf GG "%s\n",join("\t",@coln);
    for (my $i=0;$i<=$#{$datcols[0]};$i++) {
	foreach my $d (@datcols) {
	    printf GG "%s ",$d->[$i];
	}
	printf GG "\n";
    }
    close(GG);
    for (my $col=0;$col<$ncols;$col++) {
	undef $data->[$col];
    }
}

sub check_error {
    my ($stat)=@_;
    if ($stat) {
        fits_report_error(*STDERR,$stat);
    }
}
