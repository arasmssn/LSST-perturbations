#!/usr/bin/perl
use strict;
use warnings;
# script to repeatedly call CAMERA_perturbations.per
# and to package up the output into FITS and excel files
use POSIX;
use Astro::FITS::CFITSIO qw ( :longnames :constants );

my $status=0;
my $fptr=Astro::FITS::CFITSIO::open_file("output_fitstable.fits",
					 READONLY(),
					 $status);
$fptr->movnam_hdu(BINARY_TBL,"L1_S1_distortion",0,$status);
check_error($status);
my $arr=[];
# $fptr->read_col_flt(1,387,1,1,0,$arr,$status);
check_error($status);
exit;

sub check_error {
    my ($stat)=@_;
    if ($stat) {
        fits_report_error(*STDERR,$stat);
    }
}
