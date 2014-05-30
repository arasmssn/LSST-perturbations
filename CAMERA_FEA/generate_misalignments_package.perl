#!/usr/bin/perl
use strict;
use warnings;
# script to repeatedly call CAMERA_perturbations.per
# and to package up the output into FITS and excel files
use Astro::FITS::CFITSIO qw ( :longnames :constants );

my $status=0;
# $maxzern may be 28 or 45
my $maxzern=45;
my $compensate_camera=0;

my $fptr=Astro::FITS::CFITSIO::create_file("!output_fitstable.fits",$status);
$fptr->create_img(BYTE_IMG,2,[0,0],$status);
check_error($status);
fingerprint($fptr);

#
my (@names,@inputs,@input_units);
my ($tfields,$ttype,$tform,$tunit);
my %solidbodies;
my %distbodies;
my %zernbodies;
my %compbodies;

@inputs=("TCS_elevation","CAM_rotation","Soak_Temp");
@input_units=(("degree")x2,("degrees C"));


foreach my $object ("L1","L2","F","L3","FP") {
    my $label=$object."_solidbody";
    my $file=$object."solidbody_g.dat";
    @names=(@inputs,
	    $object."_tx",$object."_ty",$object."_tz",
	    $object."_rx",$object."_ry",$object."_rz");
    $tfields=$#names-$[+1;
    $ttype=[@names];
    $tform=[qw(1E)x$tfields];
    $tunit=[@input_units,("mm")x3,("radian")x3];

    $solidbodies{$label}=$file;
    $fptr->create_tbl(BINARY_TBL,0,$tfields,$ttype,$tform,$tunit,$label,$status);
}

# node-by-node representation of distortions

foreach my $object ("L1","L2","L3") {
    my $label=$object."_S1_distortion";
    my $file=$object."_S1_perturbation_g.tnt";
    @names=(@inputs);
    push(@names,$object."_S1_X");
    push(@names,$object."_S1_Y");
    push(@names,$object."_S1_Z_distortion");
    $tfields=$#names-$[+1;
    $ttype=[@names];
    $tform=[qw(1E)x($#inputs-$[+1),qw(1PE)x3];
    $tunit=[@input_units,("mm")x$tfields];
    $distbodies{$label}=$file;
    $fptr->create_tbl(BINARY_TBL,0,$tfields,$ttype,$tform,$tunit,$label,$status);
}

foreach my $object ("L1","L2","L3") {
    my $label=$object."_S2_distortion";
    my $file=$object."_S2_perturbation_g.tnt";
    @names=(@inputs);
    push(@names,$object."_S2_X");
    push(@names,$object."_S2_Y");
    push(@names,$object."_S2_Z_distortion");
    $tfields=$#names-$[+1;
    $ttype=[@names];
    $tform=[qw(1PE)x$tfields];
    $tunit=[@input_units,("mm")x$tfields];
    $distbodies{$label}=$file;
    $fptr->create_tbl(BINARY_TBL,0,$tfields,$ttype,$tform,$tunit,$label,$status);
}

# zernike expansion representation of distortions

foreach my $object ("L1","L2","L3") {
    my $label=$object."_S1_zernexp";
    my $file=$object."_S1_zern_g.dat";
    @names=(@inputs);
    for (my $i=0;$i<$maxzern;$i++) {
	push(@names,$object."_S1_Z".$i);
    }
    $tfields=$#names-$[+1;
    $ttype=[@names];
    $tform=[qw(1E)x$tfields];
    $tunit=[@input_units,("mm") x $maxzern];
    $zernbodies{$label}=$file;
    $fptr->create_tbl(BINARY_TBL,0,$tfields,$ttype,$tform,$tunit,$label,$status);
}


foreach my $object ("L1","L2","L3") {
    my $label=$object."_S2_zernexp";
    my $file=$object."_S2_zern_g.dat";
    @names=(@inputs);
    for (my $i=0;$i<$maxzern;$i++) {
	push(@names,$object."_S2_Z".$i);
    }
    $tfields=$#names-$[+1;
    $ttype=[@names];
    $tform=[qw(1E)x$tfields];
    $tunit=[@input_units,("mm") x $maxzern];
    $zernbodies{$label}=$file;
    $fptr->create_tbl(BINARY_TBL,0,$tfields,$ttype,$tform,$tunit,$label,$status);
}

if ($compensate_camera==1) {
    @names=(@inputs);
    push(@names,"CAM_comp_tx","CAM_comp_ty","CAM_comp_tz","CAM_comp_rx","CAM_comp_ry","CAM_comp_rz","Best_FOM");

    $tfields=$#names-$[+1;
    $ttype=[@names];
    $tform=[qw(1E)x$tfields];
    $tunit=[@input_units,("mm")x3,("radian")x3,"arcsec FWHM"];

    my $label="CAMcomp";
    my $file="CAMmisalign.dat";
    $compbodies{$label}=$file;
    $fptr->create_tbl(BINARY_TBL,0,$tfields,$ttype,$tform,$tunit,$label,$status);
}

# setup for table definitions finished. move on to input settings.

my $settings=[];

push(@{$settings},[ 0, 0,-99]);
push(@{$settings},[90, 0,-99]);
push(@{$settings},[90,90,-99]);

my @soaktemps=(-10,-5,0,5,10,15,20,25);
foreach my $soaktemp (@soaktemps) {
    push(@{$settings},[-99,-99,$soaktemp]);
}

if (0) {
    for (my $theta=10;$theta<=45;$theta+=5) {
	for (my $phi=-90;$phi<=90;$phi+=30) {
	    foreach my $soaktemp (@soaktemps) {
		push(@{$settings},[$theta,$phi,$soaktemp]);
	    }
	}
    }
}


my $row=0;
my $col;
my $sb=[];

foreach my $setlist (@{$settings}) {
    $row++;
    my ($theta,$phi,$soaktemp)=@{$setlist};
    my $cmd;
    my $gravargs="";
    my $tempargs="";

    if ($soaktemp != -99) {
	$tempargs=" -T $soaktemp ";
    }

    if (!(($theta == -99) || ($phi == -99))) {
	$gravargs=" -G -theta $theta -phi $phi ";
    }

    $cmd="CAMERA_perturbations.perl -z3up $gravargs $tempargs";

    `rm {L,F}*dat`;
    print "cmd=$cmd\n";
    `$cmd`;

    # now populate the tables with these values
    my $hduver=1;
    foreach my $label (keys %solidbodies) {
	$fptr->movnam_hdu(BINARY_TBL,$label,$hduver,$status);
	check_error($status);
	$col=0;
	$col++;    $fptr->write_col_flt($col,$row,1,1,[$theta],$status);
	$col++;    $fptr->write_col_flt($col,$row,1,1,[$phi],$status);
	$col++;    $fptr->write_col_flt($col,$row,1,1,[$soaktemp],$status);
	parse_solidbody($sb,$solidbodies{$label});
	for (my $ix=0;$ix<=$#{$sb};$ix++) {
	    $col++;    $fptr->write_col_flt($col,$row,1,1,[$sb->[$ix]],$status);
	}
    }

    foreach my $label (keys %distbodies) {
	$fptr->movnam_hdu(BINARY_TBL,$label,$hduver,$status);
	check_error($status);
	$col=0;
	$col++;    $fptr->write_col_flt($col,$row,1,1,[$theta],$status);
	$col++;    $fptr->write_col_flt($col,$row,1,1,[$phi],$status);
	$col++;    $fptr->write_col_flt($col,$row,1,1,[$soaktemp],$status);

	my $distortion={};

	parse_distortion($distortion,$distbodies{$label});

	$col++;    $fptr->write_col_flt($col,$row,1,
					$#{$distortion->{"x"}}-$[+1,
					$distortion->{"x"},$status);
	$col++;    $fptr->write_col_flt($col,$row,1,
					$#{$distortion->{"y"}}-$[+1,
					$distortion->{"y"},$status);
	$col++;    $fptr->write_col_flt($col,$row,1,
					$#{$distortion->{"z"}}-$[+1,
					$distortion->{"z"},$status);
    }

    foreach my $label (keys %zernbodies) {
	$fptr->movnam_hdu(BINARY_TBL,$label,$hduver,$status);
	check_error($status);
	$col=0;
	$col++;    $fptr->write_col_flt($col,$row,1,1,[$theta],$status);
	$col++;    $fptr->write_col_flt($col,$row,1,1,[$phi],$status);
	$col++;    $fptr->write_col_flt($col,$row,1,1,[$soaktemp],$status);
	parse_zernike_distortion($sb,$zernbodies{$label});
	for (my $ix=0;$ix<=$#{$sb};$ix++) {
	    next if ($ix >= $maxzern);
	    $col++;    $fptr->write_col_flt($col,$row,1,1,[$sb->[$ix]],$status);
	}
    }

    if ($compensate_camera==1) {
	# here compute camera compensation
	# this is done in raytrace coordinates for historical reasons. 
	# should change this, although solid body compensation is simple
	# and transformation from RT to CCS should be transparent
	print "cmd(rt)=$cmd -rt\n";
	`$cmd -rt`;
	`rm -rf HPOD_soln_run*`;
	`HPOD_soln.perl \*CAMmisalign_{t{z,x,y},r{x,y}}`;

	foreach my $label (keys %compbodies) {
	    $fptr->movnam_hdu(BINARY_TBL,$label,$hduver,$status);
	    check_error($status);
	    $col=0;
	    $col++;    $fptr->write_col_flt($col,$row,1,1,[$theta],$status);
	    $col++;    $fptr->write_col_flt($col,$row,1,1,[$phi],$status);
	    $col++;    $fptr->write_col_flt($col,$row,1,1,[$soaktemp],$status);
	    parse_solidbody($sb,$compbodies{$label});
	    # transform to CCS coordinate system
	    my @tmp=($sb->[1],$sb->[0],-$sb->[2],
		     $sb->[4],$sb->[3],-$sb->[5]);
	    @{$sb}=@tmp;
	    # and write
	    for (my $ix=0;$ix<=$#{$sb};$ix++) {
		$col++;    $fptr->write_col_flt($col,$row,1,1,[$sb->[$ix]],$status);
	    }
	}
	# and insert the FOM
	my $best_fom=`awk 'NR>2{print \$NF}' *tnt | sort -n | awk 'NR==1'`;
	chomp $best_fom;
	$col++;    $fptr->write_col_flt($col,$row,1,1,[$best_fom],$status);
    }
}

sub parse_solidbody {
    my ($handle,$filename)=@_;
    my $contents=`cat $filename`;
    chomp $contents;
    my @t=($contents =~ /-t\s+(\S+)\s+(\S+)\s+(\S+)/);
    my @r=($contents =~ /-r\s+(\S+)\s+(\S+)\s+(\S+)/);
    @{$handle}=(@t,@r);
}

sub parse_distortion {
    my ($handle,$filename)=@_;
    %{$handle}=();
    $handle->{"x"}=[];    $handle->{"y"}=[];    $handle->{"z"}=[];
    open(F,$filename) || die "can't open $filename";
    my $line;
    $line=<F>;
    $line=<F>;
    while ($line=<F>) {
	chomp($line);
	my @contents=split(" ",$line);
	push($handle->{"x"},$contents[0]);
	push($handle->{"y"},$contents[1]);
	push($handle->{"z"},$contents[2]);
    }
}

sub parse_zernike_distortion {
    my ($handle,$filename)=@_;
    my $contents=`cat $filename`;
    chomp $contents;
    my @t=($contents =~ /-Z\S*\s+(\S+)/);
    my @z=split(':',$t[0]);
    @{$handle}=(@z);
}

sub check_error {
    my ($stat)=@_;
    if ($stat) {
        fits_report_error(*STDERR,$stat);
    }
}

sub fingerprint {
    my ($fp)=@_;
    my $cmt;
    $cmt="generated by andy rasmussen arasmus\@slac.stanford.edu";
    $fp->write_comment($cmt,$status);
    $cmt="using generate_misalignments_package.perl in $ENV{PWD}";
    $fp->write_comment($cmt,$status);
    $fp->write_date($status);
}
