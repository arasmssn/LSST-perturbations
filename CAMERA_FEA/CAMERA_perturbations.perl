#!/usr/bin/perl -w
# this script collects up data products of FEA analysis done over various
# runs (applying correct interpretations) and computes a set of misalignments
# and perturbations appropriate for use with (1) raytrace and (2) imsim/phosim
# the inputs are (1) the transformation specs prepared on an optic basis
# from FEA analysis, (2) the zernike expansion of the optic, perpared on a 
# surface basis.

# the zernike expansion is expected to be in a form compatible with raytrace,
# i.e. "-z z0:z1:z2:z3.."
# and the solidbody transform specs are also expected to be compatible with
# raytrace, i.e. "-t tx ty tz -r rx ry rz" - where ordering of -t and -r should
# be flexible. For the time being, the "-d" switch won't be permitted.

use warnings;
use strict;
use Scalar::Util qw(looks_like_number);
use POSIX;

my $ThisScript="CAMERA_perturbations.perl";

my $DISTORTION_DB_PATH="/home/arasmus/camera_perturbations/DB"; # this can be gotten from an environmnetal variable

my $DBKEY="CAMERA_PERTURBATIONS_DB";
if (!defined($ENV{$DBKEY})) {
    printf STDERR "\n";
    printf STDERR "Please define the environment variable $DBKEY\n";
    printf STDERR "\n";
    printf STDERR "e.g.  for bash:\n";
    printf STDERR "[user\@host]\$ $DBKEY=\"my_database_directory\"\n";
    printf STDERR "[user\@host]\$ export $DBKEY\n";
    printf STDERR "      for csh/tcsh:\n";
    printf STDERR "[user\@host]\$ setenv $DBKEY \"my_database_directory\"\n";
    printf STDERR "\n$0 needs to locate either the raw MS Excel files\n";
    printf STDERR "to perform a reanalysis and to store an intermediate analysis file --\n";
    printf STDERR "or an instance of a ready-made intermediate analysis file\n";
    printf STDERR "(optionally specified by the --intermediate_tnt=<filename> switch)\n\n";
    printf STDERR "This script looks in the location specified by \$$DBKEY for these.\n\n";
    
    printf STDERR "For more information, use --help or --man switches.\n";
    printf STDERR "exiting..\n\n";
    exit(1);
}

$DISTORTION_DB_PATH=$ENV{$DBKEY};

my $system="ccs";
my @bodies=("L1","L2","L3","F","FP");
my @G_conditions=("GX","GY","GZ");
my @T_conditions=("Cold","Nom","Hot");
my %taskdir;
@taskdir{@G_conditions,@T_conditions}=(("gravity")x3,("thermal")x3);
my @sides=("S1","S2");

my %nominal_env=("GX" => "1,0,0","GY" => "0,1,0","GZ" => "0,0,1",
		 "Cold" =>  -3,	 "Nom"  =>  +7,	 "Hot"  => +17);
	  
# outlier environments applied in rogue FEA calculations:
my %exception_env=(join($;,"L3","GX") => "+1, 0, 0",
		   join($;,"F","GX")  => "-1, 0, 0",
		   join($;,"FP","GX") => " 0, 1, 0",
		   join($;,"FP","GY") => " 1, 0, 0" ); 

# read in commandline arguments to specify soak temperature and also 
# gravity load direction (from theta(zenith distance) and phi (camera rotation)
# some interpretation of phi is necessary because document-1614 is ambiguous.

my $enable_G=0;
my $enable_T=0;
my $T_camera=20;

if ($#ARGV<0) {
    printf "no args!\nexiting..\n";
    printf "usage:\n";
    printf "CAMERA_perturbations.perl \t[-G] [-T <T_camera>] \n".
	"\t\t\t\t\t[-theta <zenith dist>] [-phi <cam_rot>] [-z3up] [-zall] [-rt]\n".
	"\t\t\t\t\t[-no_precomp] [-zlim <zern. rad. order limit>]\n";
    printf "interpolates available FEA output for conditions specified.\n";
    printf "\n\n";
    printf "-G enables interpolation of the gravity load calculations\n";
    printf "-T enables interpolation of the thermal distortion calculations\n";
    printf "<T_camera> specifies the outer skin temperature for the camera body\n";
    printf "<zenith dist> specifies the gravity vector direction relative to Zccs\n";
    printf "<cam_rot> specifies the camera rotation angle [-90,+90].\n";
    printf "-no_precomp specifies that perturbations relative to zero-G and passive room temperature assembly be used.\n";
    printf "-zlim limits the radial order to <zern.rad.order.limit>.\n";
    printf "-rt performs the calculation for raytrace coordinate system:\n";
    printf "\tconnected to CCS via transformation\n";
    printf "/ 0  1  0 \\\n";
    printf "| 1  0  0 |\n";
    printf "\\ 0  0 -1 /\n";
    exit(1);
}


my ($theta,$phi)=(0,0);
my $z3up=1;
my $zall=0;

my $deg=atan2(1,1)/45.0;

#  need to work the following into commandline options at some point.
#  for now hard-code these. e.g. 

($theta,$phi)=(0,0);

my $precomp_figures=0;

my @precomp_gvector=(sin($theta*$deg)*sin(-$phi*$deg),
		     sin($theta*$deg)*cos(-$phi*$deg),
		     cos($theta*$deg));
my $precomp_Tcond=20;

# using switch -no_precomp has the effect of undefining these two parameters:
# undef(@precomp_gvector);
# undef($precomp_Tcond);

my $max_zern_rad_order=9;

while ($_=$ARGV[0],shift) {
    if (/^-/) {
	/^-G/     && ($enable_G=1,next);
	/^-z3up/  && ($z3up=1,next);
	/^-zall/  && ($z3up=0,next);
	/^-rt/    && ($system="rt",next);
	/^-theta/ && ($theta=$ARGV[0],shift,next);
	/^-phi/   && ($phi=$ARGV[0],shift,next);
	/^-T/     && ($enable_T=1,$T_camera=$ARGV[0],shift,next);
	/^-zlim/  && ($max_zern_rad_order=$ARGV[0],shift,next);
	/^-no_precomp/ && 
	    (undef(@precomp_gvector),undef($precomp_Tcond),next);
	die "didn't expect arg: $_\n";
    } else {
	die "didn't expect arg: $_\n";
    }
}


if ($max_zern_rad_order>9) {
    printf STDERR "Error - Maximum radial order for zernikes should ".
	"not exceed the radial order by which the coefficient expansion ".
	"was originally characterized. Currently this is 7. exiting..\n";
    exit;
}

my ($my_files,$my_env)=read_database();
my %my_files=%{$my_files};
my %my_env=%{$my_env};

if (!looks_like_number($T_camera)) {
    printf "argument supplied following -T ".
	"doesn't look like a temperature: $T_camera\nexiting..\n";
    exit(1);
}

# go on to set up the interpolation

my @gvector=(sin($theta*$deg)*sin(-$phi*$deg),
	     sin($theta*$deg)*cos(-$phi*$deg),
	     cos($theta*$deg));


my @scale;


my $db={};

foreach my $body (@bodies) {
    printf "%s\n",$body;
    if ($body eq "FP") {
	$db->{$body}={ "contribs" => {"solidbody" => []},
		       "weights"  => {"solidbody" => []} };
    } else {
	$db->{$body}={ "contribs" => {"solidbody" => [],
				      "S1"        => [],
				      "S2"        => [],
				      "S1"."tnt"  => [],
				      "S2"."tnt"  => []   },
		       "weights"  => {"solidbody" => [],
				      "S1"        => [],
				      "S2"        => [],
				      "S1"."tnt"  => [],
				      "S2"."tnt"  => []   }};
    }
    if ($enable_G) {
	foreach my $condition (@G_conditions) {
	    my @cond_vector=split(',',$my_env{$body,$condition});
	    my $scale=dot_prod(@gvector,@cond_vector);

	    my $precomp_scale;

	    if (@precomp_gvector) {
		$precomp_scale=dot_prod(@precomp_gvector,@cond_vector);
	    } else {
		$precomp_scale=0;
	    }

	    if ($scale*$scale > 1e-8 || $precomp_scale*$precomp_scale > 1e-8) {
		if (defined($my_files{$body,$condition})) {
		    push($db->{$body}->{"contribs"}->{"solidbody"},
			 $my_files{$body,$condition});
		    push($db->{$body}->{"weights"}->{"solidbody"},
			 $scale-$precomp_scale);
		}
		foreach my $side (@sides) {
		    # for distortions (zernike expansion format)
		    if (defined($my_files{$body,$condition,$side})) {
			push($db->{$body}->{"contribs"}->{$side},
			     $my_files{$body,$condition,$side});

			# weighting of deformations do not include
			# precompensation..
			if ($precomp_figures==1) {
			    push($db->{$body}->{"weights"}->{$side},
				 $scale-$precomp_scale);
			} else {
			    push($db->{$body}->{"weights"}->{$side},$scale);
			}
		    }

		    # for distortions (tnt format)
		    if (defined($my_files{$body,$condition,$side,"tnt"})) {
			push($db->{$body}->{"contribs"}->{$side."tnt"},
			     $my_files{$body,$condition,$side,"tnt"});

			# weighting of deformations do not include
			# precompensation..
			if ($precomp_figures==1) {
			    push($db->{$body}->{"weights"}->{$side."tnt"},
				 $scale-$precomp_scale);
			} else {
			    push($db->{$body}->{"weights"}->{$side."tnt"},
				 $scale);
			}
		    }

		}
	    }
	}
    }

    if ($enable_T) {

	my @sorted_Tconds = sort {pow($my_env{$body,$a}-$T_camera,2) <=> 
				   pow($my_env{$body,$b}-$T_camera,2)} 
	                       @T_conditions;

	$scale[0]=($T_camera-$my_env{$body,$sorted_Tconds[1]})/
	    ($my_env{$body,$sorted_Tconds[0]}-$my_env{$body,$sorted_Tconds[1]});
	$scale[1]=1-$scale[0];

	my @precomp_sorted_Tconds;
	my @precomp_scale;
	
	if (defined($precomp_Tcond)) {
	    @precomp_sorted_Tconds = sort 
	    {pow($my_env{$body,$a}-$precomp_Tcond,2) <=> 
		 pow($my_env{$body,$b}-$precomp_Tcond,2)} @T_conditions;
	    $precomp_scale[0]=
		($precomp_Tcond-$my_env{$body,$precomp_sorted_Tconds[1]})/
		($my_env{$body,$precomp_sorted_Tconds[0]}-
		 $my_env{$body,$precomp_sorted_Tconds[1]});
	    $precomp_scale[1]=1-$precomp_scale[0];
	} else {
	    @precomp_scale[0,1]=(0,0);
	}

	for (my $i=0;$i<=1;$i++) {
	    next if ($scale[$i]*$scale[$i]<1e-8);

	    my $condition=$sorted_Tconds[$i];

	    if (defined($my_files{$body,$condition})) {

		push($db->{$body}->{"contribs"}->{"solidbody"},
		     $my_files{$body,$condition});
		push($db->{$body}->{"weights"}->{"solidbody"},$scale[$i]);
	    }

	    if (defined($precomp_Tcond) && 
		defined($my_files{$body,$precomp_sorted_Tconds[$i]})) {

		my $precomp_condition=$precomp_sorted_Tconds[$i];
		push($db->{$body}->{"contribs"}->{"solidbody"},
		     $my_files{$body,$precomp_condition});
		push($db->{$body}->{"weights"}->{"solidbody"},
		     -1*$precomp_scale[$i]);
	    }

	    foreach my $side (@sides) {
		if (defined($my_files{$body,$condition,$side})) {
		    push($db->{$body}->{"contribs"}->{$side},
			 $my_files{$body,$condition,$side});
		    if ($precomp_figures==1) {
			push($db->{$body}->{"weights"}->{$side},
			     $scale[$i]-$precomp_scale[$i]);
		    } else {
			push($db->{$body}->{"weights"}->{$side},$scale[$i]);
		    }
		}
	    }
	}
    }
}

my %outputfile;

foreach my $body (@bodies) {
    $outputfile{$body,"solidbody"}=$body."solidbody_g.dat";
    $outputfile{$body,"S1"}=$body."_S1_zern_g.dat";
    $outputfile{$body,"S2"}=$body."_S2_zern_g.dat";
    $outputfile{$body,"S1"."tnt"}=$body."_S1_perturbation_g.tnt";
    $outputfile{$body,"S2"."tnt"}=$body."_S2_perturbation_g.tnt";
}

my @outfiles=();

if (0) {
    # this part is incomplete because it lacks a translation between
    # solid body tilts and z1,2. Need radius of each optic and the sense
    # relationship, tilt vs. zernike.
    foreach my $body (@bodies) {
	printf "\nFOR BODY $body:\n";
	my %return_string=();
	foreach my $key ( sort keys $db->{$body}->{"contribs"} ) {
	    if ($key eq "S1" || $key eq "S2") {
		$return_string{$key}=	
		    interp_zernike_expansions($db->{$body}->{"contribs"}->{$key},
				       $db->{$body}->{"weights"}->{$key});
	    } elsif ($key eq "solidbody") {
		$return_string{$key}=	
		    interp_solidbody($db->{$body}->{"contribs"}->{$key},
				     $db->{$body}->{"weights"}->{$key});
	    } 
	}

	my $have_zernikes=(!$z3up && 
			   (defined($db->{$body}->{"contribs"}->{"S1"})||
			    defined($db->{$body}->{"contribs"}->{"S2"})));

	printf STDERR "HAVE ZERNIKES = $have_zernikes\n";
	if ($have_zernikes==1) {
	    # to avoid double-counting, replace solidbody terms
	    # tz, rx, ry with average of terms from S1 & S2 for z0..2 
	    # and reduce their amplitudes by these solid body offsets.
	    if (! (defined($db->{$body}->{"contribs"}->{"S1"}) && 
		   defined($db->{$body}->{"contribs"}->{"S2"}))) {
		die("only one side's zernike expansion is defined.\n".
		    "exiting..\n");
	    }
	    # modify the %return_strings
	    my (@tmp,@zl);
	    $tmp[0]=[split(' ',$return_string{"S1"})];
	    $tmp[1]=[split(' ',$return_string{"S2"})];
	    $zl[0]=[split(':',$tmp[0]->[1])];
	    $zl[1]=[split(':',$tmp[1]->[1])];
	    # compute average piston, tip & tilt from @zl[0,1]
	    my ($piston,$tip,$tilt);
	    $piston=($zl[0]->[0] + $zl[1]->[0])/2.0;
	    $tip =($zl[0]->[1] + $zl[1]->[1])/2.0;
	    $tilt=($zl[0]->[2] + $zl[1]->[2])/2.0;
	    printf STDERR "will subtract off piston $piston from $zl[0]->[0] and $zl[1]->[0]\n";
	    $zl[0]->[0] -= $piston;	    $zl[1]->[0] -= $piston;
	    $zl[0]->[1] -= $tip;	    $zl[1]->[1] -= $tip;
	    $zl[0]->[2] -= $tilt;	    $zl[1]->[2] -= $tilt;
	    $tmp[0]->[1]=join(':',@{$zl[0]});  
	    $tmp[1]->[1]=join(':',@{$zl[1]});
	    $return_string{"S1"}=join(' ',@{$tmp[0]});
	    $return_string{"S2"}=join(' ',@{$tmp[1]});
	} else {
	    # either neither zernike expansion is defined, or $z3up=1.
	    # overcounting may happen
	}

	foreach my $key ( sort keys $db->{$body}->{"contribs"} ) {
	    if ($key eq "S1") {
		# replace -Z with -Z1 and output
		$return_string{$key} =~ s/-Z/-Z1/g;
	    } elsif ($key eq "S2") {
		# replace -Z with -Z2 and output
		$return_string{$key} =~ s/-Z/-Z2/g;
	    } else {
		# shouldn't reach here.
	    }
	    push(@outfiles,$outputfile{$body,$key});
	    open(F,">$outputfile{$body,$key}") || die;
	    printf F $return_string{$key};
	    close(F);
	}
    }

} else {
    # this approach simply zeros out zernikes 0,1 & 2 
    # if the -z3up switch is given, leaving those coefficients 
    # alone if -z3up is not given. This is according to the assumption that 
    # the solid body transformation terms Tz, Rx & Ry are simply related to
    # zernike coefficients zO, z1 and z2 (according to ANSI Z80.28). 
    # This assumption may break down in the regime where sag at the apex is 
    # comparable to or larger than the mean solid body axial displacement,
    # or where the density of FEA calculation nodes are nonuniform across the
    # surfaces.
    # If in doubt, do not use the zernike expansion representation of the surface
    # error - which was fit directly using the z displacement distribution across 
    # the surfaces. Use instead the corresponding *perturbation_g.tnt files, together
    # with the *solidbody.dat file. Those *perturbation_g.tnt files were generated after 
    # having removed the solid body terms from the displacement map, so this is arguably
    # the more straightforward option.
    foreach my $body (@bodies) {
	printf "\nFOR BODY $body:\n";
	foreach my $key ( sort keys $db->{$body}->{"contribs"} ) {
	    my $nontrivial_zernikes=1;
	    my $return_string;
	    if ($key eq "S1" || $key eq "S2") {
		$return_string=
		    interp_zernike_expansions($db->{$body}->{"contribs"}->{$key},
				       $db->{$body}->{"weights"}->{$key});
		{
		    # check the return
		    my @ret=split(' ',$return_string);
		    my @zl=split(':',$ret[1]);
		    foreach my $zl (@zl) {
			goto nontrivial if ($zl!=0);
		    }
		    $nontrivial_zernikes=0;
		  nontrivial:
		}

		if ($z3up==1) {
		    # replace z0..2 with zeros since they will be included in
		    # solidbody file
		    my @lst=split(' ',$return_string);
		    my @zrn=split(':',$lst[1]);
		    @zrn[0,1,2]=(0,0,0);
		    print STDERR "adjusting: $return_string\n";
		    $lst[1]=join(':',@zrn);
		    $return_string=join(' ',@lst);
		    print STDERR "to: $return_string\n";
		}

		if ($key eq "S1") {
		    # replace -Z with -Z1 and output
		    $return_string =~ s/-Z/-Z1/g;
		} elsif ($key eq "S2") {
		    # replace -Z with -Z2 and output
		    $return_string =~ s/-Z/-Z2/g;
		} else {
		    # shouldn't reach here.
		}
	    } elsif ($key eq "S1"."tnt" || $key eq "S2"."tnt") {
		$return_string=
		    interp_perturbations($db->{$body}->{"contribs"}->{$key},
					 $db->{$body}->{"weights"}->{$key});
	    } elsif ($key eq "solidbody") {
		$return_string=
		    interp_solidbody($db->{$body}->{"contribs"}->{$key},
				     $db->{$body}->{"weights"}->{$key},
				     (!$z3up && 
				      (defined($db->{$body}->{"contribs"}->{"S1"})||
				       defined($db->{$body}->{"contribs"}->{"S2"}))));
	    } else {
		printf STDERR "unexpected key: $key\nexiting..\n";
		exit(1);
	    }

	    next if (length($return_string) == 0);
	    next if (($key eq "S1" || $key eq "S2") 
		     && ($nontrivial_zernikes==0));
	    push(@outfiles,$outputfile{$body,$key});
	    open(F,">$outputfile{$body,$key}") || die;
	    printf F $return_string;
	    close(F);
	}
    }
}




printf STDOUT "\ndone! new files here (evaluated in $system system) are:\n%s\n\n",join("\n",@outfiles);

exit;
foreach my $body (@bodies) {
    foreach my $condition (@G_conditions,@T_conditions) {
	if (defined($my_files{$body,$condition})) {
	    printf "have solidbody file for $body $condition\n";
	    printf "environment: $my_env{$body,$condition}\n"
	}
	foreach my $side (@sides) {
	    if (defined($my_files{$body,$condition,$side})) {
		printf "have distortion file for $body $condition $side\n";
	    }
	}
    }
}

sub interp_perturbations {
    my ($contribs,$weights)=@_;
    my $line;
    my $ret="";
    my ($x_arr,$y_arr,$z_arr)=([],[],[]);
    for (my $i=0;$i<=$#{$contribs};$i++) {
	my $infile=`basename $contribs->[$i]`;
	chomp $infile;
	printf "will weight %s by %f\n",$infile,$weights->[$i];
	open(Q,$contribs->[$i]) || die "can't open file $contribs->[$i]";
	printf STDERR "GOT file $contribs->[$i]\n";
	if ($i==0) {
	    $line=<Q>;	    $ret .= $line;
	    $line=<Q>;	    $ret .= $line;
	    while ($line=<Q>) {
		chomp $line;
		my @ent=split(' ',$line);
		push($x_arr,$ent[0]);
		push($y_arr,$ent[1]);
		push($z_arr,$ent[2]*$weights->[$i]);
	    }
	} else {
	    $line=<Q>;
	    $line=<Q>;
	    my $j=-1;
	    while ($line=<Q>) {
		chomp $line;
		$j++;
		my @ent=split(' ',$line);
		# be diligent - check for match between x & y values
		if (($ent[0] != $x_arr->[$j]) || ($ent[1] != $y_arr->[$j])) {
		    printf STDERR "mismatch between node coordinates!\n";
		    printf STDERR "for perturbation entry $j between files:\n";
		    printf STDERR " %s and %s:\n",@{$contribs}[0,$i];
		    printf STDERR 
			"((x,y)[0],(x,y)[$i]) = ((%s,%s),(%s,%s))\n",
			@ent[0,1],$x_arr->[$j],$y_arr->[$j];
		    printf STDERR "exiting..\n";
		    exit(1);
		}
		$z_arr->[$j] += $ent[2]*$weights->[$i];
	    }
	}
	close(Q);
    }
    for (my $j=0;$j<=$#{$z_arr};$j++) {
	$ret .= sprintf("%s %s %g\n",
			1000*$x_arr->[$j],
			1000*$y_arr->[$j],
			1000*$z_arr->[$j]);
    }
    return($ret);
}

sub interp_zernike_expansions {
    my ($contribs,$weights)=@_;
    my (@contribs,@weights,$i,@contents,$line,@list,@zlist,@Zlist);

    @Zlist=();
    @contribs=@{$contribs};
    @weights=@{$weights};
    for (my $i=0;$i<=$#contribs;$i++) {
	my $infile=`basename $contribs[$i]`;
	my $zlist;
	chomp $infile;
	printf STDERR "will weight %s by %f\n",$infile,$weights[$i];
	@zlist=();
	open(F,$contribs[$i]) || die;
	printf STDERR "GOT file $contribs[$i]\n";
	@contents=<F>;
	close(F);
	foreach $line (@contents) {
	    @list=split(' ',$line);
	    while ($_ = shift @list) {
		/-Z/ && ($zlist=shift @list,next);
		printf STDERR "unexpected list in zernike expansion\n";
		printf STDERR "file = $contribs[$i] offending string $_\n";
		printf STDERR "exiting..\n";
		exit(1);
	    }
	}
	@zlist=split(':',$zlist);
	printf STDERR "GOT one with $#zlist entries.\n";
	for (my $j=0;$j<=$#zlist;$j++) {
	    if (defined($Zlist[$j])) {
		$Zlist[$j] += $zlist[$j]*$weights[$i];
	    } else {
		$Zlist[$j]   = $zlist[$j]*$weights[$i];
	    }
	}
    }

    printf STDERR "GOT grand total of $#Zlist zernike coefficient entries.\n";
#    printf STDERR "GOT1 %s\n",join('::',@Zlist);

    # originally:
    my $nzern=1+2+3+4+5+6+7-1;
    # now responding to possible commmandline specification:
    {
	$nzern=0;
	for (my $i=1;$i<=$max_zern_rad_order;$i++) {
	    $nzern += $i;
	}
	$nzern -= 1;
    }
    {
	my @tmplist=();
	for (my $i=0;$i<=$nzern;$i++) {
	    if (defined($Zlist[$i])) {
		push(@tmplist,$Zlist[$i]);
	    } else {
		push(@tmplist,0);
	    }
	}
	@Zlist=@tmplist;
    }
#    printf STDERR "GOT and here is nzern: $nzern\n";
#    printf STDERR "GOT2 %s\n",join('::',@Zlist);
    sprintf("-Z %s\n",join(':',@Zlist[0..$nzern]));
}

sub interp_solidbody {
    # ignore $have_zernikes for now
    my ($contribs,$weights,$have_zernikes)=@_;
    my (@contribs,@weights,$i,@contents,$line,@list,@tx,@ty,@tz,@rx,@ry,@rz);
    my ($degrees);
    my ($tx,$ty,$tz,$rx,$ry,$rz);

    $degrees=0;

    @tx=();    @ty=();    @tz=();
    @rx=();    @ry=();    @rz=();

    @contribs=@{$contribs};
    @weights=@{$weights};
    for ($i=0;$i<=$#contribs;$i++) {
	my $infile=`basename $contribs[$i]`;
	chomp $infile;
	printf STDOUT ("will weight %s by %f\n",$infile,$weights[$i]);
#	printf STDOUT ("will weight %s by %f (have zernikes? %d)\n",
#		       $infile,$weights[$i],$have_zernikes);
	open(F,$contribs[$i]) || die;
	@contents=<F>;
	($tx,$ty,$tz,$rx,$ry,$rz)=(0,0,0,0,0,0);
	foreach $line (@contents) {
	    @list=split(' ',$line);
	    while ($_ = shift @list) {
		/-t/ && ($tx = shift @list,
			 $ty = shift @list,
			 $tz = shift @list,
			 next);
		/-r/ && ($rx = shift @list,
			 $ry = shift @list,
			 $rz = shift @list,
			 next);
		/-T/ && next;
		/--/ && next;
		/-d/ && ($degrees=1,next);
		printf STDERR "unexpected list in solidbody tranformation\n";
		printf STDERR "file = $contribs[$i] offending string $_\n";
		printf STDERR "exiting..\n";
		exit(1);
	    }
	}
	close(F);

#	($tz,$rx,$ry)=(0,0,0) if ($have_zernikes);

	if ($degrees) {
	    $deg=atan2(1,1)/45.0;
	} else {
	    $deg=1;
	}
	push(@tx,$tx*$weights[$i]);	push(@rx,$rx*$deg*$weights[$i]);
	push(@ty,$ty*$weights[$i]);	push(@ry,$ry*$deg*$weights[$i]);
	push(@tz,$tz*$weights[$i]);	push(@rz,$rz*$deg*$weights[$i]);
    }
    # now collect up the equivalent misalignment
    $tx=sum(@tx);
    $ty=sum(@ty);
    $tz=sum(@tz);
    $rx=sum(@rx);
    $ry=sum(@ry);
    $rz=sum(@rz);
    sprintf("-T -t %g %g %g -r %g %g %g --\n",$tx,$ty,$tz,$rx,$ry,$rz);
}

sub sum {
    my(@list)=@_;
    my ($i,$sum);
    $sum=0;
    for ($i=0;$i<=$#list;$i++) {
	$sum+=$list[$i];
    }
    $sum;
}

sub dot_prod {
    my ($x1,$y1,$z1,$x2,$y2,$z2)=@_;
    $x1*$x2+$y1*$y2+$z1*$z2;
}

sub read_database {
    my (@files,%files,%environ);
    @files=();
    %files=();
    my ($match,$side);

    foreach my $body (@bodies) {
	foreach my $condition (@G_conditions,@T_conditions) {
	    $match=$DISTORTION_DB_PATH."/${body}*/$taskdir{$condition}/${body}_*${condition}*${system}_coords_full_solidbody.dat";
	    @files=`ls $match 2>/dev/null`;
	    chomp @files;
	    if ($#files<0) {
		# there is no such file in the database
		printf STDERR "no such file in the database. looking for: %s\n",$match;
		exit(1);
	    } elsif ($#files>0) {
		# ambiguous - which file to use.
		printf STDERR "ambiguity in camera distortion database - which file to use for body $body condition $condition: %s\n",join(' ',@files);
		printf STDERR "please see code for $ThisScript\nexiting..\n";
		exit(1);
	    } else {
		$files{$body,$condition}=$files[0];
		$environ{$body,$condition}=
		    (defined($exception_env{$body,$condition}))?
		     $exception_env{$body,$condition}:$nominal_env{$condition};
	    }

	    foreach $side (@sides) {
		# locate the fit parameters
		$match=$DISTORTION_DB_PATH."/${body}*/$taskdir{$condition}/${body}_*${side}_${condition}*${system}_coords_fit.dat";
		@files=`ls $match 2>/dev/null`;
		chomp @files;
		if ($#files<0) {
		    # there is no such file in the database
		} elsif ($#files>0) {
		    # ambiguous - which file to use.
		    printf STDERR "ambiguity in camera distortion database - which file to use for body $body condition $condition side $side: %s\n",join(' ',@files);
			printf STDERR "please see code for $ThisScript\nexiting..\n";
		    exit(1);
		} else {
		    # $#files==0. no ambiguity. assign.
		    $files{$body,$condition,$side}=$files[0];
		}

		# locate the distortion tnt files
		$match=$DISTORTION_DB_PATH."/${body}*/$taskdir{$condition}/${body}_*${side}_${condition}*${system}_coords_perturbation.tnt";
		@files=`ls $match 2>/dev/null`;
		chomp @files;
		if ($#files<0) {
		    # there is no such file in the database
		} elsif ($#files>0) {
		    # ambiguous - which file to use.
		    printf STDERR "ambiguity in camera distortion database - which file to use for body $body condition $condition side $side: %s\n",join(' ',@files);
			printf STDERR "please see code for $ThisScript\nexiting..\n";
		    exit(1);
		} else {
		    # $#files==0. no ambiguity. assign.
		    $files{$body,$condition,$side,"tnt"}=$files[0];
		}
	    }
	}
    }
    my ($ret_files,$ret_env);
    $ret_files={};
    $ret_env={};
    %{$ret_files}=%files;
    %{$ret_env}=%environ;
    return($ret_files,$ret_env);
}
