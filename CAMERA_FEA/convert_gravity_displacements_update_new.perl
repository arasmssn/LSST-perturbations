#!/usr/bin/perl -w
use warnings;
use POSIX;
use Math::Amoeba qw(MinimiseND);
use Scalar::Util qw(looks_like_number);

# this script is for converting J.Ku's updated FEA of 120620
# this script is for converting J.Ku's updated FEA of 120911:
# re-using surface indices 6 thru 9 (J.Ku re-used node labels 100000 thru 116368)

use POSIX;

# %optic is only for naming output files..
%optic=("0" => L3,
	"1" => L3,
	"2" => L2,
	"3" => L2,
	"4" => L1,
	"5" => L1,
	"6" => L2_new,
	"7" => L2_new,
	"8" => L1_new,
	"9" => L1_new);
# update
@optic{"6","7","8","9"}=("L1_new","L1_new","L2_new","L2_new");

%surface=("0" => 2,
	  "1" => 1,
	  "2" => 2,
	  "3" => 1,
	  "4" => 2,
	  "5" => 1,
	  "6" => 2,
	  "7" => 1,
	  "8" => 2,
	  "9" => 1);

@surface{"6","7","8","9"}=(1,2,1,2);

# this $nzern is used for evaluating zernike fields used later in fitting. large values for now..
$nzern=(1+2+3+4+5+6+7+8+9);
# (28 zernikes)
#$nzern=(1+2+3+4+5+6+7+8+9);

$file="L1-L2-Deformations/GRID\ COORDINATES.csv";

@nodelist=([],[],[],[],[],[]);
%nodecoords=();

@surfs=(6,7,8,9);

# rather than work out the zernike mixing matrix for coordinate transformations
# just specify whether this is for raytrace coordinates (or CCS/MCS).
# the tranformation between the two systems goes as:
# /x\   / 0 1  0 \ /x\
# |y| = | 1 0  0 | |y|
# \z/rt \ 0 0 -1 / \z/ccs
#

foreach $raytrace_coords ( 0,1 ) {

    $rt_string=($raytrace_coords==1)?"rt_coords":"ccs_coords";

    %ro=(0=>0.361,
	 1=>0.361,
	 2=>0.551,
	 3=>0.551,
	 4=>0.775,
	 5=>0.775,
	 6=>0.551,
	 7=>0.551,
	 8=>0.775,
	 9=>0.775);
# update
    @ro{6,7,8,9}=(0.775,0.775,0.551,0.551);

    open(F,$file) || die "can't open file $file";

    $nskip=0;

    while ($line=<F>) {
	chomp($line);
	next if ( $line !~ /^GRID/ );
	@this=split(',',$line);
	($node,$x,$y,$z)=@this[1,3,4,5];

	if (0) {
	    if ((!looks_like_number($x)) || 
		(!looks_like_number($y)) ||
		(!looks_like_number($z))) {
		$nskip++;
		next;
	    }
	} else {
	    $x=0 if (!looks_like_number($x));
	    $y=0 if (!looks_like_number($y));
	    $z=0 if (!looks_like_number($z));
	}

	if (($node-10001)*(10481-$node)>=0) {
	    $surf=0;
	} elsif (($node-10501)*(10981-$node)>=0) {
	    $surf=1;
	} elsif (($node-11001)*(11781-$node)>=0) {
	    $surf=2;
	} elsif (($node-12001)*(12781-$node)>=0) {
	    $surf=3;
	} elsif (($node-13001)*(13961-$node)>=0) {
	    $surf=4;
	} elsif (($node-14001)*(14961-$node)>=0) {
	    $surf=5;
	    # the numbers that follow reflect the new node configuration
	} elsif (($node-100000)*(102988-$node)>=0) {
	    $surf=6;
	} elsif (($node-105000)*(107988-$node)>=0) {
	    $surf=7;
	} elsif (($node-110000)*(111368-$node)>=0) {
	    $surf=8;
	} elsif (($node-115000)*(116368-$node)>=0) {
	    $surf=9;
	} else {
	    # dont use this node. 
	    next;
	}
	

	next if (sqrt($x*$x+$y*$y) > $ro{$surf});

	push(@{$nodelist[$surf]},$node);
	if ($raytrace_coords==1) {
	    $nodecoords{$node}=join(' ',$y,$x,$z);
	} else {
	    $nodecoords{$node}=join(' ',$x,$y,$z);
	}
    }
    close(F);

    %files=("X" => "L1-L2-Deformations/GX-\*.csv",
	    "Y" => "L1-L2-Deformations/GY-\*.csv",
	    "Z" => "L1-L2-Deformations/GZ-\*.csv");

    %distortions=("X" => {},
		  "Y" => {},
		  "Z" => {});

    foreach $accel ("X","Y","Z") {
	$distort=$distortions{$accel};
	foreach $file (`ls $files{$accel}`) {
	    open(F,$file) || die;
	    while ($line=<F>) {
		chomp($line);
		@this=split(',',$line);
		next if (($#this != 7) || ($this[0] !~ /^[0-9]/));
		($node,$tx,$ty,$tz,$rx,$ry,$rz)=@this[0,2,3,4,5,6,7];
		if ($raytrace_coords==1) {
		    ${$distort}{$node}=join(' ',$ty,$tx,-$tz,$ry,$rx,-$rz);
		} else {
		    ${$distort}{$node}=join(' ',$tx,$ty,$tz,$rx,$ry,$rz);
		}
	    }
	    close(F);
	}
    }

    open(G,">node_coordinates.tnt") || die;
    printf G "node coords\n%s\n",join("\t","surf","x","y","z");

    foreach $surf ( @surfs ) {
	foreach $node (@{$nodelist[$surf]}) {
	    printf G "%s\n",join(' ',$surf,$nodecoords{$node});
	}
    }
    close(G);

    open(G,">GX_node_distortion_${rt_string}_.tnt") || die;
    open(H,">GY_node_distortion_${rt_string}_.tnt") || die;
    open(J,">GZ_node_distortion_${rt_string}_.tnt") || die;

    printf G "node distortion\n%s\n",join("\t","surf","x","y","z","tx","ty","tz","rx","ry","rz");
    printf H "node distortion\n%s\n",join("\t","surf","x","y","z","tx","ty","tz","rx","ry","rz");
    printf J "node distortion\n%s\n",join("\t","surf","x","y","z","tx","ty","tz","rx","ry","rz");

    foreach $surf (@surfs) {
	foreach $node (@{$nodelist[$surf]}) {
	    printf "NODE $node\n";
	    printf G "%s\n",join(' ',$surf,$nodecoords{$node},${$distortions{"X"}}{$node});
	    printf H "%s\n",join(' ',$surf,$nodecoords{$node},${$distortions{"Y"}}{$node});
	    printf J "%s\n",join(' ',$surf,$nodecoords{$node},${$distortions{"Z"}}{$node});
	}
    }
    close(G);
    close(H);
    close(J);

# now write out the mean translations and rotations that are contained in the
# distortion files.

    %lintran=();
    %stx=();%sty=();%stz=();
    %srx=();%sry=();%srz=();
    %n=();

    foreach $surf (@surfs) {
	foreach $gv ("X","Y","Z") {
	    foreach $node (@{$nodelist[$surf]}) {
		($tx,$ty,$tz,$rx,$ry,$rz)=split(' ',${$distortions{$gv}}{$node});
		foreach $var (\$stx{$optic{$surf},$gv},\$sty{$optic{$surf},$gv},
			      \$stz{$optic{$surf},$gv},\$srx{$optic{$surf},$gv},
			      \$sry{$optic{$surf},$gv},\$srz{$optic{$surf},$gv},
			      \$n{$optic{$surf},$gv}) {
		    ${$var}=0.0 if (!defined(${$var}));
		}
		$stx{$optic{$surf},$gv}+=$tx;
		$sty{$optic{$surf},$gv}+=$ty;
		$stz{$optic{$surf},$gv}+=$tz;
		$srx{$optic{$surf},$gv}+=$rx;
		$sry{$optic{$surf},$gv}+=$ry;
		$srz{$optic{$surf},$gv}+=$rz;
		$n{$optic{$surf},$gv}++;
	    }
	}
    }

    foreach $optic ( values %optic ) {
	$optics{$optic}=1;
    }

    foreach $optic ( keys %optics ) {
	foreach $gv ("X","Y","Z") {
	    next if (!defined($n{$optic,$gv}) || $n{$optic,$gv}==0);
	    $lintran{$optic,$gv}=
		sprintf("-t %6.4f %6.4f %6.4f",
			$stx{$optic,$gv}/$n{$optic,$gv}*1e3,
			$sty{$optic,$gv}/$n{$optic,$gv}*1e3,
			$stz{$optic,$gv}/$n{$optic,$gv}*1e3);
	    $rottran{$optic,$gv}=
		sprintf("-r %g %g %g",
			$srx{$optic,$gv}/$n{$optic,$gv},
			$sry{$optic,$gv}/$n{$optic,$gv},
			$srz{$optic,$gv}/$n{$optic,$gv});
	    printf "lintran $lintran{$optic,$gv}\n";
	    $outfile=join('_',$optic,"G".$gv,$rt_string,"solidbody.dat");
	    open(F,">$outfile") || die;
	    printf F "$lintran{$optic,$gv}";
	    close(F);
	    $outfile=join('_',$optic,"G".$gv,$rt_string,"full","solidbody.dat");
	    open(F,">$outfile") || die;
	    printf F "$lintran{$optic,$gv} $rottran{$optic,$gv}";
	    close(F);
	}
    }

# now for a series of files that might be used in fitting or in regression
# analysis.
# <x> <y> <deltaZ> where <x>,<y> have been divided through by max(sqrt(x^2+y^2))
# one for each surface.

    foreach $surf (@surfs) {
	foreach $node (@{$nodelist[$surf]}) {
	    $zern{$node}=[];
	}
    }

    foreach $ax ("X","Y","Z") {
	foreach $surf (@surfs) {
	    $outfile=join('_',
			  $optic{$surf},"S".$surface{$surf},
			  "G".$ax,$rt_string.".tnt");
	    open(G,">$outfile") || die;
#	printf G "$outfile\n%s\n",join("\t","x/rm","y/rm","dz");
	    printf G "$outfile\n%s\n",join("\t","x/rm","y/rm","dx","dy","dz","z1","z2","z3","z4","z5","z6","z7","z8","z9");
#	printf G "$outfile\n%s\n",join("\t","x/rm","y/rm","dx","dy","dz");
	    foreach $node (@{$nodelist[$surf]}) {
		@coords=split(' ',$nodecoords{$node});
		($x,$y)=@coords[0,1];
		# in the rare case where either $x or $y is "#NAME?" then 
		# the following 2 lines converts it to zero (inspection of
		# M2 surface 2 seems to handle this one OK, there is one instance
		# of "#NAME?" there..
		$x /= $ro{$surf};
		$y /= $ro{$surf};
		$norm_xy{$node}=join(' ',$x,$y);
		@dists=split(' ',${$distortions{$ax}}{$node});
		$dz=$dists[2];
		$deform_z{$node}=$dists[2];
		# %zerns are organized by $j index and $node
		@{$zern{$node}}=zernike($x,$y,0 .. ($nzern-1));
#	    printf "zern: %s\n",join(' ',@{$zern{$node}});
#	    printf G "%s\n",join(' ',$x,$y,$dz);
		printf G "%s\n",join(' ',$x,$y,@dists[0,1,2],@{$zern{$node}}[1,2,3,4,5,6,7,8,9]);
#	    printf G "%s\n",join(' ',$x,$y,@dists[0,1,2]);
	    }
	    close(G);
	    # here do some decomposition of the $deform_z{$node} using the basis
	    # @{$zern{$node}}[0..($nzern-1)]
	    # setup the initial guesses by doing an unweighted sum:
#	printf STDERR "number of nodes: %d\n",$#{$nodelist[$surf]};
#	printf STDERR "number of zernikes computed: %d\n",$nzern;
	    @guess=();
	    @amp=();
	    foreach $node (@{$nodelist[$surf]}) {
		$dz=$deform_z{$node};
		for ($ix=0;$ix<$nzern;$ix++) {
		    $amp[$ix]=0 if (!defined($amp[$ix]));
		    $amp[$ix] += $dz*${$zern{$node}}[$ix];
		}
	    }
	    @tol=();
	    for ($ix=0;$ix<$nzern;$ix++) {
		{
		    my ($n,$m);
		    $n=ceil((-3+sqrt(9+8*$ix))/2);
		    $m=2*$ix-$n*($n+2);
		    $norm=sqrt(2*($n+1))/(($m==2)?2:1);
		}
		if ($ix<=2) {
		    $guess[$ix] = $amp[$ix]/($norm*($#{$nodelist[$surf]}-$[+1));
		} else {
		    $guess[$ix] = $amp[$ix]/($norm*($#{$nodelist[$surf]}-$[+1));
#		    $guess[$ix] = 0;
		}
	    }

	    for ($ix=0;$ix<$nzern;$ix++) {
		if (defined($guess[$ix])) {
		    push(@tol,$guess[$ix]);
		} else {
		    push(@tol,$amp[$ix]/($#{$nodelist[$surf]}-$[+1));
		}
	    }

	    # the figure of merit function will run through the node list
	    # @{$nodelist[$surf]};
	    $thesurf=$surf;
	    $ziter=0;

	    for ($iteration=0;$iteration<=5;$iteration++) {

		if ($iteration==0) {
		    $upper_ix=1+2+3+4-1;
		} elsif ($iteration==1) {
		    $upper_ix=1+2+3+4+5-1;
		} elsif ($iteration==2) {
		    $upper_ix=1+2+3+4+5+6-1;
		} elsif ($iteration==3) {
		    $upper_ix=1+2+3+4+5+6+7-1;
		} elsif ($iteration==4) {
		    $upper_ix=1+2+3+4+5+6+7+8-1;
		} elsif ($iteration==5) {
		    $upper_ix=1+2+3+4+5+6+7+8+9-1;
		} # can add more, just be sure to set $nzern above appropriately for sampling the zernike expansion.

		printf "\na new iteration ($upper_ix)..\n";
		
		@g=@guess[0..$upper_ix];
		@t=@tol[0..$upper_ix];
		($p,$y)=MinimiseND(\@g,\@t,\&zernikefit,1e-2,10000);
		printf "\n";
		# store the results into the "guess" in case this will repeat
		# over a larger number of coefficients.
		@guess[0..$#{$p}]=@{$p};
	    }
	    # do one last time to compute $model{$node}
	    $y=zernikefit(@{$p});
	    # and output the solution as a zernike expansion:
	    {
		$outfile=join('_',
			      $optic{$surf},"S".$surface{$surf},"G".$ax,
			      $rt_string,"fit.dat");
#	    $outfile=$optic{$surf}."_S".$surface{$surf}."_G".$ax."_fit.dat";
		my @z_exp=();
		my $zrn;
		foreach $zrn (@{$p}) {
		    push(@z_exp,$zrn*1e3); # report in mm, not m
		}
		open(G,">$outfile") || die;
		printf G "-Z %s\n",join(':',@z_exp);
		close(G);
	    }
#	printf STDERR "FOM(%s)=%g\n",join(',',@{$p}),$y;
	    # now write out a tnt file comparing model to data
#	$outfile=join('_',$optic{$surf}."_S".$surface{$surf}."_G".$ax."_fit.tnt";
	    $outfile=join('_',
			  $optic{$surf},"S".$surface{$surf},
			  "G".$ax,$rt_string,"fit.tnt");
	    open(G,">$outfile") || die;
	    printf G ("$outfile\n%s\n",
		      join("\t",
			   "x/rm","y/rm",
			   "dz","model","residual","z3+"));
	    foreach $node (@{$nodelist[$surf]}) {
		printf G ("%g %g %g %g %g %g\n",
			  split(' ',$norm_xy{$node}),
			  $deform_z{$node},
			  $model{$node},
			  $deform_z{$node}-$model{$node},
			  $z3plusmodel{$node});
	    }
	    close(G)
	}
    }
}

sub zernikefit {
    local(@zernike_amps)=@_;
    # use $thesurf
    my $zernikefit=0;
    foreach my $node (@{$nodelist[$thesurf]}) {
	my $i=-1;
	$model{$node}=0;
	$z3plusmodel{$node}=0;
	foreach my $za (@zernike_amps) {
	    $i++;
	    $model{$node} += ($za*${$zern{$node}}[$i]);
	    $z3plusmodel{$node} += ($za*${$zern{$node}}[$i]) if ($i>2);
	}
	$zernikefit += pow(($deform_z{$node}-$model{$node})/1e-6,2);
    }
    $zernikefit /= ($#{$nodelist[$thesurf]}-$[+1) ; # divide by number of nodes.. 
    printf "\riter=%d zernikefit=$zernikefit um^2           ",$ziter;
    $ziter++;
    return($zernikefit);
}

sub zernike {
    # returns z_{$j}(p,phi) where p = sqrt(x^2+y^2)/R, phi=atan2(y,x)
    local($x,$y,@js)=@_;
    my @zs=();
    my $phi=atan2($y,$x);
    my $rho=sqrt($x*$x+$y*$y);
    foreach $j (@js) {
	my $n=ceil((-3+sqrt(9+8*$j))/2);
	my $m=2*$j-$n*($n+2);
	my $s_max=($n-sqrt($m*$m))/2;
	my ($s,$Rnm,$Nnm,$amp);
	$Nnm=sqrt(2*($n+1)/(($m==0)?2:1));
	$Rnm=0;
	for ($s=0;$s<=$s_max;$s++) {
	    $amp=(pow(-1,$s)*factrl($n-$s))/
		(factrl($s)*factrl(abs($m)+$s_max-$s)*factrl($s_max-$s));
	    $Rnm+=$amp*pow($rho,$n-2*$s);
	}
	if ($m>=0) {
	    push(@zs,$Nnm*$Rnm*cos($m*$phi));
	} else {
	    push(@zs,-$Nnm*$Rnm*sin($m*$phi));
	}
    }
    return(@zs);
}

sub factrl {
    local($n)=@_;
    if ($n>1) {
	return($n*factrl($n-1));
    } else {
	return(1);
    }
}
