#!/usr/bin/perl
use strict;
use warnings;
use POSIX;

# script to take in named file containing zernike coefficient expansion
# and then output a perturbation.tnt file. The X,Y coordinates are 
# auto-generated.

my @zamp=();

while (<>) {
    chomp;
    my ($dummy,$zlist)=split(' ');
    push(@zamp,split(':',$zlist));
}

my ($x,$y,$nzern)=([],[],$#zamp-$[+1);

# populate $x,$y
my $i;
my $noise=1e-3; # 1mm out of 10mm spacing

for (my $xv=-1.0;$xv<=1.0;$xv+=2e-2) {
    for (my $yv=-1.0;$yv<=1.0;$yv+=2e-2) {
	if (pow($xv,2)+pow($yv,2)<1.0) {
	    push($x,$xv+$noise*(1-2*rand(1)));	
	    push($y,$yv+$noise*(1-2*rand(1)));
	}
    }
}


printf "zern2perturbation DAT\n%s\n",join("\t","rho . x","rho . y","zv");

my $ro=775.0;

for ($i=0;$i<scalar($x);$i++) {
    next if (!defined($x->[$i]) || !defined($y->[$i]));
    my @zvlist=zernike($x->[$i],$y->[$i],0..$nzern-1);
    my $ztot=0;
    for (my $j=0;$j<$nzern-1;$j++) {
	$ztot += $zvlist[$j]*$zamp[$j];
    }
    printf "%f %f %g\n",$ro*$x->[$i],$ro*$y->[$i],$ztot;
}

sub zernike {
    # returns z_{$j}(p,phi) where p = sqrt(x^2+y^2)/R, phi=atan2(y,x)
    my ($x,$y,@js)=@_;
    my @zs=();
    my $phi=atan2($y,$x);
    my $rho=sqrt($x*$x+$y*$y);
    foreach my $j (@js) {
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
    my ($n)=@_;
    if ($n>1) {
	return($n*factrl($n-1));
    } else {
	return(1);
    }
}

