#!/usr/bin/perl -w
use warnings;
use strict;
use POSIX;
use Math::Amoeba qw(MinimiseND);
use Spreadsheet::Read;
use Text::CSV_XS;
use Scalar::Util qw(looks_like_number);

my $update_120911=1;
my %model=();
my %z3plusmodel=();
my $thesurf;
my @nodelist;
my %distortions=();
my %zern;
my $ziter;
my %deform_z;
my @surfs;
my %optic=();
my $makedirname;

# this script is for converting J.Ku's updated FEA of 120620
# this script is for converting J.Ku's updated FEA of 120911:
# re-using surface indices 6 thru 9 (J.Ku re-used node labels 100000 thru 116368)

# unpack the zip file
my $ThisScript="extract_displacements.perl";

my $DISTORTION_DB_PATH="/home/arasmus/camera_perturbations/DB"; # this can be gotten from an environmnetal variable

my $DBKEY="LSST_PERTURBATIONS_DB";

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

my ($nodeid_array,$x_array,$y_array,$z_array);

# EDIT FOLLOWING LINES APPROPRIATELY
my %jobs=("gravity" => [("New-L1-L2-Figure-Calcs-r1.zip")],
	  "thermal" => [("Therm-Figure-R2.zip")]);

$jobs{"gravity"}->[0] = "L1-L2\*Deformations-2012-Sep-11.zip" if ($update_120911 == 1);

foreach my $task ("gravity","thermal") {
# foreach my $task ("gravity") {
    foreach my $db ( @{$jobs{$task}} ) {

	printf STDERR "db=$db\n";
	my @zipfile;
	@zipfile=glob(join('/',$DISTORTION_DB_PATH,"$db"));

	if (scalar(@zipfile)>1) {
	    # source zip file is ambiguous. use the first one
	    printf STDERR "more than one zipfile is available: %s\nambiguous.\n",join(',',@zipfile);
	    printf STDERR "proceeding with the first one (%s), giving it benefit of doubt\n",$zipfile[0];
	} elsif (scalar(@zipfile)<1) {
	    printf STDERR "no matching zipfiles expected (e.g., %s). terminating..\n",join('/',$DISTORTION_DB_PATH,"L1-L2*zip");
	}
	
# proceed. make a new directory matching zipfile name for dumping contents
	
	my @makedirname=split('/',$zipfile[0]);
	$makedirname[$#makedirname] =~ s/.zip//g;
	$makedirname[$#makedirname] =~ s/\ /-/g;
	$makedirname=join('/',$DISTORTION_DB_PATH,$makedirname[$#makedirname]);
	mkdir $makedirname if (! -d $makedirname);
	`cd "$makedirname";unzip -o "$zipfile[0]"`;

#    print STDERR `ls -lR "$makedirname"`."\n";
# now read in node coordinates
	
	my @xlsfile;
	$makedirname =~ s/\ /*/g;
	if ($db =~ /L1-L2/) {
# EDIT FOLLOWING LINES APPROPRIATELY
	    @xlsfile=glob "$makedirname/Point*csv";
	    @xlsfile=glob "$makedirname/*xlsx" if ($update_120911==1);
	    
	    printf "xlsfile = %s\n",join(' ',@xlsfile);
	    my $xls=ReadData($xlsfile[0]) || die "can't open file $xlsfile[0]\n";
	    
# EDIT FOLLOWING LINES APPROPRIATELY
	    if ($update_120911==1) {
		($nodeid_array,$x_array,$y_array,$z_array)=
		    ([ @{$xls->[1]{cell}[2]}[2..11357] ],
		     [ @{$xls->[1]{cell}[4]}[2..11357] ],
		     [ @{$xls->[1]{cell}[5]}[2..11357] ],
		     [ @{$xls->[1]{cell}[6]}[2..11357] ]);
	    } else {
		($nodeid_array,$x_array,$y_array,$z_array)=
		    ([ @{$xls->[1]{cell}[2]}[2..11285] ],
		     [ @{$xls->[1]{cell}[4]}[2..11285] ],
		     [ @{$xls->[1]{cell}[5]}[2..11285] ],
		     [ @{$xls->[1]{cell}[6]}[2..11285] ]);
	    }
	} elsif ($db =~ /Displacements/) {
	    @xlsfile=glob "$makedirname/*/Grid*";
	    
	    printf "xlsfile = %s\n",join(' ',@xlsfile);
	    my $xls=ReadData($xlsfile[0]) || die "can't open file $xlsfile[0]\n";
	    
	    ($nodeid_array,$x_array,$y_array,$z_array)=
		([ @{$xls->[1]{cell}[2]}[1..4446] ],
		 [ @{$xls->[1]{cell}[4]}[1..4446] ],
		 [ @{$xls->[1]{cell}[5]}[1..4446] ],
		 [ @{$xls->[1]{cell}[6]}[1..4446] ]);
	} elsif ($task eq "thermal") {
	    # grid files came with the gravity package, not with thermal..
	    # do nothing
	} else {
	    die "task not recognized: $task";
	}

# %optic is only for naming output files..
	%optic=("0" => "L3",
		"1" => "L3",
		"2" => "L2",
		"3" => "L2",
		"4" => "L1",
		"5" => "L1",
		"6" => "L2_120620",
		"7" => "L2_120620",
		"8" => "L1_120620",
		"9" => "L1_120620");
# update
	@optic{"6","7","8","9"}=("L1_new","L1_new","L2_new","L2_new")
	    if ($update_120911 == 1);

	my %surface=("0" => 2,
		     "1" => 1,
		     "2" => 2,
		     "3" => 1,
		     "4" => 2,
		     "5" => 1,
		     "6" => 2,
		     "7" => 1,
		     "8" => 2,
		     "9" => 1);

	@surface{"6","7","8","9"}=(1,2,1,2)  if ($update_120911 == 1);

# this $nzern is used for evaluating zernike fields used later in fitting. large values for now..
	my $nzern=(1+2+3+4+5+6+7+8+9);
# (28 zernikes)
#$nzern=(1+2+3+4+5+6+7+8+9);

# EDIT FOLLOWING LINES APPROPRIATELY
	@surfs=(6,7,8,9);

# rather than work out the zernike mixing matrix for coordinate transformations
# just specify whether this is for raytrace coordinates (or CCS/MCS).
# the tranformation between the two systems goes as:
# /x\   / 0 1  0 \ /x\
# |y| = | 1 0  0 | |y|
# \z/rt \ 0 0 -1 / \z/ccs
#

#	foreach my $raytrace_coords ( 0,1 ) {
	foreach my $raytrace_coords ( 0 ) {

	    my %sbtran=();
	    my %nodecoords=();
	    
	    my $rt_string=($raytrace_coords==1)?"rt_coords":"ccs_coords";

	    my %ro=(0=>0.361,
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
	    @ro{6,7,8,9}=(0.775,0.775,0.551,0.551) if ($update_120911 == 1);

	    my $nskip=0;

	    # clear out the nodelist if necessary
	    for my $surf (@surfs) {
		$nodelist[$surf]=[] if (!defined($nodelist[$surf]));
		@{$nodelist[$surf]}=();
	    }

	    for (my $i=0;$i<=$#{$nodeid_array};$i++) {

		my ($node,$x,$y,$z)=($nodeid_array->[$i],
				     $x_array->[$i],
				     $y_array->[$i],
				     $z_array->[$i]);

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

		my $surf;

# EDIT FOLLOWING LINES APPROPRIATELY
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
		    # the numbers that follow reflect a partitioning that works for both
		    # r1 (120620) and r2 (120911) grid node definitions #
		    # (nb. different surface IDs) - and there is no corresponding 
		    # deformation data provided on the "spoke" grid point array 
		    # described in 120911 documentation!
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
		    $nodecoords{$node}=join(' ',$y,$x,-$z);
		} else {
		    $nodecoords{$node}=join(' ',$x,$y,$z);
		}
	    }
	    close(F);

	    my %files;
	    my @conditions;
	    if ($db =~ /L1-L2/) {
		@conditions=("X","Y","Z");
		if ($update_120911==1) {
		    @files{@conditions}=("$makedirname/\*/GX-\*.csv",
					 "$makedirname/\*/GY-\*.csv",
					 "$makedirname/\*/GZ-\*.csv");
		} else {
		    @files{@conditions}=("$makedirname/GX-\*.csv",
					 "$makedirname/GY-\*.csv",
					 "$makedirname/GZ-\*.csv");
		}
	    } elsif ($db =~ /Displacements/) {
		@conditions=("X","Y","Z");
		@files{@conditions}=("$makedirname/*/X-\*.csv",
				     "$makedirname/*/Y-\*.csv",
				     "$makedirname/*/Z-\*.csv");
	    } elsif ($task eq "thermal") {
		@conditions=("Cold","Hot","Nom");
		@files{@conditions}=("$makedirname/*/Cold-\*.csv",
				     "$makedirname/*/Hot-\*.csv",
				     "$makedirname/*/Nom-\*.csv");
	    } else {
		die "don't recognize task $task.\n";
	    }
	    
	    foreach my $cond (@conditions) {
		$distortions{$cond}={} if (!defined($distortions{$cond}));
		my $distort=$distortions{$cond};
		foreach my $file (glob $files{$cond}) {
		    open(F,$file) || die "can't get file $file\n";
		    while (my $line=<F>) {
			chomp($line);
			my @this=split(',',$line);
			next if (($#this != 7) || ($this[0] !~ /^[0-9]/));
			my ($node,$tx,$ty,$tz,$rx,$ry,$rz)=@this[0,2,3,4,5,6,7];
			if ($raytrace_coords==1) {
			    $distort->{$node}=join(' ',$ty,$tx,-$tz,$ry,$rx,-$rz);
			} else {
			    $distort->{$node}=join(' ',$tx,$ty,$tz,$rx,$ry,$rz);
			}
		    }
		    close(F);
		}
	    }

	    $makedirname =~ s/\*/\ /g;

	    printf STDERR "trying to deposit into $makedirname\n";
	    open(G,">$makedirname/node_coordinates.tnt") || die;
	    printf G "node coords\n%s\n",join("\t","surf","x","y","z");

	    foreach my $surf ( @surfs ) {
		foreach my $node (@{$nodelist[$surf]}) {
		    printf G "%s\n",join(' ',$surf,$nodecoords{$node});
		}
	    }
	    close(G);
	    
	    if ($task eq "gravity") {
		open(G,">$makedirname/GX_node_distortion_${rt_string}_.tnt") || die;
		open(H,">$makedirname/GY_node_distortion_${rt_string}_.tnt") || die;
		open(J,">$makedirname/GZ_node_distortion_${rt_string}_.tnt") || die;
	    } elsif ($task eq "thermal") {
		open(G,">$makedirname/Cold_node_distortion_${rt_string}_.tnt") || die;
		open(H,">$makedirname/Hot_node_distortion_${rt_string}_.tnt") || die;
		open(J,">$makedirname/Nom_node_distortion_${rt_string}_.tnt") || die;
	    }

	    printf G "node distortion\n%s\n",
	    join("\t","surf","x","y","z","tx","ty","tz","rx","ry","rz");
	    printf H "node distortion\n%s\n",
	    join("\t","surf","x","y","z","tx","ty","tz","rx","ry","rz");
	    printf J "node distortion\n%s\n",
	    join("\t","surf","x","y","z","tx","ty","tz","rx","ry","rz");

	    foreach my $surf (@surfs) {
		foreach my $node (@{$nodelist[$surf]}) {
		    if ($task eq "gravity") {

			printf G "%s\n",join(' ',$surf,$nodecoords{$node},$distortions{"X"}->{$node})
			    if (defined($distortions{"X"}->{$node}));
			printf H "%s\n",join(' ',$surf,$nodecoords{$node},$distortions{"Y"}->{$node})
			    if (defined($distortions{"Y"}->{$node}));
			printf J "%s\n",join(' ',$surf,$nodecoords{$node},$distortions{"Z"}->{$node})
			    if (defined($distortions{"Z"}->{$node}));
		    } elsif ($task eq "thermal") {
			printf G "%s\n",join(' ',$surf,$nodecoords{$node},$distortions{"Cold"}->{$node})
			    if (defined($distortions{"Cold"}->{$node}));
			printf H "%s\n",join(' ',$surf,$nodecoords{$node},$distortions{"Hot"}->{$node})
			    if (defined($distortions{"Hot"}->{$node}));
			printf J "%s\n",join(' ',$surf,$nodecoords{$node},$distortions{"Nom"}->{$node})
			    if (defined($distortions{"Nom"}->{$node}));
		    } else {
			die "task $task not recognized.\n";
		    }
		}
	    }
	    close(G);
	    close(H);
	    close(J);

# now write out the mean translations and rotations that are contained in the
# distortion files.

	    my %lintran=();
	    my %stx=();my %sty=();my %stz=();
	    my %srx=();my %sry=();my %srz=();
	    my %n=();

	    foreach my $surf (@surfs) {
		foreach my $cond (@conditions) {

		    undef $stx{$surf,$cond};
		    undef $sty{$surf,$cond};
		    undef $stz{$surf,$cond};
		    undef $srx{$surf,$cond};
		    undef $sry{$surf,$cond};
		    undef $srz{$surf,$cond};
		    undef $n{$surf,$cond};

		    undef $stx{$optic{$surf},$cond};
		    undef $sty{$optic{$surf},$cond};
		    undef $stz{$optic{$surf},$cond};
		    undef $srx{$optic{$surf},$cond};
		    undef $sry{$optic{$surf},$cond};
		    undef $srz{$optic{$surf},$cond};
		    undef $n{$optic{$surf},$cond};

		    next if (!defined($nodelist[$surf]) || !defined($distortions{$cond}));
		    foreach my $node (@{$nodelist[$surf]}) {
			my ($tx,$ty,$tz,$rx,$ry,$rz)=
			    split(' ',$distortions{$cond}->{$node});
			foreach my $var (\$stx{$optic{$surf},$cond},
					 \$sty{$optic{$surf},$cond},
					 \$stz{$optic{$surf},$cond},
					 \$srx{$optic{$surf},$cond},
					 \$sry{$optic{$surf},$cond},
					 \$srz{$optic{$surf},$cond},
					 \$n{$optic{$surf},$cond},
					 \$stx{$surf,$cond},
					 \$sty{$surf,$cond},
					 \$stz{$surf,$cond},
					 \$srx{$surf,$cond},
					 \$sry{$surf,$cond},
					 \$srz{$surf,$cond},
					 \$n{$surf,$cond}) {
			    ${$var}=[] if (!defined(${$var}));
			}
			push($stx{$optic{$surf},$cond},$tx);
			push($sty{$optic{$surf},$cond},$ty);
			push($stz{$optic{$surf},$cond},$tz);
			push($srx{$optic{$surf},$cond},$rx);
			push($sry{$optic{$surf},$cond},$ry);
			push($srz{$optic{$surf},$cond},$rz);
			push(  $n{$optic{$surf},$cond},  1);

			push($stx{$surf,$cond},$tx);
			push($sty{$surf,$cond},$ty);
			push($stz{$surf,$cond},$tz);
			push($srx{$surf,$cond},$rx);
			push($sry{$surf,$cond},$ry);
			push($srz{$surf,$cond},$rz);
			push(  $n{$surf,$cond},  1);
		    }
		}
	    }

	    my %optics=();
	    my %rottran=();

	    foreach my $optic ( @optic{@surfs} ) {
		$optics{$optic}=1;
	    }

	    foreach my $optic ( keys %optics ) {
		printf STDERR "doing optic $optic\n";
		foreach my $cond (@conditions) {

		    next if (!defined($n{$optic,$cond}) || $n{$optic,$cond}==0);
		    $sbtran{$optic,$cond,"tx"}=sum($stx{$optic,$cond})/sum($n{$optic,$cond});
		    $sbtran{$optic,$cond,"ty"}=sum($sty{$optic,$cond})/sum($n{$optic,$cond});
		    $sbtran{$optic,$cond,"tz"}=sum($stz{$optic,$cond})/sum($n{$optic,$cond});
		    $sbtran{$optic,$cond,"rx"}=sum($srx{$optic,$cond})/sum($n{$optic,$cond});
		    $sbtran{$optic,$cond,"ry"}=sum($sry{$optic,$cond})/sum($n{$optic,$cond});
		    $sbtran{$optic,$cond,"rz"}=sum($srz{$optic,$cond})/sum($n{$optic,$cond});

		    $lintran{$optic,$cond}=
			sprintf("-t %6.4f %6.4f %6.4f",
				$sbtran{$optic,$cond,"tx"}*1e3,
				$sbtran{$optic,$cond,"ty"}*1e3,
				$sbtran{$optic,$cond,"tz"}*1e3);
		    $rottran{$optic,$cond}=
			sprintf("-r %g %g %g",
				$sbtran{$optic,$cond,"rx"},
				$sbtran{$optic,$cond,"ry"},
				$sbtran{$optic,$cond,"rz"});
#	    printf "lintran $lintran{$optic,$cond}\n";
#	    printf "rottran $rottran{$optic,$cond}\n";
		    my $outfile;
		    if ($task eq "gravity") {
			$outfile=join('_',$optic,"G".$cond,$rt_string,"solidbody.dat");
		    } elsif ($task eq "thermal") {
			$outfile=join('_',$optic,$cond,$rt_string,"solidbody.dat");
		    }
		    open(F,">$makedirname/$outfile") || die;
		    printf F "$lintran{$optic,$cond}";
		    close(F);
		    if ($task eq "gravity") {
			$outfile=join('_',$optic,"G".$cond,$rt_string,"full","solidbody.dat");
		    } elsif ($task eq "thermal") {
			$outfile=join('_',$optic,$cond,$rt_string,"full","solidbody.dat");
		    }
		    open(F,">$makedirname/$outfile") || die;
		    printf F "$lintran{$optic,$cond} $rottran{$optic,$cond}";
		    close(F);
		    printf STDERR "raytrace_coords $raytrace_coords optic $optic cond $cond $lintran{$optic,$cond} $rottran{$optic,$cond}\n";
		}
	    }

	    # also output figure maps that take into account the solidbody transform


	    foreach my $surf (@surfs) {
		foreach my $cond (@conditions) {
#		next if (!defined($sbtran{$optic,$cond,"tz"}));
		    $sbtran{$surf,$cond,"tx"}=
			sum($stx{$surf,$cond})/sum($n{$surf,$cond});
		    $sbtran{$surf,$cond,"ty"}=
			sum($sty{$surf,$cond})/sum($n{$surf,$cond});
		    $sbtran{$surf,$cond,"tz"}=
			sum($stz{$surf,$cond})/sum($n{$surf,$cond});
		    $sbtran{$surf,$cond,"rx"}=
			sum($srx{$surf,$cond})/sum($n{$surf,$cond});
		    $sbtran{$surf,$cond,"ry"}=
			sum($sry{$surf,$cond})/sum($n{$surf,$cond});
		    $sbtran{$surf,$cond,"rz"}=
			sum($srz{$surf,$cond})/sum($n{$surf,$cond});

		    my ($a,$b,$c)=(-$sbtran{$surf,$cond,"ry"},
				   +$sbtran{$surf,$cond,"rx"},
				   +$sbtran{$surf,$cond,"tz"});

		    my $outfile;
		    if ($task eq "gravity") {
			$outfile=join('_',
				      $optic{$surf},"S".$surface{$surf},
				      "G".$cond,$rt_string,"perturb_nosb.tnt");
		    } elsif ($task eq "thermal") {
			$outfile=join('_',
				      $optic{$surf},"S".$surface{$surf},
				      $cond,$rt_string,"perturb_nosb.tnt");
		    }
		    open(F,">$makedirname/$outfile") || die;
		    printf F "DAT\n%s\n",join("\t","x","y","dz","sbdz","dz-sbdz");
		    foreach my $node (@{$nodelist[$surf]}) {
			next if (!defined($distortions{$cond}->{$node}));
			my @coords=split(' ',$nodecoords{$node});
			my ($x,$y)=@coords[0,1];
			my @dists=split(' ',$distortions{$cond}->{$node});
			my $dz=$dists[2];
			my $sbdz=$a*$x+$b*$y+$c;
			printf F "%s\n",join(' ',$x,$y,$dz,$sbdz,$dz-$sbdz);
		    }
		    close(F);
		    # output a simple perturbation.tnt file containing perturbations that exclude piston/tip/tilt.

		    if ($task eq "gravity") {
			$outfile=join('_',
				      $optic{$surf},"S".$surface{$surf},
				      "G".$cond,$rt_string,"perturbation.tnt");
		    } elsif ($task eq "thermal") {
			$outfile=join('_',
				      $optic{$surf},"S".$surface{$surf},
				      $cond,$rt_string,"perturbation.tnt");
		    }
		    open(F,">$makedirname/$outfile") || die;
		    printf F "DAT\n%s\n",join("\t","x","y","dz-sbdz");
		    foreach my $node (@{$nodelist[$surf]}) {
			next if (!defined($distortions{$cond}->{$node}));
			my @coords=split(' ',$nodecoords{$node});
			my ($x,$y)=@coords[0,1];
			my @dists=split(' ',$distortions{$cond}->{$node});
			my $dz=$dists[2];
			my $sbdz=$a*$x+$b*$y+$c;
			printf F "%s\n",join(' ',$x,$y,$dz-$sbdz);
		    }
		    close(F);
		}
	    }

# now for a series of files that might be used in fitting or in regression
# analysis.
# <x> <y> <deltaZ> where <x>,<y> have been divided through by max(sqrt(x^2+y^2))
# one for each surface.

	    %zern=();

	    foreach my $surf (@surfs) {
		foreach my $node (@{$nodelist[$surf]}) {
		    $zern{$node}=[];
		}
	    }

	    my %norm_xy=();
	    %deform_z=();

	    foreach my $cond (@conditions) {
		foreach my $surf (@surfs) {
		    my $outfile;
		    if ($task eq "gravity") {
			$outfile=join('_',$optic{$surf},"S".$surface{$surf},"G".$cond,$rt_string.".tnt");
		    } else {
			$outfile=join('_',$optic{$surf},"S".$surface{$surf},$cond,$rt_string.".tnt");
		    }
		    open(G,">$makedirname/$outfile") || die;
#	printf G "$outfile\n%s\n",join("\t","x/rm","y/rm","dz");
		    printf G "$outfile\n%s\n",join("\t","x/rm","y/rm","dx","dy","dz","z1","z2","z3","z4","z5","z6","z7","z8","z9");
#	printf G "$outfile\n%s\n",join("\t","x/rm","y/rm","dx","dy","dz");
		    foreach my $node (@{$nodelist[$surf]}) {
			my @coords=split(' ',$nodecoords{$node});
			my ($x,$y)=@coords[0,1];
			# in the rare case where either $x or $y is "#NAME?" then 
			# the following 2 lines converts it to zero (inspection of
			# M2 surface 2 seems to handle this one OK, there is one instance
			# of "#NAME?" there..
			$x /= $ro{$surf};
			$y /= $ro{$surf};
			$norm_xy{$node}=join(' ',$x,$y);
			my @dists=split(' ',$distortions{$cond}->{$node});
			my $dz=$dists[2];
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
		    my @guess=();
		    my @amp=();
		    my @znorm=();
		    {
			foreach my $node (@{$nodelist[$surf]}) {
			    my $dz=$deform_z{$node};
			    for (my $ix=0;$ix<$nzern;$ix++) {
				$amp[$ix]=0 if (!defined($amp[$ix]));
				$amp[$ix] += $dz*$zern{$node}->[$ix];

				$znorm[$ix]=0 if (!defined($znorm[$ix]));
				$znorm[$ix]+=$zern{$node}->[$ix]*$zern{$node}->[$ix];

			    }
			}
		    }
		    my @tol=();
		    for (my $ix=0;$ix<$nzern;$ix++) {
			my $norm;
			{
			    my ($n,$m);
			    $n=ceil((-3+sqrt(9+8*$ix))/2);
			    $m=2*$ix-$n*($n+2);
			    $norm=sqrt(2*($n+1)/(($m==0)?2:1));
			    $norm=1;
			}
			if ($ix<=2) {
			    $guess[$ix] = $amp[$ix]/($norm*($#{$nodelist[$surf]}-$[+1));
			} else {
			    $guess[$ix] = $amp[$ix]/($norm*($#{$nodelist[$surf]}-$[+1));
#		    $guess[$ix] = 0;
			}
			$znorm[$ix] /= ($#{$nodelist[$surf]}-$[+1);
		    }
#           could close loop here if no fitting is needed.
#	    exit;
		    for (my $ix=0;$ix<$nzern;$ix++) {
			if (defined($guess[$ix])) {
			    push(@tol,0.1*$guess[$ix]);
			} else {
			    push(@tol,0.1*$amp[$ix]/($#{$nodelist[$surf]}-$[+1));
			}
		    }

		    # the figure of merit function will run through the node list
		    # @{$nodelist[$surf]};
		    $thesurf=$surf;
		    $ziter=0;
		    my $upper_ix;
		    my ($p,$y);

		for (my $iteration=0;$iteration<=5;$iteration++) {

			if ($iteration==0) {
			    $upper_ix=1+2+3+4-1;
			} elsif ($iteration==1) {
			    $upper_ix=1+2+3+4+5-1;
			    @guess[1+2+3+4..1+2+3+4+5-1] = (0)x5;
			} elsif ($iteration==2) {
			    $upper_ix=1+2+3+4+5+6-1;
			    @guess[1+2+3+4+5..1+2+3+4+5+6-1] = (0)x6;
			} elsif ($iteration==3) {
			    $upper_ix=1+2+3+4+5+6+7-1;
			    @guess[1+2+3+4+5+6..1+2+3+4+5+6+7-1] = (0)x7;
			} elsif ($iteration==4) {
			    $upper_ix=1+2+3+4+5+6+7+8-1;
			    @guess[1+2+3+4+5+6+7..1+2+3+4+5+6+7+8-1] = (0)x8;
			} elsif ($iteration==5) {
			    $upper_ix=1+2+3+4+5+6+7+8+9-1;
			    @guess[1+2+3+4+5+6+7+8..1+2+3+4+5+6+7+8+9-1] = (0)x9;
			} # can add more, just be sure to set $nzern above appropriately for sampling the zernike expansion.

			printf "\na new iteration ($upper_ix)..\n";
			
			my @g=@guess[0..$upper_ix];
			my @t=@tol[0..$upper_ix];
			($p,$y)=MinimiseND(\@g,\@t,\&zernikefit,1e-3,10000);
			printf "\n";
			# store the results into the "guess" in case this will repeat
			# over a larger number of coefficients.
			@guess[0..$#{$p}]=@{$p};
		    }
		    # do one last time to compute $model{$node}
		    $y=zernikefit(@{$p});
		    # and output the solution as a zernike expansion:
		    {
			if ($task eq "gravity") {
			    $outfile=join('_',$optic{$surf},"S".$surface{$surf},"G".$cond,$rt_string,"fit.dat");
			} else {
			    $outfile=join('_',$optic{$surf},"S".$surface{$surf},$cond,$rt_string,"fit.dat");
			}

			my @z_exp=();
			my $zrn;
			foreach $zrn (@{$p}) {
			    push(@z_exp,$zrn*1e3); # report in mm, not m
			}
			open(G,">$makedirname/$outfile") || die;
			printf G "-Z %s\n",join(':',@z_exp);
			close(G);
		    }
		    if ($task eq "gravity") {
			$outfile=join('_',$optic{$surf},"S".$surface{$surf},"G".$cond,$rt_string,"fit.tnt");
		    } else {
			$outfile=join('_',$optic{$surf},"S".$surface{$surf},$cond,$rt_string,"fit.tnt");
		    }
		    open(G,">$makedirname/$outfile") || die;
		    printf G ("$outfile\n%s\n",
			      join("\t",
				   "x/rm","y/rm",
				   "dz","model","residual","z3+"));
		    foreach my $node (@{$nodelist[$surf]}) {
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
    }
    foreach my $surf (@surfs) {
	my $dir_to_make=join('/',$DISTORTION_DB_PATH,$optic{$surf});
	mkdir $dir_to_make if (! -d $dir_to_make);
	`rm $dir_to_make/$task` if (-l "$dir_to_make/$task");
	`ln -sf $makedirname $dir_to_make/$task`;
    }
}

sub zernikefit {
    my (@zernike_amps)=@_;

    my $zernikefit=0;
    foreach my $node (@{$nodelist[$thesurf]}) {
	my $i=-1;
	$model{$node}=0;
	$z3plusmodel{$node}=0;
	foreach my $za (@zernike_amps) {
	    $i++;
	    $model{$node} += ($za*$zern{$node}->[$i]);
	    $z3plusmodel{$node} += ($za*$zern{$node}->[$i]) if ($i>2);
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

sub sum {
    my ($list)=@_;
    my $sum=0;
    return(0) if (!defined($list) || (scalar($list)==0));
    for my $entry ( sort {pow($a,2)<=>pow($b,2)} @{$list} ) {
	$sum += $entry;
    }
    $sum;
}
