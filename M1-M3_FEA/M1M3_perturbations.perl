#!/usr/bin/perl 
use strict;
use warnings;
use POSIX;
use Math::Amoeba qw(MinimiseND);
use Spreadsheet::Read;
use Getopt::Long;
use Pod::Usage;

my %d=();
my @bm_list=
    ("bm_Tip","bm_Tilt","bm_Power",
     "bm_Astig1","bm_Astig2","bm_Focus3",
     "bm_Tref4","bm_Tref5","bm_Coma6","bm_Coma7",
     "bm_Quad8","bm_Quad9","bm_SecAstig10","bm_SecAstig11","bm_Spherical12",
     "bm_Penta13","bm_Penta14","bm_15","bm_16","bm_17","bm_18",
     "bm_19","bm_20","bm_21","bm_22","bm_23","bm_24","bm_25","bm_26","bm_27");

@d{@bm_list}=(0,0,0,
	      0,0,0,
	      0,0,0,0,
	      0,0,0,0,0,
	      0,0,0,0,0,0,0,
	      0,0,0,0,0,0,0,0);

my ($reanalyse,$Focus3,
    $Tb,$Tx,$Ty,$Tz,$Tr,
    $Tb_ao,$Tx_ao,$Ty_ao,$Tz_ao,$Tr_ao,
    $theta,$precomp_scale,$ao_calib,
    $intermediate_tnt,$output_tnt)=
    (0,1,
     0,0,0,0,0,
     0,0,0,0,0,
     30,1.0,1.0,
     "M1M3_intermediate_FEA.tnt","M1M3_perturbation.tnt");

my $help=0;
my $man=0;

GetOptions('help|?'   => \$help,
	   'man'      => \$man,
	   'reanalyse'=> \$reanalyse,
	   'Focus3=i' => \$Focus3,
	   'Tb=f'    => \$Tb,
	   'Tx=f'    => \$Tx,
	   'Ty=f'    => \$Ty,
	   'Tz=f'    => \$Tz,
	   'Tr=f'    => \$Tr,
	   'Tb_ao=f' => \$Tb_ao,
	   'Tx_ao=f' => \$Tx_ao,
	   'Ty_ao=f' => \$Ty_ao,
	   'Tz_ao=f' => \$Tz_ao,
	   'Tr_ao=f' => \$Tr_ao,
	   'theta=f' => \$theta,
	   'precomp_scale=f' => \$precomp_scale,
	   'ao_calib=f'      => \$ao_calib,
	   ('d_Tip=f'        => \$d{"bm_Tip"},
	    'd_Tilt=f'       => \$d{"bm_Tilt"},
	    'd_Power=f'      => \$d{"bm_Power"},
	    'd_Astig1=f'     => \$d{"bm_Astig1"},
	    'd_Astig2=f'     => \$d{"bm_Astig2"},
	    'd_Focus3=f'     => \$d{"bm_Focus3"},
	    'd_Tref4=f'      => \$d{"bm_Tref4"},
	    'd_Tref5=f'      => \$d{"bm_Tref5"},
	    'd_Coma6=f'      => \$d{"bm_Coma6"},
	    'd_Coma7=f'      => \$d{"bm_Coma7"},
	    'd_Quad8=f'      => \$d{"bm_Quad8"},
	    'd_Quad9=f'      => \$d{"bm_Quad9"},
	    'd_SecAstig10=f' => \$d{"bm_SecAstig10"},
	    'd_SecAstig11=f' => \$d{"bm_SecAstig11"},
	    'd_Spherical12=f'=> \$d{"bm_Spherical12"},
	    'd_Penta13=f'    => \$d{"bm_Penta13"},
	    'd_Penta14=f'    => \$d{"bm_Penta14"},
	    'd_15=f'         => \$d{"bm_15"},
	    'd_16=f'         => \$d{"bm_16"},
	    'd_17=f'         => \$d{"bm_17"},
	    'd_18=f'         => \$d{"bm_18"},
	    'd_19=f'         => \$d{"bm_19"},
	    'd_20=f'         => \$d{"bm_20"},
	    'd_21=f'         => \$d{"bm_21"},
	    'd_22=f'         => \$d{"bm_22"},
	    'd_23=f'         => \$d{"bm_23"},
	    'd_24=f'         => \$d{"bm_24"},
	    'd_25=f'         => \$d{"bm_25"},
	    'd_26=f'         => \$d{"bm_26"},
	    'd_27=f'         => \$d{"bm_27"}),
	   'intermediate_tnt=s'     => \$intermediate_tnt,
	   'output_tnt=s'    => \$output_tnt) || die;;

pod2usage(-exitval=>0 , -verbose=>1) if ($help);
pod2usage(-exitval=>0 , -verbose=>2) if ($man);

my $DBKEY="M1M3_PERTURBATIONS_DB";
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

my @bm_distortion=();
foreach my $bendingmode (@bm_list) {
    my $lbl=$bendingmode;
    $lbl =~ s/bm_/d_/;
    push(@bm_distortion,"$lbl=$d{$bendingmode}") 
	if ($d{$bendingmode} != 0);
}

printf STDERR "%s\n",join("\n",
			  "reanalyse=$reanalyse",
			  "Focus3=$Focus3",
			  "Tb=$Tb",
			  "Tx=$Tx",
			  "Ty=$Ty",
			  "Tz=$Tz",
			  "Tr=$Tr",
			  "Tb_ao=$Tb_ao",
			  "Tx_ao=$Tx_ao",
			  "Ty_ao=$Ty_ao",
			  "Tz_ao=$Tz_ao",
			  "Tr_ao=$Tr_ao",
			  "theta=$theta",
			  "ao_calib=$ao_calib",
			  @bm_distortion,
			  "precomp_scale=$precomp_scale",
			  "intermediate_tnt=$intermediate_tnt",
			  "output_tnt=$output_tnt");

my $file;
my %globals;
my $deg=atan2(1,1)/45.0;

my @actuator_forces=("Actuator","Xloc","Yloc","Gy","Gz","Tb","Tr","Tx","Ty","Tz");

my %Actuators;

my %BendingModes;

my @modes=("node","Xloc","Yloc",
	   "Tip","Tilt","Power",
	   "Astig1","Astig2","Focus3",
	   "Tref4","Tref5","Coma6","Coma7",
	   "Quad8","Quad9","SecAstig10","SecAstig11","Spherical12",
	   "Penta13","Penta14",15,16,17,18,
	   19,20,21,22,23,24,25,26,27);

my %BendingModeActuatorForces;

my @fields=("Actuator","Xloc","Yloc",1..27);

my %rawdata;

my @use_rawfields=("node","Xloc","Yloc","z0","z90",
		   "Tb","Tx","Ty","Tz","Tr");

if ($reanalyse==1) {

    $file=join('/',$ENV{$DBKEY},"Raw_Data_Thermal2.xls");
    my $xls2;
    $xls2=ReadData($file) || die "can't find file $file\n";

    $file=join('/',$ENV{$DBKEY},"Actuator_Force_Summary_2006.xlsx");
    my $xls=ReadData($file) || die "can't find file $file\n";


    my $cix=1;
    foreach my $actuator_force (@actuator_forces) {
	$Actuators{$actuator_force} = [ @{$xls->[1]{cell}[$cix]}[4..163] ];
	$cix++;
    }

    push(@actuator_forces,"R");
    $Actuators{"R"}=[];
    for (my $i=0;$i<=$#{$Actuators{$actuator_forces[0]}};$i++) {
	$Actuators{"R"}->[$i] = sqrt(pow($Actuators{"Xloc"}->[$i],2)+
				     pow($Actuators{"Yloc"}->[$i],2));
    }



    my @rawfields=("surfID","node","Xloc","Yloc","tests","z0_in","z90_in",
		   "dummy","Tb_in","Tx_in","Ty_in","Tz_in","Tr_in",
		   "dummy","z0","z90","dummy",
		   "Tb","Tx","Ty","Tz","Tr");
    my @rawfields_tofit=("z0","z90","Tb","Tx","Ty","Tz","Tr");
    $cix=1;

    foreach my $field (@rawfields) {
	$rawdata{$field}=[ @{$xls2->[1]{cell}[$cix]}[4..5247] ];
	# scale input data as appropriate, to uniform units
	if (($field eq "Xloc") || ($field eq "Yloc")) {
	    for (my $i=0;$i<=$#{$rawdata{$field}};$i++) {
		$rawdata{$field}->[$i] *= 25.4; # convert to mm
	    }
	}
	foreach my $thisfield ( @rawfields_tofit ) {
	    if ($thisfield eq $field){
		for (my $i=0;$i<=$#{$rawdata{$field}};$i++) {
		    $rawdata{$field}->[$i] *= 1e-3; # convert to mm
		}
	    }
	}
	# done with scaling input data
	$cix++;
    }

    push(@rawfields,"R");
    push(@use_rawfields,"R");

    $rawdata{"R"}=[];
    for (my $i=0;$i<=$#{$rawdata{$rawfields[0]}};$i++) {
	$rawdata{"R"}->[$i] = sqrt(pow($rawdata{"Xloc"}->[$i],2)+
				   pow($rawdata{"Yloc"}->[$i],2));
    }

    my @fitmodes=("Tip","Tilt","Power",
		  "Astig1","Astig2","Focus3",
		  "Tref4","Tref5","Coma6","Coma7",
		  "Quad8","Quad9","SecAstig10","SecAstig11","Spherical12",
		  "Penta13","Penta14",15,16,17,18,
		  19,20,21,22,23,24,25,26,27);
    

    $cix=1;
    foreach my $mode (@modes) {
	$BendingModes{$mode}=[ @{$xls->[2]{cell}[$cix]}[4..5247] ];
	# scale input data as appropriate, to uniform units
	if (($mode eq "Xloc") || ($mode eq "Yloc")) {
	    for (my $i=0;$i<=$#{$BendingModes{$mode}};$i++) {
		$BendingModes{$mode}->[$i] *= 25.4; # convert to mm
	    }
	}

	# calibrate the bending modes based on current sampling RMS
	my $bm_norm;
	my @sxx;
	my @sx;
	my @n;
	my ($sxx,$sx,$n)=(0,0,0);

	for (my $i=0;$i<=$#{$BendingModes{$mode}};$i++) {
	    # avoid structural terms in FEA based on surfID field in rawdata:
	    next if ($rawdata{"surfID"}->[$i] eq "PM OD");
	    push(@sxx,pow($BendingModes{$mode}->[$i],2));
	    push(@sx,$BendingModes{$mode}->[$i]);
	    push(@n,1);
	}

	$sxx=sum(sort {$a<=>$b} @sxx);
	$sx =sum(sort {$a<=>$b} @sx);
	$n  =sum(sort {$a<=>$b} @n);
	
	# printf STDERR "for mode $mode:\n";
	# printf STDERR "\tn=%d\n",$n;
	# printf STDERR "\tmean=%g\n",$sx/$n;
	# printf STDERR "\trms=%g\n",sqrt($sxx/$n-pow($sx/$n,2));
	$bm_norm=sqrt($sxx/$n);


	foreach my $thismode ( @fitmodes ) {
	    if ($thismode eq $mode){
		for (my $i=0;$i<=$#{$BendingModes{$mode}};$i++) {
		    $BendingModes{$mode}->[$i] /= $bm_norm; # normalize (nm)
		    $BendingModes{$mode}->[$i] *= 1e-3;     # conv. nm to mm
		}
	    }
	}


	# done with scaling input data
	$cix++;
    }


    push(@modes,"R");
    $BendingModes{"R"}=[];
    for (my $i=0;$i<=$#{$BendingModes{$modes[0]}};$i++) {
	$BendingModes{"R"}->[$i] = sqrt(pow($BendingModes{"Xloc"}->[$i],2)+
					pow($BendingModes{"Yloc"}->[$i],2));
    }



    $cix=1;
    foreach my $field (@fields) {
	$BendingModeActuatorForces{$field}=[ @{$xls->[3]{cell}[$cix]}[4..163] ];
	$cix++;
    }



#my %BendingModeActuatorForces;

    $cix=1;
    foreach my $field (@fields) {
	$BendingModeActuatorForces{$field}=[ @{$xls->[3]{cell}[$cix]}[4..163] ];
	$cix++;
    }

    push(@fields,"R");
    $BendingModeActuatorForces{"R"}=[];
    for (my $i=0;$i<=$#{$BendingModeActuatorForces{$fields[0]}};$i++) {
	$BendingModeActuatorForces{"R"}->[$i] = 
	    sqrt(pow($BendingModeActuatorForces{"Xloc"}->[$i],2)+
		 pow($BendingModeActuatorForces{"Yloc"}->[$i],2));
    }


# spreadsheet now resides in %Actuators, %BendingModes and %BendingModeActuatorForces

# now describe, e.g. $rawdata{z0} in terms of expansion in $BendingModes{@fitmodes}
# and output the residual.

    my %test=%rawdata;
    my %hash;
    my @labels;
    my @test_fields=@use_rawfields;

    if (1) {
	# this part will fit the modes shown in "raw data" using the bending modes.
	open(F,">$ENV{$DBKEY}/$intermediate_tnt") || 
	    die "can't open file ".join('/',$ENV{$DBKEY},$intermediate_tnt);
	open(GG,">$ENV{$DBKEY}/rawdata_bymodes.tnt") || 
	    die "can't open file ".join('/',$ENV{$DBKEY},"rawdata_bymodes.tnt");
	printf GG "rawdata_bymodes\n%s\n",join("\t",@fitmodes);
	my @rawf;
	foreach my $raw_tofit (@rawfields_tofit) {
	    printf "fitting $raw_tofit\n";
	    my $surface=$rawdata{$raw_tofit};

	    # now find amplitudes of the bending modes that will minimize deviation

	    @globals{"shape","modehash","modekeys","iter","FOM"}=
		($surface,\%BendingModes,\@fitmodes,0,0);

	    delete($globals{"guess"}) if (exists($globals{"guess"}));
	    # now populate a $globals{"guess"} to speed up fitting
	    if (1) {
		$globals{"guess"}=[];
		foreach my $bm ( @fitmodes ) {
		    my $amp=0;
		    my $count=0;
		    for (my $i=0;$i<=$#{$surface};$i++) {
			$amp   += ($surface->[$i] * $BendingModes{$bm}->[$i]);
			$count += 1;
		    }
		    # normalize..
		    $amp /= $count;
		    # recall bending modes are normalized to rms = 1e-3
		    $amp /= pow(1e-3,2);
		    push($globals{"guess"},$amp);
		}
	    }

	    my @amplitudes=fit_modes_modebymode();

	    printf GG "! %s: RMS=%f um\n",$raw_tofit,$globals{"FOM"};
	    printf GG "%s\n",join(' ',@amplitudes);

	    my @mod=model(\@amplitudes);

	    push(@rawf,$raw_tofit); # echo the raw data

	    my $label=sprintf("%s_%s",$raw_tofit,"model");
	    push(@rawf,$label);
	    $test{$label} = [ @mod ];

	    $label=sprintf("%s_%s",$raw_tofit,"residuals");
	    push(@rawf,$label);
	    $test{$label} = [];
	    for (my $i=0;$i<=$#mod;$i++) {
		$test{$label}->[$i]=$surface->[$i]-$mod[$i];
	    }

	    if (1) {
		# now that fitting is done, zero out the focus(3) term since that is the 
		# approach described by Doug Neill
		$amplitudes[5]=0; # try that.

		@mod=model(\@amplitudes);

		my $label=sprintf("%s_%s",$raw_tofit,"model_Focus3=0");
		push(@rawf,$label);
		$test{$label} = [ @mod ];

		$label=sprintf("%s_%s",$raw_tofit,"residuals_noFocus3=0");
		push(@rawf,$label);
		$test{$label} = [];
		for (my $i=0;$i<=$#mod;$i++) {
		    $test{$label}->[$i]=$surface->[$i]-$mod[$i];
		}
	    }

	}

	# make a new hash with better labels
	my %output_bendingmodes=();
	my @output_bendingmode_labels=();
	foreach my $key ( @fitmodes ) {
	    my $newkey="bm_".$key;
	    $output_bendingmodes{$newkey}=$BendingModes{$key};
	    push(@output_bendingmode_labels,$newkey);
	}

	%hash=(%test,%output_bendingmodes);
	@labels=("Xloc","Yloc",@rawf,@output_bendingmode_labels);

	printf F "raw/models/residuals/bendingmodes\n%s\n",join("\t",@labels);
	for (my $i=0;$i<=$#{$hash{$labels[0]}};$i++) {
	    next if ($hash{"surfID"}->[$i] eq "PM OD");
	    foreach my $col (@labels) {
		printf F "%g ",$hash{$col}->[$i];
	    }
	    printf F "\n";
	}
	close(F);
	
	exit;
    }

    open(GG,">$ENV{$DBKEY}/test_bycomp.tnt") || 
	die "can't open file ".join('/',$ENV{$DBKEY},"test_bycomp.tnt");
    printf GG "test_bycomp\n%s\n",join("\t","zenith","coszen","sinzen",@fitmodes);

    for (my $zenith=2;$zenith<=90;$zenith+=5) {
	my ($sinzen,$coszen)=(sin($zenith*$deg),
			      cos($zenith*$deg));
	my $newlabel=sprintf("Z%02d",$zenith);
	push(@test_fields,$newlabel);
	$test{$newlabel}=[];
	my $surface=$test{$newlabel};
	for (my $i=0;$i<=$#{$rawdata{"z0"}};$i++) {
	    push($surface,
		 $sinzen*$rawdata{"z90"}->[$i]+($coszen-1)*$rawdata{"z0"}->[$i]);
	}
	# now find amplitudes of the bending modes that will minimize deviation
	@globals{"shape","modehash","modekeys","iter"}=
	    ($surface,\%BendingModes,\@fitmodes,0);

	my @amplitudes=fit_modes();
	printf GG "%s\n",join(' ',$zenith,$coszen,$sinzen,@amplitudes);
	my @mod=model(\@amplitudes);
	
	my $label=sprintf("%s_%s",$newlabel,"model");
	push(@test_fields,$label);
	$test{$label} = [ @mod ];

	$label=sprintf("%s_%s",$newlabel,"residuals");
	push(@test_fields,$label);
	$test{$label} = [];
	for (my $i=0;$i<$#mod;$i++) {
	    $test{$label}->[$i]=$surface->[$i]-$mod[$i];
	}
    }
    close(GG);
    {
	%hash=%test;
	@labels=@test_fields;
	
	open(F,">test.tnt") || die "can't open test.tnt";
	printf F "figures\n%s\n",join("\t",@labels);
	for (my $i=0;$i<=$#{$hash{$labels[0]}};$i++) {
	    foreach my $col (@labels) {
		printf F "%g ",$hash{$col}->[$i];
	    }
	    printf F "\n";
	}
	close(F);
    }
}

# here process commandlines to concoct a working perturbation model.
# suggested arguments: calibration scalar (near 1) to see how much of the
# possible compensation is removed; untracked distortion in gravity, thermal gradients etc.
# this is all done using a linear model, node by node.

# read in data that was collected and packaged by $reanalyse=1 above.

# single data file raw_plus_residuals.tnt

{
    my %m1m3_ao=();
    my $file=join('/',$ENV{$DBKEY},$intermediate_tnt);
    open(F,$file) || die "can't find data file $file\n";
    my @m1m3_fea_ao=<F>;
    close(F);
    chomp(@m1m3_fea_ao);
    shift @m1m3_fea_ao;
    my @keys=split("\t",shift @m1m3_fea_ao);
    foreach my $key (@keys) {
	$m1m3_ao{$key}=[];
    }
    while (my $line=shift(@m1m3_fea_ao)) {
	@_=split(' ',$line);
	my $i=0;
	foreach my $key (@keys) {
	    push($m1m3_ao{$key},$_[$i]);
	    $i++;
	}
    }
    # now output the interpolated values
    my @output_keys=("Xloc","Yloc","perturbation");
#    my $theta=10;
#    my ($Tb,$Tx,$Ty,$Tz,$Tr)=(0,0,0,0,0);
#    my ($Tb_ao,$Tx_ao,$Ty_ao,$Tz_ao,$Tr_ao);
    # $ao_calib is the scalar value by which the correction is dialed in - if $ao_calib is equal to 1, this is the optimal case because residuals are minimized there.
#    my $ao_calib=0.0;
    my ($costheta,$sintheta)=(cos($theta*$deg),sin($theta*$deg));
    # $Tb, $Tx, $Ty, $Tz and $Tr are "true" temperatures & gradients in effect (set above)
    # $Tb_ao, $Tx_ao, $Ty_ao, $Tz_ao and $Tr_ao are temperatures & gradients assumed by the AO system
#    ($Tb_ao,$Tx_ao,$Ty_ao,$Tz_ao,$Tr_ao)=($Tb+0.0,$Tx+0.0,$Ty+0.0,$Tz+0.0,$Tr+0.0);

    open(F,">$output_tnt") || die "can't open output file..\n";
    printf F "interpolated data\n%s\n",join("\t",@output_keys);
    for (my $i=0;$i<=$#{$m1m3_ao{$keys[0]}};$i++) {
	my @output;
	if ($Focus3==0) { # don't correct for Focus3
	    @output=($m1m3_ao{"Xloc"}->[$i],
		     $m1m3_ao{"Yloc"}->[$i],
		     (($costheta-$precomp_scale)*(
			  $m1m3_ao{"z0"}->[$i]-
			  $ao_calib*$m1m3_ao{"z0_model_Focus3=0"}->[$i])+
		      $sintheta*(
			  $m1m3_ao{"z90"}->[$i]-
			  $ao_calib*$m1m3_ao{"z90_model_Focus3=0"}->[$i]))+
		     ($Tb*$m1m3_ao{"Tb"}->[$i]-
		      $Tb_ao*$ao_calib*$m1m3_ao{"Tb_model_Focus3=0"}->[$i])+
		     ($Tx*$m1m3_ao{"Tx"}->[$i]-
		      $Tx_ao*$ao_calib*$m1m3_ao{"Tx_model_Focus3=0"}->[$i])+
		     ($Ty*$m1m3_ao{"Ty"}->[$i]-
		      $Ty_ao*$ao_calib*$m1m3_ao{"Ty_model_Focus3=0"}->[$i])+
		     ($Tz*$m1m3_ao{"Tz"}->[$i]-
		      $Tz_ao*$ao_calib*$m1m3_ao{"Tz_model_Focus3=0"}->[$i])+
		     ($Tr*$m1m3_ao{"Tr"}->[$i]-
		      $Tr_ao*$ao_calib*$m1m3_ao{"Tr_model_Focus3=0"}->[$i]));
	} else {
	    @output=($m1m3_ao{"Xloc"}->[$i],
		     $m1m3_ao{"Yloc"}->[$i],
		     (($costheta-$precomp_scale)*(
			  $m1m3_ao{"z0"}->[$i]-
			  $ao_calib*$m1m3_ao{"z0_model"}->[$i])+
		      $sintheta*(
			  $m1m3_ao{"z90"}->[$i]-
			  $ao_calib*$m1m3_ao{"z90_model"}->[$i]))+
		     ($Tb*$m1m3_ao{"Tb"}->[$i]-
		      $Tb_ao*$ao_calib*$m1m3_ao{"Tb_model"}->[$i])+
		     ($Tx*$m1m3_ao{"Tx"}->[$i]-
		      $Tx_ao*$ao_calib*$m1m3_ao{"Tx_model"}->[$i])+
		     ($Ty*$m1m3_ao{"Ty"}->[$i]-
		      $Ty_ao*$ao_calib*$m1m3_ao{"Ty_model"}->[$i])+
		     ($Tz*$m1m3_ao{"Tz"}->[$i]-
		      $Tz_ao*$ao_calib*$m1m3_ao{"Tz_model"}->[$i])+
		     ($Tr*$m1m3_ao{"Tr"}->[$i]-
		      $Tr_ao*$ao_calib*$m1m3_ao{"Tr_model"}->[$i]));
	}

	# and add in any specified bending mode amplitude:
	foreach my $bendingmode (@bm_list) {
	    next if ($d{$bendingmode} == 0);
	    $output[2] += $d{$bendingmode}*$m1m3_ao{$bendingmode}->[$i];
	}
	printf F "%s\n",join(' ',@output);
    }
    close(F);

    foreach my $bendingmode (@bm_list) {

	my @sxx;
	my @sx;
	my @n;
	my ($sxx,$sx,$n)=(0,0,0);

	for (my $i=0;$i<=$#{$m1m3_ao{$keys[0]}};$i++) {
	    push(@sxx,pow($m1m3_ao{$bendingmode}->[$i],2));
	    push(@sx,$m1m3_ao{$bendingmode}->[$i]);
	    push(@n,1);
	}

	$sxx=sum(sort {$a<=>$b} @sxx);
	$sx =sum(sort {$a<=>$b} @sx);
	$n  =sum(sort {$a<=>$b} @n);
	
	# printf STDERR "for mode $bendingmode:\n";
	# printf STDERR "\tn=%d\n",$n;
	# printf STDERR "\tmean=%g\n",$sx/$n;
	# printf STDERR "\trms=%g\n",sqrt($sxx/$n-pow($sx/$n,2));
	# printf STDERR "\tsqrt(var)=%g\n",sqrt($sxx/$n);
    }

}


# here compute the rms deviation for each bending mode as recorded

exit;

# fit_modes(@condition,$bendmodes,@modes);

my $output_ascii=1;
if ($output_ascii) {

# output them into 3 separate tnt files for inspection

    my %hash=%Actuators;
    my @labels=@actuator_forces;

    open(F,">a.tnt") || die;
    printf F "actuators\n%s\n",join("\t",@labels);
    for (my $i=0;$i<=$#{$hash{$labels[0]}};$i++) {
	foreach my $col (@labels) {
	    printf F "%g ",$hash{$col}->[$i];
	}
	printf F "\n";
    }
    close(F);

    %hash=%BendingModes;
    @labels=@modes;

    open(F,">b.tnt") || die;
    printf F "bending modes\n%s\n",join("\t",@labels);
    for (my $i=0;$i<=$#{$hash{$labels[0]}};$i++) {
	foreach my $col (@labels) {
	    printf F "%g ",$hash{$col}->[$i];
	}
	printf F "\n";
    }
    close(F);

    %hash=%BendingModeActuatorForces;
    @labels=@fields;

    open(F,">c.tnt") || die;
    printf F "bending mode actuator forces\n%s\n",join("\t",@labels);
    for (my $i=0;$i<=$#{$hash{$labels[0]}};$i++) {
	foreach my $col (@labels) {
	    printf F "%g ",$hash{$col}->[$i];
	}
	printf F "\n";
    }
    close(F);

    %hash=%rawdata;
    @labels=@use_rawfields;
#    @labels=@rawfields;

    open(F,">d.tnt") || die;
    printf F "raw data\n%s\n",join("\t",@labels);
    for (my $i=0;$i<=$#{$hash{$labels[0]}};$i++) {
    	foreach my $col (@labels) {
    	    printf F "%g ",$hash{$col}->[$i];
    	}
    	printf F "\n";
    }
    close(F);
}

sub fit_modes {
    # length of the mode list
    my @guess;
    if (!exists($globals{"guess"})) {
	@guess = (0) x @{$globals{"modekeys"}};
    } else {
	@guess = @{$globals{"guess"}};
    }
    my @tol = (0.01) x @{$globals{"modekeys"}}; # starting tolerance: 10nm
    my ($p,$y)=MinimiseND(\@guess,\@tol,\&modefit,1e-3,100000);
    @{$globals{"guess"}}=() if (exists($globals{"guess"}));
    $globals{"guess"} = [ @{$p} ];
    @{$p};
}

sub fit_modes_modebymode {
    # this routine zeroes out all but one tolerance at a time to speed up
    # fitting in the case of a nearly orthogonal basis set.
    # length of the mode list
    my @guess;
    my @tol;
    if (!exists($globals{"guess"})) {
	@guess = (0) x @{$globals{"modekeys"}};
	@tol   = (0.1) x @{$globals{"modekeys"}};
    } else {
	@guess = @{$globals{"guess"}};
	for (my $i=0;$i<=$#guess;$i++) {
	    $tol[$i] = $guess[$i]/10.0;
	}
    }
    
    # get the ordering of the @guess values
    my @par_order=reverse sort {pow($guess[$a],2) <=> 
				    pow($guess[$b],2)} (0..$#guess);

    my ($p,$y);
    my @itmax  =(1000,1000,1000,1000,10000);
    my @fit_tol=(1e-2,1e-3,1e-3,1e-4,1e-6);
    my @pars_at_a_time=(3,3,3,3,scalar(@par_order));
    my @atten_factor=(10,10,100,100,1000);
    my @shift_pars=(1,2,3,3,scalar(@par_order));

    for (my $fit_iter=0;$fit_iter<5;$fit_iter++) {
	# set parameter fitting order by magnitude of the guess
	# starting tolerance: 10nm

	my $pars_at_a_time=$pars_at_a_time[$fit_iter];

	for (my $startpar=0;
	     $startpar<=$#par_order;
	     $startpar+=$shift_pars[$fit_iter]) {

	    @guess=@{$globals{"guess"}};

	    my @fitting_tol=(0)x@{$globals{"modekeys"}};
	    for (my $i=$startpar;$i<$startpar+$pars_at_a_time;$i++) {
		$fitting_tol[$par_order[$i%scalar(@par_order)]]=
		    $guess[$par_order[$i%scalar(@par_order)]]/
		    $atten_factor[$fit_iter];
	    }
	    
	    ($p,$y)=MinimiseND(\@guess,\@fitting_tol,\&modefit,
			       $fit_tol[$fit_iter],$itmax[$fit_iter]);
	    $globals{"FOM"}=$y;
	    @{$globals{"guess"}}=() if (exists($globals{"guess"}));
	    $globals{"guess"} = [ @{$p} ];
	}
    }
    @{$p};
}

sub model {
    my ($amps)=@_;
    my ($modehash,$modekeys)=@globals{"modehash","modekeys"};
    my @model = (0) x @{$modehash->{$modekeys->[0]}};
    for (my $node=0;$node<=$#model;$node++) {
	for (my $mode=0;$mode<=$#{$amps};$mode++) {
	    $model[$node]+=$amps->[$mode]*$modehash->{$modekeys->[$mode]}->[$node];
	}
    }
    @model;
}

sub modefit {
    my @amps=@_;
    my ($shape,$modehash,$modekeys,@model);
    # get these from the global hash
    $shape=$globals{"shape"};

    @model = model(\@amps);

    my @FOM=();
    my $n=0;

    for (my $node=0;$node<=$#{$shape};$node++) {
	push(@FOM,pow($shape->[$node]-$model[$node],2));
	$n++;
    }

    my $FOM=sqrt(sum(sort {$a<=>$b} @FOM)/$n)*1e3; # report FOM in um rms
    $globals{"iter"}++;

    # print to screen if appropriate
    if ($globals{"iter"}%10 == 0) {
	my @printamps;
	for (my $i=0;$i<=$#amps;$i++) {
	    push(@printamps,sprintf("%5.3g",$amps[$i]));
	}
	printf STDERR ("\riter %d FOM %g (%s)",
		       $globals{"iter"},
		       $FOM,
		       join(' ',@printamps));
    }

    $FOM;
}

sub sum {
    my ($sum,@rest)=@_;
    while (scalar(@rest)>0) {
	$sum += pop(@rest);	
    }
    $sum;
}

__END__

=head1 NAME

M1M3_perturbations.perl - Obtain an instance mapping of M1/M3 monolith perturbations with input FEAs and bending mode definitions provided by LSST mirror team

=head1 SYNOPSIS

M1M3_perturbations.perl [options]

[options] include: 

[--help|--?] [--man] [--reanalyse] [--Focus3=0(default 1)] [--Tb=B<Tb>] [--Tb_ao=B<Tb_ao>] [--Tx=B<Tx>] [--Tx_ao=B<Tx_ao>] [--Ty=B<Ty>] [--Ty_ao=B<Ty_ao>] [--Tz=B<Tz>] [--Tz_ao=B<Tz_ao>] [--Tr=B<Tr>] [--Tr_ao=B<Tr_ao>] [--theta=B<theta>] [--precomp_scale=B<precomp_scale>] [--ao_calib=B<ao_calib>] [--d_<bendingmode1>=<amplitude1>] [--d_<bendingmode2>=<amplitude2>] .. [--d_<bendingmodeN>=<amplitudeN>] ..

Bending mode amplitudes may be specified on a mode-by-mode basis. Unity amplitude corresponds to 1 micron r.m.s. B<measured about zero>. Bending mode names (bendingmodeX) are drawn from the following list of 30 possibilities:

B<Tip>, B<Tilt>, B<Power>, B<Astig1>, B<Astig2>, B<Focus3>, B<Tref4>, B<Tref5>, B<Coma6>, B<Coma7>, B<Quad8>, B<Quad9>, B<secAstig10>, B<secAstig11>, B<Spherical12>, B<Penta13>, B<Penta14>, B<15> thru B<27>.

A manual page is available with the -man switch.

=head1 SUMMARY

computes and outputs a distorted figure for the M1/M3 monolith, consistent with the available FEA data ("Raw Data Thermal2.csv", captured August 1st, 2013, origin: Doug Neill). FEA calculations were produced for the following seven driving conditions:

    1) zenith pointing
    2) horizon pointing
    3) bulk temperature difference of 1 degree C
    4) lateral (x) differential of 1 degree C
    5) lateral (y) differential of 1 degree C
    6) axial (z) differential of 1 degree C
    7) radial (r) differential of 1 degree C

The mirror team has also provided a set of bending modes, based on FEA/dynamic calculations ("Actuator Force Summary 2006.xlsx", captured August 13, 2013, origin: Doug Neill). Their names, listed below, are reminiscent of a Zernike expansion but are specifically not equal to those descriptions. These terms are named:

    a) Tip
    b) Tilt
    c) Power
    d) Astig1
    e) Astig2
    f) Focus3
    g) Tref4
    h) Tref5
    i) Coma6
    j) Coma7
    k) Quad8
    l) Quad9
    m) SecAstig10
    n) SecAstig11
    o) Spherical12
    p) Penta13
    q) Penta14
    r, s, t..ad) 15, 16, 17..27

C<M1M3_perturbations.perl> directly utilizes a subset of the data contained in "Actuator Force Summary 2006.xlsx" and "Raw_Data_thermal2.xls". The first of these files contains three worksheets, named:

=over 4

=item B<Actuator Forces>

This is a listing of the forces applied by each of 160 actuators (positions defined) in order to cancel out each of the seven drivind conditions listed above. This data is not used by M1M3_perturbations.perl.

=item B<Bending Mode Surfaces>

This is a listing of the surface deformation values (delta z) for each of 5244 nodes that correspond to each of the 30 (nominally orthogonal) bending modes named above. Each node is labeled by a unique index and its corresponding (X,Y) coordinates. This data is used by M1M3_perturbations.perl.

=item B<Bending Mode Actuator Forces>

This is a listing of the forces applied by each of 160 actuators (positions defined) in order to excite each of the 30 bending modes named above. This data is not used by M1M3_perturbations.perl.

=back

The second file contains a single worksheet, comprised of the Aug 1, 2013 CSV data referred to above. The worksheet is named:

=over 4

=item B<Raw Data Thermal2.csv>

This is a listing of the surface deformation values for each of 5244 FEA nodes when under each of the seven driving conditions listed above. Each node is labeled by the same unique index & (X,Y) coordinates (as in B<Bending Mode Surfaces> above), together with one of four possible designations for that node (C<PM OD>, C<PMCA>, C<TM CA>, C<TMID>). This data is used by M1M3_perturbations.perl.

=back

=head1 DESCRIPTION

M1M3_perturbations.perl utilizes the above two tables in two steps. An analysis step is invoked by specifying the B<--reanalyse> switch. This is a somewhat time consuming step that terminates after generating and storing an intermediate analysis data file that is in turn used in subsequent invokations that do not include the B<--reanalyse> switch. The default intermediate analysis data file name can be overridden and specified via the B<--intermediate_tnt=my file> switch.

The intermediate analysis file (B<--reanalyse> switch) is generated by:

=over 3

=item 1.

Data from B<Bending Mode Surfaces> and B<Raw Data Thermal2.csv> are read in from the B<xls> and B<xlsx> files. 

=item 2.

All values read in are converted to millimeter units from inches (for mirror X,Y positions) and from meters (for figure distortions). 

=item 3.

Surface functions from both origins are trimmed to remove real estate outside of the M1/M3 clear aperture (nodes labeled C<PM OD> are removed from consideration). 

=item 4.

Each of the 30 bending modes are normalized to 1 micron (1E-3 mm) r.m.s. using the discrete sampling defined by these FEA nodes. These normalizations are on average different from the normalizations provided in B<Bending Mode Surfaces> by 3% (30nm) per bending mode.

=item 5.

Each of the seven driving conditions listed above are fit using the 30 normalized bending modes. Starting amplitudes for each bending mode are estimated from inner product calculations - which are approximate at best, because node densities are not uniform and bending modes are not strictly orthogonal. Fitting is completed via a simplex algorithm, and 30 bending mode amplitudes per driving condition are recorded (for information only) in an output file B<rawdata_bymodes.tnt>.

=item 6.

The intermediate analysis file is used to record (for each FEA node, identified now by only its X,Y coordinates): the perturbation due to the driving condition, the best fit model to the driving condition, and the residual perturbation (data minus model) for that driving condition. An additional pair (model, perturbation) of values are recorded for the case where no B<Focus3> bending mode is applied in the aO compensation model - leaving a dominant B<Focus3> term in the residuals.

=item 7.

The 30 normalized bending modes are recorded in the intermediate analysis file

=back

The resulting intermediate analysis file contains 4968 rows and 67 columns. Generating the intermediate analysis file typically takes about an hour on a typical desktop machine. Do not fool around with the contents of the intermediate analysis file or swap around column names, or unexpected/meaningless output will result! Five columns are provided for each driving condition: raw FEA (B<P>), best fit model (B<M>), residual (B<R>), model excluding B<Focus3> (B<M3=0>) and the residual corresponding to a model that excludes B<Focus3> (B<R3=0>). Because the current aO compensation strategy (reference?) specifies non-excitation of the B<Focus3> bending mode - traded for for axial compensation by M2 and/or Camera instead, this can be achieved by the replacement of B<P3=0> and B<R3=0> for B<P> and B<R>, respectively with the B<--Focus3=0> switch. (The default setting of B<--Focus3=1> provides B<Focus3> aO compensation along with the rest of the bending modes.)

To maximize flexibility of M1M3_perturbations.perl, environmental terms are specified on a per-driving condition basis, with the exception of gravity load. In that case, a single parameter, B<--theta> specifies relative contributions of the B<z0> and B<z90> patterns, with gravity load contributions according to B<cos(theta)*z0 + sin(theta)*z90>. This expression is modified by a precompensation scalar (specified by the B<--precomp_scale> switch) - where if B<--precomp_scale=1.0> (default), the M1/M3 monolith has been precompensated for zenith pointing to have zero figure errors due to acceleration load there. This (default) configuration alters the perturbation expression given above to:

B<[cos(theta)-precomp_scale]*z0+sin(theta)*z90>.

The bending mode expansion based figure correction subtraction is performed with a corresponding term where z0 and z90 are replaced by the corresponding model functions (B<*_m>) for each term - but scaled by an overall aO calibration term (specified by the B<--ao_calib> switch):

B<ao_calib*([cos(theta)-precomp_scale]*z0_m +sin(theta)*z90_m)>.


Figure error contributions due to the five thermal driving conditions (B<Tb>, B<Tx>, B<Ty>, B<Tz> and B<Tr>) are specified by two parameters for each: one for the actual environmental condition, and one for the aO system's compensation lookup value. Any difference between the C<actual> and C<assumed> conditions can be expressed by specifying different values for these two amplitude parameters (e.g., B<--Tb> and B<--Tb_ao>, respectively) where the irreducible residual contribution scales with the value of the C<actual> condition value. The basic algorithm for each driving condition can be expressed (using amplitudes B<A> and B<A_ao> to go along with generic terms B<P>, B<M> and B<R> defined above:

B<A*P - ao_calib*A_ao*M = A*R + (A - ao_calib*A_ao)*M>

where by definition,

B<R> (residual) = B<P> (perturbative distortion) - B<M> (correction model).

A summation of the terms outlined above provides a method to predicted figure errors of the M1/M3 monolith based on environmental parameters and imperfect bending mode removal, whether due to a finite number of modes applied, aO wavefront recovery system latency, absolute force calibration changes, etc. This approach utilizes the following set of 14 parameter switches:

B<--Focus3> (boolean), B<--Tb>, B<--Tb_ao>, B<--Tx>, B<--Tx_ao>, B<--Ty>, B<--Ty_ao>, B<--Tz>, B<--Tz_ao>, B<--Tr>, B<--Tr_ao> (degrees C), B<--theta> (zenith distance in degrees), B<--precomp_scale> (amplitude by which acceleration load at zenith pointing is corrected) and B<--ao_calib> (scalar term affecting all actuators).

An additional set of switches is included to express individual bending modes in the output perturbation map. These are simply named after the bending modes as defined above, with amplitudes specified in units of 1 micron r.m.s. Bending mode excess (or deficit) may be useful in evaluating sensitivity matrices. These are the following 30 parameter switches:

B<--d_Tip>, B<--d_Tilt>, B<--d_Power>, B<--d_Astig1>, B<--d_Astig2>, B<--d_Focus3>, B<--d_Tref4>, B<--d_Tref5>, B<--d_Coma6>, B<--d_Coma7>, B<--d_Quad8>, B<--d_Quad9>, B<--d_SecAstig10>, B<--d_SecAstig11>, B<--d_Spherical12>, B<--d_Penta13>, B<--d_Penta14>, B<--d_15>, B<--d_16>, B<--d_17>, B<--d_18>, B<--d_19>, B<--d_20>, B<--d_21>, B<--d_22>, B<--d_23>, B<--d_24>, B<--d_25>, B<--d_26>, and B<--d_27> (microns). 

The preceding bending mode amplitudes are all set to zero by default. Expression of a single bending mode (e.g., Astig1) with 0.5 micron r.m.s. and with the figure reversed in sign from the tabulated bending mode - may be generated simply by the invocation:

M1M3_perturbations.perl --theta=0 --precomp_scale=1.0 --d_Astig1=-0.5

Arbitrary ad-mixtures of those bending modes may of course be specified by appending more corresponding switches.

Any update to the underlying FEA calculations will probably require an update to the current script M1M3_perturbations.perl - particularly since the choice of nodes is likely to change - and because more driving conditions may be computed. 

=head2 EXAMPLES

The following examples are provided to clear up any confusion the preceding explanations may have caused.

=over 3

=item 1.

Suppose the M1/M3 monolith has been figured to express zero amplitude of the B<z0> driving condition when facing zenith. The B<z0> expression shows up in all non-zero zenith distances then. With the M1/M3 monolith operating under the same conditions as when it was fabricated and figured, the non-acceleration load driving conditions are identically zero, and suppose the aO system uses zero amplitude correction patterns for each fo these conditions. The resulting figure error map is retrieved for zenith distance B<theta> by:

M1M3_perturbations.perl --theta=B<theta> [--precomp_scale=1.0] --output_tnt="my_output.tnt"

with output deposited into the file named "my_output.tnt". Specifying the B<precomp_scale> is optional because unity is its default value. No aO correction for the B<Focus3> bending mode can be achieved by adding in an extra switch:

M1M3_perturbations.perl --theta=B<theta> --output_tnt="my_output.tnt" --Focus3=0.

An uncompensated figure map may be inspected with:

M1M3_perturbations.perl --theta=B<theta> [--precomp_scale=1.0] --ao_calib=0

=item 2.

Suppose now that the M1/M3 monolith has been figured to minimize expression of B<z0> driving condition for a more LSST-representative zenith distance, say 30 degrees. Then, the B<precomp_scale> parameter should be set equal to cos(30) or 0.866. Suppose also that the mirror's temperature is 10 degrees cooler than when it was figured in the mirror lab, with the aO system correctly removing as much of the B<Tb> as can be removed with the available bending mode expansion. Let there be in addition to this a lateral thermal differential of 2 degrees across the M1/M3 diameter in the y direction (in the plane that includes the zenith and the boresight) - since the lower side is closer to the observatory floor and receives more blackbody radiation from it. Assume that at the time of the last aO update, that thermal differential was 1.75 degrees, and wavefront recovery was performed perfectly. The predicted figure error map is retrieved for zenith distance B<theta> (with no B<Focus3> compensation) by:

M1M3_perturbations.perl --theta=B<theta> --precomp_scale=0.866 --Focus3=0 --Ty=2.0 --Ty_ao=1.75 --Tb=-10

with the output stored in the default filename, M1M3_perturbations.tnt.

=back

=head2 RESIDUAL SCALES

Root-mean-square (about zero) figure errors are tabulated here for single driving condition inputs, for B<precomp_scale>=0,1 and B<ao_calib>=0,1:

=over 2

=item B<theta>=0:

 (precomp_scale,ao_calib)=(0,0): 13.29um rms; 0.286 arcsec (dz/dr) rms
 (precomp_scale,ao_calib)=(1,0): 0um     rms; 0.00  arcsec (dz/dr) rms
 (precomp_scale,ao_calib)=(0,1): 0.046um rms; 0.033 arcsec (dz/dr) rms
 (precomp_scale,ao_calib)=(1,1): 0um     rms; 0.00  arcsec (dz/dr) rms

=item B<theta>=90:

 (precomp_scale,ao_calib)=(0,0): 3.21um  rms; 0.685 arcsec (dz/dr) rms
 (precomp_scale,ao_calib)=(1,0): 13.68um rms; 0.740 arcsec (dz/dr) rms
 (precomp_scale,ao_calib)=(0,1): 0.069um rms; 0.047 arcsec (dz/dr) rms
 (precomp_scale,ao_calib)=(1,1): 0.082um rms; 0.058 arcsec (dz/dr) rms

=item B<theta>=0,B<precomp_scale>=1,B<Tb>=1,B<Tb_ao>=1:

 ao_calib=(0,1): (16.94,0.0315)um rms; (0.032,0.027) arcsec (dz/dr) rms

=item B<theta>=0,B<precomp_scale>=1,B<Tx>=1,B<Tx_ao>=1:

 ao_calib=(0,1): (0.478,0.0124)um rms; (0.062,0.011) arcsec (dz/dr) rms

=item B<theta>=0,B<precomp_scale>=1,B<Ty>=1,B<Ty_ao>=1:

 ao_calib=(0,1): (0.478,0.0124)um rms; (0.062,0.011) arcsec (dz/dr) rms

=item B<theta>=0,B<precomp_scale>=1,B<Tz>=1,B<Tz_ao>=1:

 ao_calib=(0,1): (12.18,0.640)um rms;  (1.037,0.420) arcsec (dz/dr) rms

=item B<theta>=0,B<precomp_scale>=1,B<Tr>=1,B<Tr_ao>=1:

 ao_calib=(0,1): (3.82,0.0432)um rms;  (0.140,0.050) arcsec (dz/dr) rms

=item B<theta>=0,B<precomp_scale>=1,B<d_Astig1>=1:

 ao_calib=(0,1): (1.000,1.000)um rms;  (0.102,0.102) arcsec (dz/dr) rms

=back

=head1 LICENSE

This is released under the Artistic License. See L<perlartistic>.

=head1 AUTHOR

Andy Rasmussen - L<mailto:arasmus@slac.stanford.edu>

=head1 SEE ALSO

L<M2_perturbations.perl>, L<CAMERA_perturbations.perl>

=cut
