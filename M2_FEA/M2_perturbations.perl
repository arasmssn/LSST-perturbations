#!/usr/bin/perl 
use strict;
use warnings;
use POSIX;
use Math::Amoeba qw(MinimiseND);
use Spreadsheet::Read;
use Getopt::Long;
use Pod::Usage;

# precomp_scale is set to -1 because according to Jacques, the 
# precompensation step is taken with the mirror surface pointing up. 
# It would be +1 if precompensation is made with the mirror surface
# pointing down (as in operational)

my $deg=atan2(1,1)/45.0;
my %globals;
my %BendingModes;
my %rawdata;
my $file;

my $myung_update_130820=1;

my %d=();

my @bm_list=("bm_1","bm_2","bm_3","bm_4","bm_5",
	     "bm_6","bm_7","bm_8","bm_9","bm_10"    );

@d{@bm_list}=(0,0,0,0,0,0,0,0,0,0);

my ($reanalyse,$Focus3,
    $Tz,$Tz_ao,$Tr,$Tr_ao,
    $theta,$precomp_scale,$ao_calib,
    $intermediate_tnt,$output_tnt)=
    (0,1,
     0.0,0.0,0.0,0.0,
     30,+1.0,1.0,
     "M2_intermediate_FEA.tnt","M2_perturbation.tnt");

my $help=0;
my $man=0;

GetOptions('help|?'  => \$help,
	   'man'     => \$man,
	   'reanalyse'=> \$reanalyse,
	   'Focus3=i'=> \$Focus3,
	   'Tz=f'    => \$Tz,
	   'Tz_ao=f' => \$Tz_ao,
	   'Tr=f'    => \$Tr,
	   'Tr_ao=f' => \$Tr_ao,
	   'theta=f' => \$theta,
	   'ao_calib=f' => \$ao_calib,
	   'precomp_scale=f' => \$precomp_scale,
	   ('d_1=f'  => \$d{"bm_1"},
	    'd_2=f'  => \$d{"bm_2"},
	    'd_3=f'  => \$d{"bm_3"},
	    'd_4=f'  => \$d{"bm_4"},
	    'd_5=f'  => \$d{"bm_5"},
	    'd_6=f'  => \$d{"bm_6"},
	    'd_7=f'  => \$d{"bm_7"},
	    'd_8=f'  => \$d{"bm_8"},
	    'd_9=f'  => \$d{"bm_9"},
	    'd_10=f' => \$d{"bm_10"}),
	   'intermediate_tnt=s'  => \$intermediate_tnt,
	   'output_tnt=s' => \$output_tnt) || die;;

pod2usage(-exitval=>0, -verbose=>1) if ($help);
pod2usage(-exitval=>0, -verbose=>2) if ($man);

my $DBKEY="M2_PERTURBATIONS_DB";
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
    $lbl =~ s/bm_/d_/g;
    push(@bm_distortion,"$lbl=$d{$bendingmode}")
	if ($d{$bendingmode} != 0);
}

printf STDERR "%s\n",join("\n",
			  "reanalyse=$reanalyse",
			  "Focus3=$Focus3",
			  "Tz=$Tz",
			  "Tz_ao=$Tz",
			  "Tr=$Tr",
			  "Tr_ao=$Tr_ao",
			  "theta=$theta",
			  "ao_calib=$ao_calib",
			  @bm_distortion,
			  "precomp_scale=$precomp_scale",
			  "intermediate_tnt=$intermediate_tnt",
			  "output_tnt=$output_tnt");

if ($reanalyse==1) {
    if ($myung_update_130820==1) {

	$file=join('/',$ENV{$DBKEY},"M2_GRF_Raw2.xlsx");
	my $xls=ReadData($file) || die "can't read $file";
	my @m2=("Xloc","Yloc","z0","z90");
	my %M2=();
	foreach my $i (0..$#m2) {
	    $M2{$m2[$i]}=[ @{$xls->[1]{cell}[$i+1]}[1..9084] ];
	}
	# reverse data for z0 - this was calculated for M2 face up
	
	$file=join('/',$ENV{$DBKEY},"M2_GrTzTr.xls");
	my $xls2=ReadData($file) || die "can't read $file";
	@m2=("Xloc","Yloc","z0_ao","z90_ao","Tz","Tz_ao","Tr","Tr_ao","R");
	# essentially overwrite $m2{$m2[0]} & $m2{$m2[1]}
	# read from separate sheets, so don't do this in a loop
	$M2{$m2[0]}=[ @{$xls2->[1]{cell}[1]}[1..9084] ];
	$M2{$m2[1]}=[ @{$xls2->[1]{cell}[2]}[1..9084] ];
	$M2{$m2[2]}=[ @{$xls2->[1]{cell}[3]}[1..9084] ];
	$M2{$m2[3]}=[ @{$xls2->[1]{cell}[4]}[1..9084] ];
	$M2{$m2[4]}=[ @{$xls2->[2]{cell}[3]}[1..9084] ];
	$M2{$m2[5]}=[ @{$xls2->[2]{cell}[4]}[1..9084] ];
	$M2{$m2[6]}=[ @{$xls2->[3]{cell}[3]}[1..9084] ];
	$M2{$m2[7]}=[ @{$xls2->[3]{cell}[4]}[1..9084] ];
	$M2{$m2[8]}=[];
	# add the radial column
	for (my $i=0;$i<=$#{$M2{$m2[0]}};$i++) {
	    $M2{$m2[8]}->[$i] = sqrt(pow($M2{"Xloc"}->[$i],2)+
				     pow($M2{"Yloc"}->[$i],2));
	}

	# convert positions to mm units

	@m2=("Xloc","Yloc","z0","z0_ao","z90","z90_ao",
	     "Tz","Tz_ao","Tr","Tr_ao","R");
	for (my $i=0;$i<=$#{$M2{$m2[0]}};$i++) {
	    for (my $col=0;$col<$#m2;$col++) {
		$M2{$m2[$col]}->[$i] *= 1e3; # convert input into mm
		# reverse data for z0 & z0_ao ..
		# these were calculated for M2 face up
		$M2{$m2[$col]}->[$i] *= -1 if ($m2[$col] =~ /^z0/);
	    }
	}
	
	# read in the bending modes

	$file=join('/',$ENV{$DBKEY},"M2_FRQ10s.xlsx");
	my $xls3=ReadData($file) || die "can't read $file";
	my @bmlist=("bm_1","bm_2","bm_3","bm_4","bm_5",
		    "bm_6","bm_7","bm_8","bm_9","bm_10");
	foreach my $i (0..$#bmlist) {
	    $M2{$bmlist[$i]}=[ @{$xls3->[1]{cell}[$i+3]}[1..9084] ];
	}

	# fabricate bending modes "bm_Tip","bm_Tilt" and "bm_Power" 
        # for fitting purposes.

	@bmlist=("bm_Tip","bm_Tilt","bm_Power",@bmlist);

	$M2{"bm_Tip"}=[];
	$M2{"bm_Tilt"}=[];
	$M2{"bm_Power"}=[];

	foreach (my $i=0;$i<=$#{$M2{$m2[0]}};$i++) {
	    push($M2{"bm_Tip"},$M2{"Xloc"}->[$i]);
	    push($M2{"bm_Tilt"},$M2{"Yloc"}->[$i]);
	    push($M2{"bm_Power"},1);
	}
	
	# normalize the bending modes to 1um rms
	my %bm_norm=();
	foreach my $bm (@bmlist) {
	    my @amp=();
	    foreach my $sample (@{$M2{$bm}}) {
		push(@amp,pow($sample,2));
	    }
	    $bm_norm{$bm}=sqrt(sum(sort {$a <=> $b} @amp)/scalar(@{$M2{$bm}}));

	    # normalize
	    my $ref=$M2{$bm};

	    for (my $i=0;$i<=$#{$ref};$i++) {
		$ref->[$i] *= (1.0e-3/$bm_norm{$bm}) if defined($ref->[$i]);
	    }
	}

	# any fitting should be done here
	my @rawfields_tofit=("z0","z90","Tz","Tr");
	my @fitmodes=("bm_Tip","bm_Tilt","bm_Power",
		      "bm_1","bm_2","bm_3","bm_4","bm_5",
		      "bm_6","bm_7","bm_8","bm_9","bm_10");
	my @fits=();

	@BendingModes{@fitmodes}   = @M2{@fitmodes};
	@rawdata{@rawfields_tofit} = @M2{@rawfields_tofit};

	$file=join('/',$ENV{$DBKEY},$intermediate_tnt);
	open(F,">$file") || die "can't open $file for writing";

	$file=join('/',$ENV{$DBKEY},"rawdata_bymodes.tnt");
	open(GG,">$file") || die "can't open $file for writing";
	printf GG "rawdata_bymodes\n%s\n",join("\t",@rawfields_tofit);

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
			$amp +=($surface->[$i] * $BendingModes{$bm}->[$i]);
			$count+= 1;
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

	    my $label;
            # $label=sprintf("%s_%s",$raw_tofit,"raw");
            # push(@fits,$label); # echo the raw data                         
            # $M2{$label} = $surface;

            $label=sprintf("%s_%s",$raw_tofit,"model");
            push(@fits,$label);
            $M2{$label} = [ @mod ];

            $label=sprintf("%s_%s",$raw_tofit,"residuals");
            push(@fits,$label);
            $M2{$label} = [];

            for (my $i=0;$i<=$#mod;$i++) {
                $M2{$label}->[$i]=$surface->[$i]-$mod[$i];
            }

	    # now make up a map containing bending mode bm_5 only with finite
	    # value
	    
	    foreach my $i (0..$#amplitudes) {
		next if ($i==7);
		$amplitudes[$i]=0;
	    }
	    @mod=model(\@amplitudes);

            $label=sprintf("%s_%s",$raw_tofit,"Focus3only");
            push(@fits,$label);
            $M2{$label} = [ @mod ];

	}

	my %hash;
	my @labels;
	# %M2 already has entries for output bending modes
	%hash=%M2;
	@labels=@m2;
	@labels=(@m2,@bmlist);
	@labels=(@m2,@fits,@bmlist);

	printf F "M2_GrTzTr\n%s\n",join("\t",@labels);
	for (my $i=0;$i<=$#{$hash{$labels[0]}};$i++) {
	    foreach my $col (@labels) {
		printf F "%g ",$hash{$col}->[$i];
	    }
	    printf F "\n";
	}
	close(F);
	
    } else {
	$file=join('/',$ENV{$DBKEY},"M2_GR0_81010FromCho.xls");
	my $xls2=ReadData($file) || die "can't read $file";

	$file=join('/',$ENV{$DBKEY},"M2_TH2_81112_FromCho.xls");
	my $xls=ReadData($file) || die "can't read $file";

	my @gr=("Xloc","Yloc","z0_ao","z90_ao","R");
	my %GR=();
	$GR{$gr[0]}=[ @{$xls2->[1]{cell}[1]}[1..9084] ];
	$GR{$gr[1]}=[ @{$xls2->[1]{cell}[2]}[1..9084] ];
	$GR{$gr[2]}=[ @{$xls2->[1]{cell}[3]}[1..9084] ];
	$GR{$gr[3]}=[ @{$xls2->[1]{cell}[4]}[1..9084] ];
	$GR{$gr[4]}=[];
# add a radial column
	for (my $i=0;$i<=$#{$GR{$gr[0]}};$i++) {
	    $GR{"R"}->[$i] = sqrt(pow($GR{"Xloc"}->[$i],2)+pow($GR{"Yloc"}->[$i],2));
	}
	for (my $col=0;$col<5;$col++) {
	    for (my $i=0;$i<=$#{$GR{$gr[$col]}};$i++) {
		$GR{$gr[$col]}->[$i] *= 1e3; # convert input into mm
	    }
	}

	my @th=("Xloc","Yloc","Tz","Tr","R");
	my %TH=();
	$TH{$th[0]}=[ @{$xls->[1]{cell}[1]}[1..9084] ];
	$TH{$th[1]}=[ @{$xls->[1]{cell}[2]}[1..9084] ];
	$TH{$th[2]}=[ @{$xls->[1]{cell}[3]}[1..9084] ];
	$TH{$th[3]}=[ @{$xls->[1]{cell}[4]}[1..9084] ];
	$TH{$th[4]}=[];
# add a radial column
	for (my $i=0;$i<=$#{$TH{$th[0]}};$i++) {
	    $TH{"R"}->[$i] = sqrt(pow($TH{"Xloc"}->[$i],2)+pow($TH{"Yloc"}->[$i],2));
	}
	for (my $col=0;$col<5;$col++) {
	    for (my $i=0;$i<=$#{$TH{$th[$col]}};$i++) {
		$TH{$th[$col]}->[$i] *= 1e3; # convert input into mm
	    }
	}

	my %hash;
	my @labels;

	%hash=%GR;
	@labels=@gr;

	$file=join('/',$ENV{$DBKEY},"gr.tnt");
	open(F,">$file") || die "can't open $file for writing";
	printf F "M2_GR\n%s\n",join("\t",@labels);
	for (my $i=0;$i<=$#{$hash{$labels[0]}};$i++) {
	    foreach my $col (@labels) {
		printf F "%g ",$hash{$col}->[$i];
	    }
	    printf F "\n";
	}
	close(F);

	%hash=%TH;
	@labels=@th;

	$file=join('/',$ENV{$DBKEY},"th.tnt");
	open(F,">$file") || die "can't open $file for writing";
	printf F "M2_TH data\n%s\n",join("\t",@labels);
	for (my $i=0;$i<=$#{$hash{$labels[0]}};$i++) {
	    foreach my $col (@labels) {
		printf F "%g ",$hash{$col}->[$i];
	    }
	    printf F "\n";
	}
	close(F);
    }
} else {
    if ($myung_update_130820==1) {
	my %m2_ao=();
	$file=join('/',$ENV{$DBKEY},$intermediate_tnt);
	open(F,$file) || die "can't find data file $file\n";
	my @m2_fea_ao=<F>;
	close(F);
	chomp(@m2_fea_ao);
	shift(@m2_fea_ao);
	my @keys=split("\t",shift @m2_fea_ao);
	foreach my $key (@keys) {
	    $m2_ao{$key}=[];
	}
	while (my $line=shift(@m2_fea_ao)) {
	    @_=split(' ',$line);
	    my $i=0;
	    foreach my $key (@keys) {
		push($m2_ao{$key},$_[$i]);
		$i++;
	    }
	}
	# at this point ready to output an instance of the FEA calculation
	my ($costheta,$sintheta)=(cos($theta*$deg),sin($theta*$deg));
	open(G,">$output_tnt") || die "can't open $output_tnt";
	my @outputs=("Xloc","Yloc","perturbation");
	printf G "interpolated data\n%s\n",join("\t",@outputs);

	$m2_ao{$outputs[$#outputs]}=[];
	for (my $i=0;$i<=$#{$m2_ao{$keys[0]}};$i++) {
	    # no ao_calib for gravity (no raw data)
	    my $perturbation;

	    my ($raw,$model,$uncorrected)=(
		$m2_ao{"z0_ao"}->[$i]+$m2_ao{"z0"}->[$i],
		$m2_ao{"z0"}->[$i],
		(($Focus3==0)?$m2_ao{"z0_Focus3only"}->[$i]:0));
	    $perturbation=
		($costheta-$precomp_scale)*($raw-$model)+
		$costheta*($model - $ao_calib*($model-$uncorrected));
	    ($raw,$model,$uncorrected)=(
		$m2_ao{"z90_ao"}->[$i]+$m2_ao{"z90"}->[$i],
		$m2_ao{"z90"}->[$i],
		(($Focus3==0)?$m2_ao{"z90_Focus3only"}->[$i]:0));
	    $perturbation+= 
		$sintheta*($raw-$ao_calib*($model-$uncorrected));
	    ($raw,$model,$uncorrected)=(
		$m2_ao{"Tz"}->[$i],
		$m2_ao{"Tz"}->[$i]-$m2_ao{"Tz_ao"}->[$i],
		(($Focus3==0)?$m2_ao{"Tz_Focus3only"}->[$i]:0));
	    $perturbation+= 
		$Tz*$raw-$ao_calib*$Tz_ao*($model-$uncorrected);
	    ($raw,$model,$uncorrected)=(
		$m2_ao{"Tr"}->[$i],
		$m2_ao{"Tr"}->[$i]-$m2_ao{"Tr_ao"}->[$i],
		(($Focus3==0)?$m2_ao{"Tr_Focus3only"}->[$i]:0));
	    $perturbation+= 
		$Tr*$raw-$ao_calib*$Tr_ao*($model-$uncorrected);

	    # add in any bending modes specified
	    foreach my $bendingmode (@bm_list) {
		next if ($d{$bendingmode} == 0);
		$perturbation += $d{$bendingmode}*$m2_ao{$bendingmode}->[$i];
	    }
	    # done. push to output list
	    
	    push($m2_ao{$outputs[$#outputs]},$perturbation);
	}

	for (my $i=0;$i<=$#{$m2_ao{$keys[0]}};$i++) {
	    foreach my $key (@outputs) {
		printf G "%g ",$m2_ao{$key}->[$i];
	    }
	    printf G "\n";
	}
    } else {
	my %m2_ao=();
	$file=join('/',$ENV{$DBKEY},"gr.tnt");
	open(F,$file) || die "can't find data file $file\n";
	my @m2_gr_fea_ao=<F>;
	close(F);
	chomp(@m2_gr_fea_ao);
	shift(@m2_gr_fea_ao);
	my @keys=split("\t",shift @m2_gr_fea_ao);
	foreach my $key (@keys) {
	    $m2_ao{$key}=[];
	}
	while (my $line=shift(@m2_gr_fea_ao)) {
	    @_=split(' ',$line);
	    my $i=0;
	    foreach my $key (@keys) {
		push($m2_ao{$key},$_[$i]);
		$i++;
	    }
	}
	$file=join('/',$ENV{$DBKEY},"th.tnt");
	open(F,$file) || die "can't find data file $file\n";
	my @m2_th_fea_ao=<F>;
	close(F);
	chomp(@m2_th_fea_ao);
	shift(@m2_th_fea_ao);
	@keys=split("\t",shift @m2_th_fea_ao);
	foreach my $key (@keys) {
	    if (defined($m2_ao{$key})) {
		@{$m2_ao{$key}}=(); # overwrite any previous data
	    }
	    $m2_ao{$key}=[];
	}
	while (my $line=shift(@m2_th_fea_ao)) {
	    @_=split(' ',$line);
	    my $i=0;
	    foreach my $key (@keys) {
		push($m2_ao{$key},$_[$i]);
		$i++;
	    }
	}
	# at this point ready to output an instance of the FEA calculation
	my ($costheta,$sintheta)=(cos($theta*$deg),sin($theta*$deg));
	open(G,">$output_tnt") || die "can't output $output_tnt";
	my @outputs=("Xloc","Yloc","perturbation");
	printf G "data\n%s\n",join("\t",@outputs);

	$m2_ao{$outputs[$#outputs]}=[];
	for (my $i=0;$i<=$#{$m2_ao{$keys[0]}};$i++) {
	    my $perturbation=((($costheta-$precomp_scale)*$m2_ao{"z0_ao"}->[$i]+
			       $sintheta*$m2_ao{"z90_ao"}->[$i])+
			      ($Tz*$m2_ao{"Tz"}->[$i]+
			       $Tr*$m2_ao{"Tr"}->[$i]));
	    push($m2_ao{$outputs[$#outputs]},$perturbation);
	}

	for (my $i=0;$i<=$#{$m2_ao{$keys[0]}};$i++) {
	    foreach my $key (@outputs) {
		printf G "%g ",$m2_ao{$key}->[$i];
	    }
	    printf G "\n";
	}
    }
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

sub model {
    my ($amps)=@_;
    my ($modehash,$modekeys)=@globals{"modehash","modekeys"};
    my @model = (0) x @{$modehash->{$modekeys->[0]}};
    for (my $node=0;$node<=$#model;$node++) {
        for (my $mode=0;$mode<=$#{$amps};$mode++) {
	    $model[$node] +=
		$amps->[$mode]*$modehash->{$modekeys->[$mode]}->[$node];
        }
    }
    @model;
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

M2_perturbations.perl - Obtain an instance mapping of M2 perturbations with input FEAs and bending mode definitions provided by LSST mirror team

=head1 SYNOPSIS

M2_perturbations.perl [options]

[options] include: 

[--help|--?] [--man] [--reanalyse] [--Focus3=0(default 1)] [--Tz=B<Tz>] [--Tz_ao=B<Tz_ao>] [--Tr=B<Tr>] [--Tr_ao=B<Tr_ao>] [--theta=B<theta>] [--precomp_scale=B<precomp_scale>(default +1)] [--ao_calib=B<ao_calib>] [--d_<bendingmode1>=<amplitude1>] [--d_<bendingmode2>=<amplitude2>] ..  [--d_<bendingmodeN>=<amplitudeN>] ..

Bending mode amplitudes may be specified on a mode-by-mode basis. Unity amplitude corresponds to 1 micron r.m.s. B<measured about zero>. Bending mode names (bendingmodeX) are drawn from the following list of 13 possibilities:                                                                                      
B<Tip>, B<Tilt>, B<Power>, B<Astig1>, B<1>, B<2> thru B<10>, where the numbered bending modes correspond to the first 10 "free-free" modes provided by Myung Cho.

A manual page is available with the -man switch.

=head1 SUMMARY

computes and outputs a distorted figure for the M2 optic, consistent with the available FEA data ("M2_GRF_Raw2.xlsx" (document-15129) and "M2_GrTzTr.xls" (document-14865), captured October 7th, 2013 and August 20th, 2013, respectively; origin: Myung Cho). FEA calculations were produced for the following four driving conditions:

    1) face up (G0, opposite from zenith pointing)
    2) face horizontal (G90, or horizon pointing)
    3) axial (z) differential of 1 degree C
    4) radial (r) differential of 1 degree C

The mirror team has also provided a set of 10 bending modes, based on FEA/dynamic calculations ("M2_FRQ10s.xlsx" (document-15128), captured October 7, 2013, origin: Myung Cho). These are the data displayed in Fig. 9 of Cho, Liang & Neill, SPIE 7424-06, August 2009. This truncation in bending mode expansion appears somwhat arbitrary, and we therefore utilize a different strategy in delivering M2 surface errors on the commandline.

C<M2_perturbations.perl> directly utilizes the data contained in "M2_GRF_Raw2.xlsx" and "M2_GrTzTr.xls", which represent FEA calculations for 9084 nodes across M2. The first of these files contains a single, unnamed worksheet, allegedly containing B<raw> acceleration load FEA results, with columns B<X>, B<Y>, B<G0> and B<G90>. Myung kindly provided this file in response to our request, when we noticed inconsistency in the level of detail provided in the other file ("M2_GrTzTr.xls") for the two acceleration load cases. The latter file contains three worksheets, named:

=over 4

=item B<Gravity(Zen&Hor)>

Unnamed columns apparently contain B<X>, B<Y>, B<G0_ao> and B<G90_ao> for the two gravity load cases, but B<after> the aO system has optimally removed a high fidelity surface distortion model (an expansion in bending modes) for each case. Thus, the 3rd and 4th columns contain the best possible figure possible, for those two orientations, relative to an initial, zero load case. (The raw deformations for each of these orientations were provided in the first file, "M2_GRF_Raw2.xlsx", as indicated above.)

=item B<Tz&Tz_ao>

Unnamed columns apparently contain B<X>, B<Y>, B<Tz> and B<Tz_ao>, where B<Tz> is the raw surface error contribution resulting from a 1 degree axial temperature differential through M2, and B<Tz_ao> is the residual error, per degree differential, following optimal bending mode removal by the aO system. 

=item B<Tr&Tr_ao>

Unnamed columns apparently contain B<X>, B<Y>, B<Tr> and B<Tr_ao>, where B<Tr> is the raw surface error contribution resulting from a 1 degree radial temperature differential through M2, and B<Tr_ao> is the residual error, per degree differential, following optimal bending mode removal by the aO system.

=back

=head1 DESCRIPTION

M2_perturbations.perl utilizes the above two tables in two steps. An analysis step is invoked by specifying the B<--reanalyse> switch. This terminates after generating and storing an intermediate analysis data file that is in turn used in subsequent invokations that do not include the B<--reanalyse> switch. The default intermediate analysis data file name can be overridden and specified via the B<--intermediate_tnt=my file> switch.

The intermediate analysis file (B<--reanalyse> switch) is generated by:

=over 3

=item 1.

Data from B<M2_GRF_Raw2.xlsx> and B<M2_GrTzTr.xls> are read into memory.

=item 2.

All values read in are converted to millimeter units from meters.

=item 3.

FEA results for B<G0> input are inverted for zenith (B<z0>) pointing, while B<G90> are used as-is for horizon (B<z90>) pointing.

=item 4.

Three (3) additional terms are added to the list of 10 bending modes: we add solid body terms B<Tip>, B<Tilt> and B<Power> to be used in fitting. Each of the 13 modes are normalized to 1 micron (1E-3 mm) r.m.s. using the discrete sampling defined by these FEA nodes.

=item 5.

Each of the four driving conditions listed above are fit using the 13 normalized modes (10 bending modes + 3). Starting amplitudes for each bending mode are estimated from inner product calculations - which are approximate at best, because node densities are not uniform and bending modes are not strictly orthogonal. Fitting is completed via a simplex algorithm, and 13 mode amplitudes per driving condition are recorded (for information only) in an output file B<rawdata_bymodes.tnt>. Two curious features were noted following fitting: (1) the acceleration based driving conditions were fit with residuals that were a fraction of the amplitude of the residuals seen in the corresponding "ao" dataset, with no apparent print-through of the actuator distribution. We assume therefore, that this dataset was drawn from a "best fit model" rather than a raw FEA calculation, as advertised. (2) the temperature based driving conditions were fit with residuals that were several times greater than the amplitude of the residuals seen in the corresponding "ao" dataset. This is not surprising, because Cho et al. (2009) state that 24 bending modes were used in determining aO limited surface errors, while our database contains only the first 10 bending modes, which each have resonant frequencies that range between 37 and 220Hz. Higher order terms are expected to excite with lower amplitudes, so the currently available set of bending modes may be used to estimate vibration induced image blur by themselves. 

=item 6.

The intermediate analysis file is used to record (for each FEA node, identified now by only its X,Y coordinates), the "raw" and "ao" data sets, the 13-mode best fit model and the residual perturbation (data minus model). For example the zenith pointed case data sets include entries B<z0>, B<z0_ao>, B<z0_model> and B<z0_residuals>. However, see discussion in B<5> above on our concern for the validity of gravity induced perturbation raw data. We proceed below (for the time being) with the assumption that this data was provided (described) incorrectly. Any corrections to the above assumptions may be made, an expanded set of bending modes provided, and/or new FEA data may be used, to generate an updated modeling paradigm with a short turnaround time.

It is still unclear whether the aO system will correct for focus mode terms on M2, or intentionally leave this mode uncorrected (as in the case of M1/M3). The current approach within M2_perturbations.perl is to assume that all 24 bending modes will be excited to minimize the surface error. 

=item 7.

The 13 normalized modes (3 solid body & 10 bending modes) are recorded in the intermediate analysis file, indicated by column names B<bm_Tip>, B<bm_Tilt>, B<bm_Power>, B<bm_1> .. B<bm_10>.

=back

The resulting intermediate analysis file contains 9084 rows and 32 columns. Generating the intermediate analysis file typically takes several minutes on a typical desktop machine, because fitting using available bending modes is performed. Do not fool around with the contents of the intermediate analysis file or swap around column names, or unexpected/meaningless output will result! Just Four columns are provided for each driving condition: e.g., raw FEA (B<P>), aO limited surface error (B<P_ao>), best fit model (B<P_m>), and model fit residual (B<P_r>).

Because there are evidently a few problems in the comparisons (see above) between fitting residuals and Myung's reported "aO limited" surface errors, we proceed to utilize the raw (B<P>) and aO limited surface errors (B<P_ao>) and interpret them for their face value. From here on we no longer refer to the best fit model (B<P_m>) nor to the model fit residual (B<P_r>) and temporarily ignore the fact that the acceleration induced "raw" data contain no print-through contributions.

To maximize flexibility of M2_perturbations.perl, environmental terms are specified on a per-driving condition basis, with the exception of gravity load. In that case, a single parameter, B<--theta> specifies relative contributions of the B<z0> and B<z90> patterns, with gravity load contributions according to B<cos(theta)*z0 + sin(theta)*z90>. This expression is modified by a precompensation scalar (specified by the B<--precomp_scale> switch) - where if B<--precomp_scale=1.0> (default), the M2 monolith has been precompensated for zenith pointing (surface facing downward - different from the assumptions laid out in Cho et al. 2009, based on word of mouth from Jacques Sebag) to have zero figure errors due to acceleration load there. This (default) configuration alters the perturbation expression given above to:

B<[cos(theta)-precomp_scale]*z0+sin(theta)*z90>.

The option for M2 final figuring being performed in the face up orientation can be simulated by specifying B<--precomp_scale>=-1 on the command line (e.g. as assumed by Cho et al. 2009).

The bending mode expansion based figure correction subtraction is performed with a corresponding term where z0 and z90 are replaced by the corresponding model function proxy (B<P-P_ao>) for each term - but scaled by an overall aO calibration term (specified by the B<--ao_calib> switch). The correction term that is subtracted is:

B<[ao_calib*cos(theta)-precomp_scale]*(z0-z0_ao)+ao_calib*sin(theta)*(z90-z90_ao))>.


Figure error contributions due to the two thermal driving conditions (B<Tz> and B<Tr>) are specified by two parameters for each: one for the actual environmental condition, and one for the aO system's compensation lookup value. Any difference between the C<actual> and C<assumed> conditions can be expressed by specifying different values for these two amplitude parameters (e.g., B<--Tz> and B<--Tz_ao>, respectively) where the irreducible residual contribution scales with the value of the C<actual> condition value. The basic algorithm for each driving condition can be expressed (using amplitudes B<A> and B<A_ao> to go along with generic terms B<P> and B<P_ao> defined above:

B<A*P - ao_calib*A_ao*(P-P_ao) = A*P_ao + (A - ao_calib*A_ao)*(P-P_ao)>

where by our working assumption,

B<P_ao> (residual) = B<P> (perturbative distortion) - B<P_m> (correction model).

A summation of the terms outlined above provides a method to predicted figure errors of the M2 optic based on environmental parameters and imperfect bending mode removal, whether due to a finite number of modes applied, aO wavefront recovery system latency, absolute force calibration changes, etc. This approach utilizes the following set of 7 parameter switches:

B<--Tz>, B<--Tz_ao>, B<--Tr>, B<--Tr_ao>, (degrees C), B<--theta> (zenith distance in degrees), B<--precomp_scale> (amplitude by which acceleration load at zenith pointing is corrected) and B<--ao_calib> (scalar term affecting all actuators).

B<--Focus3=0> (boolean) has also been implemented, which specifies that the Focus3 term in the surface error will not be corrected by pusher forces, but that image performance will be recovered by solid body axial translation of M2, the Camera, or both. In this case, expressions written above for the correction terms (acceleration load and thermal) are modified slightly to so that the B<Focus3> term survives in the resulting surface error (an additional column for each driving condition corresponding to the B<Focus3> component exists in the intermediate file).

The available parameter space in B<precomp_scale>, B<ao_calib> and <Focus3> for each driving condition were carefully examined to verify expected behavior in the resulting surface error.

An additional set of switches is included to express individual bending modes in the output perturbation map. These are simply named after the bending modes as defined above, with amplitudes specified in units of 1 micron r.m.s. Bending mode excess (or deficit) may be useful in evaluating sensitivity matrices. These are the following 30 parameter switches:

B<--d_Tip>, B<--d_Tilt>, B<--d_Power>, B<--d_1>, B<--d_2> .. B<--d_10> (microns). 

The preceding bending mode amplitudes are all set to zero by default. Expression of a single bending mode (e.g., bending mode 1) with 0.5 micron r.m.s. and with the figure reversed in sign from the tabulated bending mode - may be generated simply by the invocation:

M2_perturbations.perl --theta=0 --precomp_scale=1.0 --d_1=-0.5

Arbitrary ad-mixtures of those bending modes may of course be specified by appending more corresponding switches.

Any update to the underlying FEA calculations will probably require an update to the current script M2_perturbations.perl - particularly since the choice of nodes is likely to change - and because more driving conditions may be computed. 

=head2 EXAMPLES

The following examples are provided to clear up any confusion the preceding explanations may have caused.

=over 3

=item 1.

Suppose the M2 optic has been figured to express zero amplitude of the B<z0> driving condition when facing zenith. The B<z0> expression shows up in all non-zero zenith distances then. With M2 operating under the same conditions as when it was fabricated and figured, the non-acceleration load driving conditions are identically zero, and suppose the aO system uses zero amplitude correction patterns for each fo these conditions. The resulting figure error map is retrieved for zenith distance B<theta> by:

M2_perturbations.perl --theta=B<theta> [--precomp_scale=1.0] --output_tnt="my_output.tnt"

with output deposited into the file named "my_output.tnt". Specifying the B<precomp_scale> is optional because unity is its default value.

An uncompensated figure map may be inspected with:

M2_perturbations.perl --theta=B<theta> [--precomp_scale=1.0] --ao_calib=0

=item 2.

Suppose now that the mirror's temperature has a 0.5 degree axial differential together with a 1 degree radial differential. Assume that at the time of the last aO update, that radial thermal differential was 0.75 degrees, and wavefront recovery was performed perfectly. The predicted figure error map is retrieved for zenith distance B<theta> by:

M2_perturbations.perl --theta=B<theta> [--precomp_scale=1] --Tz=0.5 --Tz_ao=0.5 --Tr=1.0 --Tr_ao=0.75

with the output stored in the default filename, M2_perturbations.tnt.

=back

=head2 RESIDUAL SCALES

Root-mean-square (about zero) figure errors are tabulated here for single driving condition inputs, for B<precomp_scale>=0,1 and B<ao_calib>=0,1:

=over 2

=item B<theta>=0:

 (precomp_scale,ao_calib)=(0,0): 0.078um rms; 0.034 arcsec (dz/dr) rms
 (precomp_scale,ao_calib)=(1,0): 0.078um rms; 0.022 arcsec (dz/dr) rms
 (precomp_scale,ao_calib)=(0,1): 0.008um rms; 0.027 arcsec (dz/dr) rms
 (precomp_scale,ao_calib)=(1,1): 0um     rms; 0.00  arcsec (dz/dr) rms

=item B<theta>=90:

 (precomp_scale,ao_calib)=(0,0): 0.035um rms; 0.041 arcsec (dz/dr) rms
 (precomp_scale,ao_calib)=(1,0): 0.036um rms; 0.049 arcsec (dz/dr) rms
 (precomp_scale,ao_calib)=(0,1): 0.010um rms; 0.038 arcsec (dz/dr) rms
 (precomp_scale,ao_calib)=(1,1): 0.013um rms; 0.046 arcsec (dz/dr) rms

=item B<theta>=0,B<precomp_scale>=1,B<Tz>=1,B<Tz_ao>=1:

 ao_calib=(0,1): (0.096,0.0008)um rms; (0.065,0.003) arcsec (dz/dr) rms

=item B<theta>=0,B<precomp_scale>=1,B<Tr>=1,B<Tr_ao>=1:

 ao_calib=(0,1): (0.078,0.000)um rms; (0.021,0.000) arcsec (dz/dr) rms

=back

=head1 LICENSE

This is released under the Artistic License. See L<perlartistic>.

=head1 AUTHOR

Andy Rasmussen - L<mailto:arasmus@slac.stanford.edu>

=head1 SEE ALSO

L<M2_perturbations.perl>, L<CAMERA_perturbations.perl>

=cut
