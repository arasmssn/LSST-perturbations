#!/usr/bin/perl -w
#using John Ku's summary tables (gravity distortion & thermal distortion)
# gravity distortion part:
# 	T0=23°C	Rigid Body Distortions w.r.t. CCS		
# 	Optic/DOF	GX	Gy	Gz
# L1	TX [microns]	345.06	5.72	-0.47
# 	TY [microns]	-1.74	-396.10	-7.22
# 	TZ [microns]	0.06	-3.00	-60.09
# 	RX [micro-radians]	-0.33	-36.94	0.09
# 	RY [micro-radians]	21.70	1.62	0.85
# 	RZ [micro-radians]	-93.82	-2.76	0.12
# L2	TX [microns]	374.93	6.85	0.07
# 	TY [microns]	-1.26	-396.57	-7.09
# 	TZ [microns]	0.58	-4.27	-93.26
# 	RX [micro-radians]	6.46	-30.62	-7.30
# 	RY [micro-radians]	25.95	-0.98	-9.53
# 	RZ [micro-radians]	0.00	0.00	0.00
# Filters	TX [microns]	-80.18	-364.48	9.84
# 	TY [microns]	0.37	-7.92	-66.78
# 	TZ [microns]	85.72	14.85	0.04
# 	RX [micro-radians]	26.63	-13.32	-0.44
# 	RY [micro-radians]	90.94	579.64	-12.46
# 	RZ [micro-radians]	-6.46	-69.69	42.47
# L3	TX [microns]	22.52	0.00	-0.07
# 	TY [microns]	0.04	-22.31	0.12
# 	TZ [microns]	0.01	-0.08	-11.37
# 	RX [micro-radians]	0.05	-5.43	0.02
# 	RY [micro-radians]	-6.43	0.00	0.07
# 	RZ [micro-radians]	0.05	0.00	0.00
# Focal Plane	TX [microns]	154.67	0.00	-0.07
# 	TY [microns]	0.02	-111.44	-0.13
# 	TZ [microns]	0.01	-0.81	-15.39
# 	RX [micro-radians]	0.06	-11.23	-4.73
# 	RY [micro-radians]	-7.24	0.02	0.22
# 	RZ [micro-radians]	-33.20	0.04	0.03
# ---------------------------------------------------
# thermal distortion part
# 	T0=23C	Rigid Body Distortions w.r.t. CCS		
# 	Optic/DOF	Case 1 Nominal	Case 2 Cold	Case 3 Hot
# L1	TX [microns]	-1.32	-2.04	-0.59
# 	TY [microns]	32.47	53.13	11.80
# 	TZ [microns]	403.08	655.14	151.01
# 	RX [micro-radians]	29.22	47.59	10.86
# 	RY [micro-radians]	-0.89	-1.41	-0.36
# 	RZ [micro-radians]	-1.56	-2.58	-0.55
# L2	TX [microns]	-1.83	-2.87	-0.80
# 	TY [microns]	17.98	29.53	6.42
# 	TZ [microns]	277.09	450.42	103.77
# 	RX [micro-radians]	29.16	47.49	10.84
# 	RY [micro-radians]	-1.04	-1.66	-0.42
# 	RZ [micro-radians]	0.00	0.00	0.00
# Filters	TX [microns]	-0.04	-0.10	0.02
# 	TY [microns]	-9.83	-15.76	-3.90
# 	TZ [microns]	21.09	34.54	7.63
# 	RX [micro-radians]	-88.96	-144.40	-33.52
# 	RY [micro-radians]	1.06	1.75	0.36
# 	RZ [micro-radians]	0.55	0.68	0.42
# L3	TX [microns]	0.31	0.51	0.11
# 	TY [microns]	-0.44	-1.23	0.35
# 	TZ [microns]	320.74	498.74	142.73
# 	RX [micro-radians]	4.78	3.49	6.08
# 	RY [micro-radians]	-0.20	-0.33	-0.07
# 	RZ [micro-radians]	-0.06	-0.10	-0.02
# Focal Plane	TX [microns]	4.91	5.11	4.71
# 	TY [microns]	-9.11	-9.35	-8.87
# 	TZ [microns]	106.46	264.92	-51.99
# 	RX [micro-radians]	15.26	12.79	17.74
# 	RY [micro-radians]	-0.10	-0.24	0.04
# 	RZ [micro-radians]	0.03	-0.02	0.07
# ----------------------------------------------------

# tabulated values are microns (for Tx,Ty,Tz) and microradians (Rx,Ry,Rz)
# gravity distortion part:
# 	T0=23°C	Rigid Body Distortions w.r.t. CCS		
# 	Optic/DOF	GX	Gy	Gz
$tmp={};$L1{GD,Tx}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(345.06,5.72,-0.47);
$tmp={};$L1{GD,Ty}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(-1.74,-396.10,-7.22);
$tmp={};$L1{GD,Tz}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(0.06,-3.00,-60.09);
$tmp={};$L1{GD,Rx}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(-0.33,-36.94,0.09);
$tmp={};$L1{GD,Ry}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(21.70,1.62,0.85);
$tmp={};$L1{GD,Rz}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(-93.82,-2.76,0.12);
$tmp={};$L2{GD,Tx}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(374.93,6.85,0.07);
$tmp={};$L2{GD,Ty}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(-1.26,-396.57,-7.09);
$tmp={};$L2{GD,Tz}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(0.58,-4.27,-93.26);
$tmp={};$L2{GD,Rx}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(6.46,-30.62,-7.30);
$tmp={};$L2{GD,Ry}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(25.95,-0.98,-9.53);
$tmp={};$L2{GD,Rz}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(0.00,0.00,0.00); # reference point.
$tmp={};$F{GD,Tx}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(-80.18,-364.48,9.84);
$tmp={};$F{GD,Ty}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(0.37,-7.92,-66.78);
$tmp={};$F{GD,Tz}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(85.72,14.85,0.04);
$tmp={};$F{GD,Rx}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(26.63,-13.32,-0.44);
$tmp={};$F{GD,Ry}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(90.94,579.64,-12.46);
$tmp={};$F{GD,Rz}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(-6.46,-69.69,42.47);
$tmp={};$L3{GD,Tx}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(22.52,0.00,-0.07);
$tmp={};$L3{GD,Ty}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(0.04,-22.31,0.12);
$tmp={};$L3{GD,Tz}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(0.01,-0.08,-11.37);
$tmp={};$L3{GD,Rx}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(0.05,-5.43,0.02);
$tmp={};$L3{GD,Ry}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(-6.43,0.00,0.07);
$tmp={};$L3{GD,Rz}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(0.05,0.00,0.00);
$tmp={};$FP{GD,Tx}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(154.67,0.00,-0.07);
$tmp={};$FP{GD,Ty}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(0.02,-111.44,-0.13);
$tmp={};$FP{GD,Tz}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(0.01,-0.81,-15.39);
$tmp={};$FP{GD,Rx}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(0.06,-11.23,-4.73);
$tmp={};$FP{GD,Ry}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(-7.24,0.02,0.22);
$tmp={};$FP{GD,Rz}=$tmp;
@{$tmp}{Gx,Gy,Gz}=(-33.20,0.04,0.03);
# ---------------------------------------------------
# thermal distortion part
# 	T0=23C	Rigid Body Distortions w.r.t. CCS		
# 	Optic/DOF	Case 1 Nominal	Case 2 Cold	Case 3 Hot
$tmp={};$L1{TD,Tx}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(-1.32,-2.04,-0.59);
$tmp={};$L1{TD,Ty}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(32.47,53.13,11.80);
$tmp={};$L1{TD,Tz}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(403.08,655.14,151.01);
$tmp={};$L1{TD,Rx}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(29.22,47.59,10.86);
$tmp={};$L1{TD,Ry}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(-0.89,-1.41,-0.36);
$tmp={};$L1{TD,Rz}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(-1.56,-2.58,-0.55);
$tmp={};$L2{TD,Tx}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(-1.83,-2.87,-0.80);
$tmp={};$L2{TD,Ty}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(17.98,29.53,6.42);
$tmp={};$L2{TD,Tz}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(277.09,450.42,103.77);
$tmp={};$L2{TD,Rx}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(29.16,47.49,10.84);
$tmp={};$L2{TD,Ry}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(-1.04,-1.66,-0.42);
$tmp={};$L2{TD,Rz}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(0.00,0.00,0.00);
$tmp={};$F{TD,Tx}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(-0.04,-0.10,0.02);
$tmp={};$F{TD,Ty}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(-9.83,-15.76,-3.90);
$tmp={};$F{TD,Tz}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(21.09,34.54,7.63);
$tmp={};$F{TD,Rx}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(-88.96,-144.40,-33.52);
$tmp={};$F{TD,Ry}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(1.06,1.75,0.36);
$tmp={};$F{TD,Rz}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(0.55,0.68,0.42);
$tmp={};$L3{TD,Tx}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(0.31,0.51,0.11);
$tmp={};$L3{TD,Ty}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(-0.44,-1.23,0.35);
$tmp={};$L3{TD,Tz}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(320.74,498.74,142.73);
$tmp={};$L3{TD,Rx}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(4.78,3.49,6.08);
$tmp={};$L3{TD,Ry}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(-0.20,-0.33,-0.07);
$tmp={};$L3{TD,Rz}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(-0.06,-0.10,-0.02);
$tmp={};$FP{TD,Tx}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(4.91,5.11,4.71);
$tmp={};$FP{TD,Ty}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(-9.11,-9.35,-8.87);
$tmp={};$FP{TD,Tz}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(106.46,264.92,-51.99);
$tmp={};$FP{TD,Rx}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(15.26,12.79,17.74);
$tmp={};$FP{TD,Ry}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(-0.10,-0.24,0.04);
$tmp={};$FP{TD,Rz}=$tmp;
@{$tmp}{Nom,Cold,Hot}=(0.03,-0.02,0.07);

use POSIX;

$deg=atan2(1,1)/45.0;
$theta=0;
$phi=45;

$choose_thermal="Nom";


# locate arguments on the command line

$enable_G=0;
$enable_T=0;

while ($_=$ARGV[0],shift) {
    if (/^-/) {
	/^-G/ &&     ($enable_G=1,next);
	/^-theta/ && ($theta=$ARGV[0],shift,next);
	/^-phi/ &&   ($phi=$ARGV[0],shift,next);
	/^-T/ &&     ($enable_T=1,$choose_thermal=$ARGV[0],shift,next);
    } else {
	die "didn't expect arg: $_\n";
    }
}

if ($choose_thermal eq "Nom" or
    $choose_thermal eq "Cold" or
    $choose_thermal eq "Hot") {
} else {
    die "allowed values for thermal model: -T (Nom|Cold|Hot)\n";
}

if ($enable_G==1) {
    @Gvec=(sin($theta*$deg)*cos($phi*$deg),
	   sin($theta*$deg)*sin($phi*$deg),
	   cos($theta*$deg));
} else {
    @Gvec=(0,0,0);
}

printf "theta $theta phi $phi gvec %s\n",join(':',@Gvec);
printf "remaining args: %s\n",join(' ',@ARGV);

# for each optic combine gravity and thermal distortion effects by superposing

%output_filenames=("L1" => "L1solidbody_g.dat",
		   "L2" => "L2solidbody_g.dat",
		   "F"  => "Fsolidbody_g.dat",
		   "L3" => "L3solidbody_g.dat",
		   "FP" => "FPsolidbody_g.dat");
foreach $optic ("L1","L2","F","L3","FP") {
    %distort=();

    if ($enable_T==1) {
	foreach $mode ("Tx","Ty","Tz") {
	    $distort{$mode}=1e-3*${${$optic}{"TD",$mode}}{$choose_thermal};
	}
	foreach $mode ("Rx","Ry","Rz") {
	    $distort{$mode}=1e-6*${${$optic}{"TD",$mode}}{$choose_thermal};
	}
    } else {
	foreach $mode ("Tx","Ty","Tz","Rx","Ry","Rz") {
	    $distort{$mode}=0;
	}
    }

    foreach $mode ("Tx","Ty","Tz") {
	$distort{$mode}+=1e-3*($Gvec[0]*${${$optic}{"GD",$mode}}{"Gx"}+
			       $Gvec[1]*${${$optic}{"GD",$mode}}{"Gy"}+
			       $Gvec[2]*${${$optic}{"GD",$mode}}{"Gz"});
    }
    foreach $mode ("Rx","Ry","Rz") {
	$distort{$mode}+=1e-6*($Gvec[0]*${${$optic}{"GD",$mode}}{"Gx"}+
			       $Gvec[1]*${${$optic}{"GD",$mode}}{"Gy"}+
			       $Gvec[2]*${${$optic}{"GD",$mode}}{"Gz"});
    }

    printf STDOUT ("optic $optic:\n");
    printf STDOUT ("-T -t %-7.4g %-7.4g %-7.4g -r %-7.4g %-7.4g %-7.4g --\n",
		   @distort{Tx,Ty,Tz,Rx,Ry,Rz});
    open(F,">$output_filenames{$optic}") || die;
    printf F ("-T -t %-7.5g %-7.5g %-7.5g -r %-7.5g %-7.5g %-7.5g --\n",
	      @distort{Tx,Ty,Tz,Rx,Ry,Rz});
    close(F);
}

printf STDOUT "done! new files here are:\n%s\n",join("\n",values %output_filenames);

exit;
foreach $optic ("L1","L2","F","L3","FP") {
    foreach $mode ("Tx","Ty","Tz","Rx","Ry","Rz") {
	foreach $dist ("GD","TD") {
	    if ($dist eq "GD") {
		@effects=("Gx","Gy","Gz");
	    } else {
		@effects=("Nom","Cold","Hot");
	    }
	    foreach $effect ( @effects ) {
		printf "optic $optic mode $mode dist $dist effect $effect: %g\n",${${$optic}{$dist,$mode}}{$effect};
	    }
	}
    }
}

