#!/usr/bin/perl
use strict;
use warnings;
use POSIX;

# script to populate missing elements from distortion FEA calculations
# these are primarily solid body misalignments recorded by John Ku but without
# node specific displacements & rotations shipped out by him.
# This script repeatedly calls the script named "distortion_misalign_gen.perl"
# and populates the directories that already exist (generating them if 
# necessary) and organizes them.

# first check for the database directory

my $ThisScript="patchup_distortions.perl";

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
my $raytrace_coords=0;
my $script;
$script="distortion_misalign_gen.perl -rt" if ($raytrace_coords==1);
$script="distortion_misalign_gen.perl"     if ($raytrace_coords==0);

my %runset=("F/gravity/"=>"ha");

my @commands;
my @bodies;
my @naming;
my @task;
my $run=-1;

$run++;
$commands[$run]="-G -theta 0 -phi 0";
$bodies[$run]=[ "F","FP" ];
$task[$run]="gravity";
$naming[$run]=
    ($raytrace_coords==1)?"GZ_rt_coords_full_":"GZ_ccs_coords_full_";

$run++;
$commands[$run]="-G -theta 90 -phi 0";
$bodies[$run]=[ "F","FP" ];
$task[$run]="gravity";
$naming[$run]=
    ($raytrace_coords==1)?"GY_rt_coords_full_":"GY_ccs_coords_full_";

$run++;
$commands[$run]="-G -theta 90 -phi 90";
$bodies[$run]=[ "F","FP" ];
$task[$run]="gravity";
$naming[$run]=
    ($raytrace_coords==1)?"GX_rt_coords_full_":"GX_ccs_coords_full_";

$run++;
$commands[$run]="-T Cold";
$bodies[$run]=[ "F","FP","L3" ];
$task[$run]="thermal";
$naming[$run]=
    ($raytrace_coords==1)?"Cold_rt_coords_full_":"Cold_ccs_coords_full_";

$run++;
$commands[$run]="-T Hot";
$bodies[$run]=[ "F","FP","L3" ];
$task[$run]="thermal";
$naming[$run]=
    ($raytrace_coords==1)?"Hot_rt_coords_full_":"Hot_ccs_coords_full_";

$run++;
$commands[$run]="-T Nom";
$bodies[$run]=[ "F","FP","L3" ];
$task[$run]="thermal";
$naming[$run]=
    ($raytrace_coords==1)?"Nom_rt_coords_full_":"Nom_ccs_coords_full_";

my $orig_suffix="solidbody_g.dat";

my $n=$run;
foreach $run (0..$n) {
    print `$script $commands[$run]`;
    foreach my $body (@{$bodies[$run]}) {
	my $src=$body.$orig_suffix;
	my $newfile=$src;
	$newfile =~ s/_g\.dat/\.dat/;
	$newfile =~ s/$body/${body}_${naming[$run]}/;
	my $destdir=join("/",$DISTORTION_DB_PATH,$body,$task[$run]);
	my $dest=join("/",$destdir,$newfile);
	printf STDERR "will try to move: $src to $dest\n";
	my @dirstruct=split("/",$destdir);
	for (my $i=0;$i<=$#dirstruct;$i++) {
	    my $tmpdir=join("/",@dirstruct[0..$i]);
	    mkdir $tmpdir if (! -d $tmpdir);
	}
	`mv $src $dest`;
    }
}
