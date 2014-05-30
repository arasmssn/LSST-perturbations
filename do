# install CPAN packages as necessary
. setup.bash

# setup relevant perturbation executables path

export LSST_PERTURBATIONS_DIR=`pwd`
export PATH=$PATH:$LSST_PERTURBATIONS_DIR/M1-M3_FEA
export PATH=$PATH:$LSST_PERTURBATIONS_DIR/M2_FEA
export PATH=$PATH:$LSST_PERTURBATIONS_DIR/CAMERA_FEA

# setup relevant perturbation database environmental variables

export M1M3_PERTURBATIONS_DB=$LSST_PERTURBATIONS_DIR/M1-M3_FEA/DB
export M2_PERTURBATIONS_DB=$LSST_PERTURBATIONS_DIR/M2_FEA/DB
export CAMERA_PERTURBATIONS_DB=$LSST_PERTURBATIONS_DIR/CAMERA_FEA/DB

# run an instance of M1/M3, M2 and CAMERA to produce an instance of 
# the operational perturbations (it should not be necessary to 
# regenerate intermediate files)


M1M3_perturbations.perl --theta=30 --Focus3=0 --Tb=-2 --Tb_ao=-2 --Ty=0.5 --Ty_ao=0.5 --Tr=0.5 --Tr_ao=0.5

M1M3_perturbations.perl --theta=30 --Tb=-2 --Tb_ao=-2 --Ty=0.5 --Ty_ao=0.5 --Tr=0.5 --Tr_ao=0.5

# See M1M3_example1.png for a figure produced with the outputs of these previous 2 commands

M2_perturbations.perl --theta 30 --Focus3=0 --Tz=1.0 --Tz_ao=1.0 --Tr=1.5 --Tr_ao=1.5 

M2_perturbations.perl --theta 30 --Tz=1.0 --Tz_ao=1.0 --Tr=1.5 --Tr_ao=1.5 

# See M2_example1.png for a figure produced with the outputs of these previous 2 commands

CAMERA_perturbations.perl -T 10 -G -theta 30 -phi 30

# See CAMERA_example1.png for a figure produced with the output of this command 

