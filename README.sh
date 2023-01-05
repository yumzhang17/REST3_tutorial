#/bin/bash

gmx=/home/zgjia/Software/gromacs/gromacs514GPU_plumed/bin/gmx_mpi

# generate processed.top

$gmx grompp -f rest.mdp -c prot.gro -p prot.top -pp processed.top

# **extra step** before commenting the following exit: modify the processed_ext.top
# 1) define the hot atoms that must have a "_" appended to the atom type,
# please refer to the link: https://www.plumed.org/doc-v2.7/user-doc/html/hrex.html
# 2) remove the unnecessary the lines [ nonbond_params ] if any,
# here, this modified c36mr-tip4pd force field will generate two lines with the [ nonbond_params ],
# just delete one. If it only has one, then please ignore this step.

exit

# scaled topology files using REST3 paprameters for 8 replicas
bash gennewtop.sh 1.000 1.000 processed.top topol0.top
bash gennewtop.sh 0.943 1.000 processed.top topol1.top
bash gennewtop.sh 0.889 1.003 processed.top topol2.top
bash gennewtop.sh 0.838 1.010 processed.top topol3.top
bash gennewtop.sh 0.790 1.020 processed.top topol4.top
bash gennewtop.sh 0.745 1.027 processed.top topol5.top
bash gennewtop.sh 0.702 1.035 processed.top topol6.top
bash gennewtop.sh 0.662 1.045 processed.top topol7.top

exit

# create tpr files: topol[0-7].top
$gmx grompp -maxwarn 3 -o topol0.tpr -f rest.mdp -p topol0.top -c prot.gro -r prot.gro
$gmx grompp -maxwarn 3 -o topol1.tpr -f rest.mdp -p topol1.top -c prot.gro -r prot.gro
$gmx grompp -maxwarn 3 -o topol2.tpr -f rest.mdp -p topol2.top -c prot.gro -r prot.gro
$gmx grompp -maxwarn 3 -o topol3.tpr -f rest.mdp -p topol3.top -c prot.gro -r prot.gro
$gmx grompp -maxwarn 3 -o topol4.tpr -f rest.mdp -p topol4.top -c prot.gro -r prot.gro
$gmx grompp -maxwarn 3 -o topol5.tpr -f rest.mdp -p topol5.top -c prot.gro -r prot.gro
$gmx grompp -maxwarn 3 -o topol6.tpr -f rest.mdp -p topol6.top -c prot.gro -r prot.gro
$gmx grompp -maxwarn 3 -o topol7.tpr -f rest.mdp -p topol7.top -c prot.gro -r prot.gro

# rm files
rm \#*
