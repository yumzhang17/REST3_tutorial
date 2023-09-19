# REST3_tutorial
----------------------------------------

This is a simple tutorial to illustrate how to configure and run REST3 simulations. It's based on the REST2 implementation in GROMACS patched with PLUMED 2.
(see https://www.plumed.org/doc-v2.7/user-doc/html/hrex.html
for more details about how to perform REST2 simulations in GROMACS).
Users of this tutorial are assumed to be famililar with GROMACS usage and REST2 algorithms. 

Compared with REST2, an additional parameter is introduced in REST3,
the scaling factors of solute−solvent van der Waals interactions.
This is to control the solute conformational
properties at different effective temperatures for maximizing sampling
efficiency. See the following reference for more details:

Reference: Yumeng Zhang, Xiaorong Liu, and Jianhan Chen, "Re-balancing Replica Exchange with Solute Tempering for Sampling Dynamic Protein Conformations", JCTC, 19, 1602-1614 (2023) 

All other details can be found in Jupytor notebook "rest3.ipynb".

Contact: Jianhan Chen (jianhanc@umass.edu)
