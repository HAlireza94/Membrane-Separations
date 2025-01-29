This is a program aiming to calculate Density across X axis in a simulation box.
The input file contains coordinates of water extracted from trajectory.
To do so, it is recommended to use gmx trjconv -f trj.xtc -s trj.tpr -pbc mol -o working.gro .
ALso, the output file will be written in a Density.xvg.
The simulation details also should be included manually in the main.cpp before compiling.
