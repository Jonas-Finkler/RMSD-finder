# RMSD-finder
Computes the tranlation, rotation and permutation of atoms that minimize the **RMSD** between two atomic configurations. 

The **RMSD** is defined as:

<img src="https://render.githubusercontent.com/render/math?math=\mathbf{\mathrm{RMSD}}=\sqrt{\frac{\sum_{i=1}^{N} \Vert R_i - r_i\Vert^2}{N}}">

The **mass-weighted RMSD** is defined as:

<img src="https://render.githubusercontent.com/render/math?math=\mathbf{\mathrm{RMSD}}=\sqrt{\frac{\sum_{i=1}^N  m_i \Vert R_i - r_i\Vert^2}{N}}">


The optimal translation is found trivially by superimposing the centers of mass of the two structures. 
Given a permuation, the optimal rotation can be found using an algorithm based on quaternions. 
Given a rotation, the optimal permutation is found using the Hungarian algorithm. 
To solve the joint problem a set of initial rotations is used for each of which the optimal permuation and the rotation are determined iteratively until a converged solution is found. The lowest **RMSD** obtained out of all iterative optimizations is returned by the program.

It can not be guaranteed, that a _globally_ optimal solution is found.
Wether such a globally optimal solution is found depends crucially on the set of inital rotation that is used. 
To minimize the numbers of rotations that have to be tried (and computation time) as well as the chance of not finding the truly optimal **RMSD**, lists of rotations that are evenly distributed over the space of all rotations are included in this code. 
Lists are included for N = 100 to 2000 in steps of 100 and can be used with the command line flag `-n<num_rotations>`. More lists can be generated with the program `createQuaternionLists.f90`.

In general the higher the **RMSD** between two structures and the more atoms are contained, the harder it is to find the optimal solution and hence more rotations are required.  



## Usage
`rmsdFinder -A<structure-A> -B<structure-B> [-n<num-rotations>] [-i<use-inversion>] [-o<out-structure>] [-m<use-mass-weighted-RMSD>] [-r<rotation-format>]`

`<structure-A>` is translated, rotated and permuted to minimize the RMSD to `<structure-B>`
All structures must be provided in xyz coordinates. 
       
`<num-rotations>` defines the number of rotations that are tried. The default is 2000. 
Possible values are 100, 200, ..., 1900, 2000.
A higher number increases the computation time but also the probability
that the true minimal RMSD is found.
       
`<use-inversion>` can either be **T** (true) or **F** (false). 
If true, inverted structures are also considered when searching the minimal RMSD.
       
`<out-structure>` is the file to which the transformed copy of `<structure-A>` is saved
   
`<use-mass-weighted-RMSD>` can either be **T** (true) or **F** (false).
If true, the mass weighted RMSD definition is used.
       
`<rotation-format>` specifies the format in which the rotation is printed.
Can be either **Q** for quaternion (default), **M** for rotation Matrix or **A** for angle-axis representation.

## Assignment convention
Atom number assignment(i) of structure A has been matched to atom i in structure B.


## Compiling the code

To compile the example you need a fortran compiler as well as an installation of Blas/LAPACK. 
Using CMake the code can be compiled with the following commands.
Two CMake flags are provided, that allow to compile with or without debug flags and with the intel or gnu compiler.

```bash
mkdir build
cd build
cmake -DDEBUG=OFF -DINTEL=ON .. # compile without debug flags and the intel fortran compiler
make
```


# Reference
An explanation of the algorithm can be found in the following paper. 
Please cite this paper if you use this code in an academic setting.


_The Journal of Chemical Physics_ 152.16 (2020): 164106.   
Finkler, Jonas A., and Stefan Goedecker. 
"Funnel hopping Monte Carlo: An efficient method to overcome broken ergodicity."
<https://doi.org/10.1063/5.0004106>



