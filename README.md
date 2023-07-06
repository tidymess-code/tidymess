Welcome to TIDYMESS!

TIDYMESS stands for "TIdal DYnamics of Multi-body ExtraSolar Systems". 
It is an N-body code with the following novel features: 

- the mass distribution of each body is described by its inertia tensor
- the deformation of a body has contributions from the centrifugal force and tidal force
- the deformation follows the tidal Creep model 
- the orbits, spins and inertia tensors are integrated directly and self-consistently

Other physical phenomena included:

- general relativity corrections up to 2.5 PN
- stellat magnetic braking
- mergers and collisions

Further details can be found in the code announcement paper:

Boekholt & Correia (MNRAS 2023, vol. 522, pp. 2885–2900).
 
TIDYMESS is a C++ code, developed and tested within the Ubuntu and Mac environments.
In order to get started: 

- simpy open a terminal, enter the tidymess directory, and compile the code using "make clean" and "make"
- the compilation will produce the executable (tidymess.exe) and the "build" folder inside the "integrator" folder
- if there are compilation errors due to bugs within the code, please contact us with the error message

- by typing "./tidymess.exe -h" or "./tidymess.exe --help" some documentation is presented on how to run the code 
- a simulation is defined by specifying two files: parameter file (e.g. tidymess.par) and initial condition file (e.g. tidymess.ic) 
- we provide a few examples in the "examples" folder, which can be used to get you on your way with TIDYMESS.

For further questions or comments, please feel free to contact us at tjardaboekholt@gmail.com, acm.correia@gmail.com, or 
through the TIDYMESS GitHub repository.

Boa Sorte!


Tjarda Boekholt and Alexandre Correia
Mar. 2023

