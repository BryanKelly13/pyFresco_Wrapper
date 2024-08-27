This is a python input file wrapper that allows easy change to FRESCO input cards.

Change the input_generator.inp file to a given excited state -> Excitation energy, spin, and preffered trasnfer configuration. (e.g. for a G.S Jpi of 7/2- to a excited 3+ 4.172 MeV state via a 2p3/2 transfer configuration for a (d,p) reaction),
this example is given with the files in the folder. This will then change the input card that fresco reads in accordingly, run fresco with a unique naming scheme for the .fri & .fro files, and then save the fort.16 files
(which contains the elastic and inelastic channels cross sections) as a unique name as well, here as 50Ti_4172_2p32_3+.txt.

*** WARNING ***
The input fresco file is hard coded for optical potentials for my 49Ti(d,p)50Ti reaction, therefore much of this needs to be changed in the fresco_input.py file for your 
use (optical potentials, masses of nuclei, charges, spin & parities of target nucleus, etc.)

Use your favorite optical potentials and mass evaluation tables and change accordingly.
This wrapper was written solely to aviod the spacing issues that fresco input cards rely heavily on, so that the things that change frequently are done so
in an automated manner and allows you to run calculations for many states at once.
