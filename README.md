# LonelyPlanets
AMUSE scripts for integrating planetary systems in a stellar cluster
Written by Simon Portegies Zwart (Leiden, 2023)

The idea of LonelyPlanets is based on the manuscript Cai, M.X.,
Portegies Zwart, S.F., Kouwenhoven, M.B., Spurzem, R., 2019, MNRAS
489, 4311 (2019MNRAS.489.4311C), On the survivability of planets in
young massive clusters and its implication of planet orbital
architectures in globular clusters.

Here we present an independent implementation which works on any
N-body integrators and stellar evolution code in AMUSE, and the effect
of the Galactic tidal field is included.

The LonelyPlanets script can be used in conjunction with any of the
N-body codes and stellar evolution codes implemented in AMUSE, it
includes stellar evolution, collisions, and a semi-analytic Galactic
potential.

For installation instructions, read the file INSTALL

For the licence, read the file LICENCE

Running LonelyPlanets is done in two steps.

Make a production-run environment.

#Stage 1:

Generate Initial conditions:
Initial conditions can be generated by any AMUSE script.  The simplest
would be to generate them from a simple cluster model.
Make a directory for the initial conditions
  %> mkdir ICs
  %> python make_initial_cluster.py -N 1024 -Rvir 1.0 -W 7

Generates a virialized King model distribution of stars with a mass
function, and writes the results to a file called: initial_cluster.amuse


# Run Stage 1:
In your production environment, generate the stage 1 directory
  %> mkdir S1

You can only run with an initial amuse file containing the cluster.
Input file: <input snapshot amuse file>
for example: ../ICs/initial_cluster.amuse
Log output filename: <log-outputfilename>
for example: cluster_S1.log

Run the code via
  %> python <LonelyPlanets source directory>run_cluster.py -f <input snapshot amuse file>  --Nnn 6 -t 1000 >& <log-outputfilename>

runs the script with 6 nearest neighbors for 1000Myr

# Stage 2:

In your production environment, generate the stage 2 directory
  %> mkdir S2

and run the script
  %> python <LonelyPlanets source directory>run_stageII_with_Galaxy.py --Nast 1000 -t 1000 -f <input S1 file>& <log-outputfilename>

You can only run after Stage 1, and with the star that was followed
identified. This file is located in your S1 directory.  The 
star-file can be the following:
lps_key_9974767140939965756.amuse
Here 9974767140939965756 is the key that identifies the star in Stage 1.

Input file: <input S1 file>
for example: ../S1/lps_key_9974767140939965756.amuse
Log output filename: <log-outputfilename>
for example: cluster_S2_key_9974767140939965756.log


