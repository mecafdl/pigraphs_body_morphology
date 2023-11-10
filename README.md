This is the accompanying repository for the manuscript titled Machine Learning Driven Self-Discovery of the Robot Body Morphology. The repository contains:

- Links to data sets and 
- MATLAB source code 

# Requirements

 - MATLAB's Robotics System Toolbox
 - MATLAB's Optimization Toolbox
 - Tool­boxes for opti­mization on manifolds and matrices  [MANOPT](https://www.manopt.org/)
 - Java Information Dynamics Toolkit [JIDT](https://github.com/jlizier/jidt)
 
 NOTE: <u>MATLAB 2021b</u> was used.

## Getting the datasets
All datasets are publicly available at Kaggle; these are the corresponding links:

 - Simulated robot manipulator with moving base [here](https://www.kaggle.com/datasets/fernandodazledezma/franka-proprioception-simulated-moving-base)
 - Simulated robot manipulator with fixed base [here](https://www.kaggle.com/datasets/fernandodazledezma/frankaproprioceptionsimulatedfixedbase)
 - Physical manipulator experiment (fixed base) [here](https://www.kaggle.com/datasets/fernandodazledezma/frankaproprioceptionrealfixedbase)
 - Simulated hexpod robot [here](https://www.kaggle.com/datasets/fernandodazledezma/phantomx-proprioception)
 - Simulated humanoid robot [here](https://www.kaggle.com/datasets/fernandodazledezma/humanoid-proprioception)
 
 All datasets must be stored in the  <u>`data`</u>  directory.

## Usage
Versions of the toolboxes used are included in this repository, check if you want to udate them. To reproduce the results:

1. Add directories to the search path in MATLAB's command line: `run init_path.m`

2. To see the results for the simulated robot manipulator open, `frankaArmSimulated.m`. 

3. To see the results for the physical robot manipulator open, `frankaArmPhysical.m`.

4. To see the results for the hexapod robot, open `phantomxHexapod.m`.

5. To see the results for the humanoid robot, open `poppyHumanoid.m`.

The different sections in the scripts are self-explanatory.
