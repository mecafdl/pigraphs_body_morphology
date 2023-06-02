
This is the accompanying repository for the manuscript titled `Close Your Eyes and Tell Me What You See: From Proprioceptive Sensing to Robot Morphology`. The repository contains:

- Links to data sets and 
- MATLAB source code 

# Requirements

 - MATLAB's Robotics System Toolbox
 - Tool­boxes for opti­mization on manifolds and matrices  [MANOPT](https://www.manopt.org/)
 - Java Information Dynamics Toolkit [JIDT](https://github.com/jlizier/jidt)

## Getting the datasets
All datasets are publicly available at Kaggle; these are the corresponding links:

 - Simulated robot manipulator with moving base [here](https://www.kaggle.com/datasets/fernandodazledezma/franka-proprioception-simulated-moving-base)
 - Simulated robot manipulator with fixed base [here](https://www.kaggle.com/datasets/fernandodazledezma/frankaproprioceptionsimulatedfixedbase)
 - Simulated hexpod robot [here](https://www.kaggle.com/datasets/fernandodazledezma/phantomx-proprioception)

## Usage
Versions of the toolboxes used are included in this repository, check if you want to udate them. To reproduce the results:

1. Add directories to search path, in MATLAB's command line: `run init_path.m`

2. To see the results for the robot manipulator open `franka.m`. The different sections of the code are self-explanatory

3. To see the results for the robot manipulator open `phantomx.m`. The different sections of the code are self-explanatory