[See illustration video here.](https://youtu.be/RAdbk1a0JFY)


This repository contains:

- Links to data sets and MATLAB source code to reproduce the results shown in Tools for running artificial skin simulation on the manuscript titled `Close Your Eyes and Tell Me What You See: From Proprioceptive Sensing to Robot Morphology` 


# Requirements

 - MATLAB's Robotics Systems Toolbox
 - Tool­boxes for opti­mization on manifolds and matrices  [MANOPT](https://www.manopt.org/)
 - Java Information Dynamics Toolkit [JIDT](https://github.com/jlizier/jidt)

## Getting the datasets
All datesets are publicly available at Kaggle, you will have to download by yourself.

 - Simulated robot manipulator with moving base [here](https://www.kaggle.com/datasets/fernandodazledezma/franka-proprioception-simulated-moving-base)
 - Simulated robot manipulator with fixed base [here](https://www.kaggle.com/datasets/fernandodazledezma/frankaproprioceptionsimulatedfixedbase)
 - Simulated hexpod robot [here](https://www.kaggle.com/datasets/fernandodazledezma/phantomx-proprioception)

## Usage
Versions of the toolboxes used are included in this repo, check if you want to udate them. To reproduce the results:

1. Add directories to search path, in MATLAB's command line: `run init_path.m`

2. To see the results for the robot manipulator open `franka.m`. The different sections of the code are self-explanatory

3. To see the results for the robot manipulator open `phantomx.m`. The different sections of the code are self-explanatory