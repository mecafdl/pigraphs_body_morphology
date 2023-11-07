% Change directory to active file    
cd(fileparts(which(mfilename)));

% Add paths to function and toolbox directories
addpath(genpath(fullfile(pwd,"functions")))
addpath(genpath(fullfile(pwd,"data")))
addpath(genpath(fullfile(pwd,"parameters")))
addpath(genpath(fullfile(pwd,"external_toolboxes","infodynamics-dist-1.5")))
addpath(genpath(fullfile(pwd,"external_toolboxes","mi.0.912")))
addpath(genpath(fullfile(pwd,"external_toolboxes","minimum_bounding_spheres")))
addpath(genpath(fullfile(pwd,"external_toolboxes","Manopt_5.0")))
addpath(genpath(fullfile(pwd,"external_toolboxes","npy-matlab")))
addpath(genpath(fullfile(pwd,"external_toolboxes","BCT")))

% Add JIDT jar library to the path, and disable the warnings indicating
% that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath(fullfile(pwd,"external_toolboxes",'infodynamics-dist-1.5','infodynamics.jar'));
    
% Add utilities to the path
addpath(fullfile(pwd,"external_toolboxes",'infodynamics-dist-1.5','demos','octave'));

clc
