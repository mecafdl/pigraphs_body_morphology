%%
%%  Java Information Dynamics Toolkit (JIDT)
%%  Copyright (C) 2012, Joseph T. Lizier
%%  
%%  This program is free software: you can redistribute it and/or modify
%%  it under the terms of the GNU General Public License as published by
%%  the Free Software Foundation, either version 3 of the License, or
%%  (at your option) any later version.
%%  
%%  This program is distributed in the hope that it will be useful,
%%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%  GNU General Public License for more details.
%%  
%%  You should have received a copy of the GNU General Public License
%%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%

% To recreate the plots in:
% J. T. Lizier, "Measuring the dynamics of information processing on a local scale in time and space",
%  in "Directed Information Measures in Neuroscience", edited by M. Wibral, R. Vicente, J. T. Lizier,
%  pp. 161-193, Springer, Berlin/Heidelberg, 2014; doi: 10.1007/978-3-642-54474-3_7

clear all;

% Set up simulation options:
cells = 10000;
timeSteps = 600;
neighbourhood = 3;
caStates = 2;

% Set up options for information dynamics analysis, and which segment of the CA to plot
measureParams.k=16;
options.saveImages = true;
options.saveImagesFormat = 'pdf';
options.plotOptions.scaleColoursToSubsetOfPlot = true;
scaleColoursToExtremesDefault = false;
options.plotOptions.scaleColoursToExtremes = scaleColoursToExtremesDefault;
% Turn up the contrast so that the small values aren't disproportionately visible: (0.15, 0.30 was good, except for separable which was better with 0.15, 0.35)
options.plotOptions.scalingMainComponent = 0.15;
options.plotOptions.scalingScdryComponent = 0.30;
options.plotOptions.gammaPower = 0.5;

%%%%%%%%%
% Examining rule 54:
options.plotOptions.plotRows = 35;
options.plotOptions.plotCols = 35;
options.plotOptions.plotStartRow = 20+20;
options.plotOptions.plotStartCol = 1+10;
options.seed = 3; % Set up the random number generator to give reproducible initial states for all measurements
if (exist('initialStates/DirectedMeasuresChapterDemo2013-initialStates.txt', 'file'))
	% A file specifying the initial state exists -- this
	%  ensures that Matlab and Octave use the same initial state
	%  (otherwise only Octave recreates the same initial state used in our chapter).
	%  (You can delete/move the initial state file if you want them generated from scratch.)
	options.initialState = load('initialStates/DirectedMeasuresChapterDemo2013-initialStates.txt');
elseif (isfield(options, 'initialState'))
	options = rmfield(options, 'initialState');
end
fprintf('\nStarting rule 54 ...\n');
fprintf('\nPlotting active info storage ...\n');
plotLocalInfoMeasureForCA(neighbourhood, caStates, 54, cells, timeSteps, 'active', measureParams, options);
options.plotRawCa = false;
fprintf('\nPress any key when ready for apparent transfer entropy j = 1 ...\n');
pause
% Use the full red scale for transfer and separable info, since we need to see the extreme negative values properly
options.plotOptions.scaleColoursToExtremes = true;
options.plotOptions.scalingScdryComponent = 0.35;
options.plotOptions.scalingMainComponent = 0.35; 
measureParams.j = 1;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 54, cells, timeSteps, 'transfer', measureParams, options);
fprintf('\nPress any key when ready for apparent transfer entropy j = -1 ...\n');
pause
measureParams.j = -1;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 54, cells, timeSteps, 'transfer', measureParams, options);
fprintf('\nPress any key when ready to apply to the next rule\n')
pause
options.plotOptions.scaleColoursToExtremes = scaleColoursToExtremesDefault; % return to default value 
options.plotOptions.scalingScdryComponent = 0.30; % return to previous value

%%%%%%%%%
% Examining rule 18: (50, 50, 20, 900 is not too bad)
options.plotRawCa = true;
options.plotOptions.plotRows = 50;
options.plotOptions.plotCols = 50;
options.plotOptions.plotStartRow = 20;
options.plotOptions.plotStartCol = 900;
options.seed = 3; % Set up the random number generator to give reproducible initial states for all measurements
if (exist('initialStates/DirectedMeasuresChapterDemo2013-initialStates.txt', 'file'))
	% A file specifying the initial state exists -- this
	%  ensures that Matlab and Octave use the same initial state
	%  (otherwise only Octave recreates the same initial state used in our chapter).
	%  (You can delete/move the initial state file if you want them generated from scratch.)
	options.initialState = load('initialStates/DirectedMeasuresChapterDemo2013-initialStates.txt');
elseif (isfield(options, 'initialState'))
	options = rmfield(options, 'initialState');
end
fprintf('\nStarting rule 18 ...\n');
fprintf('\nPlotting active info storage ...\n');
options.plotOptions.scalingScdryComponent = 0.45; % Make the moderately strong values easier to see:
plotLocalInfoMeasureForCA(neighbourhood, caStates, 18, cells, timeSteps, 'active', measureParams, options);
options.plotRawCa = false;
fprintf('\nPress any key when ready for apparent transfer entropy j = -1 ...\n');
pause
% Use the full red scale for transfer and separable info, since we need to see the extreme negative values properly
options.plotOptions.scaleColoursToExtremes = true;
measureParams.j = -1;
plotLocalInfoMeasureForCA(neighbourhood, caStates, 18, cells, timeSteps, 'transfer', measureParams, options);
fprintf('\nPress any key when ready for complete transfer entropy j = -1 ...\n');
pause
plotLocalInfoMeasureForCA(neighbourhood, caStates, 18, cells, timeSteps, 'transfercomplete', measureParams, options);
options.plotOptions.scaleColoursToExtremes = scaleColoursToExtremesDefault; % return to default value 
options.plotOptions.scalingScdryComponent = 0.30; % return to previous value
fprintf('\nAll done, press any key to continue ...\n');
pause

