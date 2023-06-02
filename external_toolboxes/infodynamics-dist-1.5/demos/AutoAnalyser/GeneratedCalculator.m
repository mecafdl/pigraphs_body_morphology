% Add JIDT jar library to the path, and disable warnings that it's already there:
warning('off','MATLAB:Java:DuplicateClass');
javaaddpath('/home/fernando/Dropbox/PhD_WORK/matlab_irt/external_toolboxes/infodynamics-dist-1.5/infodynamics.jar');
% Add utilities to the path
addpath('/home/fernando/Dropbox/PhD_WORK/matlab_irt/external_toolboxes/infodynamics-dist-1.5/demos/octave');

% 0. Load/prepare the data:
data = load('/home/fernando/Dropbox/PhD_WORK/matlab_irt/external_toolboxes/infodynamics-dist-1.5/demos/data/panda_data.txt');
% 1. Construct the calculator:
calc = javaObject('infodynamics.measures.continuous.kernel.MutualInfoCalculatorMultiVariateKernel');
% 2. Set any properties to non-default values:
calc.setProperty('KERNEL_WIDTH', '0.8');

% Compute for all pairs:
for s = 1:42
	for d = 1:42
		% For each source-dest pair:
		if (s == d)
			continue;
		end
		% Column indices start from 1 in Matlab:
		source = octaveToJavaDoubleArray(data(:, s));
		destination = octaveToJavaDoubleArray(data(:, d));

		% 3. Initialise the calculator for (re-)use:
		calc.initialise();
		% 4. Supply the sample data:
		calc.setObservations(source, destination);
		% 5. Compute the estimate:
		result = calc.computeAverageLocalOfObservations();

		fprintf('MI_Kernel(col_%d -> col_%d) = %.4f bits\n', ...
			s, d, result);
	end
end
