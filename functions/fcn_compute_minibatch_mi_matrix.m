% ########################################################################
% 
% 
%                    Mutual Information calculation                       
% 
% 
% #########################################################################

function MI_contracted = ...
               fcn_compute_minibatch_mi_matrix(miniBatch_centered, method, constPar)
    
    if strcmp(method,'jidt')
        % JIDT settings ===================================================        
        % Add JIDT jar library to the path, and disable the warnings indicating
        % that it's already there:
        warning('off','MATLAB:Java:DuplicateClass');
        javaaddpath(fullfile(fcn_get_matlab_files_dir,'external_toolboxes','infodynamics-dist-1.5','infodynamics.jar'));
        
        % Add utilities to the path
        addpath(fullfile(fcn_get_matlab_files_dir,'external_toolboxes','infodynamics-dist-1.5','demos','octave'));
        
        % 1. Construct the calculator:
        %calc = javaObject('infodynamics.measures.continuous.kernel.MutualInfoCalculatorMultiVariateKernel');
        calc = javaObject('infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian');
                
        % 2. Set any properties to non-default values:
        KERNEL_WIDTH = 0.8;
        calc.setProperty('KERNEL_WIDTH', num2str(KERNEL_WIDTH));
        calc.setProperty('NORMALISE', 'false');
        calc.setProperty('k','12');
%         cprintf('*yellow', '>> Computing MI using JIDT...\n')
    else
%         cprintf('*yellow', '>> Computing MI using MEX...\n')
    end
    
    % Load/prepare the data, change to row vectores required to use JIDT
    signal_indices = [1:3*constPar.noj, ...
                      3*constPar.noj+1:3*constPar.noj + 3*constPar.nob, ...
                      3*constPar.noj + 3*3*constPar.nob + 1:3*constPar.noj + 4*3*constPar.nob]; 

    miniBatch_centered = zscore(transpose(miniBatch_centered(signal_indices,:)));
        N_signals          = size(miniBatch_centered,2);
    MI_miniBatch       = zeros(N_signals,N_signals);
    % Compute for all pairs:
    for s = 1:1:size(MI_miniBatch,1)
        % Select source signal
        source = miniBatch_centered(:, s);
        for d = s:1:size(MI_miniBatch,2)
            if s==d
                continue;
            end
            % Select destination signal
            destination = miniBatch_centered(:, d);
            if strcmp(method, 'jidt')
                % 3. Initialise the calculator for (re-)use:
                calc.initialise();
                % 4. Supply the sample data:
                calc.setObservations(source, destination);
                % 5. Compute the estimate:
                result = max(0,calc.computeAverageLocalOfObservations());
            else
                result = max(0,mutualinfo(source, destination));
            end
            MI_miniBatch(s,d) = result;
            MI_miniBatch(d,s) = result;
        end
    end
    MI_contracted = fcn_contract_mi_matrix(MI_miniBatch, constPar);
end