function [hexapodSignals, proprioceptiveSignals] = fcn_loadHexapodRobotData(constPar)
    % Load datasets
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    hexapodData = load([filepath,'/../data/phantomx_proprioception.mat'], 'phantomx_proprioception');
    
    % The rows of each of the datasets is organized as follows:
    % - q_indices       = 1:18;
    % - dq_indices      = 19:36;
    % - torque_indices  = 37:54;
    % - ang_vel_indices = 55:111;
    % - lin_acc_indices = 112:168;
    
    q_indices       = 1:18;
    dq_indices      = 19:36;
    torque_indices  = 37:54;
    ang_vel_indices = 55:111;
    lin_acc_indices = 112:168;
    
    %% Create the proprioceptive signals matrix
    clc
    close all
    proprioception = [];
    
    fn = fieldnames(hexapodData.phantomx_proprioception);
    for k=1:2%:numel(fn)
        if( isnumeric(hexapodData.phantomx_proprioception.(fn{k})) )
            proprioception_aux = ...
                             [hexapodData.phantomx_proprioception.(fn{k})([q_indices,dq_indices,torque_indices,ang_vel_indices],:); ...
                              NaN*hexapodData.phantomx_proprioception.(fn{k})(ang_vel_indices,:); ...
                              NaN*hexapodData.phantomx_proprioception.(fn{k})(lin_acc_indices,:); ...
                              hexapodData.phantomx_proprioception.(fn{k})(lin_acc_indices,:)];
            proprioception = [proprioception proprioception_aux];
        end
    end
    lin_acc_indices = lin_acc_indices + 2*3*constPar.nol;
    disp('done')
    
    WITH_NOISE = 1;
    if WITH_NOISE == 1
        SNR     = 20;
        proprioceptiveSignals = zeros(size(proprioception));
        for signal = 1:size(proprioceptiveSignals,1)
            proprioceptiveSignals(signal,:) = awgn(proprioception(signal,:), SNR, 'measured');
        end
    else
        proprioceptiveSignals = proprioception;
    end
    
    %% Filtering of the signals ------------------------------------------------
    
    % Zero-phase filter (zpf)
    Fs = 100; % signal sampling frequency [Hz]
    Fc =  10; % cutoff frequency [Hz], obtained from the previous block as the average of the frequencies for all joints where the amplitude becomes negligible
    d  = designfilt('lowpassiir','FilterOrder',6, ...
        'HalfPowerFrequency',(2/Fs)*Fc,'DesignMethod','butter'); 
    %                       |_________|
    %                            |_______________ This is the Nyquist normalized frequency
    [b_zpf, a_zpf] = tf(d);
    
    % Estimated Cartesian angular acceleration
    hexapodSignals.samplingTime                  = ...
        1E-2;
    hexapodSignals.jointPosition.raw             = ...
        proprioceptiveSignals(q_indices,:);
    hexapodSignals.jointPosition.zeroPhaseFilter = ...
        transpose(filtfilt(b_zpf, a_zpf, transpose(hexapodSignals.jointPosition.raw)));
    hexapodSignals.jointVelocity.raw             = ...
        proprioceptiveSignals(dq_indices,:);
    hexapodSignals.jointVelocity.zeroPhaseFilter = ...
        transpose(filtfilt(b_zpf, a_zpf, transpose(hexapodSignals.jointVelocity.raw)));
    hexapodSignals.bodyAngularVelocity.raw = ...
        proprioceptiveSignals(ang_vel_indices,:);
    hexapodSignals.bodyAngularVelocity.zeroPhaseFilter = ...
        transpose(filtfilt(b_zpf, a_zpf, transpose(hexapodSignals.bodyAngularVelocity.raw)));
    hexapodSignals.bodyAngularAcceleration.numerical = ...
        gradient(hexapodSignals.bodyAngularVelocity.zeroPhaseFilter, hexapodSignals.samplingTime);
    hexapodSignals.bodyAngularAcceleration.zeroPhaseFilter = ...
        transpose(filtfilt(b_zpf, a_zpf, transpose(hexapodSignals.bodyAngularAcceleration.numerical)));
    hexapodSignals.bodyLinearAcceleration.raw = ...
        proprioceptiveSignals(lin_acc_indices,:);
    hexapodSignals.bodyLinearAcceleration.zeroPhaseFilter = ...
        transpose(filtfilt(b_zpf,a_zpf, transpose(hexapodSignals.bodyLinearAcceleration.raw)));
    
    cprintf('*green', '>> Filtered signals and numerical gradients ready!\n')
end
