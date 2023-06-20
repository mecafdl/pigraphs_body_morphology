%% Filtering of the signals ------------------------------------------------
function [domg_signals_zpf, acc_signals_zpf] = fcn_getFilteredSignals(omg_signals, acc_signals)
    % Zero-phase filter (zpf)
    Fs = 100; % signal sampling frequency [Hz]
    Fc =  10; % cutoff frequency [Hz], obtained from the previous block as the average of the frequencies for all joints where the amplitude becomes negligible
    % Fc = f_cut;
    d  = designfilt('lowpassiir','FilterOrder',6, ...
        'HalfPowerFrequency',(2/Fs)*Fc,'DesignMethod','butter'); 
    %                       |_________|
    %                            |_______________ This is the Nyquist normalized frequency
    [b_zpf, a_zpf] = tf(d);
    
    % Filtering of the signals
    omg_signals_zpf     = transpose(filtfilt(b_zpf, a_zpf, transpose(omg_signals)));
    domg_signals        = gradient(omg_signals_zpf,1E-2);
    domg_signals_zpf    = transpose(filtfilt(b_zpf, a_zpf, transpose(domg_signals)));
    ang_acc_imu_num_zpf = domg_signals_zpf;
    acc_signals_zpf     = transpose(filtfilt(b_zpf,a_zpf, transpose(acc_signals)));
end