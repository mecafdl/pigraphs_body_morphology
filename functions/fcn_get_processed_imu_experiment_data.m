function armIMUExperiment = fcn_get_processed_imu_experiment_data(armIMUExperiment)
    % The INITIAL content of armIMUExperiment is as follows
    % signals_imu = [q; dq; ddq; omg; vel; tau; acc];
    %
    % * NOTE: measurements of the joint acceleration and link Cartesian velo-
    %         city are not available

    fprintf('>> Processing experiment signals...\n')
    
    armIMUExperiment.samples      = numel(armIMUExperiment.time); 
    armIMUExperiment.samplingTime = armIMUExperiment.time(2)-armIMUExperiment.time(1); 
    
    % Zero-Phase Filtering of the signals =====================================
    Fs = 1E3; % signal sampling frequency [Hz]
    Fc =  2;  % cutoff frequency [Hz], obtained from the previous block as the average of the frequencies for all joints where the amplitude becomes negligible
    % Fc = f_cut;
    d  = designfilt('lowpassiir','FilterOrder',6, ...
        'HalfPowerFrequency',(2/Fs)*Fc,'DesignMethod','butter'); 
    %                       |_________|
    %                            |_______________ This is the Nyquist normalized frequency
    
    [b_zpf, a_zpf] = tf(d);
    
    for snsr = 1:7
        % Filtering of body linear acceleration
        signal                                        = armIMUExperiment.bodyLinearAcceleration.raw((snsr-1)*3+1:snsr*3,:);
        armIMUExperiment.bodyLinearAcceleration.zeroPhaseFilter((snsr-1)*3+1:snsr*3,:) = transpose(filtfilt(b_zpf, a_zpf, signal'));
        % Filtering of body angular velocity
        signal                                        = armIMUExperiment.bodyAngularVelocity.raw((snsr-1)*3+1:snsr*3,:);
        armIMUExperiment.bodyAngularVelocity.zeroPhaseFilter((snsr-1)*3+1:snsr*3,:) = transpose(filtfilt(b_zpf, a_zpf, signal'));
    end

    % Filtering of joint velocity
    armIMUExperiment.jointVelocity.zeroPhaseFilter = transpose(filtfilt(b_zpf, a_zpf, armIMUExperiment.jointVelocity.raw'));    
    
    % Numerical differentiation of the body angular velocity ==============
    armIMUExperiment.jointAcceleration.numerical       = gradient(armIMUExperiment.jointVelocity.zeroPhaseFilter,armIMUExperiment.samplingTime);
    armIMUExperiment.bodyAngularAcceleration.numerical = gradient(armIMUExperiment.bodyAngularVelocity.zeroPhaseFilter,armIMUExperiment.samplingTime);
    
% % Third order approximation (toa) of the Cartesian angular velocity derivative
% N_joints = 7;
% imu_experiment.domg_toa = zeros(3*N_joints,armIMUExperiment.samples);
% for i = 1:armIMUExperiment.samples-2
%     if i>2
%         imu_experiment.domg_toa(:,i) = ...
% 	                (imu_experiment.omg_zpf(:,(i-2)) ...
%                      - (N_joints+1)*imu_experiment.omg_zpf(:,(i-1)) ...
%                      + (N_joints+1)*imu_experiment.omg_zpf(:,(i+1)) ...
%                      - imu_experiment.omg_zpf(:,(i+2)))/(12*1E-3);
%     end
% end
    
    % Moving average filtering (maf) of the signals ===========================
    windowSize             = 100; 
    b_maf                  = (1/windowSize)*ones(1,windowSize);
    a_maf                  = 1;
    armIMUExperiment.bodyAngularVelocity.movAvgFilter    = filter(b_maf,a_maf,armIMUExperiment.bodyAngularVelocity.raw,[],2);
    armIMUExperiment.bodyLinearAcceleration.movAvgFilter = filter(b_maf,a_maf,armIMUExperiment.bodyLinearAcceleration.raw,[],2);
    
    % q_ref           = imu_experiment.q;
    % dq_ref          = imu_experiment.dq;
    % ddq_ref         = gradient(imu_experiment.dq_zpf,1E-3);
    
    % Add fictitious base link measurements
    armIMUExperiment.bodyLinearVelocity.numerical = ...
                                      [zeros(3,armIMUExperiment.samples);
                                       cumtrapz(armIMUExperiment.bodyLinearAcceleration.zeroPhaseFilter)];
    
    armIMUExperiment.bodyAngularVelocity.raw = ...
                                      [zeros(3, armIMUExperiment.samples);
                                       armIMUExperiment.bodyAngularVelocity.raw];
    armIMUExperiment.bodyLinearAcceleration.raw = ...
                                      [[0; 0; 9.8100]*ones(1,armIMUExperiment.samples);
                                       armIMUExperiment.bodyLinearAcceleration.raw];
%     armIMUExperiment.bodyAngularAcceleration.Numerical = ...
%                                       [zeros(3, armIMUExperiment.samples);
%                                        armIMUExperiment.bodyAngularAcceleration.Numerical];                                
                                  
    % Measurements without noise
    armIMUExperiment.bodyAngularVelocity.zeroPhaseFilter = ...
                                      [zeros(3, armIMUExperiment.samples);
                                       armIMUExperiment.bodyAngularVelocity.zeroPhaseFilter];

    armIMUExperiment.bodyAngularVelocity.movAvgFilter = ...
                                      [zeros(3, armIMUExperiment.samples);
                                       armIMUExperiment.bodyAngularVelocity.movAvgFilter];

    armIMUExperiment.bodyLinearAcceleration.zeroPhaseFilter = ...
                                      [[0; 0; 9.8100]*ones(1,armIMUExperiment.samples);
                                       armIMUExperiment.bodyLinearAcceleration.zeroPhaseFilter];

    armIMUExperiment.bodyLinearAcceleration.movAvgFilter = ...
                                      [[0; 0; 9.8100]*ones(1,armIMUExperiment.samples);
                                       armIMUExperiment.bodyLinearAcceleration.movAvgFilter];    

    armIMUExperiment.bodyAngularAcceleration.numerical = ...
                                      [zeros(3, armIMUExperiment.samples);
                                       armIMUExperiment.bodyAngularAcceleration.numerical];
                                  
    fprintf('>> IMU acceleration/velocities signals filtered!\n')
end