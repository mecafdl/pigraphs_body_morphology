function armSignals = fcn_loadRobotArmData()

% Load configuration parameters
% load('./data/franka_parameters.mat')
    [filepath,~,~] = fileparts(mfilename('fullpath'));
    q_indices       = 1:7;
    dq_indices      = 8:14;
    ddq_indices     = 15:21;
    torque_indices  = 22:28;
    ang_vel_indices = 29:52;
    lin_vel_indices = 53:76;
    ang_acc_indices = 77:100;
    lin_acc_indices = 101:124;   
    
    armData = load([filepath,'/../data/frankaProprioceptionFixedBaseSimulated.mat'], 'frankaProprioceptionFixedBaseSimulated');
    armSignals.fixedBase.jointPosition           = armData.frankaProprioceptionFixedBaseSimulated(q_indices, :);
    armSignals.fixedBase.jointVelocity           = armData.frankaProprioceptionFixedBaseSimulated(dq_indices, :);
    armSignals.fixedBase.jointTorque             = armData.frankaProprioceptionFixedBaseSimulated(torque_indices, :);
    armSignals.fixedBase.bodyAngularVelocity     = armData.frankaProprioceptionFixedBaseSimulated(ang_vel_indices, :);
    armSignals.fixedBase.bodyLinearVelocity      = armData.frankaProprioceptionFixedBaseSimulated(lin_vel_indices, :);
    armSignals.fixedBase.bodyAngularAcceleration = armData.frankaProprioceptionFixedBaseSimulated(ang_acc_indices, :);
    armSignals.fixedBase.bodyLinearAcceleration  = armData.frankaProprioceptionFixedBaseSimulated(lin_acc_indices, :);
    
    armData = load([filepath,'/../data/frankaProprioceptionMovingBaseSimulated.mat'], 'frankaProprioceptionMovingBaseSimulated');
    armSignals.movingBase.jointPosition           = armData.frankaProprioceptionMovingBaseSimulated(q_indices, :);
    armSignals.movingBase.jointVelocity           = armData.frankaProprioceptionMovingBaseSimulated(dq_indices, :);
    armSignals.movingBase.jointTorque             = armData.frankaProprioceptionMovingBaseSimulated(torque_indices, :);
    armSignals.movingBase.bodyAngularVelocity     = armData.frankaProprioceptionMovingBaseSimulated(ang_vel_indices, :);
    armSignals.movingBase.bodyLinearVelocity      = armData.frankaProprioceptionMovingBaseSimulated(lin_vel_indices, :);
    armSignals.movingBase.bodyAngularAcceleration = armData.frankaProprioceptionMovingBaseSimulated(ang_acc_indices, :);
    armSignals.movingBase.bodyLinearAcceleration  = armData.frankaProprioceptionMovingBaseSimulated(lin_acc_indices, :);
    
    % Measurements sampling time
    armSignals.samplingTime  = 1E-3;

end