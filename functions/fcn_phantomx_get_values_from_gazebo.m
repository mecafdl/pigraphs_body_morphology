%% **************************************************************************
%                        Get values from Gazebo                           *
% *************************************************************************
function dataGazebo = fcn_phantomx_get_values_from_gazebo(dataPath, ROTATE_FRAMES, constPar)
    clc
    
    % Load saved data matrices collected from the simulation in Gazebo
    % =========================================================================
    if numel(dataPath)<1
    states           = transpose(readNPY(fullfile(dataPath,'states.npy')));
    actions          = transpose(readNPY(fullfile(dataPath,'actions.npy')));
    else
        states       = [];
        actions      = [];
        for path_entry = 1:numel(dataPath)
            states  = [states, transpose(readNPY(fullfile(dataPath{path_entry},'states.npy')))];
            actions = [actions,transpose(readNPY(fullfile(dataPath{path_entry},'actions.npy')))];
        end
    end
    
    imu_signals = states(1:constPar.nol*6,:);
    omg_signals = zeros(3*constPar.nol, size(states,2));
    acc_signals = zeros(3*constPar.nol, size(states,2));
    
    
    rng('default')
    R_sensors       = randrot(3, constPar.nol); % random rotation matrices
    n_points        = size(states,2);
%     ROTATE_FRAMES   = 0;
    for i=1:constPar.nol
        if(ROTATE_FRAMES == 1 && i~=1)
            warning('The Cartesian angular veclocities will be rotated at random.')
            omg_signals(3*(i-1)+1:3*i,:) = R_sensors(:,:,i)*imu_signals(6*(i-1)+1:6*(i-1)+3,:);
            acc_signals(3*(i-1)+1:3*i,:) = R_sensors(:,:,i)*imu_signals(6*(i-1)+4:6*(i-1)+6,:);
        else    
            omg_signals(3*(i-1)+1:3*i,:) = imu_signals(6*(i-1)+1:6*(i-1)+3,:);
            acc_signals(3*(i-1)+1:3*i,:) = imu_signals(6*(i-1)+4:6*(i-1)+6,:);    
        end
    end
    
    % Define the joint signals
    joint_signals = [states(constPar.nol*6+1:end,:); 
                     actions]; 

    % IMU sensor signals from gazebo are re arranged to match the order of the
    % joint signals
    omg_signals_aux = NaN(size(omg_signals));
    acc_signals_aux = NaN(size(acc_signals));
    j = [0,1,2,3,7,8,9,4,5,6,10,11,12,16,17,18,13,14,15]+1;
    for i =1:19
        omg_signals_aux(3*(i-1)+1:3*i,:) = omg_signals(3*(j(i)-1)+1:3*j(i),:);
        acc_signals_aux(3*(i-1)+1:3*i,:) = acc_signals(3*(j(i)-1)+1:3*j(i),:);
    end
    
    omg_signals = omg_signals_aux;
    acc_signals = acc_signals_aux;
    
    dataGazebo  = [joint_signals; ...
                   omg_signals; ...
                   acc_signals];

end