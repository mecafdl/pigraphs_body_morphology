%% **********************************************************************
%            OFFLINE LEARNING (INTERIOR POINT OPTIMIZATION)               *
% *************************************************************************

function [gamma_hat, rho_hat,constPar] = ...
    fcn_robot_kinematics_offline_optimization(A_kin_pi, ...
                                              armSignals, ...
                                              N_samples, ...
                                              MOVING_BASE, ...
                                              IS_ACCELERATION, ...
                                              constPar)

    clc
    close all
    
    if ~MOVING_BASE
        signals = armSignals.fixedBase;
    else
        signals = armSignals.movingBase;
    end
    
    [parents, children] = find(triu(A_kin_pi) == 1);
    
    %% ROTATION MATRICES/ROTATION AXES  BETWEEN SENSORS =====================
    
    clc
    rng('default')
    % N_samples   = 100;
    test_points = sort(randperm(size(signals.jointPosition,2),N_samples));% choose k unique data points from data set
    
    % Initial point -----------------------------------------------------------
    xi_0   = rand(3,1);
    xi_0   = xi_0/norm(xi_0);
    phi_0  = deg2rad(randi(360));
    zeta_0 = rand(3,1);
    zeta_0 = zeta_0/norm(zeta_0);
    x0     = [xi_0;phi_0;zeta_0];
    
    % Optimizer settings ------------------------------------------------------
    constPar.MaxIter     = 100;
    constPar.MaxFunEvals = 1000;
    options = optimset('Display','iter-detailed','MaxIter',constPar.MaxIter,'MaxFunEvals', 1000,...
        'TolFun',1e-6,'UseParallel',false);
    
    gamma_hat = NaN(7,constPar.noj);
    for j = 1:constPar.noj
        %p = j;
        %c = j+1;
        p = parents(j);
        c = children(j);
        cprintf('*yellow', ['>> Checking: ' num2str(p-1) '-->' num2str(c-1) '\n'])
    
        [gamma_hat(:,j), ~] = fmincon(...
        @(x) fcn_get_rotation_err(signals.bodyAngularVelocity(3*(p-1)+1:3*p,test_points), ...
                                  signals.bodyAngularVelocity(3*(c-1)+1:3*c,test_points), ...
                                  signals.jointPosition(c-1,test_points), ...
                                  signals.jointVelocity(c-1,test_points), x),...
                                  x0,[],[],[],[],[],[],...
                                      @(x)fcn_axis_nonlConstraints(x),options);
        
        disp('--------------------------------------------------------------------')
    end
    % Sensor-to-Sensor rotation matrix
    constPar.p_R_c = @(q_c, x) expm(skew([x(1);x(2);x(3)])*(x(4)))*expm(skew([x(5);x(6);x(7)]).*q_c);
    cprintf('*green', '>> Rot mat between sensors and joint axes found!!!\n')
    
    %% FIND JOINT POSITION VECTORS ==========================================
    
    clc
    % test_points = sort(randperm(size(q_ref,2),3000));% choose k unique data points from data set
    
    % IS_ACCELERATION = 1;
    
    if IS_ACCELERATION == 1
        reference_signal = 'ACC';
    else
        reference_signal = 'VEL';
    end
    
    rho_hat = NaN(6,constPar.noj);
    for j = 1:constPar.noj
        %p = j;
        %c = j+1;
        p = parents(j);
        c = children(j);
    
        x0 = 1E-3*ones(6,1);
        
        if IS_ACCELERATION == 0
            cprintf('*yellow', ['>> Checking: ' num2str(p-1) '-->' num2str(c-1) ' using ', reference_signal,'\n'])
            [rho_hat(:,j), ~] = fmincon(...
            @(x) fcn_get_cartesian_velocity_error_extended(...
                     signals.bodyLinearVelocity(3*(p-1)+1:3*p,test_points), ...
                     signals.bodyLinearVelocity(3*(c-1)+1:3*c,test_points), ...
                     signals.bodyAngularVelocity(3*(p-1)+1:3*p,test_points), ...
                     signals.bodyAngularVelocity(3*(c-1)+1:3*c,test_points), ...
                     signals.jointPosition(j,test_points), ...
                     constPar.p_R_c, gamma_hat(:,j), ...
                     x), x0,[],[],[],[],[],[], [],options);    
    
        else
            cprintf('*yellow', ['>> Checking: ' num2str(p-1) '-->' num2str(c-1) ' using ', reference_signal,'\n'])
            [rho_hat(:,j), ~] = fmincon(...
            @(x) fcn_get_cartesian_acceleration_error_extended(...
                     signals.bodyLinearAcceleration(3*(p-1)+1:3*p,test_points), ...
                     signals.bodyLinearAcceleration(3*(c-1)+1:3*c,test_points), ...
                     signals.bodyAngularVelocity(3*(p-1)+1:3*p,test_points), ...
                     signals.bodyAngularVelocity(3*(c-1)+1:3*c,test_points), ...
                     signals.bodyAngularAcceleration(3*(p-1)+1:3*p,test_points), ...
                     signals.bodyAngularAcceleration(3*(c-1)+1:3*c,test_points), ...
                     signals.jointPosition(j,test_points), ...
                     constPar.p_R_c, gamma_hat(:,j), ...
                     0.1, x), x0,[],[],[],[],[],[], [],options);
        end
    disp('--------------------------------------------------------------------')
    end
    cprintf('*Green', '>> Joint center point relative to sensor CS found!!!\n')
end