function [lambda_hat_online, lambda_rams_mav, J_log] = ...
    fcn_robot_kinematics_rams_online_learning(signals, ...
                                              lambda_hat_0, ...
                                              IS_ACC, ...
                                              A_kin_pi, ...
                                              epochs, ...
                                              buffer, ...
                                              batchSize, ...
                                              constPar)

%     arguments
%         signals (:,:) double
%         A_kin_pi (:,:) double
%         constPar struct
%     end
% *************************************************************************
%                         RAMS Gradient Descent                           *
% *************************************************************************

assert(isstruct(constPar),'consrPar needs to be struct')

assert(isa(signals,'double') && all(size(signals,1) == (3*constPar.noj+ 4*3*constPar.nob)) && isreal(signals), ...
  ['Signals vector needs too be [',num2str(3*constPar.noj+ 4*3*constPar.nob) ,' x 1] double']); 

assert(isa(IS_ACC,'logical') && all(size(IS_ACC) == [1 1]) && isreal(IS_ACC), ...
  'Parameter IS_ACC has to be logical')

assert(isa(A_kin_pi,'double') && all(size(A_kin_pi) == [constPar.nob,constPar.nob]) && isreal(A_kin_pi), ...
  ['A_kin_pi matrix needs too be [',num2str(constPar.nob) ,' x ,',num2str(constPar.nob),'] double']); 


clc
close all

[parents, children] = find(triu(A_kin_pi) == 1);

ang_vel_indices = constPar.noj*3 + 1:constPar.noj*3 + 3*constPar.nob;
lin_vel_indices = ang_vel_indices(end) + 1:ang_vel_indices(end) + 3*constPar.nob;
ang_acc_indices = lin_vel_indices(end) + 1:lin_vel_indices(end) + 3*constPar.nob;
lin_acc_indices = ang_acc_indices(end) + 1:ang_acc_indices(end) + 3*constPar.nob;

% Moving average filter (maf)
windowSize  = 100; 
b_maf       = (1/windowSize)*ones(1,windowSize);
a_maf       = 1;
signals_maf = filter(b_maf,a_maf,signals,[],2);        

% Approximate angular acceleratrion with central difference
signals_maf(ang_acc_indices,:) = gradient(signals_maf(ang_vel_indices,:),1E-3);

% Define the data pool
dataPool    = transpose(signals_maf);

% Manifold settings
S2_man      = spherefactory(2);

% Shuffle the data points
% rng('default')
k  = randperm(size(dataPool,1));

% Cost function and gradient settings
GRAD_TYPE    = 'Riemannian';
% GRAD_TYPE    = 'Euclidean';

% Replay buffer size and initial values
dataBuff    = zeros(buffer.Size, size(dataPool,2));
rho_G_avg_n = zeros(buffer.Size, 1); % gradient norm per sample point

%                                  ||
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
CONTINUE  = 0;
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%                                  ||

if ~CONTINUE
    current_sample = 1;
    close all
    % Riemannian AMS Initial settings =====================================
    m         = zeros(size(lambda_hat_0,1),constPar.noj);
    v         = zeros(constPar.noj, 1);
    v_hat     = zeros(constPar.noj, 1);    
    tau       = zeros(size(lambda_hat_0,1),constPar.noj);
    eta_max   = 0.01*ones(constPar.noj,1);

    % Euclidean AMS Initial settings ======================================
    m_ams     = zeros(7,1);
    v_ams     = zeros(7,1);
    v_hat_ams = zeros(7,1);
    alpha     = 0.001;
    % =====================================================================

    % Parameter update initial value
    lambda_update = repmat(lambda_hat_0,1,constPar.noj);
    % Settings for RAMS gradient descent
    %epochs      = size(dataPool,1);
    lambda_rams = zeros(size(lambda_hat_0,1), constPar.noj, epochs);
    J_log       = zeros(epochs, 1);
    rho_G_avg   = zeros(epochs, 1);
    SHUFFLE     = 0;
    SORT_BUFF   = 1;
    SIMPLIFIED  = false;

    clc
    cprintf('*yellow', '>> Running %s GD with:\n', GRAD_TYPE)
    cprintf('*yellow', '>> - LEARN RATE  : %d\n', alpha)
    cprintf('*red',    '>> - SHUFFLE     : %d\n', SHUFFLE)
    cprintf('*yellow', '>> - EPOCHS      : %d\n', epochs)
    cprintf('*yellow', '>> - MB SIZE     : %d\n',batchSize)
    cprintf('*yellow', '>> - BUFFER SIZE : %d\n',buffer.Size)
    disp('>> Initial Point:')
    disp(lambda_update);
    tic
    pause(3)
else
    current_sample  = iter;
    cprintf('*yellow', '>> Resuming %s GD with:\n', GRAD_TYPE)
    disp('>> Current Point:')
    lambda_update = squeeze(lambda_rams(:,:,current_sample-1));
%     lambda_update = reshape(permute(lambda_rams(:,current_sample-1,:),[1,3,2]),70,1);
%     disp(reshape(lambda_update,10,7));
    tic     
end
pause(1)
t = 1;
update_grad  = 1;
for iter = current_sample:epochs       
    % Store sample in replay buffer
    buffer.storeSample(transpose(dataPool(iter,:)));

    if iter >= buffer.Size && update_grad == 1  
        t = t + 1;
        % Calculate cost function and Euclidean/Riemannian gradients ------
        % [update_grad, J, riem_gradJ_mb, wrench_hat, P_i, ...
        %     rho_G_avg_n, dJ_dtheta_mb, dJ_dP_i_agg] = ...
        %     fcn_franka_axis_learning_manifold_MANOPT(lambda_est, ...
        %     Buff, learn_kin, S1, S2)
        J_tmp         = zeros(batchSize,constPar.noj);
        dJ_lambda_tmp = zeros(numel(lambda_update(:,1)), batchSize, constPar.noj); 
        rgrad_J       = zeros(numel(lambda_update(:,1)),7);
        for jnt = 1:constPar.noj
            pr = parents(jnt);
            ch = children(jnt);
            omg_indices   = reshape(ang_vel_indices, 3, constPar.nob);
            vel_indices   = reshape(lin_vel_indices, 3, constPar.nob);
            domg_indices  = reshape(ang_acc_indices, 3, constPar.nob);
            dvel_indices  = reshape(lin_acc_indices, 3, constPar.nob);            
            %entries       = [jnt; constPar.noj +  jnt; omg_indices(:,jnt);omg_indices(:,jnt+1)];
            entries       = [(ch-1); 
                             constPar.noj + (ch-1); ...
                             omg_indices(:,pr); ...
                             omg_indices(:,ch);...
                             vel_indices(:,pr);...
                             vel_indices(:,ch); ...
                             domg_indices(:,pr); ...
                             domg_indices(:,ch);...
                             dvel_indices(:,pr);...
                             dvel_indices(:,ch);];
            % Draw samples from buffer
            miniBatch = transpose(buffer.drawSamples(batchSize));
            miniBatch = miniBatch(:,entries);
            for n = 1:size(miniBatch,1)
                q_n         = miniBatch(n,1)';
                dq_n        = miniBatch(n,2)';
                p_OMG_p_n   = miniBatch(n,3:5)';
                c_OMG_c_n   = miniBatch(n,6:8)';
                p_VEL_p_n   = miniBatch(n,9:11)';
                c_VEL_c_n   = miniBatch(n,12:14)';                
                p_dOMG_p_n  = miniBatch(n,15:17)';
                c_dOMG_c_n  = miniBatch(n,18:20)';
                p_dVEL_p_n  = miniBatch(n,21:23)';
                c_dVEL_c_n  = miniBatch(n,24:26)';                                
                
                c_OMG_c_hat = expm(-skew(lambda_update(5:7,jnt))*q_n)*expm(-skew(lambda_update(1:3,jnt))*lambda_update(4,jnt))*p_OMG_p_n ...
                                               + dq_n*lambda_update(5:7,jnt);

%                 beta   = 0.1;
%                 beta = 0.5;
                beta   = 0.01;
                beta = constPar.beta;
                if IS_ACC == 0
                    % Cost
                    J_tmp(n, jnt)          = fcn_J_morph_vel_mex(p_OMG_p_n, p_VEL_p_n, c_OMG_c_n, c_VEL_c_n, q_n, dq_n, lambda_update(:,jnt), beta);                
                    % Sample gradient
                    dJ_lambda_tmp(:,n,jnt) = gradest(@(x)  fcn_J_morph_vel_mex(p_OMG_p_n, p_VEL_p_n, c_OMG_c_n, c_VEL_c_n, q_n, dq_n, x, beta),lambda_update(:,jnt));
                    %dJ_lambda_tmp(:,n,jnt) = dJ_morph_vel_dlambda(p_OMG_p_n, p_VEL_p_n, c_OMG_c_n, c_VEL_c_n, q_c, dq_c, lambda_hat)
                else
                    % Cost
                    J_tmp(n, jnt)          = fcn_J_morph_acc_mex(p_OMG_p_n, p_dOMG_p_n, p_dVEL_p_n, c_OMG_c_n, c_dOMG_c_n, c_dVEL_c_n, q_n, dq_n, lambda_update(:,jnt), beta);
                    % Sample gradient
                    %dJ_lambda_tmp(:,n,jnt) = gradest(@(x) fcn_J_morph_acc_mex(p_OMG_p_n, p_dOMG_p_n, p_dVEL_p_n, c_OMG_c_n, c_dOMG_c_n, c_dVEL_c_n, q_n, dq_n, x, beta),lambda_update(:,jnt));                   
                    dJ_lambda_tmp(:,n,jnt) = fcn_dJ_morph_acc_dlambda_mex(p_OMG_p_n, p_dOMG_p_n, p_dVEL_p_n, c_OMG_c_n, c_dOMG_c_n, c_dVEL_c_n, q_n, dq_n, lambda_update(:,jnt), beta);
                end                                      
            end
            % Total (Euclidean) gradient per batch
            dJ_lambda_mb   = mean(dJ_lambda_tmp(:,:,jnt),2);
            % Batch Riemannian gradient                    
            rgrad_J(:,jnt) = [S2_man.egrad2rgrad(lambda_update(1:3,jnt), dJ_lambda_mb(1:3)); 
                              dJ_lambda_mb(4);
                              S2_man.egrad2rgrad(lambda_update(5:7,jnt), dJ_lambda_mb(5:7));
                              dJ_lambda_mb(8:13)];
        end
        
        
        % Calculate parameter update --------------------------------------
        switch GRAD_TYPE
            case 'Euclidean'
                % Based on Euclidean AMS gradient descent
                warning('Section undefined')
                break
             case 'Riemannian'
                % Based on RAMS gradient descent
                learnRate = eta_max(1);
                [m, v, v_hat, tau, lambda_update, rho_G_avg(iter), eta_avg] = ...
                        fcn_riemannian_amsgrad_for_kinematics(...
                                m, v, v_hat, tau, eta_max, rgrad_J, lambda_update, ...
                                update_grad, S2_man, 'RAMS', SIMPLIFIED, constPar);
        end
        % Total cost per batch
        J_log(iter)  = mean(J_tmp,'All'); 


    elseif iter < buffer.Size  
        % Total cost per batch
        J_log(iter)     = 0;
        rho_G_avg(iter) = 0;            
    end


    % Averaged gradient's Riemannian norm
    rho_G_avg(iter) = (1/constPar.noj)*rho_G_avg(iter);
    
    %if iter > buffer.Size && rho_G_avg(iter) < 1E-6
%     if iter > buffer.Size && J_log(iter) < 1E-10        
%         cprintf('*yellow', '>> Training stopped!\n')
%         update_grad = 0;
%         J_log(iter)     = J_log(iter-1);
%         rho_G_avg(iter) = rho_G_avg(iter-1);
%         break
%     end
    
    % Log intermiediate values
    %disp(iter)
%     disp(size(gamma_update,2))
    %lambda_rams(:,:,iter) = lambda_update;

    % Re sample initial point if optimization stuck
    lambda_rams(:,:,iter) = lambda_update;% + ~mod(iter,500)*rand(size(lambda_update));%*J_log(iter);

% KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
    if (iter > buffer.Size + batchSize) && (mod(iter,4*batchSize) == 0)
        if mean(J_log(iter-batchSize:iter))>10^1
            clc
            warning('Resetting ICs and gradient!')
            pause(0.5)
            lambda_update = cell2mat(arrayfun(@(i) fcn_resampleInitialPoint(),1:constPar.noj,'UniformOutput',false));
            m         = zeros(size(lambda_hat_0,1),constPar.noj);
            v         = zeros(constPar.noj, 1);
            v_hat     = zeros(constPar.noj, 1);    
            tau       = zeros(size(lambda_hat_0,1),constPar.noj);
            eta_max   = 0.01*ones(constPar.noj,1);
        elseif mean(J_log(iter-5*batchSize:iter))<1E-3
            break
        end
    end
% KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK


    % Display progress
    if iter <= buffer.Size
        warning('No yet enough samples in the replay buffer!')
    elseif mod(iter,1) == 0
        disp(['Iter: ', num2str(iter),'/',num2str(epochs),' | Cost: ', num2str(J_log(iter)), ' | Grad norm: ' num2str(rho_G_avg(iter)), ' | eta_avg: ', num2str(eta_avg)]);
    end    


end

cprintf('*yellow', '>> Computation finished\n')

% SAVE_RESULTS = false;
% 
% if SAVE_RESULTS == true
%     disp('Saving...');
%     save(['/home/diaz/Dropbox/PhD_WORK/matlab_irt/franka/learning/morphology_learning/tests/',...
%         'poppy_kin_RAMS:', GRAD_TYPE, ...
%         '_LR:', strrep(num2str(learnRate),'.','p'), ...
%         '_B:' , num2str(BUFF_SIZE), ...
%         '_X:' , num2str(batchSize), ...
%         '_D:' , strcat(date,'_',datestr(now,'HH:MM'))], 'lambda_rams','J_log','rho_G_avg')
%     pause(3) 
%     disp('Saving completed!');
% end

% if SAVE_RESULTS == true
%     disp('Saving...');
%     save(['/media/fernando/extreme_ssd/panda_experiments/ip_learning_results/tests/',...
%         'panda_GD:', GRAD_TYPE, ...
%         '_LR:', strrep(num2str(learnRate),'.','p'), ...
%         '_IC:', initialPoint, ...
%         '_M:' , num2str(MEASUREMENTS), ...
%         '_B:' , num2str(BUFF_SIZE), ...
%         '_X:' , num2str(batchSize), ...
%         '_BM:', num2str(ACTION_TORSO),...
%         '_N:' , num2str(NOISE),...
%         '_LK' , num2str(LOAD_KNOWN),...
%         '_D:' , strcat(date,'_',datestr(now,'HH:MM'))], 'theta_rams','J_log','rho_G_avg')
%     pause(3) 
%     disp('Saving completed!');
% end
%
% figure('Color','w')
% for p = 1:13
%     subplot(3,5,p)
%     plot(BUFF_SIZE:iter,permute(lambda_rams(8:13,p,BUFF_SIZE:iter),[1,3,2]))
%     ylabel(['$\boldmath{\lambda}_' num2str(p), '(t)$'],'interpreter','latex')
%     xlabel(['Time [s]'])
% end

%  Filter estimates (if needed) and average over the last n samples =======

B = 1/100*ones(1,100);
A = 1; %denominator coefficients
lambda_rams_mav = filter(B,A,permute(lambda_rams(:,:,buffer.Size:iter),[3,1,2])); %filter input x and get result in y
lambda_rams_mav = permute(lambda_rams_mav,[2,3,1]);

figure('Color','w')
for p = 1:13
    subplot(3,5,p)
    plot(1:size(lambda_rams_mav,3),squeeze(lambda_rams_mav(p,:,:)))
    ylabel(['$\boldmath{\lambda}_{' num2str(p), '}(t)$'],'interpreter','latex')
    xlabel(['Time [s]'])
end

lambda_hat_online  = mean(lambda_rams_mav(:,:,end-100:end),3);

end

function lambda_hat_0 = fcn_resampleInitialPoint()
%     xi_0         = rand(3,1);
%     xi_0         = xi_0/norm(xi_0);
%     phi_0        = wrapToPi(deg2rad(randi(360)));
%     zeta_0       = rand(3,1);
%     zeta_0       = zeta_0/norm(zeta_0);
%     gamma_hat_0  = [xi_0;phi_0;zeta_0];
%     rho_hat_0    = 1E-3*ones(6,1);
%     lambda_hat_0 = [gamma_hat_0;
%                     rho_hat_0];

    xi_0         = ones(3,1);
    xi_0         = xi_0/norm(xi_0);
    phi_0        = wrapToPi(deg2rad(randi(360)));
    zeta_0       = ones(3,1);
    zeta_0       = zeta_0/norm(zeta_0);
    gamma_hat_0  = [xi_0;phi_0;zeta_0];
    rho_hat_0    = 1E-3*ones(6,1);
    lambda_hat_0 = [gamma_hat_0;
                    rho_hat_0];
    
end