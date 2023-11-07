%% **********************************************************************
%                         RAMS Gradient Descent                           *
% *************************************************************************
function [lambda_hat_online] = fcn_poppy_morphology_online_learning_rams(proprioceptiveSignals, parents, children, constPar)
    close all
    
    ang_vel_indices = constPar.noj*3 + 1:constPar.noj*3 + 3*constPar.nob;
    lin_vel_indices = ang_vel_indices(end) + 1:ang_vel_indices(end) + 3*constPar.nob;
    ang_acc_indices = lin_vel_indices(end) + 1:lin_vel_indices(end) + 3*constPar.nob;
    lin_acc_indices = ang_acc_indices(end) + 1:ang_acc_indices(end) + 3*constPar.nob;
    
    % Moving average filter (maf)
    windowSize             = 100; 
    b_maf                  = (1/windowSize)*ones(1,windowSize);
    a_maf                  = 1;
    
    
    signals_maf                    = filter(b_maf,a_maf,proprioceptiveSignals,[],2);        
    %signals_maf(ang_acc_indices,:) = gradient(signals_maf(ang_vel_indices,:),1E-3);
    dataPool                       = transpose(proprioceptiveSignals);
    
    % Manifold settings
    S2_man  = spherefactory(2);
    
    % Shuffle the data points
    % rng('default')
    k  = randperm(size(dataPool,1));
    
    % Initial point choice
    % xi_0         = rand(3,1);
    % xi_0         = xi_0/norm(xi_0);
    % phi_0        = wrapToPi(deg2rad(randi(360)));
    % zeta_0       = rand(3,1);
    % zeta_0       = zeta_0/norm(zeta_0);
    
    xi_0         = rand(3,1);
    xi_0         = xi_0/norm(xi_0);
    phi_0        = wrapToPi(deg2rad(randi(360)));
    zeta_0       = rand(3,1);
    zeta_0       = zeta_0/norm(zeta_0);
    gamma_hat_0  = [xi_0;phi_0;zeta_0];
    rho_hat_0    = 1E-3*ones(6,1);
    lambda_hat_0 = [gamma_hat_0;
                    rho_hat_0];
    
    % Cost function and gradient settings
    GRAD_TYPE    = 'Riemannian';
    % GRAD_TYPE    = 'Euclidean';
    
    % Size of the data batch (buffer)
    MB_SIZE     = 10;
    % MB_SIZE   = 50;
    
    % Replay buffer size and initial values
    BUFF_SIZE   = 1000;
    dataBuff    = zeros(BUFF_SIZE, size(dataPool,2));
    % dataBuff    = zeros(BUFF_SIZE, 8);
    rho_G_avg_n = zeros(BUFF_SIZE, 1); % gradient norm per sample point
    
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
        %v         = zeros(size(gamma_hat_0,1), 1);
        %v_hat     = zeros(size(gamma_hat_0,1), 1);
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
        epochs      = size(dataPool,1);
        lambda_rams = zeros(size(lambda_hat_0,1), constPar.noj, epochs);
        J_log       = zeros(epochs, 1);
        rho_G_avg   = zeros(epochs, 1);
        SHUFFLE     = 0;
        SORT_BUFF   = 1;
        SIMPLIFIED  = false;
    
        clc
        cprintf('*yellow', '>> Running %s GD with:\n', GRAD_TYPE)
        cprintf('*yellow', '>> - LEARN RATE  : %d\n', alpha)
    %     cprintf('*yellow', '>> - MEASUREMENTS: %d\n', MEASUREMENTS)
        %cprintf('*yellow', '>> - ACTION_TORSO: %d\n', ACTION_TORSO)
        cprintf('*red',    '>> - SHUFFLE     : %d\n', SHUFFLE)
        cprintf('*yellow', '>> - EPOCHS      : %d\n', epochs)
        cprintf('*yellow', '>> - MB SIZE     : %d\n',MB_SIZE)
        cprintf('*yellow', '>> - BUFFER SIZE : %d\n',BUFF_SIZE)
        disp('>> Initial Point:')
        disp(lambda_update);
        tic
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
        % Shuffle samples
        if SHUFFLE == 1
            dataPoint = dataPool(k(iter),:);       
        else
            dataPoint = dataPool(iter,:);       
        end
           
        % Crate a buffer of samples sorted by Riemannian norm
        if iter <= BUFF_SIZE
            dataBuff(iter,:)               = dataPoint;
        else
            dataBuff(randi(BUFF_SIZE,1),:) = dataPoint;    
        end
        if iter == BUFF_SIZE
            cprintf('*yellow', '>> Buffer full!', GRAD_TYPE)
        end
        %if mod(iter, MB_SIZE) == 0 && iter >= BUFF_SIZE 
        if iter >= BUFF_SIZE && update_grad == 1  
            t = t + 1;
            
            if iter < BUFF_SIZE
                p = randperm(iter, MB_SIZE);
            else
                p = randperm(BUFF_SIZE, MB_SIZE);
            end
            % Calculate cost function and Euclidean/Riemannian gradients ------
            % [update_grad, J, riem_gradJ_mb, wrench_hat, P_i, ...
            %     rho_G_avg_n, dJ_dtheta_mb, dJ_dP_i_agg] = ...
            %     fcn_franka_axis_learning_manifold_MANOPT(lambda_est, ...
            %     Buff, learn_kin, S1, S2)
            J_tmp         = zeros(MB_SIZE,constPar.noj);
            dJ_lambda_tmp = zeros(numel(lambda_update(:,1)), MB_SIZE, constPar.noj); 
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
                MB            = dataBuff(p,entries);
                for n = 1:size(MB,1)
                    q_n         = MB(n,1)';
                    dq_n        = MB(n,2)';
                    p_OMG_p_n   = MB(n,3:5)';
                    c_OMG_c_n   = MB(n,6:8)';
                    p_VEL_p_n   = MB(n,9:11)';
                    c_VEL_c_n   = MB(n,12:14)';                
                    p_dOMG_p_n  = MB(n,15:17)';
                    c_dOMG_c_n  = MB(n,18:20)';
                    p_dVEL_p_n  = MB(n,21:23)';
                    c_dVEL_c_n  = MB(n,24:26)';                                
                    
                    c_OMG_c_hat = expm(-skew(lambda_update(5:7,jnt))*q_n)*expm(-skew(lambda_update(1:3,jnt))*lambda_update(4,jnt))*p_OMG_p_n ...
                                                   + dq_n*lambda_update(5:7,jnt);
                    % Sample cost
                    %J_tmp(n, jnt)          = (1/2)*transpose(c_OMG_c_n - c_OMG_c_hat)*(c_OMG_c_n - c_OMG_c_hat);                    
    %                 J_tmp(n, jnt)          = (1/2)*transpose(c_OMG_c_n - c_OMG_c_hat)*(c_OMG_c_n - c_OMG_c_hat) + ...
    %                                          fcn_Jr(p_OMG_p_n, p_VEL_p_n, c_OMG_c_n, c_VEL_c_n, lambda_update(8:13,jnt));
    
                    J_tmp(n, jnt) = fcn_J_morph_mex(q_n, dq_n, p_OMG_p_n, p_VEL_p_n, c_OMG_c_n, c_VEL_c_n, lambda_update(:,jnt));
                    % Sample gradient
    %                 dJ_lambda_tmp(:,n,jnt) = [gradest(@(x) fcn_Jpsi_no_norm_mex(q_n, dq_n,p_OMG_p_n, c_OMG_c_n,x), lambda_update(1:7,jnt)), ...
    %                                           gradest(@(x)  fcn_Jr(p_OMG_p_n, p_VEL_p_n, c_OMG_c_n, c_VEL_c_n, x), lambda_update(8:13,jnt))];
                                          
    
                    dJ_lambda_tmp(:,n,jnt) = gradest(@(x)  fcn_J_morph_mex(q_n, dq_n, p_OMG_p_n, p_VEL_p_n, c_OMG_c_n, c_VEL_c_n, x),lambda_update(:,jnt));
    
    %                 dJ_lambda_tmp(:,n,jnt) = [gradest(@(x) fcn_Jpsi_no_norm_mex(q_n, dq_n,p_OMG_p_n, c_OMG_c_n,x), lambda_update(1:7,jnt)), ...
    %                                           gradest(@(x) fcn_Jr_acc(p_dVEL_p_n, c_dVEL_c_n, p_dOMG_p_n, c_dOMG_c_n, p_OMG_p_n, c_OMG_c_n, x, 0.001), lambda_update(8:13,jnt))];
    %                                       
                                          
                                          
                                          
    
    %                 dJ_lambda_tmp(:,n,jnt) = [gradest(@(x) fcn_Jpsi_no_norm_mex(q_n, dq_n,p_OMG_p_n, c_OMG_c_n,x), lambda_update(1:7,jnt)), ...
    %                                           gradest(@(x)  fcn_Jr_mex(p_OMG_p_n, p_VEL_p_n, c_OMG_c_n, c_VEL_c_n, x), lambda_update(8:13,jnt))];                
                                          
    %                 dJ_lambda_tmp(:,n,jnt) = [gradest(@(x) fcn_Jpsi_no_norm_2020a_mex(q_n, dq_n,p_OMG_p_n, c_OMG_c_n,x), lambda_update(1:7,jnt)), ...
    %                                           gradest(@(x)  fcn_Jr(p_OMG_p_n, p_VEL_p_n, c_OMG_c_n, c_VEL_c_n, x), lambda_update(8:13,jnt))];                                      
                                          
                                          
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
        elseif iter < BUFF_SIZE  
            % Total cost per batch
            J_log(iter)     = 0;
            rho_G_avg(iter) = 0;            
        end
        % Averaged gradient's Riemannian norm
        rho_G_avg(iter) = (1/constPar.noj)*rho_G_avg(iter);
        
        %if iter > BUFF_SIZE && rho_G_avg(iter) < 1E-6
    %     if iter > BUFF_SIZE && J_log(iter) < 1E-10        
    %         cprintf('*yellow', '>> Training stopped!\n')
    %         update_grad = 0;
    %         J_log(iter)     = J_log(iter-1);
    %         rho_G_avg(iter) = rho_G_avg(iter-1);
    %         break
    %     end
        
        % Log intermiediate values
        %disp(iter)
    %     disp(size(gamma_update,2))
        lambda_rams(:,:,iter) = lambda_update;
    
        % Display progress
        % disp(['Iter: ', num2str(iter) ' Cost: ', num2str(J_log(iter)), ' Grad norm: ' num2str(rho_G_avg(iter))]);
        if iter <= BUFF_SIZE
            %disp(['Iter: ', num2str(iter) ' Cost: ', num2str(J_log(iter)), ' Grad norm: ' num2str(rho_G_avg(iter))]);
            warning('No yet enough samples in the replay buffer!')
        elseif mod(iter,1) == 0
    %     elseif J_log(iter) < min(J_log(BUFF_SIZE+1:iter-1))
            disp(['Iter: ', num2str(iter) ' | Cost: ', num2str(J_log(iter)), ' | Grad norm: ' num2str(rho_G_avg(iter)), ' | eta_avg: ', num2str(eta_avg)]);
        end    
    end
    
    cprintf('*yellow', '>> Computation finished\n')
    
    SAVE_RESULTS = false;
    
    if SAVE_RESULTS == true
        disp('Saving...');
        save(['/home/diaz/Dropbox/PhD_WORK/matlab_irt/franka/learning/morphology_learning/tests/',...
            'poppy_kin_RAMS:', GRAD_TYPE, ...
            '_LR:', strrep(num2str(learnRate),'.','p'), ...
            '_B:' , num2str(BUFF_SIZE), ...
            '_X:' , num2str(MB_SIZE), ...
            '_D:' , strcat(date,'_',datestr(now,'HH:MM'))], 'lambda_rams','J_log','rho_G_avg')
        pause(3) 
        disp('Saving completed!');
    end
    
    % if SAVE_RESULTS == true
    %     disp('Saving...');
    %     save(['/media/fernando/extreme_ssd/panda_experiments/ip_learning_results/tests/',...
    %         'panda_GD:', GRAD_TYPE, ...
    %         '_LR:', strrep(num2str(learnRate),'.','p'), ...
    %         '_IC:', initialPoint, ...
    %         '_M:' , num2str(MEASUREMENTS), ...
    %         '_B:' , num2str(BUFF_SIZE), ...
    %         '_X:' , num2str(MB_SIZE), ...
    %         '_BM:', num2str(ACTION_TORSO),...
    %         '_N:' , num2str(NOISE),...
    %         '_LK' , num2str(LOAD_KNOWN),...
    %         '_D:' , strcat(date,'_',datestr(now,'HH:MM'))], 'theta_rams','J_log','rho_G_avg')
    %     pause(3) 
    %     disp('Saving completed!');
    % end
    %%
    figure('Color','w')
    for p = 1:13
        subplot(3,5,p)
        plot(BUFF_SIZE:iter,permute(lambda_rams(8:13,p,BUFF_SIZE:iter),[1,3,2]))
        ylabel(['$\boldmath{\lambda}_' num2str(p), '(t)$'],'interpreter','latex')
        xlabel(['Time [s]'])
    end
    %%
    lambda_hat_online  = mean(lambda_rams(:,:,iter-100:iter),3);
end