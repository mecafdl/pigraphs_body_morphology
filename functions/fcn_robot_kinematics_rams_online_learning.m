% *************************************************************************
%                         RAMS Gradient Descent                           *
% *************************************************************************

function [lambda_hat_online, lambda_rams_mav, J_log] = ...
    fcn_robot_kinematics_rams_online_learning(signals, ...
                                              lambda_hat_0, ...
                                              IS_ACC, ...
                                              A_kin_pi, ...
                                              epochs, ...
                                              buffer, ...
                                              batchSize, ...
                                              samplingTime, ...
                                              constPar)

    assert(isstruct(constPar),'constPar needs to be struct')
    
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
    signals_maf(ang_acc_indices,:) = gradient(signals_maf(ang_vel_indices,:), samplingTime);
    
    % Define the data pool
    dataPool    = transpose(signals_maf);
    
    % Manifold settings
    S2_man      = spherefactory(2);
    
    % Cost function and gradient settings
    GRAD_TYPE    = 'Riemannian';
    
    % Riemannian AMS Initial settings =====================================
    m         = zeros(size(lambda_hat_0,1),constPar.noj);
    v         = zeros(constPar.noj, 1);
    v_hat     = zeros(constPar.noj, 1);    
    tau       = zeros(size(lambda_hat_0,1),constPar.noj);
    eta_max   = 0.01*ones(constPar.noj,1);
    % =====================================================================
    
    % Parameter update initial value
    lambda_update = repmat(lambda_hat_0,1,constPar.noj);
    % Settings for RAMS gradient descent
    lambda_rams = zeros(size(lambda_hat_0,1), constPar.noj, epochs);
    J_log       = zeros(epochs, 1);
    rho_G_avg   = zeros(epochs, 1);
    SIMPLIFIED  = false;
    
    clc
    cprintf('*yellow', '>> Running %s GD with:\n', GRAD_TYPE)
    cprintf('*yellow', '>> - EPOCHS      : %d\n', epochs)
    cprintf('*yellow', '>> - MB SIZE     : %d\n',batchSize)
    cprintf('*yellow', '>> - BUFFER SIZE : %d\n',buffer.Size)
    disp('>> Initial Point:')
    disp(lambda_update);
    tic
    pause(3)
    
    first_sample = 0;
    flag         = 0;
    update_grad  = 1;
    for iter = 1:epochs       
        % Store sample in replay buffer
        buffer.storeSample(transpose(dataPool(iter,:)));
    
        %if iter >= buffer.Size && update_grad == 1  
        if buffer.Status == 1 && update_grad == 1
            if flag == 0
                first_sample = first_sample + 1;
                flag         = 1; 
            end
            % Calculate cost function and Riemannian gradients ============
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
                % Draw samples from buffer ================================
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
                    
                    %c_OMG_c_hat = expm(-skew(lambda_update(5:7,jnt))*q_n)*expm(-skew(lambda_update(1:3,jnt))*lambda_update(4,jnt))*p_OMG_p_n ...
                    %                               + dq_n*lambda_update(5:7,jnt);
    
                    %beta   = 0.1;
                    %beta   = 0.5;
                    %beta   = 0.01;
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
                % Average gradient per batch ==============================
                dJ_lambda_mb   = mean(dJ_lambda_tmp(:,:,jnt),2);

                % Batch Riemannian gradient ===============================                   
                rgrad_J(:,jnt) = [S2_man.egrad2rgrad(lambda_update(1:3,jnt), dJ_lambda_mb(1:3)); 
                                  dJ_lambda_mb(4);
                                  S2_man.egrad2rgrad(lambda_update(5:7,jnt), dJ_lambda_mb(5:7));
                                  dJ_lambda_mb(8:13)];
            end

            % Calculate parameter update ==================================
            switch GRAD_TYPE
                case 'Euclidean'
                    % Based on Euclidean AMS gradient descent
                    warning('Section undefined')
                    break
                 case 'Riemannian'
                    % Based on RAMS gradient descent
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
    
    % KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
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
    % KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
    
        % Display progress
        if buffer.Status == 0%iter <= buffer.Size
            warning('No yet enough samples in the replay buffer!')
        elseif mod(iter,1) == 0
            disp(['Iter: ', num2str(iter),'/',num2str(epochs),' | Cost: ', num2str(J_log(iter)), ' | Grad norm: ' num2str(rho_G_avg(iter)), ' | eta_avg: ', num2str(eta_avg)]);
        end
    end
    
    cprintf('*yellow', '>> Computation finished\n')
    
    % Moving average filter
    windowSize      = 100;
    B               = 1/windowSize*ones(1,windowSize);
    A               = 1; %denominator coefficients
    %lambda_rams_mav = filter(B,A,permute(lambda_rams(:,:,buffer.Size:iter),[3,1,2])); %filter input x and get result in y
    lambda_rams_mav = filter(B,A,permute(lambda_rams(:,:,first_sample:iter),[3,1,2])); %filter input x and get result in y
    lambda_rams_mav = permute(lambda_rams_mav,[2,3,1]);
    
    figure('Color','w')
    for p = 1:13
        subplot(3,5,p)
        plot(1:size(lambda_rams_mav,3),squeeze(lambda_rams_mav(p,:,:)))
        ylabel(['$\boldmath{\lambda}_{' num2str(p), '}(t)$'],'interpreter','latex')
        xlabel('Samples')
        xlim([0 iter])
    end
    % The final values are the avergage of the last <windowSize> filtered points
    lambda_hat_online  = mean(lambda_rams_mav(:,:,end-windowSize:end),3);

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