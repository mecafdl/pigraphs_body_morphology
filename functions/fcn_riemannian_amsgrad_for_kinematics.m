function [m, v, v_hat, tau, lambda_update, rho_G_avg, eta_avg] = ...
    fcn_riemannian_amsgrad_for_kinematics( ...
    m, v, v_hat, tau, eta_max, rgrad_J, P, update_grad, ...
    Mfd, GD_METHOD ,SIMPLE, constPar)


    % Constants ***********************************************************
    NOJ              = constPar.noj;
    N_PARAM_PER_LINK = 13;
    
    % Space allocation for variables **************************************
    eta_adaptive   = zeros(NOJ,1);
    P_update       = zeros(N_PARAM_PER_LINK,NOJ);
    rho_G_agg      = 0;
    

    % Loop over JOINTS ****************************************************
    for jnt = 1:NOJ
       
        P_jnt     = P(:,jnt);
        G_jnt     = rgrad_J(:,jnt);
                
        
        if ~update_grad
            P_update(:,jnt) = P_jnt;
            rho_G_agg     = 0;
        else
            
            %GD_METHOD = 'RSGD';
            
            % Riemannian gradient of the i-th manifold
            
            switch GD_METHOD
                case 'RSGD'
                % Standard Riemannian Stochastic Gradient Descent
                    alpha             = 0.001;
                    rho_G_jnt(1)      = Mfd.norm(P_jnt(1:3), G_jnt(1:3));
                    rho_G_jnt(2)      = Mfd.norm(P_jnt(4), G_jnt(4));
                    rho_G_jnt(3)      = Mfd.norm(P_jnt(5:7), G_jnt(5:7));
                    rho_G_jnt(4)      = norm(P_jnt(8:13));
                    
                    
                    rho_G_agg         = rho_G_agg + mean(rho_G_jnt);                    
                    P_update(1:3,jnt) = Mfd.exp(P_jnt(1:3), -alpha*G_jnt(1:3));
                    P_update(4,jnt)   = P_jnt(4) - alpha*G_jnt(4);
                    P_update(5:7,jnt) = Mfd.exp(P_jnt(5:7), -alpha*G_jnt(5:7));                         
                    P_update(8:13,jnt)   = P_jnt(8:13) - alpha*G_jnt(8:13);
                    
                case 'RAMS'
                % Riemannian AMSGrad            
                    beta_1          = 0.9;
                    beta_2          = 0.999;%1 - (1/k);
                    epsilon         = 1E-8;

                    rho_G_jnt(1)    = Mfd.norm(P_jnt(1:3), G_jnt(1:3));
                    rho_G_jnt(2)    = Mfd.norm(P_jnt(4), G_jnt(4));
                    rho_G_jnt(3)    = Mfd.norm(P_jnt(5:7), G_jnt(5:7));
                    rho_G_jnt(4)    = norm(P_jnt(8:13));
                    rho_G_agg       = rho_G_agg + mean(rho_G_jnt);
                    
                    % - Update biased 1st moment estimate
                    m(1:3,jnt)        = beta_1*tau(1:3,jnt) + (1 - beta_1)*G_jnt(1:3);  
                    m(4,jnt)          = beta_1*tau(4,jnt)   + (1 - beta_1)*G_jnt(4);  
                    m(5:7,jnt)        = beta_1*tau(5:7,jnt) + (1 - beta_1)*G_jnt(5:7);
                    m(8:13,jnt)       = beta_1*tau(8:13,jnt)   + (1 - beta_1)*G_jnt(8:13);  

                    % - Update biased 2nd raw moment estimate 
                    v(jnt)            = beta_2*v(jnt) + (1 - beta_2)*mean(rho_G_jnt);
                    v_hat(jnt)        = max(v_hat(jnt), v(jnt));
                    eta_adaptive(jnt) = eta_max(jnt)/(sqrt(v_hat(jnt)) + epsilon);

                    if ~SIMPLE                        
                        P_update(1:3,jnt)  = Mfd.exp(P_jnt(1:3),m(1:3,jnt),-eta_adaptive(jnt));                            
                        P_update(4,jnt)    = P_jnt(4) - eta_adaptive(jnt)*m(4,jnt);
                        P_update(5:7,jnt)  = Mfd.exp(P_jnt(5:7),m(5:7,jnt),-eta_adaptive(jnt));
                        P_update(8:13,jnt) = P_jnt(8:13) - eta_adaptive(jnt)*m(8:13,jnt);
                        
                        tau(1:3,jnt)       = Mfd.transp(P_jnt(1:3), P_update(1:3,jnt), m(1:3,jnt));
                        tau(4,jnt)         = m(4,jnt);
                        tau(5:7,jnt)       = Mfd.transp(P_jnt(5:7), P_update(5:7,jnt), m(5:7,jnt));
                        tau(8:13,jnt)      = m(8:13,jnt);
                    else
                        % Poor man's vector transport: exploit the fact that all tangent spaces
                        % are the set of symmetric matrices, so that the identity is a sort of
                        % vector transport. It may perform poorly if the origin and target (X1
                        % and X2) are far apart though. This should not be the case for typical
                        % optimization algorithms, which perform small steps.
                        P_update(1:3,jnt)  = Mfd.retr(P_jnt(1:3),m(1:3,jnt),-eta_adaptive(jnt));                            
                        P_update(4,jnt)    = P_jnt(4) - eta_adaptive(jnt)*m(4,jnt);
                        P_update(5:7,jnt)  = Mfd.retr(P_jnt(5:7),m(5:7,jnt),-eta_adaptive(jnt));
                        P_update(8:13,jnt) = P_jnt(8:13) - eta_adaptive(jnt)*m(8:13,jnt);
                        
                        tau(:,jnt)         = m(:,jnt);
                    end
                        
            end          
        end
        
        
        rho_G_avg = rho_G_agg/jnt;
%         % Update inertial parameters ======================================
%         % Density weighted covariance matrix Sigma
%         Sigma_i_update = P_update(1:3,1:3,jnt);
% 
%         % Compute updated inertia matrix based on Sigma
%         I_i_update   = trace(Sigma_i_update)*eye(3) - Sigma_i_update;
% 
%         theta_i_update(:,jnt) = [P_update(4,4,jnt),...
%                                P_update(1,4,jnt),...
%                                P_update(2,4,jnt),...
%                                P_update(3,4,jnt),...
%                                I_i_update(1,1),...
%                                I_i_update(1,2),...
%                                I_i_update(1,3),...
%                                I_i_update(2,2),...
%                                I_i_update(2,3),...
%                                I_i_update(3,3)]';
    end  

    lambda_update = P_update;
    eta_avg       = 1/constPar.noj*sum(eta_adaptive);

end


%     m          = beta1.*m + (1 - beta1) .* (gradJ);       % - Update biased 1st moment estimate
%     v          = beta2.*v + (1 - beta2) .* (gradJ.^2);    % - Update biased 2nd raw moment estimate 
%     vHat       = max(vHat, v);                            % - Compute bias-corrected 2nd raw moment estimate
%     learn_rate = alpha./(sqrt(vHat) + epsilon);
%     theta_hat  = theta_hat - learn_rate.*m;                            