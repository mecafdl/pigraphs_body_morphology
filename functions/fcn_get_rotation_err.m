function cost = fcn_get_rotation_err(omg_p, omg_c, q_c, dq_c, x)
    xi   = [x(1);x(2);x(3)]/norm([x(1);x(2);x(3)]);%w_pc  = [x(1);x(2);x(3)];
    phi  = x(4);
    zeta = [x(5);x(6);x(7)]/norm([x(5);x(6);x(7)]);%w_cz  = [x(5);x(6);x(7)];

    samples  =  size(omg_c,2);

    omg_c_hat = zeros(size(omg_c));
    err       = zeros(size(omg_c));
    
    for i=1:samples
%         R_p_c0 = expm(skew(xi)*(phi));
%         R_c0_c = expm(skew(zeta)*q_c(i));
%         R_p_c  = R_p_c0*R_c0_c;

%         R_c_p = expm(transpose(skew(zeta))*q_c(i))*expm(transpose(skew(xi))*phi);
        R_c_p = expm(-skew(zeta)*q_c(i))*expm(-skew(xi)*phi);
        
        omg_c_hat(:,i) = R_c_p*omg_p(:,i) + dq_c(i)*zeta;
        err(:,i)       = omg_c(:,i) - omg_c_hat(:,i); 
    end
    
    cost = (1/samples)*(1/2)*sum(sum(err.^2));
    
end

% function cost = fcn_get_rotation_err(omg_p, omg_c, q_c, dq_c, x)
%     w_pc  = [x(1);x(2);x(3)]/norm([x(1);x(2);x(3)]);%w_pc  = [x(1);x(2);x(3)];
%     theta = x(4);
%     w_cz  = [x(5);x(6);x(7)]/norm([x(5);x(6);x(7)]);%w_cz  = [x(5);x(6);x(7)];
% 
%     samples  =  size(omg_c,2);
% 
%     omg_c_hat = zeros(size(omg_c));
%     err       = zeros(size(omg_c));
%     
%     for i=1:samples
%         R_p_c0 = expm(skew(w_pc)*(theta));
%         R_c0_c = expm(skew(w_cz)*q_c(i));
%         R_p_c  = R_p_c0*R_c0_c;
%         
%         omg_c_hat(:,i) = (R_p_c)'*omg_p(:,i) + dq_c(i)*w_cz;
%         err(:,i)       = omg_c(:,i) - omg_c_hat(:,i); 
%     end
%     
%     cost = (1/samples)*(1/2)*sum(sum(err.^2));
%     
% end