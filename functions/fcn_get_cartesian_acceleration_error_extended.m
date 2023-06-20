% function cost = fcn_get_cartesian_acceleration_err(p_v_p,c_v_c,p_w_p,c_w_c, p_R_c, x)
function cost = fcn_get_cartesian_acceleration_error_extended(p_a_p, c_a_c, p_w_p, ...
                                                   c_w_c, p_dw_p, c_dw_c, q_c, p_R_c, zeta,...
                                                   beta, x)
    samples  =  size(p_a_p,2);
    err      = zeros(3,samples);
    p_r_j    = [x(1);x(2);x(3)];
    c_r_j    = [x(4);x(5);x(6)];
    for i=1:samples
        err(:,i) = p_a_p(:,i) - skew(p_w_p(:,i))*skew(p_w_p(:,i))*p_r_j - skew(p_dw_p(:,i))*p_r_j - ...
                 p_R_c(q_c(i), zeta)*(c_a_c(:,i) - skew(c_w_c(:,i))*skew(c_w_c(:,i))*c_r_j - skew(c_dw_c(:,i))*c_r_j);
    end    
    cost = (1/samples)*(1/2)*sum(sum(err.^2)) + ...
           beta*(p_r_j')*p_r_j + ...
           beta*(c_r_j')*c_r_j;
    
end