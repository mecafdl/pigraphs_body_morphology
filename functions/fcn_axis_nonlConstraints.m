function  [c,ceq] = fcn_axis_nonlConstraints(x)
    % Nonlinear inequality constraints function 
    % c <= 0

    c = [];

    % Nonlinear equality constraints function 
    % c = 0

    ceq1 = 1 - norm(x(1:3));
    ceq2 = 1 - norm(x(5:7));
    ceq = [ceq1;ceq2];
end
