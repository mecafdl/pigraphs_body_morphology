function MI_ctrct = fcn_contract_mi_matrix(MI_full, constPar)
% Contract the MI matrix by getting the Frobenius norm of the angular
% velocities/linear acceleration 3x3 blocks and 3x1 blocks 
% * NOTE: this assumes that the entries of the MI matrix follow the 
%         following order:
%         - q, dq, tau, omg, acc     
    
    if ~(size(MI_full,1) == size(MI_full,2))
        error('Matrix is no square')
    end
    if isfield(constPar,'nob') && isfield(constPar,'noj')
        N_signals = 3*constPar.noj + 6*constPar.nob; 
    else
        error('constPar structure does not contain required nob and noj fields')
    end

%     assert(isa(MI_full,'double') && size(MI_full, 1) == N_signals ...
%         && size(MI_full,2) == N_signals && isreal(MI_full), ...
%       ['MI matrix needs to be ', num2str(N_signals),'x', num2str(N_signals), ' double']); 



    if size(MI_full,1) == N_signals
        MI_ctrct = zeros(3*constPar.noj + 2*constPar.nob);
    
        % Copy the scalar entries of the MI matrix ----------------------------
        MI_ctrct(1:3*constPar.noj,1:3*constPar.noj) = 1*MI_full(1:3*constPar.noj,1:3*constPar.noj);
    
        % Contract blocks based on vector/matrix norm -------------------------
        for i=1:3*constPar.noj
            for j=1:2*constPar.nob
                c_range = ((3*j-2):3*j) + 3*constPar.noj;
                MI_ctrct(i, j + 3*constPar.noj) = norm(MI_full(i,c_range),'fro');
                MI_ctrct(j + 3*constPar.noj, i) = norm(MI_full(i,c_range),'fro');
            end
        end
    
        for i = 1:2*constPar.nob
            r_range = ((3*i-2):3*i) + 3*constPar.noj;
            for j=1:2*constPar.nob
                if i==j
                    continue;
                else
                    c_range                            = ((3*j-2):3*j) + 3*constPar.noj;
                    MI_ctrct(i+3*constPar.noj,j+3*constPar.noj) = norm(MI_full(r_range,c_range),'fro');
                end
            end
        end    
    elseif size(MI_full) == 3*constPar.noj + 2*6*constPar.nob
        MI_ctrct = zeros(3*constPar.noj + 4*constPar.nob);
    
        % Copy the scalar entries of the MI matrix ----------------------------
        MI_ctrct(1:3*constPar.noj,1:3*constPar.noj) = 1*MI_full(1:3*constPar.noj,1:3*constPar.noj);
    
        % Contract blocks based on vector/matrix norm -------------------------
        for i=1:3*constPar.noj
            for j=1:2*constPar.nob
                c_range = ((3*j-2):3*j) + 3*constPar.noj;
                MI_ctrct(i, j + 3*constPar.noj) = norm(MI_full(i,c_range),'fro');
                MI_ctrct(j + 3*constPar.noj, i) = norm(MI_full(i,c_range),'fro');
            end
        end
    
        for i = 1:4*constPar.nob
            r_range = ((3*i-2):3*i) + 3*constPar.noj;
            for j=1:4*constPar.nob
                if i==j
                    continue;
                else
                    c_range                            = ((3*j-2):3*j) + 3*constPar.noj;
                    MI_ctrct(i+3*constPar.noj,j+3*constPar.noj) = norm(MI_full(r_range,c_range),'fro');
                end
            end
        end  

    end
    
    % Ensure matrix to be symmetric ---------------------------------------
    MI_ctrct  = (MI_ctrct + transpose(MI_ctrct))/2; 
end


% function MI_ctrct = fcn_contract_mi_matrix(MI_full, constPar)
% % Contract the MI matrix by getting the Frobenius norm of the angular
% % velocities/linear acceleration 3x3 blocks and 3x1 blocks 
% % * NOTE: this assumes that the entries of the MI matrix follow the 
% %         following order:
% %         - q, dq, tau, omg, acc     
% 
%     if isfield(constPar,'nob') && isfield(constPar,'noj')
%         N_signals = 3*constPar.noj + 6*constPar.nob; 
%     else
%         error('constPar structure does not contain required nob and noj fields')
%     end
% 
%     assert(isa(MI_full,'double') && size(MI_full, 1) == N_signals ...
%         && size(MI_full,2) == N_signals && isreal(MI_full), ...
%       ['MI matrix needs to be ', num2str(N_signals),'x', num2str(N_signals), ' double']); 
% 
%     MI_ctrct = zeros(3*constPar.noj + 2*constPar.nob);
% 
%     % Copy the scalar entries of the MI matrix ----------------------------
%     MI_ctrct(1:3*constPar.noj,1:3*constPar.noj) = 1*MI_full(1:3*constPar.noj,1:3*constPar.noj);
% 
%     % Contract blocks based on vector/matrix norm -------------------------
%     for i=1:3*constPar.noj
%         for j=1:2*constPar.nob
%             c_range = ((3*j-2):3*j) + 3*constPar.noj;
%             MI_ctrct(i, j + 3*constPar.noj) = norm(MI_full(i,c_range),'fro');
%             MI_ctrct(j + 3*constPar.noj, i) = norm(MI_full(i,c_range),'fro');
%         end
%     end
% 
%     for i = 1:2*constPar.nob
%         r_range = ((3*i-2):3*i) + 3*constPar.noj;
%         for j=1:2*constPar.nob
%             if i==j
%                 continue;
%             else
%                 c_range                            = ((3*j-2):3*j) + 3*constPar.noj;
%                 MI_ctrct(i+3*constPar.noj,j+3*constPar.noj) = norm(MI_full(r_range,c_range),'fro');
%             end
%         end
%     end    
%     
%     % Ensure matrix to be symmetric ---------------------------------------
%     MI_ctrct  = (MI_ctrct + transpose(MI_ctrct))/2; 
% end