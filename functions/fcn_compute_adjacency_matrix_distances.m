function [spectral_distance, frobenius_distance, compound_distance] = ...
    fcn_compute_adjacency_matrix_distances(W_previous, W_current)


    [L_previous_norm, ~] = fcn_get_normalized_laplacian(W_previous);
    [L_current_norm, ~] = fcn_get_normalized_laplacian(W_current);
    % Spectral distance
    eigW_previous    = sort(eig(L_previous_norm));
    eigW_previous    = eigW_previous./sum(eigW_previous);
    eigW_current     = sort(eig(L_current_norm));
    eigW_current     = eigW_current./sum(eigW_current);
    spectral_distance = norm(eigW_current - eigW_previous,'fro')/norm(eigW_previous,'fro');





%     % Spectral distance
%     eigW_previous    = sort(eig(W_previous));
%     eigW_current     = sort(eig(W_current));
%     %err_spectral      = eigW_previous - eigW_current;
%     %spectral_distance = sqrt(transpose(err_spectral)*err_spectral)^2;
%     spectral_distance = norm(eigW_current - eigW_previous,'fro')/norm(eigW_previous,'fro');
    % Frobenius distance
    %err_matrix         = W_previous(:) - W_current(:);
    %frobenius_distance = sqrt(transpose(err_matrix)*err_matrix)^2;
    frobenius_distance = norm(W_current - W_previous,'fro')/norm(W_previous,'fro');
    compound_distance  = frobenius_distance^2 + spectral_distance^2;%norm([spectral_distance;frobenius_distance]).^2;
end


function [L_norm, W_norm] = fcn_get_normalized_laplacian(W)
    D         = diag(sum(W,2));
    W_norm    = ((D)^(-1/2))*W*((D)^(-1/2));   
    L_norm    = ((D)^(-1/2))*(D - W)*((D)^(-1/2));   
end