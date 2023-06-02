% *************************************************************************
% NOTE: This files plots different graphs depending on the subgraph
%       selected from G_MI.
% *************************************************************************

function [out] = fcn_plot_robot_mi_matrix_kinematics_subgraphs(MI_mat, NORMALIZE, PRUNE, ENHANCE, constPar)

    clear cb tx fig ax leg
    clc  
    
    out = struct;
    
    N_joints  = constPar.noj;
    N_links   = constPar.nol;
    
    % Weighted adjacency matrix
    W      = MI_mat - diag(diag(MI_mat));
    % Degree matrix
    D      = diag(sum(W,2));
    
    % * NOTE: Normalizing the adjacency matrix makes its largest eigenvalue 1
    % W_norm is doing some sort of average of the node neighbours while taking 
    % into account their number of neighbours (being connected to a node 
    % connected to all nodes gives less information than if it's connected 
    % only to said node). The square roots make sure that the largest eigenvaue
    % has value 1.
    
    if NORMALIZE == 1
        W_norm = ((D)^(-1/2))*W*((D)^(-1/2));   
        L_norm = ((D)^(-1/2))*(D - W)*((D)^(-1/2));   
	    W_aux  = W_norm;
        W_aux  = (W_aux+transpose(W_aux))/2; % force matrix to be symmetric
    else
        W_aux  = W;
    end
    
    
    kinIndices = [N_joints + 1:2*N_joints, ...
               3*N_joints + 1:3*N_joints + N_links];
    
    W_aux = W_aux(kinIndices,kinIndices);
    
    
    % Plot relationships between joint and Cartesian angular velocities
    cprintf('*yellow', '>> Plotting relationships between joint and Cartesian angular velocities...\n')
    W_aux  = (W_aux + transpose(W_aux))/2; % force matrix to be symmetric
    [T,~]  = minspantree(graph(-W_aux,'upper'));
    T.Edges.Weight = abs(T.Edges.Weight);
    delta  = PRUNE*(min(T.Edges.Weight));
    disp(['The prune threshold is: ', num2str(delta)])
    W_aux(W_aux<delta) = 0;

    if ENHANCE == 1
        % Remove extra edges from omg nodes
        W_omg  = W_aux(constPar.noj+1:end,constPar.noj+1:end);
        [T,~]  = minspantree(graph(-W_omg,'upper'));
        T.Edges.Weight = abs(T.Edges.Weight);
        delta  = PRUNE*(min(T.Edges.Weight));    
        W_omg(W_omg<delta) = 0;
        W_aux(constPar.noj+1:end,constPar.noj+1:end) = W_omg;
    end
    % Normalize entries to sum up to 1
%     W_aux = W_aux / sum(W_aux(:));
    
    out.W_aux = W_aux;
       
    % Plots ---------------------------------------------------------------
    [G_kin, G_kin_mst, c_kin] = fcn_pigraph_kinematics_contracted_mi_with_bodies(1, ...
           W_aux, 1:constPar.noj, 1:constPar.nol, 'force', constPar.showClusters,constPar);    
    
    out.G_kin     = G_kin;
    out.G_kin_mst = G_kin_mst;
    out.c_kin     = c_kin;

end

