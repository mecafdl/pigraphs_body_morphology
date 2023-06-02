% *************************************************************************
% NOTE: This files plots different graphs depending on the subgraph
%       selected from G_MI.
% *************************************************************************

function [out] = fcn_plot_robot_mi_matrix_graphs(MI_mat, SIGNAL_CHOICE, NORMALIZE, PRUNE, constPar)

% close all; 
clear cb tx fig ax leg
clc  

subindices = ['0ABCDEFGHIJKLMNOPQRSTUVWXYZ'];


out = struct;


% CONTROL FLAGS ===========================================================
ENHANCE      = 0;
% PRUNE        = 1;
CONTRACT     = 1;
flag_vec     = [NORMALIZE,ENHANCE,PRUNE,CONTRACT];


N_joints  = constPar.noj;
N_links   = constPar.nol;
w_indices = [[0:N_joints]',(reshape(1:(N_joints+1)*3,3,N_joints+1))' + N_joints];

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

% switch SIGNAL_CHOICE
%     case 1 %'[omg]'
%         indices = [3*N_joints + 1:3*N_joints + N_links];
%     case 2 %'[joint_vel, omg]'
%         indices = [N_joints + 1:2*N_joints, ...
%                    3*N_joints + 1:3*N_joints + N_links];
%     case 3 %'[joint_torque, acc]'
%         indices = [2*N_joints + 1:3*N_joints, ...
%                    3*N_joints + N_links + 1:3*N_joints + 2*N_links];        
%     case 4 %'[all]'
%         indices = 1:size(W_aux,1);
%     case 5 %'[acc]'
%         indices = [30:37];
%     case 6 %'[acc]'
%         indices = [N_joints + 1:2*N_joints, ...
%                    2*N_joints + 1:3*N_joints];
%     case 7 %'[acc]'
%         indices = [1:N_joints, ...
%                    2*N_joints + 1:3*N_joints];               
% end
switch SIGNAL_CHOICE
    case 'omg'
        indices = [3*N_joints + 1:3*N_joints + N_links];
    case 'dq-omg'
        indices = [N_joints + 1:2*N_joints, ...
                   3*N_joints + 1:3*N_joints + N_links];
    case 'tau-omg'
        indices = [2*N_joints + 1:3*N_joints, ...
                   3*N_joints + 1:3*N_joints + N_links];
    case 'acc'
        indices = [3*N_joints + N_links + 1: 3*N_joints + 2*N_links];        
    case 'tau-acc'
        indices = [2*N_joints + 1:3*N_joints, ...
                   3*N_joints + N_links + 1:3*N_joints + 2*N_links];        
    case 'all'
        indices = 1:size(W_aux,1);
    case 5 %'[acc]'
        indices = [30:37];
    case 'dq-tau'
        indices = [N_joints + 1:2*N_joints, ...
                   2*N_joints + 1:3*N_joints];
    case 'q_dq'
        indices = [1:N_joints, ...
                   2*N_joints + 1:3*N_joints];               
end
W_aux = W_aux(indices,indices);

if strcmp(SIGNAL_CHOICE,'omg')
    % Plot only relationships between the Cartesian angular velocities
    cprintf('*yellow', '>> Plotting relationships between the Cartesian angular velocities...\n')
    %W_aux = W_aux / sum(W_aux(:));
    [T,~]          = minspantree(graph(-W_aux,'upper'));
    T.Edges.Weight = abs(T.Edges.Weight);
    delta          = PRUNE*(min(T.Edges.Weight));
    disp(['The prune threshold is: ', num2str(delta)])
    %figure;bar((sort(W_aux(:)))); hold on;plot(1:numel(W_aux),delta*ones(1,numel(W_aux)),'r--')
    W_aux(W_aux<delta) = 0; 
    close all
    [out.G_omg_mst] = fcn_get_robot_cart_ang_vel_graph(W_aux, 'omg', constPar);
    
elseif strcmp(SIGNAL_CHOICE,'acc')
    % Plot only relationships between the Cartesian angular velocities
    cprintf('*yellow', '>> Plotting relationships between the Cartesian angular velocities...\n')
    %W_aux = W_aux / sum(W_aux(:));
    [T,~]          = minspantree(graph(-W_aux,'upper'));
    T.Edges.Weight = abs(T.Edges.Weight);
    delta          = PRUNE*(min(T.Edges.Weight));
    disp(['The prune threshold is: ', num2str(delta)])
    %figure;bar((sort(W_aux(:)))); hold on;plot(1:numel(W_aux),delta*ones(1,numel(W_aux)),'r--')
    W_aux(W_aux<delta) = 0; 
    close all
    [out.G_omg_mst] = fcn_get_robot_cart_ang_vel_graph(W_aux, 'acc', constPar);
    % Plots
    % ---------------------------------------------------------------------
    % h1 =  findall(groot,'Type','figure');
    % close(h1(1))
    % clc
    % SAVE_FIG = 1;
    % if SAVE_FIG == 1
    %     export_fig(fullfile(fileparts(matlab.desktop.editor.getActiveFilename),'figures_gazebo',['pxmarkIV_gazebo_fb_pigraph_' num2str(index)]),'-jpg')
    %     close(gcf);
    % end    

elseif strcmp(SIGNAL_CHOICE,'dq-omg')
    % Plot relationships between joint and Cartesian angular velocities
    cprintf('*yellow', '>> Plotting relationships between joint and Cartesian angular velocities...\n')
    W_aux  = (W_aux + transpose(W_aux))/2; % force matrix to be symmetric
%     out.W_aux = W_aux;
    [T,~]  = minspantree(graph(-W_aux,'upper'));
    T.Edges.Weight = abs(T.Edges.Weight);
    delta  = PRUNE*(min(T.Edges.Weight));
    disp(['The prune threshold is: ', num2str(delta)])
    W_aux(W_aux<delta) = 0;

    % Remove extra edges from omg nodes
    W_omg  = W_aux(constPar.noj+1:end,constPar.noj+1:end);
    [T,~]  = minspantree(graph(-W_omg,'upper'));
    T.Edges.Weight = abs(T.Edges.Weight);
    delta  = PRUNE*(min(T.Edges.Weight));    
    W_omg(W_omg<delta) = 0;
    W_aux(constPar.noj+1:end,constPar.noj+1:end) = W_omg;
    
    % Normalize entries to sum up to 1
    W_aux = W_aux / sum(W_aux(:));

    out.W_aux = W_aux;
       
    % Plots ---------------------------------------------------------------
    PLOT_HEATMAP_FLAG  = 0;
    if(PLOT_HEATMAP_FLAG == 1)
        subindices = ['0ABCDEFGHIJKLMNOPQRSTUVWXYZ'];
%         figure('Color','w','Name','dq -> omg')
%             h                      = heatmap(W_aux(1:constPar.noj,constPar.noj+1:end));
%             h.Colormap             = flipud(autumn);
%             ax                     = gca;
%             ax.XDisplayLabels      = arrayfun(@(index){['\omega_{', subindices(index), '}']}, 1:constPar.nob);
%             ax.YDisplayLabels      = arrayfun(@(index){['dq_{', num2str(index), '}']}, 1:constPar.noj);
%             axp                    = struct(gca);       %you will get a warning
%             axp.Axes.XAxisLocation = 'top';
%             title('\bf W_{(\omega, dq)}')
        figure('Color','w')
            h                      = heatmap(W_aux);
            h.Colormap             = flipud(autumn);
            ax                     = gca;
            ax.XDisplayLabels      = [arrayfun(@(index){['dq_{', num2str(index), '}']}, 1:constPar.noj), arrayfun(@(index){['\omega_{', subindices(index), '}']}, 1:constPar.nob)];
            ax.YDisplayLabels      = [arrayfun(@(index){['dq_{', num2str(index), '}']}, 1:constPar.noj), arrayfun(@(index){['\omega_{', subindices(index), '}']}, 1:constPar.nob)];
            axp                    = struct(gca);       %you will get a warning
            axp.Axes.XAxisLocation = 'top';
            title('\bf W_{kin}')
    
    end    
    % Plots
%     fcn_franka_pigraph_kinematics_contracted_mi(1, W_aux, 1:constPar.noj, 1:constPar.nol, 'Force', false ,constPar)
%     [G_kin, G_kin_mst] = fcn_pigraph_kinematics_contracted_mi_generic(1, ...
%            W_aux, 1:constPar.noj, 1:constPar.nol, 'Force', true ,constPar);

    [G_kin, G_kin_mst, c_kin] = fcn_pigraph_kinematics_contracted_mi_with_bodies(1, ...
           W_aux, 1:constPar.noj, 1:constPar.nol, 'force', constPar.showClusters,constPar);    
    
    out.G_kin     = G_kin;
    out.G_kin_mst = G_kin_mst;
    out.c_kin     = c_kin;
    %         qIndex = 3:7;
%         wIndex = (4:8)+7;
%         W_aux  = W_ctrct([qIndex,wIndex],[qIndex,wIndex]);
%         W_aux  = W_aux./sum(W_aux(:));
%         fcn_franka_pigraph_kinematics_contracted_mi_generic(1, ...
%             W_aux, qIndex, [4:8], 'Force', false ,constPar)
    SAVE_FIG = 0;
    if SAVE_FIG == 1
        h1 =  findall(groot,'Type','figure');
        close(h1(1))
        clc
        export_fig(fullfile(fileparts(matlab.desktop.editor.getActiveFilename),'figures','panda_pi_graph_kinematics_simulated'),'-pdf')
        close(gcf);
    end


elseif strcmp(SIGNAL_CHOICE,'tau-omg')
    % Plot relationships between joint torque and Cartesian angular velocities
    cprintf('*yellow', '>> Plotting relationships between joint and Cartesian linear accelerations...\n')
    W_aux  = (W_aux + transpose(W_aux))/2; % force matrix to be symmetric
    [T,~]  = minspantree(graph(-W_aux,'upper'));
    T.Edges.Weight = abs(T.Edges.Weight);
    delta  = PRUNE*(min(T.Edges.Weight));
    disp(['The prune threshold is: ', num2str(delta)])
    W_aux(W_aux<delta) = 0;    
    
    % Normalize entries to sum up to 1
    W_aux = W_aux / sum(W_aux(:));
       
    % Plots
    PLOT_HEATMAP_FLAG = 1;
    if(PLOT_HEATMAP_FLAG)
        figure('Color','w','Name','MI_heatmap')
            h                      = heatmap(W_aux(1:constPar.noj,constPar.noj+1:end));
            h.Colormap             = flipud(autumn);
            ax                     = gca;
            %ax.XDisplayLabels      = [arrayfun(@(index){['\tau_', num2str(index)]}, 1:constPar.noj),arrayfun(@(index){['a_', subindices(index)]}, 1:constPar.nob)];
            %ax.YDisplayLabels      = [arrayfun(@(index){['\tau_', num2str(index)]}, 1:constPar.noj),arrayfun(@(index){['a_', subindices(index)]}, 1:constPar.nob)];
            ax.XDisplayLabels      = [arrayfun(@(index){['\omega_{', subindices(index), '}']}, 1:constPar.nob)];
            ax.YDisplayLabels      = [arrayfun(@(index){['\tau_{', num2str(index), '}']}, 1:constPar.noj)];

            axp                    = struct(gca);       %you will get a warning
            axp.Axes.XAxisLocation = 'top';
            title('\bf W^{\tau, a}')
    end    
    
%     fcn_poppy_pigraph_dynamics_contracted_mi_generic(1, ...
%         W_aux, 1:constPar.noj, 1:constPar.nol, 'Force', false ,constPar)

    [G_dyn, G_dyn_mst] = fcn_pigraph_contracted_mi_scalar2vector_signals(1, ...
           W_aux, 1:constPar.noj, 1:constPar.nol, 'tau', 'omg','Force', false ,constPar);

%  [G_cntr, G_mst] = fcn_pigraph_contracted_mi_scalar2vector_signals(PLOT_MST, ...
%             W_cntr, qIndex, wIndex, scalar_signal, vector_signal, myLayout_type, SHOW_CLUSTERS, constPar)    
    
    out.G_dyn     = G_dyn;
    out.G_dyn_mst = G_dyn_mst;    
    
    SAVE_FIG = 0;
    if SAVE_FIG == 1
        h1 =  findall(groot,'Type','figure');
        close(h1(1))
        clc
        export_fig(fullfile(fileparts(matlab.desktop.editor.getActiveFilename),'figures','panda_pi_graph_dynamics_simulated'),'-pdf')
        close(gcf);
    end  

elseif strcmp(SIGNAL_CHOICE,'tau-acc')
    % Plot relationships between joint torque and Cartesian angular velocities
    cprintf('*yellow', '>> Plotting relationships between joint and Cartesian linear accelerations...\n')
    W_aux  = (W_aux + transpose(W_aux))/2; % force matrix to be symmetric
    [T,~]  = minspantree(graph(-W_aux,'upper'));
    T.Edges.Weight = abs(T.Edges.Weight);
    delta  = PRUNE*(min(T.Edges.Weight));
    disp(['The prune threshold is: ', num2str(delta)])
    W_aux(W_aux<delta) = 0;    
    
    % Normalize entries to sum up to 1
    W_aux = W_aux / sum(W_aux(:));
       
    % Plots
    PLOT_HEATMAP_FLAG = 1;
    if(PLOT_HEATMAP_FLAG)
        figure('Color','w','Name','MI_heatmap')
            h                      = heatmap(W_aux(1:constPar.noj,constPar.noj+1:end));
            h.Colormap             = flipud(autumn);
            ax                     = gca;
            %ax.XDisplayLabels      = [arrayfun(@(index){['\tau_', num2str(index)]}, 1:constPar.noj),arrayfun(@(index){['a_', subindices(index)]}, 1:constPar.nob)];
            %ax.YDisplayLabels      = [arrayfun(@(index){['\tau_', num2str(index)]}, 1:constPar.noj),arrayfun(@(index){['a_', subindices(index)]}, 1:constPar.nob)];
            ax.XDisplayLabels      = [arrayfun(@(index){['a_{', subindices(index), '}']}, 1:constPar.nob)];
            ax.YDisplayLabels      = [arrayfun(@(index){['\tau_{', num2str(index), '}']}, 1:constPar.noj)];

            axp                    = struct(gca);       %you will get a warning
            axp.Axes.XAxisLocation = 'top';
            title('\bf W^{\tau, a}')
    end    
    
%     fcn_poppy_pigraph_dynamics_contracted_mi_generic(1, ...
%         W_aux, 1:constPar.noj, 1:constPar.nol, 'Force', false ,constPar)

    [G_dyn, G_dyn_mst] = fcn_pigraph_dynamics_contracted_mi_generic(1, ...
           W_aux, 1:constPar.noj, 1:constPar.nol, 'Force', false ,constPar);
    
    out.G_dyn     = G_dyn;
    out.G_dyn_mst = G_dyn_mst;    
    
    SAVE_FIG = 0;
    if SAVE_FIG == 1
        h1 =  findall(groot,'Type','figure');
        close(h1(1))
        clc
        export_fig(fullfile(fileparts(matlab.desktop.editor.getActiveFilename),'figures','panda_pi_graph_dynamics_simulated'),'-pdf')
        close(gcf);
    end  
    
elseif strcmp(SIGNAL_CHOICE,'all')
    W_ctrct = W_aux;
    W_ctrct = (W_ctrct + transpose(W_ctrct))/2; % force matrix to be symmetric
    [T,~]   = minspantree(graph(-W_ctrct,'upper'));
    T.Edges.Weight = abs(T.Edges.Weight);
    delta  = PRUNE*(min(T.Edges.Weight));
    disp(['The prune threshold is: ', num2str(delta)])
    W_ctrct(W_ctrct<delta) = 0;
    
    % Normalize entries to sum up to 1
    W_ctrct = W_ctrct / sum(W_ctrct(:));

    % Plots
    PLOT_HEATMAP_FLAG = 0;
    if(PLOT_HEATMAP_FLAG)
        figure('Color','w','Name','MI_heatmap')
            h                      = heatmap(W_aux(1:constPar.noj,constPar.noj+1:end));
            h.Colormap             = flipud(autumn);
            ax                     = gca;
            %ax.XDisplayLabels      = [arrayfun(@(index){['\tau_', num2str(index)]}, 1:constPar.noj),arrayfun(@(index){['a_', subindices(index)]}, 1:constPar.nob)];
            %ax.YDisplayLabels      = [arrayfun(@(index){['\tau_', num2str(index)]}, 1:constPar.noj),arrayfun(@(index){['a_', subindices(index)]}, 1:constPar.nob)];
            ax.XDisplayLabels      = [arrayfun(@(index){['a_{', subindices(index), '}']}, 1:constPar.nob)];
            ax.YDisplayLabels      = [arrayfun(@(index){['\tau_{', num2str(index), '}']}, 1:constPar.noj)];

            axp                    = struct(gca);       %you will get a warning
            axp.Axes.XAxisLocation = 'top';
            title('\bf W^{\tau, a}')
    end    


%         W_dpi  = fcn_phantomx_plot_contracted_MI(W_ctrct, 0, constPar);   
%         close all
    wIndex = [0, 7, 8, 9, 4, 5, 6, 1, 2, 3, 16, 17, 18, 13, 14, 15, 10, 11, 12];
    qIndex = [7 8 9 1 2 3 4 5 6 16 17 18 10 11 12 13 14 15];
%         fcn_phantomx_pigraph_contracted_mi(1, W_ctrct, qIndex, wIndex, 'Force')
    [G_mst, idx] = fcn_robot_contracted_mi_pigraph(1, W_ctrct, 1:constPar.noj, 1:constPar.nol, 'Force', false ,constPar);
    out.G_pi = G_mst;
    return
%             fcn_poppy_pigraph_contracted_mi(PLOT_MST, W_cntr, qIndex, wIndex, myLayout_type, show_clusters, constPar)
%         G_K = subgraph(T, find(idx == 1));
%         G_D = subgraph(T, find(idx == 2));        

%         idx = fcn_franka_pigraph_contracted_mi_no_joint_angle(1, W_ctrct(constPar.noj+1:end,constPar.noj+1:end), 1:constPar.noj, 1:constPar.nol, 'Force', false ,constPar);
%         G_K = subgraph(graph(W_ctrct(constPar.noj+1:end,constPar.noj+1:end)), find(idx == 1));
%         G_D = subgraph(graph(W_ctrct(constPar.noj+1:end,constPar.noj+1:end)), find(idx == 2));
    
    G_K = subgraph(graph(W_ctrct), [8:14,21+1:21+8]);
    G_D = subgraph(graph(W_ctrct), [15:21,29+1:29+8]);        
    
%         figure
%         fcn_franka_pigraph_kinematics_contracted_mi(1, adjacency(G_K,'weighted'), 1:constPar.noj, 1:constPar.nol, 'Force', true, constPar)
%         fcn_franka_pigraph_kinematics_contracted_mi(1, W_ctrct([8:14,22:29],[8:14,22:29]), 1:constPar.noj, 1:constPar.nol, 'Force', false ,constPar,{'\tau','\dot{v}'})
%         
fcn_franka_pigraph_kinematics_contracted_mi_generic(1, ...
    adjacency(G_K,'weighted'), 1:constPar.noj, 1:constPar.nol, 'Force', true ,constPar)        
A_D = adjacency(G_D,'weighted');
fcn_franka_pigraph_dynamics_contracted_mi_generic(1, ...
    A_D, 1:constPar.noj, 1:constPar.nol, 'Force', false, constPar)            
%     fcn_franka_pigraph_dynamics_contracted_mi_generic(1, ...
%         A_D(constPar.noj+1:end,constPar.noj+1:end), 1:constPar.noj, 1:constPar.nol, 'Force', false ,constPar)        
    
    SAVE_FIG = 0;
    if SAVE_FIG == 1
        h1 =  findall(groot,'Type','figure');
        close(h1(1))
        clc            
        export_fig(fullfile(fileparts(matlab.desktop.editor.getActiveFilename),'figures',['phantomx_pigraph_' num2str(index)]),'-pdf')
        close(gcf);
    end        
elseif SIGNAL_CHOICE == 5
    % Plot relationships between joint velocity, torque and Cartesian 
    % angular velocities
    if ENHANCE == 1
        % Ignore relationships between joint angular velocities -----------
        W_aux(1:N_joints,1:N_joints) = zeros(N_joints);
        % Artificially relating sensor x,y,z components -------------------
        for i=1:N_joints+1
            range = ((3*i-2):3*i) + N_joints; 
            W_aux(range,range) = ones(3) - eye(3);
            % W_aux(range,range) = 1.1*max(W_aux,[],'all')*(ones(3) - eye(3));
        end
    end
    % Matrix contraction --------------------------------------------------
    if CONTRACT == 1
        W_ctrct = zeros(2*N_joints + N_links);
        W_ctrct(1:2*N_joints,1:2*N_joints) = 1*W_aux(1:2*N_joints,1:2*N_joints);
        % Clean inter-dq/inter-tau connections
        %  W_ctrct(1:N_joints,1:N_joints) = zeros(N_joints,N_joints);
        W_ctrct(N_joints+1:2*N_joints,N_joints+1:2*N_joints) = zeros(N_joints,N_joints);
        % Contract blocks based on vector/matrix norms
        for i=1:2*N_joints
            for j=1:N_links
                c_range = ((3*j-2):3*j) + 2*N_joints;
                W_ctrct(i,j+2*N_joints) = norm(W_aux(i,c_range),'fro');
                W_ctrct(j+2*N_joints,i) = W_ctrct(i,j+2*N_joints);
            end
        end
        for i = 1:N_links
            r_range = ((3*i-2):3*i) + 2*N_joints;
            for j=1:N_links
                if i==j
                    continue;
                else
                    c_range                            = ((3*j-2):3*j) + 2*N_joints;
                    W_ctrct(i+2*N_joints,j+2*N_joints) = norm(W_aux(r_range,c_range),'fro');
                end
            end
        end    
        W_ctrct  = (W_ctrct + transpose(W_ctrct))/2; % force matrix to be symmetric
        [T,~]  = minspantree(graph(-W_ctrct,'upper'));
        T.Edges.Weight = abs(T.Edges.Weight);
        delta  = PRUNE*(min(T.Edges.Weight));
        disp(['The prune threshold is: ', num2str(delta)])
        W_ctrct(W_ctrct<delta) = 0; % prune matrix
        W_ctrct = W_ctrct / sum(W_ctrct(:)); % normalize matrix
        wIndex = 1:7;
        qIndex = 1:7;
        fcn_franka_pigraph_contracted_mi_ang_var(1, W_ctrct, qIndex, wIndex, 'Layered', constPar)
        SAVE_FIG = 0;
        if SAVE_FIG == 1
            h1 =  findall(groot,'Type','figure');
            close(h1(1))
            clc
            export_fig(fullfile(fileparts(matlab.desktop.editor.getActiveFilename),'figures',['phantomx_pigraph_' num2str(index)]),'-pdf')
            close(gcf);
        end
    end
elseif strcmp(SIGNAL_CHOICE,'dq-tau')
    W_aux  = (W_aux + transpose(W_aux))/2; % force matrix to be symmetric
%     [T,~]  = minspantree(graph(-W_aux,'upper'));
%     T.Edges.Weight = abs(T.Edges.Weight);
%     delta  = PRUNE*(min(T.Edges.Weight));
%     disp(['The prune threshold is: ', num2str(delta)])
%     W_aux(W_aux<delta) = 0;
%     W_aux = W_aux / sum(W_aux(:));
    
    subindices = ['0ABCDEFGHIJKLMNOPQRSTUVWXYZ'];
    PLOT_FLAG = 1;
    if(PLOT_FLAG)
        figure('Color','w','Name','MI_heatmap')
            h                      = heatmap(W_aux(1:constPar.noj,constPar.noj+1:end));
            h.Colormap             = flipud(autumn);
            ax                     = gca;
            %ax.XDisplayLabels      = [arrayfun(@(index){['\tau_', num2str(index)]}, 1:constPar.noj),arrayfun(@(index){['a_', subindices(index)]}, 1:constPar.nob)];
            %ax.YDisplayLabels      = [arrayfun(@(index){['\tau_', num2str(index)]}, 1:constPar.noj),arrayfun(@(index){['a_', subindices(index)]}, 1:constPar.nob)];
            ax.XDisplayLabels      = [arrayfun(@(index){['dq_{', num2str(index), '}']}, 1:constPar.noj)];
            ax.YDisplayLabels      = [arrayfun(@(index){['\tau_{', num2str(index), '}']}, 1:constPar.noj)];

            axp                    = struct(gca);       %you will get a warning
            axp.Axes.XAxisLocation = 'top';
            title('\bf W^{\tau, dq}')
    end    
    
elseif SIGNAL_CHOICE == 7
    W_aux  = (W_aux + transpose(W_aux))/2; % force matrix to be symmetric
%     [T,~]  = minspantree(graph(-W_aux,'upper'));
%     T.Edges.Weight = abs(T.Edges.Weight);
%     delta  = PRUNE*(min(T.Edges.Weight));
%     disp(['The prune threshold is: ', num2str(delta)])
%     W_aux(W_aux<delta) = 0;
%     W_aux = W_aux / sum(W_aux(:));
    
    subindices = ['0ABCDEFGHIJKLMNOPQRSTUVWXYZ'];
    PLOT_FLAG = 1;
    if(PLOT_FLAG)
        figure('Color','w','Name','MI_heatmap')
            h                      = heatmap(W_aux(1:constPar.noj,constPar.noj+1:end));
            h.Colormap             = flipud(autumn);
            ax                     = gca;
            %ax.XDisplayLabels      = [arrayfun(@(index){['\tau_', num2str(index)]}, 1:constPar.noj),arrayfun(@(index){['a_', subindices(index)]}, 1:constPar.nob)];
            %ax.YDisplayLabels      = [arrayfun(@(index){['\tau_', num2str(index)]}, 1:constPar.noj),arrayfun(@(index){['a_', subindices(index)]}, 1:constPar.nob)];
            ax.XDisplayLabels      = [arrayfun(@(index){['q_{', num2str(index), '}']}, 1:constPar.noj)];
            ax.YDisplayLabels      = [arrayfun(@(index){['\tau_{', num2str(index), '}']}, 1:constPar.noj)];

            axp                    = struct(gca);       %you will get a warning
            axp.Axes.XAxisLocation = 'top';
            title('\bf W^{\tau, q}')
    end    
    
end
% end
