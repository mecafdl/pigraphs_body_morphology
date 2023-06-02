
%% **************************************************************************
%                Computation of the MI matrix using JIDT                  *
% *************************************************************************
clc
close all
% path_root    = fullfile('/','media','diaz','extreme_ssd','science_paper_material','final_results','phantomx');
load('./phantomx_parameters.mat')
path_root    = fullfile(pwd,'data','phantomx');
clear  spectral_distance matrix_distance sum_squared_distances g_ssd_mavg edge_error


%%

phantomx_proprioception.dataset1 = fcn_phantomx_get_values_from_gazebo({fullfile(path_root, strcat('pxmark4_memory_b1'))}, 0, constPar);
phantomx_proprioception.dataset2 = fcn_phantomx_get_values_from_gazebo({fullfile(path_root, strcat('pxmark4_memory_b2'))}, 0, constPar);
phantomx_proprioception.dataset3 = fcn_phantomx_get_values_from_gazebo({fullfile(path_root, strcat('pxmark4_memory_b3'))}, 0, constPar);
phantomx_proprioception.dataset4 = fcn_phantomx_get_values_from_gazebo({fullfile(path_root, strcat('pxmark4_memory_b4'))}, 0, constPar);
phantomx_proprioception.dataset5 = fcn_phantomx_get_values_from_gazebo({fullfile(path_root, strcat('pxmark4_memory_b5'))}, 0, constPar);

%%
buffer_index = perms(1:5);
buffer_index = buffer_index(randperm(size(buffer_index,1)),:);

total_loops = 1;
n_datasets  = 1;
spectral_distance_mat      = zeros(total_loops,n_datasets*1000-1);
matrix_distance_mat        = zeros(total_loops,n_datasets*1000-1);
sum_squared_distances_mat  = zeros(total_loops,n_datasets*1000-1);
g_ssd_mavg_mat             = zeros(total_loops,n_datasets*1000-1);
edge_error_mat             = zeros(total_loops,n_datasets*1000);
totalMutInf_mat            = zeros(total_loops,n_datasets*1000);

BATCH_SIZE  = 100;
BUFFER_SIZE = 10000;
CONTRACT    = 1;
METHOD      = 'standard';

figure
for loop = 1:total_loops
    memory_path  = arrayfun(@(i) ...
        fullfile(path_root, strcat('pxmark4_memory_b', num2str(buffer_index(loop,i)))), 1:n_datasets, ...
        'UniformOutput', false);

    dataGazebo   = fcn_phantomx_get_values_from_gazebo(memory_path, 0, constPar);
    % Adding noise to signals ---------------------------------------------
    WITH_NOISE = 1;
    if WITH_NOISE == 1
        SNR     = 20;
        msrmnts = zeros(size(dataGazebo));
        for signal = 1:size(msrmnts,1)
            msrmnts(signal,:) = awgn(dataGazebo(signal,:), SNR, 'measured');
        end
    else
        msrmnts = dataGazebo;
    end
    cprintf('*yellow', ['>> Measurements selected, WITH_NOISE: ' ,num2str(WITH_NOISE) , ' \n'])
    
    disp('>> Computing MI matrix using JIDT...\n')
    % Computation of the pairwise MI matrix -----------------------------------
    % The following function computes the MI matrix incrementally by
    % maintaining a buffer of pseudo-uniformly distributed points and sampling 
    % at random a mini batch from it. 
    [MI_cmean, ...
     ~, ...
     ~, ...
     ~, ...
     ~, ...
     ~, ...
     ~, ...
     totalMutInf, ...
     edge_error] = ...
                   fcn_compute_mi_matrix_for_streaming_signals(msrmnts, ...
                                                               BATCH_SIZE, ...
                                                               BUFFER_SIZE,...
                                                               CONTRACT, ...
                                                               METHOD, constPar);
    edge_error_mat(loop,:)     = edge_error;
    totalMutInf_mat(loop,:)    = totalMutInf;
    disp('>> MI matrix ready!\n')
    plot(1:n_datasets*1000,edge_error_mat(loop,:))
    hold on
end
%% **********************************************************************
%                COMPUTATION OF THE MI MATRIX 
% *************************************************************************
clc
close all

load('./data/franka_proprioception.mat')
% This matrix contains
% - joint position
% - joint velocity
% - joint torque
% - link angular velocty
% - link linear acceleration


load('franka_parameters.mat')

%%
close all
clc

buffer_index = perms(1:5);
buffer_index = buffer_index(randperm(size(buffer_index,1)),:);

BATCH_SIZE  = 100;
BUFFER_SIZE = 10000;
CONTRACT    = 1;
METHOD      = 'standard';

total_loops = 5;
n_datasets  = 5;
spectral_distance_mat      = zeros(total_loops,n_datasets*1200-1);
matrix_distance_mat        = zeros(total_loops,n_datasets*1200-1);
sum_squared_distances_mat  = zeros(total_loops,n_datasets*1200-1);
g_ssd_mavg_mat             = zeros(total_loops,n_datasets*1200-1);
edge_error_mat             = zeros(total_loops,n_datasets*1200);
totalMutInf_mat            = zeros(total_loops,n_datasets*1200);

figure
for loop = 1:total_loops
%     memory_path  = arrayfun(@(i) ...
%         fullfile(path_root, strcat('pxmark4_memory_b', num2str(buffer_index(loop,i)))), 1:n_datasets, ...
%         'UniformOutput', false);
    dataSet   = repmat(panda_proprioception,1,n_datasets);
    %dataGazebo = dataGazebo(:,randperm(size(dataGazebo,2)));
    % Adding noise to signals ---------------------------------------------
    WITH_NOISE = 1;
    if WITH_NOISE == 1
        SNR     = 20;
        msrmnts = zeros(size(dataSet));
        for signal = 1:size(msrmnts,1)
            msrmnts(signal,:) = awgn(dataSet(signal,:), SNR, 'measured');
        end
    else
        msrmnts = dataSet;
    end
    disp(['>> Measurements selected, WITH_NOISE: ' ,num2str(WITH_NOISE) , ' \n'])
    
    disp('>> Computing MI matrix using JIDT...\n')
    % Computation of the pairwise MI matrix -----------------------------------
    % The following function computes the MI matrix incrementally by
    % maintaining a buffer of pseudo-uniformly distributed points and sampling 
    % at random a mini batch from it. 
    [MI_cmean, ...
     ~, ...
     ~, ...
     ~, ...
     ~, ...
     ~, ...
     ~, ...
     totalMutInf, ...
     edge_error] = ...
                   fcn_compute_mi_matrix_for_streaming_signals(msrmnts, ...
                                                               BATCH_SIZE, ...
                                                               BUFFER_SIZE,...
                                                               CONTRACT, ...
                                                               METHOD, constPar);
    edge_error_mat(loop,:)     = edge_error;
    totalMutInf_mat(loop,:)    = totalMutInf;
    disp('>> MI matrix ready!\n')
    plot(1:n_datasets*1200,edge_error_mat(loop,:))
    hold on
end
%% Plot proprioceptive kinematics graph
clc
close all
constPar.showClusters         = 1;
constPar.useNodeLabels        = 1;
constPar.plotting.interpreter = 'latex';
NORMALIZE    = 0;
PRUNE        = 1;
ENHANCE      = 1;
mi_graphs    = fcn_plot_robot_mi_matrix_kinematics_subgraphs((MI_cmean(:,:,end)) - 0*0.9999, ...
                                 NORMALIZE, PRUNE, ENHANCE, constPar);
A_kin_mst    = full(mi_graphs.G_kin_mst.adjacency);

%% ************************************************************************
%               ONLINE LEARNING (RAMS GRADIENT DESCENT)                   *
% *************************************************************************

%run scrpt_phantomx_morphology_online_learning_rams
clc
% signals_RAMS = signals([q_indices, dq_indices, torque_indices, ang_vel_indices, lin_vel_indices, ang_acc_indices, lin_acc_indices],:);

signals_RAMS = [franka_proprioception(1:3*constPar.noj,:); ...
                franka_proprioception([3*constPar.noj+1:3*constPar.noj + 3*constPar.nob],:); ...
                0*franka_proprioception([3*constPar.noj+1:3*constPar.noj + 3*constPar.nob],:); ...
                0*franka_proprioception([3*constPar.noj+1:3*constPar.noj + 3*constPar.nob],:); ...
                franka_proprioception([3*constPar.noj+3*constPar.nob+1:3*constPar.noj+6*constPar.nob],:)];
A_kin_pi     = mi_graphs.G_kin_mst.adjacency;
A_kin_pi     = A_kin_pi(constPar.noj+1:end,constPar.noj+1:end);
[parents, children] = find(triu(A_kin_pi) == 1);
MOVING_BASE =1;
if MOVING_BASE == 0
    [lambda_hat_online_fix, lambda_rams_mav_fix] = fcn_robot_morphology_rams_online_learning(signals_RAMS, true, A_kin_pi, 5000, constPar);
else
    [lambda_hat_online_float, lambda_rams_mav_float] = fcn_robot_morphology_rams_online_learning(signals_RAMS, true, A_kin_pi, 120000, constPar);
end
%% Online
clc
% constPar.displayJointConfiguration =  [0 0 0 0 0 deg2rad(120) 0]';
constPar.displayJointConfiguration =  q_panda(:,randi(numel(time)));
h = figure('color','w');
% fcn_poppy_display_learned_morphology(mean(lambda_rams(:,:,end-1000:end),3), p_R_c, poppy, constPar, mi_graphs, parents, children, h)
fcn_franka_display_learned_morphology(lambda_hat_online_float, p_R_c, panda, constPar, mi_graphs, parents, children, MOVING_BASE, h)
