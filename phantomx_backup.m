
%% **************************************************************************
%                Computation of the MI matrix using JIDT                  *
% *************************************************************************
clc
close all
load('./data/phantomx_parameters.mat')
path_root = fullfile(pwd,'data','phantomx');
% clear  spectral_distance matrix_distance sum_squared_distances g_ssd_mavg edge_error

%% Load datasets

phantomx_proprioception.dataset1 = fcn_phantomx_get_values_from_gazebo({fullfile(path_root, strcat('pxmark4_memory_b1'))}, 0, constPar);
phantomx_proprioception.dataset2 = fcn_phantomx_get_values_from_gazebo({fullfile(path_root, strcat('pxmark4_memory_b2'))}, 0, constPar);
phantomx_proprioception.dataset3 = fcn_phantomx_get_values_from_gazebo({fullfile(path_root, strcat('pxmark4_memory_b3'))}, 0, constPar);
phantomx_proprioception.dataset4 = fcn_phantomx_get_values_from_gazebo({fullfile(path_root, strcat('pxmark4_memory_b4'))}, 0, constPar);
phantomx_proprioception.dataset5 = fcn_phantomx_get_values_from_gazebo({fullfile(path_root, strcat('pxmark4_memory_b5'))}, 0, constPar);

%% Load parameters
load('./data/phantomx_proprioception.mat')

% The rows of each of the datasets is organized as follows:
% - q_indices       = 1:18;
% - dq_indices      = 19:36;
% - torque_indices  = 37:54;
% - ang_vel_indices = 55:111;
% - lin_acc_indices = 112:168;

q_indices       = 1:18;
dq_indices      = 19:36;
torque_indices  = 37:54;
ang_vel_indices = 55:111;
lin_acc_indices = 112:168;

%% Create the proprioceptive signals matrix
clc
close all
proprioception = [];

fn = fieldnames(phantomx_proprioception);
for k=1:numel(fn)
    if( isnumeric(phantomx_proprioception.(fn{k})) )
        proprioception_aux = [phantomx_proprioception.(fn{k})([q_indices,dq_indices,torque_indices,ang_vel_indices],:); ...
                          0*phantomx_proprioception.(fn{k})(ang_vel_indices,:); ...
                          0*phantomx_proprioception.(fn{k})(ang_vel_indices,:); ...
                          phantomx_proprioception.(fn{k})(lin_acc_indices,:)];
        proprioception = [proprioception proprioception_aux];
    end
end
disp('done')

%% **************************************************************************
%                       COMPUTATION OF THE MI MATRIX                      *
% *************************************************************************
close all
clc

BATCH_SIZE  = 100;
BUFFER_SIZE = 10000;

% Initialize the replay buffer
buffer = replayBuffer(BUFFER_SIZE, size(proprioception,1));

n_mi_matrices = round(length(proprioception)/BATCH_SIZE);

total_loops                = 1;
n_datasets                 = 1;

for loop = 1:total_loops
    dataSet   = repmat(proprioception,1,n_datasets);
    % Adding noise to signals ---------------------------------------------
    WITH_NOISE = 1;
    if WITH_NOISE == 1
        SNR     = 20;
        propioceptiveSignals = zeros(size(dataSet));
        for signal = 1:size(propioceptiveSignals,1)
            propioceptiveSignals(signal,:) = awgn(dataSet(signal,:), SNR, 'measured');
        end
    else
        propioceptiveSignals = dataSet;
    end
    cprintf('*yellow', ['>> Measurements selected, WITH_NOISE: ' ,num2str(WITH_NOISE) , ' \n'])
    cprintf('*yellow', '>> Computing MI matrix using JIDT...\n')
    pause(1)
    % Computation of the pairwise MI matrix -----------------------------------
    % The following function computes the MI matrix incrementally by
    % maintaining a buffer of pseudo-uniformly distributed points and sampling 
    % at random a mini batch from it. 

    N_points                = size(propioceptiveSignals,2);
    N_signals               = size(propioceptiveSignals,1);
    
    num_batches             = floor(N_points/BATCH_SIZE);
    % Full MI matrix
    MI_batch                = zeros(N_signals,N_signals,num_batches);
    N_nodes                 = 3*constPar.noj + 2*constPar.nob;
    MI_cmean                = zeros(N_nodes,N_nodes,num_batches);
    edge_error              = NaN(1,num_batches);
    totalMutalInformation     = NaN(1,num_batches);
    diffTotalMutalInformation = NaN(1,num_batches);
        
    window                  = 100; % sliding window samples

    batch = 0;
    for sample = 1:length(propioceptiveSignals)
        % Store sample point in buffer
        buffer.storeSample(propioceptiveSignals(:,sample)); 
        if mod(sample,BATCH_SIZE) == 0 && buffer.LogIndex >= BATCH_SIZE
            batch = batch + 1;
            disp(sample)
            fprintf('BATCH: %d/%d | BUFFER_STATUS: %d %% | SAMPLE: %d/%d\n', ...
                batch, ...
                num_batches, ...
                buffer.Capacity, ...
                sample, ...
                N_points);
            MI_bacth_contracted = ...
                   fcn_compute_minibatch_mi_matrix(buffer.drawSamples(BATCH_SIZE), 'mex', constPar);
            if batch == 1
                MI_cmean(:,:,batch) = MI_bacth_contracted;
            else
                % Cumulative mean
                MI_cmean(:,:,batch) = (batch-1)/batch*MI_cmean(:,:,batch-1) + MI_bacth_contracted/batch;
            end
            % Get the total information in the matrix
            D_cmean                        = diag(sum(MI_cmean(:,:,batch),2));
            totalMutalInformation(batch) = 0.5*trace(D_cmean);

            if batch == 1
                diffTotalMutalInformation(batch) = 1E6;
            else
                diffTotalMutalInformation(batch) = totalMutalInformation(batch) - totalMutalInformation(batch-1);
            end
            disp(['TMI: ',num2str(totalMutalInformation(batch)),' | Grad_MTI: ',num2str(diffTotalMutalInformation(batch)),' | edge error: ', num2str(edge_error(batch))])      
        
            % Stopping criteria
            kin_indices = [constPar.noj + 1:2*constPar.noj, ...
                           3*constPar.noj + 1:3*constPar.noj + constPar.nob];
            MI_kin                 = MI_cmean(kin_indices,kin_indices,batch);
            [G_kin_mst,~]          = minspantree(graph(-MI_kin,'upper'));
            G_kin_mst.Edges.Weight = abs(G_kin_mst.Edges.Weight);
            Dist_mat = distances(G_kin_mst,'Method','unweighted');
            if all(min(Dist_mat(1:constPar.noj,1:constPar.noj)+10*eye(constPar.noj),[],2)>=3) && batch>window
                if mean(abs(diffTotalMutalInformation(batch-(window-1):batch))) < 5E-3
                %if all(abs(diffTotalMutalInformation(batch-(window-1):batch)) < 5E-3)
                %if edge_error(batch) == 0
                    break;
                end
            end
        end
    end
end
cprintf('*green', '>> MI matrix ready!\n')
figure('color','w')
semilogy(abs(diffTotalMutalInformation(1:batch)))
hold on;
area(1:batch,5E-3*ones(size(diffTotalMutalInformation(1:batch))),'FaceColor','r','FaceAlpha',0.1,'EdgeColor','none');
plot(5E-3*ones(size(diffTotalMutalInformation(1:batch))),'r--')
ylabel('|\Delta TMI|')
xlabel('Batches')
xlim([1 batch])

%% Plot proprioceptive kinematics graph
clc

constPar.showClusters         = 1;
constPar.plotting.interpreter = 'latex'; 
constPar.useNodeLabels        = 1;

GRAPH_CHOICE = 'dq-omg';
NORMALIZE    = 0;
PRUNE        = 1; 

W_MI         = squeeze(MI_cmean(:,:,batch));    
mi_graphs    = fcn_plot_robot_mi_matrix_graphs(exp(W_MI) - 1*0.9999, ...
                                 GRAPH_CHOICE, NORMALIZE, PRUNE, constPar);

%% **************************************************************************
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
