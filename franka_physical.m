%% **************************************************************************
%         MORPHOLOGY LEARNING OF THE FRANKA EMIKA ROBOT MANIPULATOR
% *************************************************************************

clearvars
clc
close all

% Load configuration parameters
load('./data/franka_parameters.mat')

% Load the dataset
load('./data/franka_proprioception_physical_fixed_base.mat')

q_indices       = 1:7;
dq_indices      = 8:14;
ddq_indices     = 15:21;
torque_indices  = 22:28;
ang_vel_indices = 29:52;
lin_vel_indices = 53:76;
ang_acc_indices = 77:100;
lin_acc_indices = 101:124;   

proprioception = signals([q_indices, ...
                          dq_indices, ...
                          torque_indices, ...
                          ang_vel_indices, ...
                          lin_vel_indices, ...
                          ang_acc_indices, ...
                          lin_acc_indices],:);
%% **************************************************************************
%                       COMPUTATION OF THE MI MATRIX                      *
% *************************************************************************
close all
clc

BATCH_SIZE  = 100;
BUFFER_SIZE = 10000;

% Initialize the replay buffer
buffer = replayBuffer(BUFFER_SIZE, size(proprioception,1));

n_mi_matrices = round(length(signals)/BATCH_SIZE);

total_loops                = 1;
n_datasets                 = 5;

for loop = 1:total_loops
%     dataSet   = repmat(proprioception,1,n_datasets);
%     % Adding noise to signals ---------------------------------------------
%     WITH_NOISE = 1;
%     if WITH_NOISE == 1
%         SNR     = 20;
%         proprioceptiveSignals = zeros(size(dataSet));
%         for signal = 1:size(proprioceptiveSignals,1)
%             proprioceptiveSignals(signal,:) = awgn(dataSet(signal,:), SNR, 'measured');
%         end
%         %return
%         %propioceptiveSignals    = dataSet + WITH_NOISE*normrnd(0,0.1/3,size(dataSet,1),size(dataSet,2));
%     else
%         proprioceptiveSignals = dataSet;
%     end
%     cprintf('*yellow', ['>> Measurements selected, WITH_NOISE: ' ,num2str(WITH_NOISE) , ' \n'])
    
    proprioceptiveSignals   = repmat(proprioception,1,n_datasets);
    cprintf('*yellow', '>> Computing MI matrix using JIDT...\n')
    pause(1)
    % Computation of the pairwise MI matrix -----------------------------------
    % The following function computes the MI matrix incrementally by
    % maintaining a buffer of pseudo-uniformly distributed points and sampling 
    % at random a mini batch from it. 

    N_points                = size(proprioceptiveSignals,2);
    N_signals               = size(proprioceptiveSignals,1);
    
    num_batches             = floor(N_points/BATCH_SIZE);
    % Full MI matrix
    N_nodes                 = 3*constPar.noj + 2*constPar.nob;
    MI_batch_contracted     = zeros(N_nodes,N_nodes,num_batches);
    MI_cmean                = zeros(N_nodes,N_nodes,num_batches);
    edge_error              = NaN(1,num_batches);
    totalMutalInformation     = NaN(1,num_batches);
    diffTotalMutalInformation = NaN(1,num_batches);
        
    window                  = 100; % sliding window samples

    batch = 0;
    for sample = 1:length(proprioceptiveSignals)
        % Store sample point in buffer
        buffer.storeSample(proprioceptiveSignals(:,sample)); 
        if mod(sample,BATCH_SIZE) == 0 && buffer.LogIndex >= BATCH_SIZE
            batch = batch + 1;
            disp(sample)
            fprintf('BATCH: %d/%d | BUFFER_STATUS: %d %% | SAMPLE: %d/%d\n', ...
                batch, ...
                num_batches, ...
                buffer.Capacity, ...
                sample, ...
                N_points);
            MI_batch_contracted(:,:,batch) = ...
                   fcn_compute_minibatch_mi_matrix(buffer.drawSamples(BATCH_SIZE), 'mex', constPar);
            if batch == 1
                MI_cmean(:,:,batch) = MI_batch_contracted(:,:,batch);
            else
                % Cumulative mean
                MI_cmean(:,:,batch) = (batch-1)/batch*MI_cmean(:,:,batch-1) + MI_batch_contracted(:,:,batch)/batch;
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
 
%% Display the MI matrix
W_MI                          = squeeze(MI_cmean(:,:,batch));    
constPar.showClusters         = 1;
constPar.plotting.interpreter = 'latex'; 
constPar.useNodeLabels        = 1;
clc

% Plot proprioceptive kinematics graph
GRAPH_CHOICE = 'dq-omg';
NORMALIZE    = 0;
PRUNE        = 1; 
mi_graphs    = fcn_plot_robot_mi_matrix_graphs(exp(W_MI) - 1*0.9999, ...
                                 GRAPH_CHOICE, NORMALIZE, PRUNE, constPar);
% Extract the kinematics subgraph
A_kin_mst    = mi_graphs.G_kin_mst.adjacency;
A_kin_pi     = A_kin_mst(constPar.noj+1:end,constPar.noj+1:end);
% * NOTE: the adjacency matrix is not totally correct since in the fixed 
%         based case the base IMU does not provide information
% Artificial correction
A_kin_pi(1,4) = 0;
A_kin_pi(4,1) = 0;
A_kin_pi(1,2) = 1;
A_kin_pi(2,1) = 1;

%% ************************************************************************
%               ONLINE LEARNING (RAMS GRADIENT DESCENT)                   *
% *************************************************************************

clc
% buffer.flush();
IS_ACC     = true;
max_epochs = 5000;
total_runs = 5;
J_log      = zeros(max_epochs,total_runs);

% Penalty factor for the joint center point vectors
constPar.beta = 0.1;
samplingTime  = 1E-3;
% rng('default')
disp('Computing...')
for loop = 1:total_runs
    xi_0         = ones(3,1);
    xi_0         = xi_0/norm(xi_0);
    phi_0        = wrapToPi(deg2rad(randi(360)));
    zeta_0       = ones(3,1);
    zeta_0       = zeta_0/norm(zeta_0);
    gamma_hat_0  = [xi_0;phi_0;zeta_0];
    rho_hat_0    = 1E-3*ones(6,1);
    lambda_hat_0 = [gamma_hat_0;
                    rho_hat_0];

    [lambda_hat_online_aux(:,:,loop), ~, J_log(:,loop)] = ...
        fcn_robot_kinematics_rams_online_learning(proprioception, ...
                                                  lambda_hat_0, ...
                                                  IS_ACC, ...
                                                  A_kin_pi, ...
                                                  max_epochs, ...
                                                  buffer, ...
                                                  BATCH_SIZE, ...
                                                  samplingTime, ...
                                                  constPar);
end

[J_min, J_index] = min(J_log(end,:));
lambda_hat_online = squeeze(lambda_hat_online_aux(:,:,J_index));
%% Display the learned kinematic description
clc
close all
[parents, children]                = find(triu(A_kin_pi) == 1);
constPar.displayJointConfiguration =  [0 0 0 0 0 deg2rad(120) 0]';
h = figure('color','w','units','normalized','outerposition',[0 0 1 1]);
fcn_franka_display_learned_morphology(lambda_hat_online, ...
    constPar.p_R_c, ...
    constPar.panda, ...
    constPar, ...
    mi_graphs, ...
    parents, ...
    children, 0, h);
