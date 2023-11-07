%% ************************************************************************
% HUMANOID MORPHOLOGY LEARNING (PROTO KINEMATICS) BASED ON                *
% IMU MEASUREMENTS                                                        *
% *************************************************************************
% This file assumes the mechanical arrangement of the robot to be known and
% uses measurements from IMUs (one per link) plus eqs. from differential
% kinematics to learn the kinematics of the robot.

%% **********************************************************************
%            MORPHOLOGY LEARNING OF THE POPPY HUMANOID ROBOT              *
% *************************************************************************
clearvars
clc
close all

% Change to the directory hosting the current file
cd(fileparts(matlab.desktop.editor.getActiveFilename));

% Load configuration parameters
load('./parameters/poppy_parameters.mat','constPar')
path_root = fullfile(pwd,'data','phantomx');

load('./data/humanoidSignals.mat','humanoidSignals')
% [hexapodSignals, proprioceptiveSignals]  = fcn_loadHexapodRobotData(constPar);

%% **********************************************************************
%                       COMPUTATION OF THE MI MATRIX                      *
% *************************************************************************
close all
clc

proprioceptiveSignals = [humanoidSignals.jointPosition.raw; ...
                         humanoidSignals.jointVelocity.raw; ...
                         humanoidSignals.jointTorque.raw; ...
                         humanoidSignals.bodyAngularVelocity.raw; ...
                         humanoidSignals.bodyLinearVelocity.raw;
                         humanoidSignals.bodyAngularAcceleration.raw; ...
                         humanoidSignals.bodyLinearAcceleration.raw];

BATCH_SIZE  = 100;
BUFFER_SIZE = 1000;

% Initialize the replay buffer
buffer = replayBuffer(BUFFER_SIZE, size(proprioceptiveSignals,1));

n_mi_matrices = round(length(proprioceptiveSignals)/BATCH_SIZE);

total_loops                = 50;
n_datasets                 = 1;

for loop = 1:total_loops
    cprintf('*yellow', '>> Computing MI matrix using JIDT...\n')
    pause(1)
    % Computation of the pairwise MI matrix -------------------------------
    % The following function computes the MI matrix incrementally by
    % maintaining a buffer of pseudo-uniformly distributed points and sampling 
    % at random a mini batch from it. 

    N_points                = size(proprioceptiveSignals,2);
    N_signals               = size(proprioceptiveSignals,1);
    
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

%% Plot proprioceptive kinematics graph ---------------------------------
clc

constPar.showClusters         = 1;
constPar.plotting.interpreter = 'latex'; 
constPar.useNodeLabels        = 1;

GRAPH_CHOICE = 'dq-omg';
NORMALIZE    = 0;
PRUNE        = 1; 

W_MI         = squeeze(MI_cmean(:,:,batch));    
mi_graphs    = fcn_plot_robot_mi_matrix_graphs(exp(W_MI) - 0.9999, ...
                                 GRAPH_CHOICE, NORMALIZE, PRUNE, constPar);

A_kin               = full(adjacency(mi_graphs.G_kin_mst));
A_dq_omg            = A_kin(1:constPar.noj,constPar.noj+1:end);
A_omg               = A_kin(constPar.noj+1:end,constPar.noj+1:end);    
[parents, children] = find(triu(A_omg) == 1);

%%
A_kin_ref = A_kin;
A_kin_ref(1:25,27:end) = eye(25);
A_kin_ref(27:end,1:25) = eye(25);
figure('Color','w');
imagesc(A_kin_ref);
axis square

%% **********************************************************************
%               ONLINE LEARNING (RAMS GRADIENT DESCENT)                   *
% *************************************************************************

[lambda_hat_online] = fcn_poppy_morphology_online_learning_rams(proprioceptiveSignals, parents, children, constPar);

%% DISPLAY LEARNED KINEMATICS ===========================================
% 
clc

% Offline
constPar.displayConfiguration = constPar.poppy.homeConfiguration;

h = figure('color','w');
fcn_poppy_display_learned_morphology(lambda_hat_online, constPar.p_R_c, constPar.poppy, constPar, mi_graphs, parents, children, h)
