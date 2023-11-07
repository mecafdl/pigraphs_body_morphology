%% **********************************************************************
%         MORPHOLOGY LEARNING OF THE FRANKA EMIKA ROBOT MANIPULATOR       *
% *************************************************************************

clearvars
clc
close all

% Change to the directory hosting the current file
cd(fileparts(matlab.desktop.editor.getActiveFilename));

% Addpahts
addpath(genpath('./data'))
addpath(genpath('./functions'))
addpath(genpath('./external_toolboxes'))

% Load configuration parameters
load('./parameters/franka_parameters.mat')

% Load the dataset
load('./data/armIMUExperiment.mat')
armIMUExperiment = fcn_get_processed_imu_experiment_data(armIMUExperiment);
signals = [armIMUExperiment.jointPosition.raw; ...
           armIMUExperiment.jointVelocity.raw; ...
           armIMUExperiment.jointAcceleration.numerical; ...
           armIMUExperiment.jointTorque.raw; ...
           armIMUExperiment.bodyAngularVelocity.zeroPhaseFilter; ...
           zeros(size(armIMUExperiment.bodyLinearAcceleration.raw)); ...
           armIMUExperiment.bodyAngularAcceleration.numerical; ...
           armIMUExperiment.bodyLinearAcceleration.zeroPhaseFilter ...
           ];

cprintf('*yellow', '>> Real recorded IMU measurements loaded!\n')

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

% Artificially added noise
WITH_NOISE     = 1;
proprioception = proprioception + WITH_NOISE*normrnd(0,0.1/3,size(proprioception,1),size(proprioception,2));

%% **********************************************************************
%                      COMPUTATION OF THE MI MATRIX                       *    
% *************************************************************************
close all
clc

BATCH_SIZE  = 1E2;
BUFFER_SIZE = 1E4;

% Initialize the replay buffer
buffer = replayBuffer(BUFFER_SIZE, size(proprioception,1));

n_mi_matrices = round(length(signals)/BATCH_SIZE);
total_loops   = 1;
n_datasets    = 5;

for loop = 1:total_loops   
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
warning('the adjacency matrix is not totally correct since, in the fixed based case, the base IMU does not provide information.')


% Artificial correction
A_kin_pi(1,5) = 0;
A_kin_pi(5,1) = 0;
A_kin_pi(1,2) = 1;
A_kin_pi(2,1) = 1;

% Extract body-to-body relationships from graph
[parents, children] = find(triu(A_kin_pi) == 1);

%% **********************************************************************
%            ROTATION MATRICES/ROTATION AXES  BETWEEN SENSORS             *
% *************************************************************************
clc

% rng('default')
test_points = randperm(size(armIMUExperiment.jointPosition.raw,2),3000);% choose k unique data points from data set

% Initial point -----------------------------------------------------------
xi_0   = rand(3,1);
xi_0   = xi_0/norm(xi_0);
phi_0  = deg2rad(randi(360));
zeta_0 = rand(3,1);
zeta_0 = zeta_0/norm(zeta_0);
x0     = [xi_0;phi_0;zeta_0];

% Optimizer settings ------------------------------------------------------
constPar.MaxIter     = 100;
constPar.MaxFunEvals = 1000;
options = optimset('Display','iter-detailed','MaxIter',constPar.MaxIter,'MaxFunEvals', 1000,...
    'TolFun',1e-6,'UseParallel',false);

zeta_hat = NaN(7,constPar.noj);
for j = 1:constPar.noj
    p = parents(j);
    c = children(j);
    cprintf('*yellow', ['>> Checking: ' num2str(p-1) '-->' num2str(c-1) '\n'])
    [zeta_hat(:,j), ~] = fmincon(...
    @(x) fcn_get_rotation_err(armIMUExperiment.bodyAngularVelocity.raw(3*(p-1)+1:3*p,test_points), ...
                              armIMUExperiment.bodyAngularVelocity.raw(3*(c-1)+1:3*c,test_points), ...
                              armIMUExperiment.jointPosition.raw(c-1,test_points), ...
                              armIMUExperiment.jointVelocity.raw(c-1,test_points), x),...
                              x0,[],[],[],[],[],[],...
                                  @(x)fcn_axis_nonlConstraints(x),options);
    disp('--------------------------------------------------------------------')
end
% Sensor-to-Sensor rotation matrix
p_R_c = @(q_c, x) expm(skew([x(1);x(2);x(3)])*(x(4)))*expm(skew([x(5);x(6);x(7)]).*q_c);
cprintf('*green', '>> Rot mat between sensors and joint axes found!!!\n')

%% **********************************************************************
%                      FIND JOINT POSITION VECTORS                        *
% *************************************************************************
clc
close all

test_points = randperm(size(armIMUExperiment.jointPosition.raw,2),3000);% choose k unique data points from data set

% Zerp-phase filter (zpf)
Fs = 1E3; % signal sampling frequency [Hz]
Fc =  2;  % cutoff frequency [Hz], obtained from the previous block as the average of the frequencies for all joints where the amplitude becomes negligible
d  = designfilt('lowpassiir','FilterOrder',6, ...
    'HalfPowerFrequency',(2/Fs)*Fc,'DesignMethod','butter'); 
%                       |_________|
%                            |_______________ This is the Nyquist normalized frequency

[b_zpf, a_zpf] = tf(d);

% Pre allocate
r_hat_sensor = zeros(6,constPar.noj);
r_hat_model  = NaN(6,constPar.noj);

% Initial point
rng('default')
x0 = 1E-3*ones(6,1);

for j = 1:constPar.noj
    p = parents(j);
    c = children(j);
    cprintf('*yellow', ['>> Checking: ' num2str(p-1) '-->' num2str(c-1) ' using ', 'ACC','\n'])

    [r_hat_sensor(:,j), ~] = fmincon(...
    @(x) fcn_get_cartesian_acceleration_error_extended(...
             armIMUExperiment.bodyLinearAcceleration.zeroPhaseFilter(3*(p-1)+1:3*p,test_points), ...
             armIMUExperiment.bodyLinearAcceleration.zeroPhaseFilter(3*(c-1)+1:3*c,test_points), ...
             armIMUExperiment.bodyAngularVelocity.zeroPhaseFilter(3*(p-1)+1:3*p,test_points), ...
             armIMUExperiment.bodyAngularVelocity.zeroPhaseFilter(3*(c-1)+1:3*c,test_points), ...
             armIMUExperiment.bodyAngularAcceleration.numerical(3*(p-1)+1:3*p,test_points), ...
             armIMUExperiment.bodyAngularAcceleration.numerical(3*(c-1)+1:3*c,test_points), ...
             armIMUExperiment.jointPosition.raw(j,test_points), ...
             p_R_c, zeta_hat(:,j), ...
             0.01, x), x0,[],[],[],[],[],[], [],options);      
end

cprintf('*Green', '>> Joint center point relative to sensor CS found!!!\n')

%% DISPLAY RESULTS ======================================================

lambda_hat_offline_fix = [zeta_hat;
                          r_hat_sensor];

clc
close all
h = figure('color','w');
fcn_franka_display_learned_morphology(lambda_hat_offline_fix, p_R_c, constPar.panda, constPar, mi_graphs, parents, children, 0, h)
