%% **************************************************************************
%             MORPHOLOGY LEARNING OF THE PHANTOMX HEXAPOD ROBOT
% *************************************************************************

clearvars
clc
close all

% Load configuration parameters
load('./data/phantomx_parameters.mat')
path_root = fullfile(pwd,'data','phantomx');

% Load datasets
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
for k=1%:numel(fn)
    if( isnumeric(phantomx_proprioception.(fn{k})) )
        proprioception_aux = ...
                         [phantomx_proprioception.(fn{k})([q_indices,dq_indices,torque_indices,ang_vel_indices],:); ...
                          NaN*phantomx_proprioception.(fn{k})(ang_vel_indices,:); ...
                          NaN*phantomx_proprioception.(fn{k})(lin_acc_indices,:); ...
                          phantomx_proprioception.(fn{k})(lin_acc_indices,:)];
        proprioception = [proprioception proprioception_aux];
    end
end
lin_acc_indices = lin_acc_indices + 2*3*constPar.nol;
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

W_MI         = readNPY([getuserdir, '/catkin_ws/src/interbotix_xshexapod_gazebo/src/memory/mi_cntr_bar.npy']);
% W_MI         = squeeze(MI_cmean(:,:,batch));    
mi_graphs    = fcn_plot_robot_mi_matrix_graphs(exp(W_MI) - 1*0.9999, ...
                                 GRAPH_CHOICE, NORMALIZE, PRUNE, constPar);

A_kin    = full(adjacency(mi_graphs.G_kin_mst));
A_dq_omg = A_kin(1:constPar.noj,constPar.noj+1:end);
A_omg    = A_kin(constPar.noj+1:end,constPar.noj+1:end);    
%% Filtering of the signals ------------------------------------------------

% Zero-phase filter (zpf)
Fs = 100; % signal sampling frequency [Hz]
Fc =  10; % cutoff frequency [Hz], obtained from the previous block as the average of the frequencies for all joints where the amplitude becomes negligible
d  = designfilt('lowpassiir','FilterOrder',6, ...
    'HalfPowerFrequency',(2/Fs)*Fc,'DesignMethod','butter'); 
%                       |_________|
%                            |_______________ This is the Nyquist normalized frequency
[b_zpf, a_zpf] = tf(d);

% Estimated Cartesian angular acceleration
sampling_time       = 1E-2;
q_signals           = propioceptiveSignals(q_indices,:);
q_signals_zpf       = transpose(filtfilt(b_zpf, a_zpf, transpose(q_signals)));
dq_signals          = propioceptiveSignals(dq_indices,:);
dq_signals_zpf      = transpose(filtfilt(b_zpf, a_zpf, transpose(dq_signals)));
omg_signals         = propioceptiveSignals(ang_vel_indices,:);
omg_signals_zpf     = transpose(filtfilt(b_zpf, a_zpf, transpose(omg_signals)));
domg_signals        = gradient(omg_signals_zpf, sampling_time);
domg_signals_zpf    = transpose(filtfilt(b_zpf, a_zpf, transpose(domg_signals)));
ang_acc_imu_num_zpf = domg_signals_zpf;
acc_signals         = propioceptiveSignals(lin_acc_indices,:);
acc_signals_zpf     = transpose(filtfilt(b_zpf,a_zpf, transpose(acc_signals)));

cprintf('*green', '>> Filtered signals and numerical gradients ready!\n')

%% ##########################################################################
%
% 
% 
%            OFFLINE LEARNING (INTERIOR POINT OPTIMIZATION)                
% 
% 
% 
% #########################################################################


%% ************************************************************************
%            ROTATION MATRICES/ROTATION AXES BETWEEN SENSORS              *
% *************************************************************************
clc
close all
% rng('default')


n_iter     = 1;
gamma_iter = NaN(constPar.nol,constPar.nol,7,n_iter);
J_rot_iter = zeros(constPar.nol,constPar.nol,n_iter);

for iter = 1:n_iter

test_points = sort(randperm(size(propioceptiveSignals,2),1000));% choose k unique data points from data set

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

gamma_hat = NaN(constPar.nol,constPar.nol,7);
J_rot     = NaN(constPar.nol);


for j = 1:constPar.nol
    p = j;
    for k = 1:constPar.nol
        c = k;
        if A_omg(p,c) == 0 || p==c
            continue;
        else
            cprintf('*yellow', ['>> Checking: ' constPar.pxmark4.Bodies{constPar.imuBodyIndices(p)}.Name, ' -> ', constPar.pxmark4.Bodies{constPar.imuBodyIndices(c)}.Name '\n'])

            omg_parent = omg_signals_zpf(3*(p-1)+1:3*p,test_points);
            omg_child  = omg_signals_zpf(3*(c-1)+1:3*c,test_points);
            
            if any(A_dq_omg(:,c) == 1)
                q_child  = q_signals_zpf(find(A_dq_omg(:,c)),test_points);
                dq_child = dq_signals_zpf(find(A_dq_omg(:,c)),test_points);
            else
                q_child  = zeros(1,numel(test_points));
                dq_child = zeros(1,numel(test_points));
            end
             
            
            [theta_tmp, J_rot(p,c)] = fmincon(...
                        @(x) fcn_get_rotation_err(omg_parent, ...
                                                  omg_child, ...
                                                   q_child, ...
                                                  dq_child, ...
                                                  x),...
                        x0,[],[],[],[],[],[],...
                        @(x)fcn_axis_nonlConstraints(x),options);
        end
        
        
        gamma_hat(p,c,:) = theta_tmp;                   
%         if J(p,c) < 1E-4
%             cprintf('*green', ['>> Conection found!\n'])
%             zeta_hat(p,c,:) = zeta_tmp;
%         end
        disp('--------------------------------------------------------------------')
    end
end
gamma_iter(:,:,:,iter) = gamma_hat;
gamma_hat_off          = gamma_hat;
figure;heatmap(J_rot)

J_rot_iter(:,:,iter) = J_rot;

end
% Get directed adjacency matrix -------------------------------------------

J_up  = triu(J_rot);
J_low = tril(J_rot);
max_up  = max(J_up,[],'all');
min_up  = min(J_up(J_up>0),[],'all');
max_low = max(J_low,[],'all');
min_low = min(J_low(J_low>0),[],'all');

if(max_up < min_low)
    A_dir = double(J_up>0);
else
    A_dir = double(J_low>0);
end


suffixes =  {'b', ...
            'lbc',...
            'lbf',...
            'lbt',...
            'lfc',...
            'lff',...
            'lft',...
            'lmc',...
            'lmf',...
            'lmt',...
            'rbc',...
            'rbf',...
            'rbt',...
            'rfc',...
            'rff',...
            'rft',...
            'rmc',...
            'rmf',...
            'rmt'};
G                  = digraph(A_dir);
close all; figure('color','w')
G_dir              = plot(G);
G_dir.NodeLabel    = suffixes;
G_dir.Marker       = '^';
G_dir.MarkerSize   = 30;
G_dir.NodeFontSize = 15;

p_R_c = @(q_c, x) expm(skew([x(1);x(2);x(3)])*(x(4)))*expm(skew([x(5);x(6);x(7)]).*q_c);
cprintf('*green', '>> Rotation matrices between sensors defined and joint axes found!!!\n')

%% ************************************************************************
%                      FIND JOINT POSITION VECTORS                        *
% *************************************************************************
clc

% test_points = sort(randperm(size(q_ref,2),3000));% choose k unique data points from data set
test_points = sort(randperm(size(q_signals,2),1000));% choose k unique data points from data set

rho_hat = zeros(6,constPar.noj);
r_hat_vel =rho_hat;
r_hat_alt =rho_hat;
[parents, children] = find(A_dir == 1);

for j = 1:numel(parents)
    p = parents(j);
    c = children(j);

    cprintf('*yellow', ['>> Checking: ' constPar.pxmark4.Bodies{constPar.imuBodyIndices(p)}.Name, ' -> ', constPar.pxmark4.Bodies{constPar.imuBodyIndices(c)}.Name '\n'])
    
    x0 = 1E-3*ones(6,1);
    [rho_hat(:,j), ~] = fmincon(...
    @(x) fcn_get_cartesian_acceleration_error_extended(...
              acc_signals_zpf(3*(p-1)+1:3*p, test_points), ...
              acc_signals_zpf(3*(c-1)+1:3*c, test_points), ...
              omg_signals_zpf(3*(p-1)+1:3*p, test_points), ...
              omg_signals_zpf(3*(c-1)+1:3*c, test_points), ...
             domg_signals_zpf(3*(p-1)+1:3*p, test_points), ...
             domg_signals_zpf(3*(c-1)+1:3*c, test_points), ...
             q_signals_zpf(c-1,test_points), ...
             p_R_c, squeeze(gamma_hat_off(p,c,:)), ...
             0.0, x), x0, [], [], [], [], [], [], [], options);     
disp('--------------------------------------------------------------------')
end
cprintf('*green', '>> Joint center point vectors found!!!\n')

%% **************************************************************************
%                         FORM PARAMETER VECTOR                           *
% *************************************************************************

lambda_hat_offline = NaN(13,constPar.noj);
for j = 1:numel(parents)
    p = parents(j);
    c = children(j);
    lambda_hat_offline(:,j) = [squeeze(gamma_hat_off(p,c,:));
                               rho_hat(:,j)];
end

cprintf('*green', '>> Parameters vector ready!!!\n')

%% ************************************************************************
%                      DISPLAY LEARNED KINEMATICS                         *
% *************************************************************************
clc

% Offline
h = figure('color','w');
fcn_phantomx_display_learned_morphology(lambda_hat_offline, p_R_c, constPar.pxmark4, constPar, mi_graphs, parents, children, h)