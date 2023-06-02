% ########################################################################
% 
% 
%                    Mutual Information calculation                       
% 
% 
% #########################################################################

function [MI_cmean, ...
          MI_batch, ...
          MI_kin_mst, ...
          spectralDistance, ...
          frobeniusDistance, ...
          sumSquaredDistance, ...
          grad_ssd_mavg, ...
          total_mutal_information, ...
          edge_error] = ...
               fcn_compute_mi_matrix_for_streaming_signals(signals, batchSize, buffer_size, CONTRACT, method, constPar)
    
    clc
    
    kin_indices = [constPar.noj + 1:2*constPar.noj, ...
                   3*constPar.noj + 1:3*constPar.noj + constPar.nob];

    % JIDT settings =======================================================
    % * NOTE: Kernel and Gaussian versions were tested
    
    % Construct the calculator:
    % calc = javaObject('infodynamics.measures.continuous.kernel.MutualInfoCalculatorMultiVariateKernel');
    calc = javaObject('infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian');
    % calc = javaObject('infodynamics.measures.continuous.kozachenko.EntropyCalculatorMultiVariateKozachenko');
    % calc = javaObject('infodynamics.measures.continuous.kraskov.MutualInfoCalculatorMultiVariateKraskov2');
    
    
    % Set any properties to non-default values:
    KERNEL_WIDTH = 0.8;
    calc.setProperty('KERNEL_WIDTH', num2str(KERNEL_WIDTH));
    calc.setProperty('NORMALISE', 'false');
    calc.setProperty('k','12');
    % =====================================================================
    
    data = transpose(signals);
    N_points                = size(data,1);
    N_signals               = size(data,2);
    
    num_batches             = floor(N_points/batchSize);
    % Full MI matrix
    MI_batch                = zeros(N_signals,N_signals,num_batches);

    if ~CONTRACT
        N_nodes             = N_signals;
        MI_cmean            = zeros(N_nodes,N_nodes,num_batches);
    else
        % Contracted MI matrix
        N_nodes             = 3*constPar.noj + 2*constPar.nob;
        MI_cmean            = zeros(N_nodes,N_nodes,num_batches);
%         MI_kin_mst              = [];
    end
    edge_error              = NaN(1,num_batches);
    total_mutal_information = NaN(1,num_batches);
    spectralDistance        = NaN(1,num_batches-1);
    frobeniusDistance       = NaN(1,num_batches-1);
    sumSquaredDistance      = NaN(1,num_batches-1);
    sumSquaredDistance_mavg = NaN(1,num_batches-1);
    grad_ssd_mavg           = NaN(1,num_batches-1);
    
    window                  = 100; % sliding window samples
    
    % Buffer properties
    buffer = zeros(1,buffer_size);
    
    disp(['>> Computing MI using ' method '...\n'])
    disp(['With ', num2str(num_batches), ' batches of ', num2str(batchSize), ' samples.'])
    pause(3)
    tic
    %rng('default')
    ssd_prev      = 0;
    batch         = 0;
    sampleCounter = 0;
    SampleSet     = 1:N_points;%[1:N_points,1:N_points,1:N_points,1:N_points];
    for sample = SampleSet
        sampleCounter = sampleCounter + 1;
        
        % Store data in buffer --------------------------------------------
        if sampleCounter<=buffer_size    
            % Store sample in order as it comes
            buffer(sample) = sample;
            buffer_full    = 0;
        else
            % Reservoir sampling, see: https://sisudata.com/blog/simple-streaming-reservoir-sampling
            k = randi(sampleCounter);
            if k < buffer_size
                buffer(k) = sample;
            end
            buffer_full = 1;
        end
    
        % Estimate MI every <batchSize> samples ---------------------------
        if sampleCounter>=batchSize && mod(sampleCounter, batchSize) == 0
            % Sample a mini batch from buffer
            if buffer_full == 0
                sample_points   = randperm(sample,batchSize);
            elseif buffer_full == 1
                % Randomly sample from data without replacement
                sample_points   = datasample(buffer,batchSize,'Replace',false);
            end
            batch = batch + 1;
            % Normalize selected samples
            miniBatch      = data(sample_points, :);
            miniBatch_mean = mean(miniBatch,1);
            miniBatch_std  = std(miniBatch,[],1); miniBatch_std(miniBatch_std == 0) = 1; % to avoid division by zero when centering w.r.t. the standard deviation
            miniBatch      = (miniBatch - miniBatch_mean)./miniBatch_std;             
%             disp(['BATCH: ',num2str(batch), '/',num2str(num_batches),' | BUFFER_STATUS: ', num2str(min(round(100*sampleCounter/buffer_size),100)),' % | SAMPLE: ', num2str(sampleCounter),'/',num2str(numel(SampleSet))])
            fprintf('BATCH: %d/%d | BUFFER_STATUS: %d %% | SAMPLE: %d/%d\n',batch,num_batches,min(round(100*sampleCounter/buffer_size),100),sampleCounter,numel(SampleSet));
            % Compute for all pairs:
            for s = 1:1:size(MI_batch,1)
                % Select source signal
                source      = miniBatch(:, s);
                for d = s:1:size(MI_batch,2)
                    if s==d
                        continue;
                    end
                    % Select destination signal
                    destination      = miniBatch(:, d);
                    switch method
                        case 'jidt'
                            % Initialise the calculator for (re-)use:
                            calc.initialise();
                            % Supply the sample data:
                            calc.setObservations(source, destination);
                            % Compute the estimate:
                            result = calc.computeAverageLocalOfObservations();
                        case 'standard'
                            result = mutualinfo(source, destination);
                    end

                    if num_batches == 1
                                fprintf('MI_Kernel(row_%d -> col_%d) = %.4f bits, BATCH: %d\n', ...
                                    s, d, result, batch);
                    end
                    MI_batch(s,d,batch) = result;
                    MI_batch(d,s,batch) = result;
                end
            end
            
            if ~CONTRACT
                MI_b = MI_batch(:,:,batch);
            else
                MI_b = fcn_contract_mi_matrix(MI_batch(:,:,batch), constPar);
            end


            if batch == 1
                MI_cmean(:,:,batch) = MI_b;
            else
                % Cumulative mean
                MI_cmean(:,:,batch) = (batch-1)/batch*MI_cmean(:,:,batch-1) + MI_b/batch;
                
                [spectralDistance(batch-1), ...
                 frobeniusDistance(batch-1), ...
                 sumSquaredDistance(batch-1)] = ...
                        fcn_compute_adjacency_matrix_distances(MI_cmean(:,:,batch-1), MI_cmean(:,:,batch));
                if  batch < window + 1
                    % If no. of badges smaller than window compute
                    % comulative average
                    sumSquaredDistance_mavg(batch-1) = 1/window*sum(sumSquaredDistance(1:batch-1));
                else
                    % Moving average of the sumSquaredDistance with a
                    % sliding window of 100 steps
                    sumSquaredDistance_mavg(batch-1) = 1/window*sum(sumSquaredDistance(batch-window:batch-1));                    
                end
                grad_ssd_mavg(batch-1)           = abs(sumSquaredDistance_mavg(batch-1) - ssd_prev);
                ssd_prev                         = sumSquaredDistance_mavg(batch-1);                
            end
            MI_kin                 = MI_cmean(kin_indices,kin_indices,batch);
            [G_kin_mst,~]          = minspantree(graph(-MI_kin,'upper'));
            G_kin_mst.Edges.Weight = abs(G_kin_mst.Edges.Weight);
%             edge_error(batch)      = 0.5*sum(triu(abs(full(G_kin_mst.adjacency) - constPar.A_kin_ref)),'all');
            
            % Get the total information in the matrix
            D_cmean                        = diag(sum(MI_cmean(:,:,batch),2));
            total_mutal_information(batch) = 0.5*trace(D_cmean);

            if batch == 1
                disp(['TMI: ',num2str(total_mutal_information(batch)),' | Grad_MTI: ', 'NA',' | edge error: ', num2str(edge_error(batch))])      
            else
                dtotal_mutal_information = total_mutal_information(batch) - total_mutal_information(batch-1);
                disp(['TMI: ',num2str(total_mutal_information(batch)),' | Grad_MTI: ',num2str(dtotal_mutal_information),' | edge error: ', num2str(edge_error(batch))])      
            end
            % Uncomment if function is to exit when convergence is achieved
            %Dist_mat = distances(G_kin_mst,'Method','unweighted');
            %if all(min(Dist_mat(1:constPar.noj,1:constPar.noj)+10*eye(constPar.noj),[],2)>=3) && batch>1
            %    if edge_error(batch) == 0
            %        break;
            %    end
            %end

        end
    end
    % Vectorize MI MST
    At         = (full(G_kin_mst.adjacency('weighted'))).';
    m          = triu(true(size(At)),1);
    MI_kin_mst = transpose(At(m).');        
    disp('>> Calculation completed!\n')
    toc
end
