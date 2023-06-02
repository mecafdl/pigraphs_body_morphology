classdef replayBuffer < handle
   
    properties
        Size     = [];
        Memory   = [];
        Status   = 0;
        Capacity = 0;
        LogIndex = 0;
    end

    methods
        function obj = replayBuffer(bufferSize, dataPointDimension)
            obj.Size   = bufferSize;
            obj.Memory = zeros(dataPointDimension, bufferSize);
        end

        function storeSample(obj, dataPoint)
            obj.LogIndex = obj.LogIndex + 1; 
            if obj.LogIndex <= obj.Size   
                obj.Memory(:,obj.LogIndex) = dataPoint;
            else
                % Change buffer status to full
                obj.Status = 1;
                k = randi(obj.LogIndex);
                if k < obj.Size
                    obj.Memory(:,k) = dataPoint;
                end
            end
            % Compute the percentage of buffer occupancy
            obj.Capacity = floor(min(100,100*(obj.LogIndex)/obj.Size)); 
        end
        
        function sample = drawSamples(obj, batchSize)
            sample = [];
            if obj.Status ~=1
               if (obj.LogIndex >= batchSize)                              
                    sampleIndices = randperm(obj.LogIndex, batchSize);
                    sample        = obj.Memory(:,sampleIndices);
               else
                    warning('Not enough samples in storage.')
               end
            else
                sampleIndices = randperm(obj.Size, batchSize);
                sample        = obj.Memory(:,sampleIndices);
            end
        end

        function flush(obj)
            obj.Memory   = zeros(size(obj.Memory));
            obj.Status   = 0;
            obj.Capacity = 0;
            obj.LogIndex = 0;    
        end        

    end
end