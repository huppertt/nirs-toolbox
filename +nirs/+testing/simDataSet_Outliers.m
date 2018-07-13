function [data, truth] = simDataSet_Outliers( noise, ngroup, stimFunc, beta, channels, testingfcn )
%% simDataSet - simulates a dataset with outliers for ROC testing
% 
% Args:
%     noise - a list of nirs.core.Data objects to use as baseline noise
%     ngrou - number of groups to split data into
%     stimFunc - a function which takes in a time vector and returns a stim design
%     beta     - nconditions x ngroups giving the magnitude of responses of 
%                each condition for each group
%     channels - n x 2 array specifying the SD pairs to add activity to

    if nargin < 1 || isempty(noise)
        for i = 1:30
            noise(i,1) = nirs.testing.simARNoise();
            lst=randi(size(noise(i,1).data,2),4,1); lst=unique(lst);
            
            % this adds extra random noise to a few channels
             noise(i,1).data(:,lst)=noise(i,1).data(:,lst)+randn(size( noise(i,1).data(:,lst)))*10;
        end
          
    end
    
    if(~isa(noise,'nirs.core.Data'))
        id=noise; clear noise;
         for i = 1:id
            noise(i,1) = nirs.testing.simARNoise();
            
             % this adds extra random noise to a few channels
             lst=randi(size(noise(i,1).data,2),4,1); lst=unique(lst);
             noise(i,1).data(:,lst)=noise(i,1).data(:,lst)+randn(size( noise(i,1).data(:,lst)))*10;
        end
    end
    
    if nargin < 2 || isempty(ngroup)
        ngroup = 1;
    end
    
    if nargin < 3 || isempty(stimFunc)
        stimFunc = @(t) nirs.testing.randStimDesign(t, 2, 7, 1);
    end
    
    if(iscell(stimFunc))
        s=stimFunc{1};
    else
        
        s = stimFunc(noise(1).time);
    end
    if nargin < 4 || isempty(beta)
        beta = 7*ones( length(s.keys), ngroup )/sqrt(length(noise));
    end
    
    if size(beta,1) == length(s.keys)
        % oxy; deoxy
        b = [beta; -beta/3];
    end
        
    if nargin < 5 || isempty(channels)
        sd = unique([noise(1).probe.link.source noise(1).probe.link.detector], 'rows');
        channels = sd(1:round(end/2),:);
    end
    
    if nargin < 6 || isempty(testingfcn)
        testingfcn=@nirs.testing.simData;
    end
    
    data = noise;
    
    % group index
	gidx = nirs.testing.randSplit(length(noise), ngroup);
    
    % loop through
    for i = 1:length(noise)
        data(i).demographics('group') = ['G' num2str(gidx(i))];
        data(i).demographics('subject') = ['S' num2str(i)];
        if(iscell(stimFunc))
            [data(i), truth] = feval(testingfcn,data(i), stimFunc{i}, b(:, gidx(i)), channels);
        else
            [data(i), truth] = feval(testingfcn,data(i), stimFunc(data(i).time), b(:, gidx(i)), channels);
        end
    end
    
end