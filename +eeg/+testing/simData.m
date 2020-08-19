function [data, truth] = simData( noise, stim, beta, channels, basis )
%SIMDATA Simulates EEG data by adding a task to baseline noise.
% 
% Args:
%     noise  -  a raw data file
%     stim   -  dictionary containing stimulus objects using the stim name as key
%     beta   -  vector of magnitude of response for each stim conditon
%     basis  -  dictionary containing basis objects using stim condition as key
%     
% Example:
%     noise = eeg.testing.simARNoise();
%     stim  = nirs.testing.randStimDesign(noise.time, .1, 3, 1);
%     beta  = [3]';
%     
%     sd = noise.probe.link.electrode;
%     channels = sd(1:round(end/2),:);
% 
%     [data, truth] = simData( noise, stim, beta, channels )
    
    if nargin < 1 || isempty(noise)
        noise = eeg.testing.simARNoise();
    end
    
    if nargin < 2 || isempty(stim)
        stim = nirs.testing.randStimDesign(noise.time, .1, 3, 1);
    end
    
    if nargin < 3
        beta = 3*ones( length(stim.keys), 1 );
    end
        
    if nargin < 5 || isempty(basis)
        % default to canonical basis
        basis = Dictionary({'default'}, {eeg.design.basis.ERP()});
    end
    
    if nargin < 4 || isempty(channels)
        % default to first half of channels
        sd = noise.probe.link.electrode;
        channels = sd(1:round(end/2));
    end
    
    % loop through and add
    data = noise.sorted();
    
    link = data.probe.link;
    Y    = data.data;
    truth = zeros(size(Y,2), 1);
    
    beta=beta*mad(Y(:));
    X = nirs.design.createDesignMatrix( stim, data.time, basis);
    lst=ismember(data.probe.link.electrode,channels);
    
    Y(:,lst)=Y(:,lst)+X*(beta*ones(1,length(find(lst))));
    
    truth(lst) = 1;        
    
    data.data = Y;
    data.stimulus = stim;
end

