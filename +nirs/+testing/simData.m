function [data, truth] = simData( noise, stim, beta, channels, basis )
%SIMDATA Simulates NIRS data by adding a task to baseline noise.
% 
% Args:
%     noise  -  a raw data file
%     stim   -  dictionary containing stimulus objects using the stim name as key
%     beta   -  vector of magnitude of response for each stim conditon
%     basis  -  dictionary containing basis objects using stim condition as key
%     
% Example:
%     noise = nirs.testing.simARNoise();
%     stim  = nirs.testing.randStimDesign(noise.time, 2, 7, 3);
%     beta  = [3 2 1]';
%     
%     sd = unique([noise.probe.link.source noise.probe.link.detector], 'rows');
%     channels = sd(1:round(end/2),:);
    
    if nargin < 1
        noise = nirs.testing.simARNoise();
    end
    
    if nargin < 2
        stim = nirs.testing.randStimDesign(noise.time, 2, 7, 1);
    end
    
    if nargin < 3
        beta = ones( length(stim.keys), 1 );
    end
    
    if nargin < 5
        % default to canonical basis
        basis = Dictionary({'default'}, {nirs.design.basis.Canonical()});
    end
    
    if nargin < 4
        % default to first half of channels
        sd = unique([noise.probe.link.source noise.probe.link.detector], 'rows');
        channels = sd(1:round(end/2),:);
    end
    
    % loop through and add
    data = noise.sorted();
    
    link = data.probe.link;
    Y    = data.data;
    truth = zeros(size(Y,2), 1);
    
    for i = 1:size(channels, 1)
        lst = link.source == channels(:,1) & link.detector == channels(:,2);
        
        truth(lst) = 1;
        %%% TODO 
        
    end
    

end

