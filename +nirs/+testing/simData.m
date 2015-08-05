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
        beta = 10*ones( length(stim.keys), 1 );
    end
    
    if length(beta) == length(stim.keys)
        % oxy; deoxy
        b = [beta; -beta/3];
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
    
    % optical density
    m = mean(Y);
    Y = bsxfun(@plus, -log(Y), log(m));
        
    for i = 1:size(channels, 1)
        lst = find(link.source == channels(i,1) ...
            & link.detector == channels(i,2));
        
        % extincion coefs
        lambda = link.type(lst);
        e = nirs.media.getspectra(lambda);
        e = e(:,1:2);
        
        
        % sd distance
        l = data.probe.distances(lst);
                
        % design mat
        Xhbo = nirs.design.createDesignMatrix( stim, data.time, basis, 'hbo' );
        Xhbr = nirs.design.createDesignMatrix( stim, data.time, basis, 'hbr' );
        
        % add to channels according to MBLL
        for j = 1:length(lst)
            Yact = [Xhbo*e(j,1)*l(j) Xhbr*e(j,1)*l(j)] * b * 5/50 * 1e-6;
            Y(:,lst(j)) = Y(:,lst(j)) + Yact;
        end        

        truth(lst) = 1;        
    end
    
    Y = exp( -bsxfun(@minus, Y, log(m)) );
    
    data.data = Y;
    data.stimulus = stim;
end

