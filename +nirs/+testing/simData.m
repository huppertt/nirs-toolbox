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
%
%     [data, truth] = simData( noise, stim, beta, channels )

if nargin < 1 || isempty(noise)
    noise = nirs.testing.simARNoise();
end

if nargin < 2 || isempty(stim)
    stim = nirs.testing.randStimDesign(noise.time, 2, 7, 1);
end

if nargin < 3 || isempty(beta)
    beta = 7*ones( length(stim.keys), 1 );
elseif(isstr(beta))
    snr = str2num(beta(strfind(beta,'SNR:')+4:end));
    beta=snr*sqrt(var(noise.data(:)));
end

if length(beta) == length(stim.keys)
    % oxy; deoxy
    b = [beta; -beta/2];
else
    b = beta;
end

if nargin < 5 || isempty(basis)
    % default to canonical basis
    basis = Dictionary({'default'}, {nirs.design.basis.Canonical()});
end

if nargin < 4 || isempty(channels)
    % default to first half of channels
    sd = [noise.probe.link.source noise.probe.link.detector];
    
    % make sure there are no short distance here
    if(ismember('ShortSeperation',noise.probe.link.Properties.VariableNames))
        lstSS=find(noise.probe.link.ShortSeperation);
        sd(lstSS,:)=[];
    else
        lstSS=[];
    end
    sd=unique(sd,'rows');
    channels = sd(1:round(end/2),:);
end

% loop through and add
data = noise.sorted();

link = data.probe.link;
Y    = data.data;
truth = zeros(size(Y,2), 1);

if(~(iscellstr(link.type) && any(ismember(link.type,{'hbo','hbr'}))))
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
            Yact = [Xhbo*e(j,1)*l(j) Xhbr*e(j,2)*l(j)] * b * 5/50 * 1e-6;
            Y(:,lst(j)) = Y(:,lst(j)) + Yact;
        end
        
        truth(lst) = 1;
    end
    
    Y = exp( -bsxfun(@minus, Y, log(m)) );
else
    for i = 1:size(channels, 1)
        lstHbO = find(link.source == channels(i,1) ...
            & link.detector == channels(i,2) & ismember(link.type,'hbo'));
        lstHbR = find(link.source == channels(i,1) ...
            & link.detector == channels(i,2) & ismember(link.type,'hbr'));
        Xhbo = nirs.design.createDesignMatrix( stim, data.time, basis, 'hbo' );
        Xhbr = nirs.design.createDesignMatrix( stim, data.time, basis, 'hbr' );
        Y(:,lstHbO)=Y(:,lstHbO)+Xhbo*b(1)*1e-6;
        Y(:,lstHbR)=Y(:,lstHbR)+Xhbr*b(2)*1e-6;
        truth(lstHbO) = 1;
        truth(lstHbR) = 1;
    end
    
    
end

truth(lstSS)=NaN;  % use NaN to mask out the short seperation channels

data.data = Y;
data.stimulus = stim;
end

