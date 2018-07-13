function [data, truth] = simPhysioNoise( data , cardiac_amp , resp_amp , mayer_amp )

if nargin<1 || ~exist('data','var') || isempty(data)
   [data, truth]=nirs.testing.simData; 
end
if nargin<2 || isempty(cardiac_amp)
    cardiac_amp = .25;
end
if nargin<3 || isempty(resp_amp)
    resp_amp = .25;
end
if nargin<3 || isempty(mayer_amp)
    mayer_amp = .25;
end

time = data.time;
[nsamp,nchan] = size(data.data);

% Generate cardiac oscillations
cardiac_freq = 1 + .1*randn;
cardiac_phase = cumsum( .01*2*pi*randn(nsamp,1) , 1 );
cardiac_data = sin( 2*pi*cardiac_freq*time + cardiac_phase );

% Generate respiratory oscillations
resp_freq = .25 + .025*randn;
resp_phase = cumsum( .01*2*pi*randn(nsamp,1) , 1 );
resp_data = sin( 2*pi*resp_freq*time + resp_phase );

% Generate Mayer waves
mayer_freq = .1 + .01*randn;
mayer_phase = cumsum( .01*2*pi*randn(nsamp,1) , 1 );
mayer_data = sin( 2*pi*mayer_freq*time + mayer_phase );

% Add noise to hb concentrations
link = data.probe.link;
Y    = data.data;

if(~(iscellstr(link.type) && any(ismember(link.type,{'hbo','hbr'}))))

    % optical density
    m = mean(Y);
    Y = bsxfun(@plus, -log(Y), log(m));
    
    % sort channels
    channels = nirs.util.uniquerows(link(:,1:2));    
    
    for j = 1:height(channels)
        
        lst = find(link.source==channels.source(j) & link.detector==channels.detector(j));
        
        assert( length(lst) > 1 )

        lambda = link.type(lst);
        ext = nirs.media.getspectra( lambda );

        clist = [1 2]; % hbo and hbr; need to fix this

        % extinction coefficients
        E = ext(:,clist);
        
        % sd distance
        L = data.probe.distances(lst);
        L=max(L,1);  % avoid issues with the short (0) seperation values

        % mbll model
        PPF = .1;
        EL = bsxfun( @times, E, L*PPF );
        iEL = pinv(EL);

        % calculates chromophore concentration (uM)
        hb = (Y(:,lst)*iEL') * 1e6;
        
        % add physiological noises
        for k = 1:size(hb,2)
            
            sigma = std(hb(:,k),0,1);
            hb(:,k) = hb(:,k) + cardiac_amp * sigma * cardiac_data ...
                              + resp_amp * sigma * resp_data ...
                              + mayer_amp * sigma * mayer_data;
        end
        
        % Convert back to OD
        Y(:,lst) = hb./1e6 * EL';
        
    end
    
    % Convert OD to intensity
    Y = exp( -bsxfun(@minus, Y, log(m)) );
    
else
    
    % add physiological noise
    for k = 1:size(Y,2)

        sigma = std(Y(:,k),0,1);
        Y(:,k) = Y(:,k) + cardiac_amp * sigma * cardiac_data ...
                        + resp_amp * sigma * resp_data ...
                        + mayer_amp * sigma * mayer_data;
    end

end

data.data = Y;

if(nargout>1)
    if(~exist('truth','var'))
        truth=[];
    end

end