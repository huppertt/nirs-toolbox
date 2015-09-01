function [data,truth] = simData_connectivity(noise)

if nargin < 1 || isempty(noise)
    noise = nirs.testing.simARNoise([],[],[],0);
end

pmax=4;

data = noise;
Y =  data.data;
link = data.probe.link;

channels = unique([noise.probe.link.source noise.probe.link.detector], 'rows');

nchan=size(channels,1);
truth = (rand(nchan)>.5);
truth=max(truth,eye(nchan));

ehbo = randn(length(noise.time),nchan);
ehbr = randn(length(noise.time),nchan);
hbo = zeros(length(noise.time),nchan);
hbr = zeros(length(noise.time),nchan);


for i=1:nchan
    for j=1:nchan
        if(truth(i,j))
            a = flipud( cumsum( rand(pmax, 1) ) );
            a = a / sum(a) * 0.99;
            hbo(:,i)=hbo(:,i) + filter(1, [1; -a], ehbo(:,i));
            hbr(:,i)=hbr(:,i) + filter(1, [1; -a], ehbr(:,i));
        end
    end
end



% optical density
m = mean(Y);
Y = bsxfun(@plus, -log(Y), log(m));

mm = var(Y,[],1);

for i = 1:size(channels, 1)
    lst = find(link.source == channels(i,1) ...
        & link.detector == channels(i,2));
    
    % extincion coefs
    lambda = link.type(lst);
    e = nirs.media.getspectra(lambda);
    e = e(:,1:2);
    
    % sd distance
    l = data.probe.distances(lst);
    
    % add to channels according to MBLL
    for j = 1:length(lst)
        Yact = (e(j,1)*l(j)*hbo(:,i) + e(j,2)*l(j)*hbr(:,i))*1E-6;
        Y(:,lst(j)) = Y(:,lst(j)) + Yact;
    end
end

Y = Y./(ones(length(noise.time),1)*(sqrt(var(Y,[],1))./mm));

Y = exp( -bsxfun(@minus, Y, log(m)) );

data.data = Y;
data.stimulus('none') = nirs.design.StimulusEvents();
end