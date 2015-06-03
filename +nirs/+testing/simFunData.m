function [data, truth] = simFunData( data, beta, stimLength, stimSpace, basis )
    if nargin < 5
        basis = Dictionary({'default'}, {nirs.design.basis.Canonical()});
    end
    
    if nargin < 4, stimSpace    = 10;   end
    if nargin < 3, stimLength   = 2;    end
    if nargin < 2, beta         = 5;    end
    
    if nargin < 1
        data = nirs.testing.simARNoise();
    end
    
    [data.probe.link, i] = sortrows(data.probe.link, {'source', 'detector', 'type'});
    data.data = data.data(:,i);
    
    % sd pairs
    [SD, ~, iSD] = unique( data.probe.link(:,1:2), 'rows' );

	nSD = max(iSD);
    
    % num channels to add stim to
    nAct = floor( nSD / 2 );
    if mod(nSD, 2)
        nAct = nAct + (rand() < 0.5);
    end
    
    % choose channels
	iAct = randperm(nSD, nAct);
    
    % design mat
    [X, stim] = nirs.testing.randDesignMat( data.time, stimLength, stimSpace, basis );
    data.stimulus = stim;
    
    % add stimulus
    truth = zeros( size( data, 1), 1 );
    for i = 1:length(iAct)
        lst = iSD == iAct(i);
        
        lambda = data.probe.link.type(lst);
        e = nirs.media.getspectra( lambda );
        
        xhbo =  X*beta;
        xhbr = -X*beta;
        
        tmp =  [xhbo xhbr] * e(:,1:2)' * 1e-6 * diag(data.probe.distances(lst)) * 5/50;
        
        data.data(:,lst) = data.data(:,lst) + tmp;
        
        truth = truth + lst;
    end
    
end

