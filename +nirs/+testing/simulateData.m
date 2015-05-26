function data = simulateData( data, beta, stimLength, stimSpace, basis )
    if nargin < 5
        basis = Dictionary({'default'}, {nirs.design.basis.Canonical()});
    end
    
	nChan = size( data(1).data,2 );

    % choose channels
	iChan = randperm(nChan,nChan/2);
    
end

