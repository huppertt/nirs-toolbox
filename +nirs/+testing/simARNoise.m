function data = simARNoise( probe, t, P, sigma)
    
    if (nargin < 4 || isempty(sigma)), sigma=.33; end;
    if (nargin < 3  || isempty(P)), P = 10; end
    if (nargin < 2 || isempty(t)), t = (0:1/10:300)'; end
    if (nargin < 1 || isempty(probe)), probe = defaultProbe(); end
    
    nchan = size(probe.link,1);
    
    % noise mean and spatial covariance
    mu = zeros(nchan,1);
    S = toeplitz( [1 sigma*ones(1,nchan-1)] );
    
    e = mvnrnd( mu, S, length(t) );

    % add temporal covariance
    for i = 1:size(e,2)
        a = randAR( P );
        e(:,i) = filter(1, [1; -a], e(:,i));
    end
    
    % output
    data = nirs.core.Data();
    data.data   = 100 * exp( - e * 5e-3 );
    data.probe  = probe;
    data.time   = t;
  
end

function a = randAR( P )
    % random Pth order AR coef    
    a = flipud( cumsum( rand(P, 1) ) );
    a = a / sum(a) * 0.99;
end

function probe = defaultProbe()

    
    
    srcPos(:,1) = (-80:20:80)';
    srcPos(:,2:3) = 0;
    
    detPos(:,1) = (-70:20:70)';
    detPos(:,2) = 25;
    detPos(:,3) = 0;
    
    probe = nirs.core.Probe(srcPos,detPos);
    
    link = [1	1	690
        2	1	690
        2	2	690
        3	2	690
        3	3	690
        4	3	690
        4	4	690
        5	4	690
        5	5	690
        6	5	690
        6	6	690
        7	6	690
        7	7	690
        8	7	690
        8	8	690
        9	8	690
        1	1	830
        2	1	830
        2	2	830
        3	2	830
        3	3	830
        4	3	830
        4	4	830
        5	4	830
        5	5	830
        6	5	830
        6	6	830
        7	6	830
        7	7	830
        8	7	830
        8	8	830
        9	8	830];
    
    link = sortrows(link);
    
    probe.link = table(link(:,1), link(:,2), link(:,3), ...
        'VariableNames', {'source', 'detector', 'type'});
    
end