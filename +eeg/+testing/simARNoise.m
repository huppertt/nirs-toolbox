function data = simARNoise( probe, t, P, sigma)
    
    if (nargin < 4 || isempty(sigma)), sigma=.33; end;
    if (nargin < 3  || isempty(P)), P = 10; end
    if (nargin < 2 || isempty(t)), t = (0:1/200:300)'; end
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
    data = eeg.core.Data();
    data.data   = e * 5e-3;
    data.probe  = probe;
    data.time   = t;
  
end

function a = randAR( P )
    % random Pth order AR coef    
    a = flipud( cumsum( rand(P, 1) ) );
    a = a / sum(a) * 0.99;
end

function probe = defaultProbe()
    labels={'Fp1'    'Fp2'    'F7'    'F3'    'Fz'    'F4'    'F8'    'FC5',...
            'FC1'    'FC2'    'FC6'    'T7'    'C3'    'Cz'    'C4'    'T8',...
            'TP9'    'CP5'    'CP1'    'CP2'    'CP6'    'TP10'    'P7'    'P3',...
            'Pz'    'P4'    'P8'    'POz'    'O1'    'Oz'    'O2'};
    probe = eeg.core.Probe(labels);
    
end