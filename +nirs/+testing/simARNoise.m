function data = simARNoise( probe, t, P, sigma, N_files )
    
    if (nargin < 5 || isempty(N_files)), N_files=1; end
    if (nargin < 4 || isempty(sigma)), sigma=.33; end;
    if (nargin < 3  || isempty(P)), P = 10; end
    if (nargin < 2 || isempty(t)), t = (0:1/10:300)'; end
    if (nargin < 1 || isempty(probe)), probe = defaultProbe(); end
    
    if(strcmp(class(probe),'double'))
        probe = defaultProbe(probe);
    end
    
    if (N_files>1)
        data(1:N_files,1) = nirs.core.Data;
        for i = 1:N_files
            data(i) = nirs.testing.simARNoise( probe, t, P, sigma, 1 );
        end
        return;
    end
    
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

function probe = defaultProbe(lambda)

if(nargin==0)
    lambda=[690 830];
end
    
    
    srcPos(:,1) = (-80:20:80)';
    srcPos(:,2:3) = 0;
    
    detPos(:,1) = (-70:20:70)';
    detPos(:,2) = 25;
    detPos(:,3) = 0;
    
    probe = nirs.core.Probe(srcPos,detPos);
    
    link = [1	1	
        2	1	
        2	2	
        3	2	
        3	3	
        4	3	
        4	4	
        5	4	
        5	5	
        6	5	
        6	6	
        7	6	
        7	7	
        8	7	
        8	8	
        9	8	];
    link=[repmat(link,length(lambda),1) reshape(repmat(lambda(:)',size(link,1),1),1,[])'];
    
    link = sortrows(link);
    
    probe.link = table(link(:,1), link(:,2), link(:,3), ...
        'VariableNames', {'source', 'detector', 'type'});
    
end

function [r,T]=mvnrnd(mu,sigma,cases)


[T,err] = chol(sigma);
r = randn(cases,size(T,1)) * T + ones(cases,1)*mu';

end

