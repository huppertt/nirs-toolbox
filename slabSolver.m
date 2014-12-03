function optProp = slabSolver( data, varargin )
%SLABSOLVER Summary of this function goes here
%   Detailed explanation goes here
    
    warning('off','stats:statrobustfit:IterationLimit')
    
    if ~isa( data, 'nirs.Data' )
        error( 'Data must be a Data object.' )
    end

    if nargin > 1 && varargin{1} >= 1
        ri = varargin{1};
    else
        ri = 1.3;
    end

    data.data = double( data.data );
    
    optProp = nirs.OpticalProperties();
    optProp.ri = ri;

    iLambda = sort( unique(data.probe.link(:,3)) );
    for i = 1:length( iLambda )

        lst = data.probe.link(:,3) == iLambda(i);

        % should replace median with better estimator
%         ac = abs( median(data.data(:,lst),1) )';
%         phs = median( angle( data.data(:,lst) ),1 )';
        ac = abs( data.data(:,lst) );
        phs = unwrap( angle(data.data(:,lst)),[],2 );
        d = data.probe.distances(lst);

        X = [d ones(size(d))];
        
        for j = 1:size(ac,1)
            y = log( ac(j,:)' .* d.^2 );
% %             alpha = robustfit(X,y,[],[],'off');
            alpha = X \ y;
            alpha = alpha(1);


            y = -phs(j,:)';
%             phi = robustfit(X,y,[],[],'off');
            phi = X \ y;
            phi = phi(1);

            mua(i,j) = (data.modFreq*2*pi)/optProp.c/2 * (phi/alpha - alpha/phi);
            mus(i,j) = (alpha^2 - phi^2)/3/mua(i,j) - mua(i,j);
        end
 
    end
    
%     mua = trimmean(mua,20,2);
%     mus = trimmean(mus,20,2);
    mua = mean(mua,2);
    mus = mean(mus,2);
    
    optProp = nirs.OpticalProperties( mua,mus,ri,data.probe.lambda(iLambda) );

end

