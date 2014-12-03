function optProp = nlSlabSolver( data, varargin )
%SLABSOLVER Summary of this function goes here
%   Detailed explanation goes here
    
    warning('off','stats:statrobustfit:IterationLimit')
    
    if ~isa( data, 'nirs.Data' )
        error( 'Data must be a Data object.' )
    end

    if nargin > 1 && varargin{1} >= 1
        ri = varargin{1};
    else
        ri = 1.45;
    end

    data.data = double( data.data );
    
    optProp = nirs.OpticalProperties();
    optProp.ri = ri;

    iLambda = sort( unique(data.probe.link(:,3)) );
    for i = 1:length( iLambda )

        lst = data.probe.link(:,3) == iLambda(i);

        ac = abs( data.data(:,lst) );
        phs = unwrap( angle(data.data(:,lst)),[],2 );
        d = data.probe.distances(lst);
        
        if size(ac,1) > 1
            w_ac = 1./std( log(ac .* repmat((d.^2).',[size(ac,1) 1])) ).';
            w_phs = 1./std( phs ).';
        else
            w_ac = ones(size(ac,2),1);
            w_phs = ones(size(phs,2),1);
        end
        
        X = [d ones(size(d))];
        
        for j = 1:size(ac,1)
            y = log( ac(j,:).' .* d.^2 );
%              alpha = robustfit(X,y,[],[],'off');
            alpha = (diag(w_ac)*X) \ (w_ac.*y);
            alpha = alpha(1);


            y = -phs(j,:)';
%             phi = robustfit(X,y,[],[],'off');
            phi = (diag(w_phs)*X) \ (w_phs.*y);
            phi = phi(1);

            thisMua = (data.modFreq*2*pi)/optProp.c/2 * (phi/alpha - alpha/phi);
            thisMus = (alpha^2 - phi^2)/3/thisMua - thisMua;
            
            lastMua = Inf;
            
            iter = 0;
            while abs( (lastMua - thisMua)/thisMua ) > 1e-3 && iter < 100
                lastMua = thisMua;
                
                DD = 1/3/(thisMua + thisMus);
                x =  (data.modFreq*2*pi)/optProp.c/thisMua;
                
                % alpha
%                 zb = 1.0; % mm
%                 z0 = 0.6; % mm
%                 z = 1; % mm
%                 
%                 Vp = sqrt(sqrt(1+x^2)+1);
                
                G = 1 + d * sqrt(2*thisMua/DD) + d.^2 * thisMua/DD * sqrt(1+x^2);
                
%                 Fac = (zb + z0) * (z + 3*DD*(...
%                     1-((zb+z0)^2+3*z^2)./(2*d.^2) ...
%                     .*(2 + (1+d.*sqrt(thisMua/2/DD)*Vp./G + d*sqrt(2*thisMua/DD)*Vp))...
%                     ));
                
%                 y = log( ac(j,:)' .* d.^3 ./ ...
%                 sqrt( G ) ./ Fac ...
%                 );
                y = log( ac(j,:)' .* d.^3 ./ ...
                    sqrt( G ));
            
                alpha = (diag(w_ac)*X) \ (w_ac.*y);
                alpha = alpha(1);
                
                % phi
                y = -phs(j,:)' + atan( d*sqrt(thisMua/2/DD)*sqrt(sqrt(1+x^2)-1) ./ ...
                    (d*sqrt(thisMua/2/DD)*sqrt(sqrt(1+x^2)+1)) ...
                    );
                
                phi = (diag(w_phs)*X) \ (w_phs.*y);
                phi = phi(1);
                
                thisMua = (data.modFreq*2*pi)/optProp.c/2 * (phi/alpha - alpha/phi);
                thisMus = (alpha^2 - phi^2)/3/thisMua - thisMua;
                
                iter = iter + 1;
            end
            
            mua(i,j) = thisMua;
            mus(i,j) = thisMus;
        end
 
    end
    
%     mua_s = std(mua,[],2);
%     mua = mean(mua,2);
%     mus_s = std(mus,[],2);
%     mus = mean(mus,2);
    
    optProp = nirs.OpticalProperties( mua.',mus.',ri,data.probe.lambda(iLambda));
%     optProp.mua_s = mua_s.';
%     optProp.mus_s = mus_s.';


end