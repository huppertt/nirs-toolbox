function prop = slabSolver( data, ri )
%SLABSOLVER Summary of this function goes here
%   Detailed explanation goes here
    
    warning('off','stats:statrobustfit:IterationLimit')

    if nargin < 2, ri = 1.45; end
    
    data.data = double( data.data );
    
    if length( data.probe.lambda ) == 1
        prop = nirs2.OpticalProperties();
        prop.ri = ri;
    else
        prop = nirs2.SpectralProperties( 0.7, 60, data.probe.lambda );
        prop.ri = ri;
    end

    iLambda = sort( unique(data.probe.link(:,3)) );
    for i = 1:length( iLambda )

        lst = data.probe.link(:,3) == iLambda(i);

        ac = abs( data.data(:,lst) );
        phs = angle( data.data(:,lst) );
%         phs = unwrap( angle(data.data(:,lst)),[],2 );
        d = data.probe.distances(lst);
        
        if size(ac,1) > 1
%             w_ac = 1./mad(log(ac),1,1).'; 
% %             %1./std( log(ac .* repmat((d.^2).',[size(ac,1) 1])) ).';
%             w_phs = 1./mad( phs,1,1 ).';
            ac = mean(ac,1);
            phs = mean(phs,1);
            w_ac = ones(size(ac,2),1);
            w_phs = ones(size(phs,2),1);
        else
            w_ac = ones(size(ac,2),1);
            w_phs = ones(size(phs,2),1);
        end
        
        X = [d ones(size(d))];
        
        for j = 1:size(ac,1)
            y = log( ac(j,:).' .* d.^2 );
%          	alpha = robustfit(diag(w_ac)*X,w_ac.*y,[],[],'off');
            alpha = (diag(w_ac)*X) \ (w_ac.*y);
            alpha = alpha(1);


            y = -phs(j,:)';
%             phi = robustfit(diag(w_phs)*X,w_phs.*y,[],[],'off');
            phi = (diag(w_phs)*X) \ (w_phs.*y);
            phi = phi(1);

            thisMua = (data.Fm*2*pi * 1e6)/prop.v/2 * (phi/alpha - alpha/phi);
            thisMus = (alpha^2 - phi^2)/3/thisMua - thisMua;
            
            lastMua = Inf;
            
            iter = 0;
            while abs( (lastMua - thisMua)/thisMua ) > 1e-6 && iter < 100
                lastMua = thisMua;
                
                DD = 1/3/(thisMua + thisMus);
                x =  (data.Fm*2*pi*1e6)/prop.v/thisMua;
       
                Vp = sqrt(sqrt(1+x^2)+1);
                Vm = sqrt(sqrt(1+x^2)-1);
                
                G = 1 + d * sqrt(2*thisMua/DD) * Vp + d.^2 * thisMua/DD * sqrt(1+x^2);

                y = log( ac(j,:)' .* d.^3 ./ sqrt( G ) );
            
                alpha = (diag(w_ac)*X) \ (w_ac.*y);
%                 alpha = robustfit(diag(w_ac)*X,w_ac.*y,[],[],'off');
                alpha = alpha(1);
                
                % phi
                y = -phs(j,:)' + atan( d*sqrt(thisMua/2/DD)*Vm ./ ...
                    (1 + d*sqrt(thisMua/2/DD)*Vp) ...
                    );
                
                phi = (diag(w_phs)*X) \ (w_phs.*y);
%                 phi = robustfit(diag(w_phs)*X,w_phs.*y,[],[],'off');
                phi = phi(1);
                
                thisMua = (data.Fm*2*pi*1e6)/prop.v/2 * (phi/alpha - alpha/phi);
                thisMus = (alpha^2 - phi^2)/3/thisMua - thisMua;
                
                iter = iter + 1;
            end
            
            mua(i,j) = thisMua;
            mus(i,j) = thisMus;
        end
 
    end
    
    prop.mua = mua;
    prop.mus = mus;


end