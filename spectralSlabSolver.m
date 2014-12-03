function tissue = spectralSlabSolver( meas, varargin )
%SPECTRALSLABSOLVER Summary of this function goes here
%   Detailed explanation goes here

    if nargin > 1;
        flag = logical( varargin{1} );
    else
        flag = logical( [1 1 0 0 0] ); %[HbO HbR Water Lipid CytC]
    end

%     prop = nirs.slabSolver( meas );
    prop = nirs.nlSlabSolver( meas );
    
%     if any( prop.mua < 0 )
%         c = [NaN NaN];
%         mus = NaN;
%     else
        ext = nirs.getSpectra( meas.probe.lambda );

        mua = prop.mua.' - 0.7 * repmat(ext(:,3),[1,size(prop.mua.',2)]);
        
        if size(mua,2) > 1
            w_mua = 1./std(prop.mua,1).';
        else
            w_mua = ones(size(mua,1),1);
        end
        
        ext = ext(:,flag);
%         ext(:,3) = [0.0159 0.0177 0.0164 0.0191]';

        c = (diag(w_mua)*ext) \ (diag(w_mua)*mua);
        
        mus = prop.mus;
%     end
    
    tissue = nirs.OpticalTissue( c(1,:).', c(2,:).', mus, prop.ri );

end

