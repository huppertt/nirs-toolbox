function out = linear_interp( pos, vol, seg, n )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    out = zeros(size(pos,1),1);
    
    for idx = 1:size(pos,1)
        
        % coordinates to interpolate
        x = pos(idx,1);
        y = pos(idx,2);
        z = pos(idx,3);

        % indices of the sub-volume to fit a quad surface to
        [X, Y, Z] = meshgrid( ...
            max(floor(x-n),1):min(ceil(x+n),size(vol,1)),...
            max(floor(y-n),1):min(ceil(y+n),size(vol,2)),...
            max(floor(z-n),1):min(ceil(z+n),size(vol,3)) );

        % sub-volume
        ii = sub2ind( size(vol), X, Y, Z );
        V = vol( ii );
        S = seg( ii );

        % mask out voxels that are zero
        mask = ( S > 0 & ~isinf(V) & ~isnan(V) );
        
        % weight matrix
        D =  (x-X).^2 + (y-Y).^2 + (z-Z).^2 ;
        W = exp( -D/4 );
        W = sqrt( W(mask) );


        % regression matrix
        A = [X(mask) Y(mask) Z(mask) ones(size(Z(mask)))];
        
        A = bsxfun(@times,W,A);
        V = W.*V(mask);

        % coefficients
        b = A \ V;


        % interpolation
        out(idx) = [x y z 1] * b;
        
    end


end

