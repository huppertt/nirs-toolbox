function out = median_interp( pos, vol, seg, n )
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
        V = vol( sub2ind( size(vol), X, Y, Z ) );
        S = seg( sub2ind( size(seg), X, Y, Z ) );

        % mask out voxels that are zero
        mask = ( S > 0 & ~isinf(V) & ~isnan(V) );


%         % regression matrix
%         A = [X(mask) Y(mask) Z(mask) ones(size(Z(mask)))];
% 
%         % coefficients
%         b = pinv(A) * V(mask);


        % interpolation
%         out(idx) = [x y z 1] * b;
        out(idx) = median(V(mask));
        
    end


end

