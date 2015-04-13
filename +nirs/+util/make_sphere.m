function img = make_sphere( radii, dim  )
%MAKE_SPHERE Summary of this function goes here
%   Detailed explanation goes here

    assert( isvector(radii) && isscalar(dim) )
    
    R = sort( radii,'descend' ) / dim; % radii in grid units
    
    vol = uint8( zeros( (2*R(1) + 3)*[1 1 1] ) );
    
    c = ceil(size(vol)/2);
    [x,y,z] = meshgrid(...
        1:size(vol,1),...
        1:size(vol,2),...
        1:size(vol,3)...
        );

    r = sqrt( (x(:)-c(1)).^2 ...
        + (y(:)-c(2)).^2 ...
        + (z(:)-c(3)).^2 ...
        );

    for i = 1:length(R)
        vol( r <= R(i) ) = i;
    end 

    img = nirs.Image(vol,dim*[1 1 1],c);

end

