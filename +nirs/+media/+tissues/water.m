function prop = water( lambda )
%CSF Summary of this function goes here
%   Detailed explanation goes here

    ext = nirs.media.getspectra( lambda );
    
    mua = ext(:,3);
    mus = 0.01*ones(size(lambda));
    
	prop = nirs.media.OpticalProp( mua,mus,lambda,1.33 );
    
end

