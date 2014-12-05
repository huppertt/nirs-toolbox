function prop = water( lambda )
%CSF Summary of this function goes here
%   Detailed explanation goes here

    ext = nirs.utilities.getSpectra( lambda );
    
    mua = ext(:,3);
    mus = 0.01*ones(size(lambda));
    
	prop = nirs.OpticalProperties( mua,mus,lambda,1.33 );
    
end

