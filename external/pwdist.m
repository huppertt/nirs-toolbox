function DST = pwdist( X,Y )

    if nargin == 1
        Y = X;
    end
    
    DST = sqrt(bsxfun(@plus,dot(X,X,2)',dot(Y,Y,2))-2*Y*X');
    
end

