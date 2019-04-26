function C = grangers( x, y, Pmax )

if(nargin==2)
    % assume that the fcn was grangers(X,Pmax) where X is a matrix
    n=size(x,2);
    C=zeros(n,n); 
    Pmax=y;
    for i=1:n
        for j=1:n
            if(i~=j)
                C(i,j)=nirs.math.grangers(x(:,i),x(:,j),Pmax);
            end
        end
    end
    return

end




    n = length(y);
    
    X = nirs.math.lagmatrix(x, 1:Pmax);
    Y = nirs.math.lagmatrix(y, 1:Pmax);
    
    [b1, res1] = nirs.math.stepwise([ones(n, 1) Y], y);
    [b2, res2] = nirs.math.stepwise( ...
        [ones(n, 1) Y(:, 1:length(b1)-1) X], y );

    C = max(log( mad(res1,0)^2 / mad(res2,0)^2 ), 0);
end

