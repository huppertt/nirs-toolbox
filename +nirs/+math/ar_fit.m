function [coef, res, yhat] = ar_fit( y,Pmax )
    % Copyright (c) 2014, Jeffrey W Barker (jwb52@pitt.edu)
    % All rights reserved.
    % 
    % Redistribution and use in source and binary forms, with or without
    % modification, are permitted provided that the following conditions
    % are met:
    % 
    % 1. Redistributions of source code must retain the above copyright
    % notice, this list of conditions and the following disclaimer.
    % 
    % 2. Redistributions in binary form must reproduce the above copyright
    % notice, this list of conditions and the following disclaimer in the
    % documentation and/or other materials provided with the distribution.
    % 
    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    % "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    % LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    % A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    % HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
    % INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
    % BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
    % OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
    % AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
    % LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
    % WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    % POSSIBILITY OF SUCH DAMAGE.
    
    % y(t) = coef(1) + coef(2)*y(t-1) + coef(3)*y(t-2) ... 
    %
    % y is a column vector with the time-series we are fitting
    %
    % Pmax is the max model order (not counting the constant)
    
    % we will minimize forward and backwards errors
    yf = y;
    yb = flipud(y);
    
    % the truncated portion of y we will fit
    yf = yf(1:end-Pmax,1);
    yb = yb(1:end-Pmax,1);
    
    % form the regression matrix from time shifted versions of y 
    % and a column of ones 
    Xf = convmtx(yf,Pmax+1);
    Xf = Xf(1:length(yf),2:end);
    Xf = [ones(length(yf),1) Xf];
    
    Xb = convmtx(yb,Pmax+1);
    Xb = Xb(1:length(yf),2:end);
    Xb = [ones(length(yf),1) Xb];
    
    yt = [yf; yb];
    X  = [Xf; Xb];
    
    % qr factorization will speed up stepwise regression significantly
    [Q,R] = qr(X,0); % zero is very important for performance
    invR = pinv(R);  % save this for performance reasons
    
    BIC = zeros(Pmax+1,1);
    for i = 1:length(BIC)
        % get residual for each fit for orders 0 to Pmax
        b = invR(1:i,1:i) * Q(:,1:i)'*yt;
        r = yt - X(:,1:i)*b;
        
        % calculate information criterion
        n = length(r);
        LL = -n/2*log( 2*pi*var(r) ) - n/2; % log-likelihood
        
        % uncomment the line you want to use for information criterion
        % I find that BIC produces more consistent results, but some
        % studies suggest AIC is superior
        BIC(i) = -LL + i/2*log(n);
        %AIC(i) = -LL + i + i*(i+1)/(n-i-1);
    end
    
    % uncomment the line you want to use for information criterion
	[~,Popt] = min( BIC ); % optimal model order
	%[~,Popt] = min( AIC ); % optimal model order

    Popt = Popt - 1;
    
    % finally, our output
    coef = invR(1:Popt+1,1:Popt+1) * Q(:,1:Popt+1)'*yt;
    res = filter([1; -coef(2:end)], 1, y-coef(1));
    
    yhat = y - res;
    
end