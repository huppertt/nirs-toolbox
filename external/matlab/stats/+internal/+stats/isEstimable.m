function [F,Q,pXTX,FUN,condXTX,VQ] = isEstimable(C,varargin)
%ISESTIMABLE Check contrast estimability.
%   F = ISESTIMABLE(C,'DesignMatrix',X) checks the contrasts in rows of M 
%   by P matrix C for estimability when using a linear regression model with 
%   the N by P design matrix X. The output F is a M by 1 logical vector with 
%   true in position K if row K of C contains an estimable contrast.
%
%   F = ISESTIMABLE(C,'param1','val1',...) accepts parameter name/value 
%   pairs to control the calculation of F. Parameters are:
%
%      'NormalMatrix'       Uses the P by P normal matrix Y. Y is related to 
%                           the N by P design matrix X via Y = X'*X. 
%
%      'pinvNormalMatrix'   Uses the P by P pseudo inverse Z of the normal 
%                           matrix Y. Z is related to the N by P design matrix 
%                           X via Z = pinv(Y) = pinv(X'*X). 
%
%      'TolSVD'             Uses the real positive scalar tolerance TolSVD 
%                           for checking 0 singular values instead of the 
%                           default value eps(class(T)) where T is either X 
%                           or Y or Z.
%
%      'TolC'               Uses the real positive scalar tolerance tolC for 
%                           checking contrast estimability instead of the 
%                           default value eps(class(C))^(1/4).
%
%   [F,Q] = ISESTIMABLE(C,...) also returns the rank of = X'*X in Q.
%
%   [F,Q,pXTX] = ISESTIMABLE(C,...) also returns the pseudo inverse of Y = X'*X.
%
%   [F,Q,pXTX,FUN] = ISESTIMABLE(C,...) also returns a function handle FUN 
%   that accepts a M by P contrast matrix C and returns a M by 1 logical vector
%   F with true in position K if row K of C is estimable using the currently
%   supplied Design matrix, Normal matrix or pseudo inverse of Normal matrix.
%
%   [F,Q,pXTX,FUN,condXTX] = ISESTIMABLE(C,...) also returns the condition 
%   number of X'*X in condXTX.
%
%   [F,Q,pXTX,FUN,condXTX,VQ] = ISESTIMABLE(C,...) also returns the left singular 
%   vectors of X from its SVD corresponding to non-zero singular values when 
%   using the option 'DesignMatrix' and [] when using 'NormalMatrix' or 
%   'pinvNormalMatrix'.
%
%   You must supply either the 'DesignMatrix', 'NormalMatrix' or 
%   'pinvNormalMatrix'. 
%
%   Example:
%      % Create a 6 by 3 design matrix.
%      X = [ones(6,1),[ones(3,1);zeros(3,1)],[zeros(3,1);ones(3,1)]]
%      % Create a 2 by 3 contrast matrix.
%      C = [1 1 1;2 1 1]
%      % Check if the contrasts in rows of C are estimable.
%      [F,Q,pXTX,FUN,condXTX] = isEstimable(C,'DesignMatrix',X,'TolSVD',eps^(3/4))
%      % For fixed X, you can reuse FUN to test different contrasts.
%      FUN([0,1,-1;[0,1,1]])
%

%   Copyright 2012 The MathWorks, Inc.


    % Check number of input arguments.
    if ( nargin < 3 )
        error(message('stats:internal:isEstimable:TooFewInputs'));
    end

    % Provide default values for optional parameters.
    dfltX = []; dfltY = []; dfltZ = []; dfltTolSVD = []; dfltTolC = [];
    
    % Set defaults for optional parameters.
    names = {'DesignMatrix', 'NormalMatrix', 'pinvNormalMatrix',     'TolSVD',     'TolC'};
    dflts = {    dfltX     ,      dfltY    ,       dfltZ       ,   dfltTolSVD,   dfltTolC};
    
    % Process optional param name/value pairs.
    [X,Y,Z,TolSVD,TolC]= internal.stats.parseArgs(names,dflts,varargin{:});

    % Define indicator strings for what has been supplied.
    usingX = 'X'; % DesignMatrix given.
    usingY = 'Y'; % NormalMatrix given.
    usingZ = 'Z'; % Pseudo inverse of NormalMatrix given.
    
    % Validate parameters.
    [C,X,Y,Z,TolSVD,TolC,using] = ValidateParameters(C,X,Y,Z,TolSVD,TolC,usingX,usingY,usingZ);
    
    % Get size of C.
    [m,p] = size(C);
    
    % Do either of the following:
    % (1) If X or Y is given, compute SVD: [U,S,~] = svd(X'*X) where U and S 
    %     are p by p matrices such that U*S*U' = Y = X'*X. 
    % (2) If Z = pinv(Y) is given then compute SVD of Z and invert non-zero 
    %     singular values to get SVD of Y.
    %
    % Note on check for zero singular values. Even if true singular value
    % is zero, when using Sd = svd(T), floating point arithmetic can introduce 
    % an error of around: eps * norm(T,'fro') = eps * sqrt(sum(Sd.^2)) . Use 
    % TolSVD instead of eps to account for intermediate roundoff errors: 
    % idx = ( Sd <=  TolSVD * sqrt(sum(Sd.^2)) ).
    switch using
        case usingX
            % If X is given, compute SVD of X for better accuracy.            
            n = size(X,1);           
            % U is p by p, Sx is p by p, V is n by p, for n > p 
            % U is p by p, Sx is n by p, V is n by n, for n <= p
            [V,Sx,U] = svd(X,0); 
            % Get singular values.
            if ( n > p )
                % Sx is p by p (n > p) and so diag(Sx) will be p by 1.
                if ( p == 1 )
                    Sd = abs(Sx(1,1));
                else
                    Sd = abs(diag(Sx));
                end
            else
                % Sx is n by p (n <= p) and so diag(Sx) will be n by 1.
                Sd = zeros(p,1);
                if ( n == 1 )
                    Sd(1:n) = abs(Sx(1,1));
                else
                    Sd(1:n) = abs(diag(Sx));
                end
            end
            % Square Sd, since we want SVD of X'*X.
            Sd = Sd.^2;   
            % Check for zero singular values.
            idx = ( Sd <=  2 * TolSVD * sqrt(sum(Sd.^2)) );                    
        case usingY
            [U,S,~] = svd(Y);
            % Get singular values.
            Sd = abs(diag(S));
            % Check for zero singular values.
            idx = ( Sd <=  TolSVD * sqrt(sum(Sd.^2)) );
        case usingZ
            [U,S,~] = svd(Z);
            Sd = abs(diag(S));
            % Invert non-zero singular values.
            idx = ( Sd <=  TolSVD * sqrt(sum(Sd.^2)) );
            Sd(~idx) = 1./Sd(~idx);
            % Check for zero singular values.
            idx = ( Sd <=  TolSVD * sqrt(sum(Sd.^2)) );
        otherwise
            error(message('stats:internal:isEstimable:XYZEmpty'));    
    end
                
    if any(idx)
        % Some singular values are zero.
        
        % C*U is of size m by p.
        CU = abs(C*U); 
        
        % Any contrast that satisfies C*U(:,idx) = 0 is estimable. Even if
        % the true value of max(CU(:,idx),[],2) is 0, floating point arithmetic 
        % can introduce an error as large as eps * sqrt(sum(abs(C).^2,2)).
        % Use TolC instead of eps to account for intermediate roundoff
        % errors.
        F = ( max(CU(:,idx),[],2) <= TolC * sqrt(sum(abs(C).^2,2)) );
    else
        % All singular values are non-zero. 
        
        % All contrasts are estimable.
        F = true(m,1);
    end
    
    if ( nargout > 1 )
        % Output the rank of X'*X.
        Q = length(find(~idx));
    end
    
    if ( nargout > 2 )
        % Output the pseudo inverse.
        Un = U(:,~idx);
        pXTX = bsxfun(@(x,y) x.*y, Un, (1./Sd(~idx))') * Un';
    end
     
    if ( nargout > 3 )
        % Create convenience function FUN. We could do this:
        %         if any(idx)
        %             FUN = @(Cnew) max(abs(Cnew*U(:,idx)),[],2) <= TolC * sqrt(sum(abs(Cnew).^2,2));
        %         else
        %             % All contrasts are estimable.
        %             FUN = @(Cnew) true(size(Cnew,1),1);
        %         end
        % Instead, we will call the subfunction makeFUN to reduce the
        % number of variables stored in FUN.
        FUN = makeFUN(idx,U(:,idx),TolC);
    end
    
    if ( nargout > 4 )
        % Compute condition number of X'*X.
        condXTX = max(Sd)/min(Sd);
    end
    
    if ( nargout > 5 )
        if strcmpi(using,usingX)
            % Get left singular vectors for non-zero singular values.
            VQ = V(:,~idx);
        else
            VQ = [];
        end
    end
    
end % End of isEstimable.

%==== Subfunction makeFUN ====
function FUN = makeFUN(idx,Uidx,TolC)
% Create convenience function FUN. This is a subfunction just to 
% create the function handle FUN. This is done to avoid storing
% unnecessary information from the workspace of isEstimable into
% the function handle.

        if any(idx)
            FUN = @(Cnew) max(abs(Cnew*Uidx),[],2) <= TolC * sqrt(sum(abs(Cnew).^2,2));
        else
            % All contrasts are estimable.
            FUN = @(Cnew) true(size(Cnew,1),1);
        end

end % End of makeFUN.

%==== Subfunction ValidateParameters ====
function [C,X,Y,Z,TolSVD,TolC,using] = ValidateParameters(C,X,Y,Z,TolSVD,TolC,usingX,usingY,usingZ)    

    % Make sure all of X, Y and Z are not empty.
    if ( isempty(X) && isempty(Y) && isempty(Z) )
        error(message('stats:internal:isEstimable:XYZEmpty'));
    end    
  
    % Prefer to use X, Y and Z in order.
    if ~isempty(X)
        % DesignMatrix given.
        using = usingX;
    elseif ~isempty(Y)
        % NormalMatrix given.
        using = usingY;
    else
        % Pseudo inverse of NormalMatrix given.
        using = usingZ;
    end
    
    % Make sure C is 2-D numeric matrix.
    [~, p] = size(C);
    if ( ~isnumeric(C) || ~ismatrix(C) )
        error(message('stats:internal:isEstimable:BadC'));
    end
    
    % Make sure we have sensible X, Y, Z.
    switch using        
        case usingX
            % X is given. Make sure X is a 2-D numeric matrix.
            if ( ~isnumeric(X) || ~ismatrix(X) )
                error(message('stats:internal:isEstimable:BadX'));
            end
            [~,p2] = size(X);       
        case usingY
            % Y is given. Make sure Y is a 2-D numeric matrix.
            if ( ~isnumeric(Y) || ~ismatrix(Y) )
                error(message('stats:internal:isEstimable:BadY'));
            end
            [~,p2] = size(Y);
        case usingZ
            % Z is given. Make sure Z is a 2-D numeric matrix.
            if ( ~isnumeric(Z) || ~ismatrix(Z) )
                error(message('stats:internal:isEstimable:BadZ'));
            end
            [~,p2] = size(Z);         
    end
                
    % Make sure C and (X or Y or Z) are of the right size.
    if ( p2 ~= p || p2 < 1 || p < 1 )        
        error(message('stats:internal:isEstimable:SizeMismatch'));
    end        
    
    % If TolSVD is not empty (user supplied), make sure it is sensible.
    if ~isempty(TolSVD)
       if ( ~isnumeric(TolSVD) || ~isscalar(TolSVD) || ~isreal(TolSVD) || (TolSVD < 0) )
            error(message('stats:internal:isEstimable:BadTolSVD'));
       end
    else
       % Set default TolSVD. 
       switch using        
           case usingX
               TolSVD = eps(class(X));
           case usingY
               TolSVD = eps(class(Y));
           case usingZ
               TolSVD = eps(class(Z));
       end 
    end
    
    % If TolC is not empty (user supplied), make sure it is sensible.
    if ~isempty(TolC)
       if ( ~isnumeric(TolC) || ~isscalar(TolC) || ~isreal(TolC) || (TolC < 0) )
            error(message('stats:internal:isEstimable:BadTolC'));
       end
    else
       % Set default TolC. 
       TolC = eps(class(C))^(1/4);
    end
    
    % If using a Normal matrix or its pseudo inverse, make sure the matrices 
    % are square and symmetric.
    badY = false; badZ = false;
    switch using
        case usingY
             emptyY = isempty(Y);
            squareY = ~emptyY && ( size(Y,1) == size(Y,2) );
              symmY = ~emptyY && squareY && ( max(max(abs(Y - Y'))) <= sqrt(eps(class(Y))) * max(max(abs(Y)+abs(Y'))) );   
               badY = ~(emptyY || symmY);
        case usingZ
             emptyZ = isempty(Z);
            squareZ = ~emptyZ && ( size(Z,1) == size(Z,2) );
              symmZ = ~emptyZ && squareZ && ( max(max(abs(Z - Z'))) <= sqrt(eps(class(Z))) * max(max(abs(Z)+abs(Z'))) );
               badZ = ~(emptyZ || symmZ);
    end
    if ( badY || badZ )
        error(message('stats:internal:isEstimable:MatrixNotSquareOrSymmetric'));
    end    
    
end % End of ValidateParameters.

    
