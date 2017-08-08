function [sch,rankJW] = getscheffeparam(varargin)
%GETSCHEFFEPARAM Compute the Scheffe parameter for a simultaneous confidence interval.
%   [SCH] = GETSCHEFFEPARAM('UnweightedJacobian',J) computes the Scheffe 
%   parameter for a simultaneous prediction interval. Suppose X is an N by R 
%   design matrix containing observations in rows and variables in columns. 
%   Then, J is a scalar multiple of an N by P Jacobian matrix of a nonlinear 
%   regression function of P variables at the N observations in X evaluated 
%   at the P by 1 least squares or robust regression solution vector BETA. 
%   The assumption is that a 'constant' error model will be used for both model 
%   fitting and making predictions.
%
%   [SCH] = GETSCHEFFEPARAM('param1','val1',...) accepts parameter name/value 
%   pairs to control the calculation of SCH. Parameters are:
%
%      'ErrorVarianceFit'   A N by 1 vector VF containing the estimated 
%                           variance of the error model EF used during model 
%                           fitting at the N observations in X. The form of 
%                           VF(i) is SIG2*G(i). Calculation of SCH depends 
%                           only on G and is insensitive to the value of SIG2. 
%                           The default value of VF is ones(N,1). If you supply
%                           VF, all its elements must be > 0. The default 
%                           assumption is that the error model used for making
%                           predictions is the same as the one used for model 
%                           fitting. For a regression with weight vector W,
%                           the sensible value for 'ErrorVarianceFit' is 1./W.
%
%     'ErrorVariancePred'   A N by 1 vector VP containing the variance of the 
%                           error model EP to be used for making predictions 
%                           at the N observations in X. The form of VP(i) is 
%                           SIG2_bar*G_bar(i). Calculation of SCH depends only 
%                           on G_bar and is insensitive to the value of SIG2_bar. 
%                           The default assumption is that the error model used 
%                           for making predictions is the same as the one used 
%                           for model fitting. Hence, the default value of VP 
%                           is VF. If you supply VP, all its elements must be > 0.
%
%      'WeightedJacobian'   Any scalar multiple of the N by P matrix JW which 
%                           is related to the unweighted Jacobian J via:
%                               JW(i,:) = J(i,:)/sqrt(VF(i)) 
%                           It is assumed that the co-variance matrix of the 
%                           estimated BETA is given by a multiple of the pseudo 
%                           inverse of (JW'*JW). This will be true when BETA has
%                           been computed using generalized least squares (such 
%                           as in nlinfit). 
%
%                'Intopt'   Either 'observation' (the default) for simultaneous 
%                           prediction intervals or 'curve' for simultaneous 
%                           confidence intervals on function values.
%
%                'TolSVD'   Uses the specified real positive scalar tolerance 
%                           for determining the rank of JW using SVD. The default 
%                           value is either eps(class(J)) or eps(class(JW)).
%
%                  'TolE'   Uses the specified real positive scalar tolerance for
%                           checking the condition VQ*VQ'*D = D. (see below). 
%                           The default value is sqrt(eps(class(J))) or 
%                           sqrt(eps(class(JW))).
%
%                    'VQ'   An N by Q matrix containing the left singular vectors
%                           of JW corresponding to non-zero singular values. 
%                           In this case an SVD of JW will not be computed and 
%                           rank of JW will be set equal to the number of columns 
%                           in VQ.
%
%   Suppose rank(JW) = Q and let JW = VQ * SQ * UQ' be the economy size SVD of 
%   JW where VQ is a N by Q matrix containing the Q left singular vectors of 
%   JW corresponding to non-zero singular values of JW. Also, suppose D is an 
%   N by 1 vector with D(i) = sqrt( VP(i)/VF(i) ). The Scheffe parameter is 
%   set as:
%             SCH = rank(JW)    , if VQ*VQ'*D = D or if Intopt = 'curve' and 
%                 = rank(JW) + 1, otherwise    
%
%   You must specify either J or JW. If you supply J instead of JW, be sure to 
%   supply the correct VF used during model fitting so that JW is properly 
%   calculated. 
%
%   [SCH,rankJW] = GETSCHEFFEPARAM(...) also returns the rank of JW.
%

%   References:
%    [1] Thomas P. Lane and William H. DuMouchel. Simultaneous confidence 
%        intervals in multiple regression. The American Statistician, 
%        48(4):315-321, 1994.

%   Copyright 2011-2013 The MathWorks, Inc.



    % Check number of input arguments.
    if ( nargin < 2 )
        error(message('stats:internal:getscheffeparam:TooFewInputs'));
    end
        
    % Allowed values of Intopt.
    IntoptObs = 'observation'; IntoptCurve = 'curve';
    allowedIntopt = {IntoptObs,IntoptCurve};
    
    % Provide default values for optional parameters.
    dfltJ  = []; dfltVF = []; dfltVP = []; dfltJW = []; 
    dfltIntopt = IntoptObs; dfltTolSVD = []; dfltTolE = []; dfltVQ = [];
         
    % Set defaults for optional parameters.
    names = { 'UnweightedJacobian', 'ErrorVarianceFit', 'ErrorVariancePred', 'WeightedJacobian',     'Intopt',       'TolSVD',       'TolE',      'VQ'  };
    dflts = {        dfltJ        ,       dfltVF      ,       dfltVP       ,       dfltJW      ,  dfltIntopt ,     dfltTolSVD,     dfltTolE,     dfltVQ };
    
    % Process optional param name/value pairs.
    [J,VF,VP,JW,Intopt,TolSVD,TolE,VQ] = internal.stats.parseArgs(names,dflts,varargin{:});

    % Validate parameters.
    [J,VF,VP,JW,Intopt,TolSVD,TolE,VQ,usingJ] = ValidateParameters(J,VF,VP,JW,Intopt,TolSVD,TolE,VQ,allowedIntopt);   
        
    %==== Compute JW if required and get VQ and rankJW ====
    if isempty(VQ)
        % VQ is not given.
        
        % Get weighted Jacobian based on error model used during fitting.
        if usingJ
           % J is given and we need JW(i,:) = J(i,:)/sqrt(VF(i)). This gives 
           % the true JW upto a scalar multiple.
           JW =  bsxfun(@(x,y) x.*y, J, 1./sqrt(VF)); 
        end
        % Compute SVD of JW = V*S*U'.
        [V,S,~] = svd(JW,'econ');
        % Check for zero singular values.
        Sd = abs(diag(S));
        idx = ( Sd <=  TolSVD * sqrt(sum(Sd.^2)) );
        % Get the rank of JW = length(find(~idx));
        rankJW = sum(~idx);    
        % Keep only columns of V corresponding to non-zero singular values.
        VQ = V(:,~idx);    
    else
        % VQ is given.
        rankJW = size(VQ,2);       
    end
    
    %==== If required, check if VQ*(VQ'*d) = d, set Scheffe parameter ====
    if strcmpi(Intopt,IntoptCurve)
       % For simultaneous CI for function value, just output the rank of JW.
       sch = rankJW;
    else
        % We would like to check if VQ*VQ'*d = d where d is an n by 1 vector
        % given by d(i) = sqrt( VP(i) ./ VF(i) ).
        d = sqrt(VP./VF);

        % Get dpred = VQ*(VQ'*d).
        dpred = VQ*(VQ'*d);

        % Compare dpred and d.    
        if ( max(abs(dpred - d)) <= TolE * norm(d) )
           % Can get away with rankJW.
           sch = rankJW;
        else
           % Must use rankJW + 1.
           sch = rankJW + 1;
        end
    end
    
end % End of getscheffeparam.
    
%==== Subfunction ValidateParameters ====
function [J,VF,VP,JW,Intopt,TolSVD,TolE,VQ,usingJ] = ValidateParameters(J,VF,VP,JW,Intopt,TolSVD,TolE,VQ,allowedIntopt)   
       
    % Prefer to use J and JW in that order.
    if ~isempty(J)
        % Use J.
        usingJ = true;
        % Validate J = n by p matrix.
        [n,p] = size(J);
        if ( ~isnumeric(J) || ~isreal(J) || ~ismatrix(J) )
            error(message('stats:internal:getscheffeparam:BadJ'));
        end
        % Store class of J.
        inputClass = class(J);
        % Make JW empty.
        JW = [];    
    elseif ~isempty(JW)
        % Use JW.
        usingJ = false;
        % Validate JW = n by p matrix.
        [n,p] = size(JW);
        if ( ~isnumeric(JW) || ~isreal(JW) || ~ismatrix(JW) )
            error(message('stats:internal:getscheffeparam:BadJW'));
        end
        % Store class of JW.
        inputClass = class(JW);
        % Make J empty.
        J = [];
    else
        % Both J and JW cannot be empty.
        error(message('stats:internal:getscheffeparam:SupplyJorJW'));    
    end
    
    % Validate VF after making it a column vector.   
    if isempty(VF)
        % Use the default value.
        VF = ones(n,1);
    else
        VF = VF(:);
        % Make sure VF has the right size.
        if ( size(VF,1) ~= n || ~isnumeric(VF) || ~isreal(VF) || ~ismatrix(VF) || any(VF <= 0) )
            error(message('stats:internal:getscheffeparam:BadVF'));
        end
    end
    
    % Validate VP after making it a column vector.
    if isempty(VP)
        % By default, we assume VP = VF.
        VP = VF;
    else
        VP = VP(:);
        % Make sure VP has the right size.
        if ( size(VP,1) ~= n || ~isnumeric(VP) || ~isreal(VP) || ~ismatrix(VP) || any(VP <= 0) )
            error(message('stats:internal:getscheffeparam:BadVP'));
        end
    end
    
    % Validate Intopt.     
    if ( isempty(Intopt) || ~ischar(Intopt) || ~any(strcmpi(Intopt,allowedIntopt)) )
        error(message('stats:internal:getscheffeparam:BadIntopt',allowedIntopt{1}, allowedIntopt{2}));
    end
           
    % Validate TolSVD.
    if isempty(TolSVD)
        % Set default tolerance.
        TolSVD = eps(inputClass);     
    else
        % Make sure TolSVD is sensible.
        if ( ~isnumeric(TolSVD) || ~isscalar(TolSVD) || ~isreal(TolSVD) || (TolSVD <= 0) )
            error(message('stats:internal:getscheffeparam:BadTolSVD'));
        end         
    end   

    % Validate TolE.
    if isempty(TolE)
        % Set default tolerance.
        TolE = sqrt(eps(inputClass));     
    else
        % Make sure TolE is sensible.
        if ( ~isnumeric(TolE) || ~isscalar(TolE) || ~isreal(TolE) || (TolE <= 0) )
            error(message('stats:internal:getscheffeparam:BadTolE'));
        end         
    end   
    
    % Validate VQ = n by q matrix.
    [nQ,pQ] = size(VQ); 
    if ~isempty(VQ) && ( nQ ~= n || pQ > p || ~isnumeric(VQ) || ~isreal(VQ) || ~ismatrix(VQ) )
        error(message('stats:internal:getscheffeparam:BadVQ'));
    end
    
    % Get rid of any NaNs in the input.
    if usingJ
       badRows = any(isnan(J),2);
    else
       badRows = any(isnan(JW),2);
    end
    badRows = ( badRows | any(isnan(VF),2) | any(isnan(VP),2) );    
    if usingJ
       J(badRows,:) = [];
       n = size(J,1);
    else
       JW(badRows,:) = [];
       n = size(JW,1);
    end
    VF(badRows) = [];
    VP(badRows) = [];
    if ~isempty(VQ)
        VQ(badRows,:) = [];
    end
    
    if ( n < 1 )
        error(message('stats:internal:getscheffeparam:BadN'));
    end
end % End of ValidateParameters.
