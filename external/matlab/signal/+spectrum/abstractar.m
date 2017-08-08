classdef (CaseInsensitiveProperties=true) abstractar < spectrum.abstractspectrum 
  %spectrum.abstractar class
  %   spectrum.abstractar extends spectrum.abstractspectrum.
  %
  %    spectrum.abstractar properties:
  %       EstimationMethod - Property is of type 'string' (read only)
  %       Order - Property is of type 'udouble user-defined'
  %
  %    spectrum.abstractar methods:
  %       ar_thisloadobj -   Load this object.
  %       confinterval -  Confidence Interval for AR spectrum estimation methods.
  %       initialize -   Set common AR properties.
  %       thisloadobj -   Load this object.
  
  
  properties (AbortSet, SetObservable, GetObservable)
    %ORDER Property is of type 'udouble user-defined'
    Order = 4;
  end
  
  
  methods
    function set.Order(obj,value)
      % User-defined DataType = 'udouble user-defined'
      
      validateattributes(value,{'numeric'},{'real','finite','integer','>',0})
      value = double(value);
      
      obj.Order = value;
    end
    
    function ar_thisloadobj(this, s)
      %AR_THISLOADOBJ   Load this object.
      
      this.Order = s.Order;
      
    end
    
    
    function CI = confinterval(this,x,Pxx,~,CL,~)
      %CONFINTERVAL  Confidence Interval for AR spectrum estimation methods.
      %   CI = CONFINTERVAL(THIS,X,PXX,W,CL,FS) calculates the confidence
      %   interval CI for spectrum estimate PXX based on confidence level CL. THIS is a
      %   spectrum object and W is the frequency vector. X is the data used for
      %   computing the spectrum estimate PXX.
      %
      %   Reference : Steven M.Kay, "Modern spectral Estimation",
      %   Prentice Hall, 1988, Chapter 6, pp 194-195
      
      alfa = 1-CL;
      normval = norminverse(1-alfa/2,0,1);
      
      N = length(x);
      p = this.Order;
      
      if( N/(2*p) > normval^2)
        beta = sqrt(2*p/N)* normval;
        CI = Pxx*[1-beta 1+beta];
      else
        warning(message('signal:spectrum:abstractar:confinterval:InsufficientData'));
        CI = [];
      end
      
      
    end
    
    function initialize(this, order)
      %INITIALIZE   Set common AR properties.
      %    If ORDER is not specified use default values.
      
      if nargin < 2,
        order = 4;
      end
      this.Order = order;
      
      
    end
    
    function thisloadobj(this, s)
      %THISLOADOBJ   Load this object.
      
      ar_thisloadobj(this, s);
      
    end
    
  end  %% public methods
  
end  % classdef

function val = setorder(this,val)
% Allow integers only but allow 0.

flr_val = floor(val);
if ((val~=0) & flr_val==0) | ( (val~=0) & rem(val,flr_val) ),
  error(message('signal:spectrum:abstractar:schema:MustBeInteger'));
end
end  % setorder

%--------------------------------------------------------------------------
function [x] = norminverse(p,mu,sigma)
%NORMINVERSE Inverse of the normal cumulative distribution function (cdf).
%   X = NORMINVERSE(P,MU,SIGMA) returns the inverse cdf for the normal
%   distribution with mean MU and standard deviation SIGMA, evaluated at
%   the value in P.
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%
%   References:
%      [1] Abramowitz, M. and Stegun, I.A. (1964) Handbook of Mathematical
%          Functions, Dover, New York, 1046pp., sections 7.1, 26.2.
%      [2] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley, 170pp.

if nargin<1
  error(message('signal:spectrum:abstractar:confinterval:Nargchk'));
end

if nargin < 2
  mu = 0;
end

if nargin < 3
  sigma = 1;
end

if(sigma <=0)
  error(message('signal:spectrum:abstractar:confinterval:InvalidValue', 'sigma'));
end

if(p < 0 || 1 < p)
  error(message('signal:spectrum:abstractar:confinterval:InvalidValue', 'P'));
end

x0 = -sqrt(2).*erfcinv(2*p);
x = sigma*x0 + mu;

end



% [EOF]
