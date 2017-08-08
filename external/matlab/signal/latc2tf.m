function [num,den]=latc2tf(K,V)
%LATC2TF Lattice filter to transfer function conversion.
%   [NUM,DEN] = LATC2TF(K,V) finds the transfer function numerator
%   NUM and denominator DEN from the IIR lattice coefficients K and
%   ladder coefficients V.
%
%   [NUM,DEN] = LATC2TF(K,'allpole') assumes that K is associated with an
%   all-pole IIR lattice filter.
%
%   [NUM,DEN] = LATC2TF(K,'allpass') assumes that K is associated with an
%   allpass IIR lattice filter.
%
%   NUM = LATC2TF(K) and NUM = LATC2TF(K,'fir') assumes that K is
%   associated with an FIR lattice filter structure and that the upper
%   output of the structure is being used (the one corresponding to the
%   first output of LATCFILT for the FIR case).
%
%   NUM = LATC2TF(K,'min'), where abs(K) <= 1, assumes that K is
%   associated with a minimum-phase FIR lattice filter structure.
%
%   NUM = LATC2TF(K,'max'), where abs(K) <= 1, assumes that K is
%   associated with a maximum-phase FIR lattice filter structure and
%   that the lower output of the structure is being used (the one
%   corresponding to the second output of LATCFILT for the FIR case).
%
%   % Example:
%   %   Generate a max-phase lattice filter.
%
%   k = [1/6 1/1.4];                % Lattice coefficients
%   bmax = latc2tf(k,'max');        % Convert to transfer function form
%   max_flag = ismaxphase(bmax)     % Determine whether filter is max phase
%
%   See also LATCFILT and TF2LATC.

% Reference:[1] J.G. Proakis, D.G. Manolakis, Digital Signal Processing,
%            3rd ed., Prentice Hall, N.J., 1996, Chapter 7.
%           [2] S. K. Mitra, Digital Signal Processing, A Computer
%           Based Approach, McGraw-Hill, N.Y., 1998, Chapter 6.
%

%   Copyright 1988-2013 The MathWorks, Inc.

narginchk(1,2);

% Handle empty cases immediately:
% if isempty(K) | ...
%         ( (nargin>1) && ~ischar(varargin{1}) && isempty(varargin{1}) ),
%     handle_emptycase(K,varargin{:});
% end

% Parse input args:
switch nargin,
    case 1,
        % FIR
        [num,den] = latc2fir(K,'fir');
    case 2,
        
        if ischar(V),
            
            switch(lower(V))
                case 'allpass',
                    [num,den] = latc2allpass(K);
                    
                case {'allpole','iir'}, % 'iir' is allowed for backwards compatibility
                    [num,den] = latc2iir(K,1);
                    
                case 'fir',
                    [num,den] = latc2fir(K,'fir');
                    
                case 'min',
                    [num,den] = latc2fir(K,'min');
                    
                case 'max'
                    [num,den] = latc2fir(K,'max');
                    
                otherwise
                    error(message('signal:latc2tf:InvalidEnum'));
            end
            
        else
            % Cast to enforce Precision Rules
            if any([signal.internal.sigcheckfloattype(K,'single','latc2tf','K')...
                signal.internal.sigcheckfloattype(V,'single','latc2tf','V')])
              K = single(K);
              V = single(V);
            end
            [num,den] = latc2iir(K,V);    % IIR
        end
end
%-------------------------------------------------------------------------
function [num,den] = latc2fir(K,ph)
% Convert min. phase or max. phase lattice FIR to TF

if isempty(K),
    num = 1;
    den = 1;
else
    num = rc2poly(K);
    % Cast to enforce Precision Rules
    if isa(K,'single')
      den = ones(1,1,'single');
    else
      den = ones(1,1);
    end        
    
    if ~strcmpi(ph,'fir') && max(abs(K)) > 1,
        warning(message('signal:latc2tf:InvalidParam'));
    end
    
    % If max. phase, reverse num
    if strcmpi(ph,'max'),
        num = conj(num(end:-1:1));
    end
end

%-------------------------------------------------------------------------
function [num,den] = latc2allpass(K)
% Convert allpass lattice to TF

if isempty(K),
    num = 1;
    den = 1;
else
    
    den = rc2poly(K);
    num = conj(den(end:-1:1));
end

%-------------------------------------------------------------------------
function [num,den] = latc2iir(K,V)
% Convert allpole or arma (general iir) lattice to TF

if isempty(K) && isempty(V),
    num = [];
    den = [];
elseif isempty(K) && length(V) == 1,
    num = V;
    den = 1;
else
    % Solve for IIR lattice or lattice-ladder coefficients:
    K=K(:); V=V(:);
    
    % Make sure V is length(K)+1:
    ordiff = length(V)-length(K)-1;
    if ordiff>0,
        K = [K; zeros(ordiff,1)];
        % error('signal:latc2tf:InvalidDimensions','length(V) must be <= 1+length(K).');
    elseif ordiff<0,
        V = [V; zeros(-ordiff,1)];
    end
    
    % We still use rc2poly to compute the den
    den = rc2poly(K);
    
    % To compute the num coefficients we solve the equations (see [2] pp. 384):
    % num(end)   = V(1)
    % num(end-1) = V(2) + conj(den(1))*V(1)
    % num(end-2) = V(3) + conj(den(1)#)*V(2) + conj(den(2))*V(1)
    % num(end-3) = V(4) + conj(den(1)##)*V(3) + conj(den(2)#)*V(2) + conj(den(3))*V(1)
    % etc.
    % where den(m)# denotes the mth coefficient of the reduced (using levdown)
    % order polynomial; den(m)## denotes the mth coefficient of the 2 step reduced (using levdown twice)
    % order polynomial; etc.
    % Note that these equations are the same used for finding V in tf2latc, except
    % we are solving for den instead of for V. In the present case, no recursive
    % solution is needed.
    
    % We will use a matrix with the denominators of lower orders
    % in each column, this matrix is the same as the one used in tf2latc
    [r,tempmatrix] = rlevinson(den,1); %#ok
    num = tempmatrix*V;
    num = num.'; % it is a polynomial, make it a row
end


% EOF

