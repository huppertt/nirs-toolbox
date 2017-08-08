function varargout = scale(this,pnorm,varargin)
%SCALE  Second-order section scaling.
%
%   See help dfilt/scale.

%   Copyright 2003-2009 The MathWorks, Inc.

if nargout>0,
    Hd = copy(this);
else
    Hd = this;
end

hfdesign = getfdesign(Hd);
hfmethod = getfmethod(Hd);

if nargin < 2,
    pnorm = 'Linf';
else
    if ~any(strcmp(pnorm,{'l1', 'L1', 'l2', 'L2', 'linf', 'Linf'}));
        error(message('signal:dfilt:abstractsos:scale:invalidnorm'));
    end
end

opts = uddpvparse('fdopts.sosscaling', {'scaleopts', Hd}, varargin{:});

% Cache current arithmetic setting
arith = Hd.Arithmetic;

% Set arithmetic temporarily to double
Hd.Arithmetic = 'double';

% Reorder prior to scaling
reorder(Hd,opts.sosReorder);

% Perform scaling
performscale(Hd,pnorm,opts);

% Set arithmetic back to what it was
Hd.Arithmetic = arith;

% Reflect norm and Opts changes in the fmethod object if it exists
if ~isempty(hfmethod)
    if strcmpi(opts.sosReorder,'auto')
        opts.sosReorder = getsosreorder(hfmethod);
    end
    set(hfmethod,'SOSScaleNorm',pnorm)
    set(hfmethod,'SOSScaleOpts',opts);
end

setfdesign(Hd, hfdesign);
setfmethod(Hd, hfmethod);

if nargout>0,
    varargout{1} = Hd;
end

%--------------------------------------------------------------------------
function performscale(Hd,pnorm,opts)

L = nsections(Hd);

% Attach pnorm to object
if ~isprop(opts,'pnorm'), adddynprop(opts,'pnorm','mxArray'); end
set(opts,'pnorm',pnorm);

% Attach max state value to object
if ~isprop(opts,'smax'), adddynprop(opts,'smax','mxArray'); end
set(opts,'smax',1);

% Compute unconstrained scaling
super_unconstrainedscale(Hd,opts,L);

% Compute offset when compensating for adjustment factors
if shiftsecondary(Hd),
    offset = 1;
else
    offset = 0;
end

% Constrain coeffs. Always do this before scale values!
constraincoeffs(Hd,opts,L,offset);

% Constrain scale values
constrainsv(Hd,opts,L,offset);

% Remove dynamic properties from opts
rmdynprop(opts,'pnorm','smax');

%--------------------------------------------------------------------------
function constrainsv(Hd,opts,L,offset)

if ~strcmpi(opts.ScaleValueConstraint,'unit'),
    % Extract sos matrix for readability
    sosm = Hd.sosMatrix;
    sv = Hd.ScaleValues;
    cmax = opts.MaxNumerator;
    svmax = opts.MaxScaleValue;

    for n = 1:L,
        if abs(sv(n)) > svmax,
            adjfact = svmax/sv(n);
            sv(n) = sign(sv(n))*svmax;
            % Try to incorporate into numerator when it hasn't been
            % constrained, but don't exceed cmax
            if (n-offset > 0) && strcmpi(opts.NumeratorConstraint,'none') && ...
                    (max(abs(sosm(n-offset,1:3)))/adjfact <= cmax),
                sosm(n-offset,1:3) = sosm(n-offset,1:3)/adjfact;
            else
                % Pass to next scale value
                sv(n+1) = sv(n+1)/adjfact;
            end
        end
    end

    if strcmpi(opts.ScaleValueConstraint,'po2'),
        % Adjust to powers of two one by one to minimize nonoptimallity
        newsv = zeros(size(sv(1:L)));
        
        % First scale value should always be 'floored' to ensure overflow
        % avoidance
        newsv(1) = 2.^(floor(log2(sv(1))));
        
        for n = 2:L,
            % Try to make next scale value larger to compensate 
            newsv(n) = 2.^(ceil(log2(sv(n))));
            
            % If the cumulative product exceeds original, cut back
            if prod(newsv(1:n)) > prod(sv(1:n)),
                newsv(n) = 2.^(floor(log2(sv(n))));
            end
        end
        
        % Compute adjustment factor to maintain overall gain
        adjfact = prod(newsv./sv(1:L));
        
        % Set new scale values
        sv(1:L) = newsv;
        
        % Compensate for adjustment
        sv(L+1) = sv(L+1)/adjfact;
    end

    Hd.sosMatrix = sosm;
    Hd.ScaleValues = sv;
end

%--------------------------------------------------------------------------
function constraincoeffs(Hd,opts,L,offset)

% Extract sos matrix for readability
sosm = Hd.sosMatrix;
sv = Hd.ScaleValues;

cmax = opts.MaxNumerator;

for n = 1:L,
    if strcmpi(opts.NumeratorConstraint,'unit'),
        % Get leading coefficient
        blead = sosm(n,1);

        % Find the adjustment factor
        adjfact = 1/blead;
        
        % Adjust other coefficients by adjustment factor
        sosm(n,1:3) = [1,sosm(n,2:3)*adjfact];

        % Saturate coefficients larger than cmax
        indx = find(abs(sosm(n,2:3)) > cmax); % Don't change leading coeff. It's a wire
        sosm(n,indx+1) = sign(sosm(n,indx+1)).*cmax;

        % Incorporate adjustment factor in scale value or next num
        [sosm,sv] = incorporate_adjfact(opts,adjfact,sosm,sv,n,L,offset);
    
    else
        % Determine maximum num coeff (in magnitude)
        nmax = max(abs(sosm(n,1:3)));
        
        % If we exceed Cmax, constrain. Always do so for 'normalize'
        if nmax > cmax || strcmpi(opts.NumeratorConstraint,'normalize'),
            adjfact = cmax/nmax;
            sosm(n,1:3) = sosm(n,1:3)*adjfact;
            [sosm,sv] = incorporate_adjfact(opts,adjfact,sosm,sv,n,L,offset);
        end
        
        % Enforce leading power of 2 if requested
        if strcmpi(opts.NumeratorConstraint,'po2'), 
            % Get leading coefficient
            blead = sosm(n,1);

            % Compute new leading num coeff
            newblead = 2.^(min(floor(log2(abs(blead))),floor(log2(cmax))));

            % Compute adjustment factor
            adjfact = newblead/blead;

            % Compensate for adjustment factor
            sosm(n,1:3) = [newblead,sosm(n,2:3)*adjfact];
            [sosm,sv] = incorporate_adjfact(opts,adjfact,sosm,sv,n,L,offset);
        end
    end
end
Hd.sosMatrix = sosm;
Hd.ScaleValues = sv;

%--------------------------------------------------------------------------
function [sosm,sv] = incorporate_adjfact(opts,adjfact,sosm,sv,n,L,offset)
% Incorporate the adjustment factor.

if strcmpi(opts.ScaleValueConstraint,'unit'),
    if n < L,
        % Try to incorporate in next section num
        sosm(n+1,1:3) = sosm(n+1,1:3)/adjfact;
    else
        % Incorporate in last scale value
        sv(n+1) = sv(n+1)/adjfact;
    end
else
    % For df2 and df1t, make sure to not increase norms
    % Remove this if statement if we always want to maintain optimal
    % scaling for the input to the states
    if adjfact > 1 && offset == 1,
        offset = 0; % Include adjustment in previous scale value
    end
    sv(n+offset) = sv(n+offset)/adjfact;
end


% [EOF]
