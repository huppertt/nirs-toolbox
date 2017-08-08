function varargout = ifir(this, varargin)
%IFIR   Design an interpolated FIR filter.

%   Author(s): J. Schickler
%   Copyright 2005-2014 The MathWorks, Inc.

% Make sure that IFIR is valid for this set of specifications.
if ~isdesignmethod(this, 'ifir')
    error(message('signal:fdesign:interpolator:ifir:invalidMethod', 'IFIR', this.SpecificationType));
end

% Parse the inputs for the multirate filter structure.
[filtstruct, varargin] = parsestruct(this, 'firinterp', 'ifir', varargin{:});

filtstruct = ['mfilt.', filtstruct];

ifactor = get(this, 'InterpolationFactor');

% If System object is requested, first design mfilt then later convert to
% System object if possible
sysObjFlag = false;
sysObjIdx = find(cellfun(@(x)strcmpi('SystemObject',x),varargin),1);
if ~isempty(sysObjIdx)
    if varargin{sysObjIdx+1} == 1
        sysObjFlag = true;
        varargin(sysObjIdx:sysObjIdx+1) = [];  % remove the property so that result is mfilt at first
    end
end

[Hd, upfactor] = ifir(this.CurrentFDesign, varargin{:}, 'rcf', ifactor);

fm = getfmethod(Hd);

if isa(Hd, 'dfilt.cascade')
    % Integrate noble identities in Cascase
    [Hm,final_upfactor] = cascadenoble(Hd,upfactor,ifactor,filtstruct);
    if sysObjFlag
        Hm = sysobj(Hm);  % return dsp.FilterCascade
    end
elseif isa(Hd, 'dfilt.parallel'),
    if sysObjFlag
        % Error out because parallel structure is not supported by System objects
        error(message('signal:fdesign:decimator:ifir:ParallelNotSupportedBySystemObjects'));
    end
    % Find cascade branch within parallel
    if isa(Hd.stage(1),'dfilt.cascade'),
        s1 = 1;
        s2 = 2;
    else
        s1 = 2;
        s2 = 1;
    end
    [Hm1,final_upfactor] = cascadenoble(Hd.stage(s1),upfactor,ifactor,filtstruct);
    Hvec = interpdelaynoble(this,Hd.stage(s2),ifactor);
    Hvec = interpmergedelays(this,Hvec);
    Hm2 = cascade(Hvec);
    Hm = parallel(Hm1,Hm2);
else

    % If we are returned a single section, the design was simple enough
    % that we can just use a single interpolator.
    Hm = feval(filtstruct, ifactor, ifactor*Hd.Numerator);
    final_upfactor = 1;
    
    if sysObjFlag
        Hm = sysobj(Hm);  % return dsp.FilterCascade
    end
end

if isa(Hm, 'dsp.private.FilterAnalysis') % if returning a System object
    setMetaData(Hm,this,fm);
else
    Hm.setfdesign(this);
    Hm.setfmethod(fm);
end

if nargout
    varargout = {Hm, final_upfactor};
else
    if this.NormalizedFrequency,
        inputs = {'NormalizedFrequency', 'On'};
    else
        inputs = {'Fs', this.Fs};
    end

    fvtool(Hm, inputs{:});
end
%--------------------------------------------------------------------------
function [Hm,final_upfactor] = cascadenoble(Hd,upfactor,ifactor,filtstruct)

%  /\ LK -- H1(z) -- H2(z^L)
%
%  /\ LK -- H2(z^L) -- H1(z)
%
%  /\ K -- H2(z) -- /\L -- H1(z)
%
%  H2_interp_K(z) -- H1_interp_L(z)
%
%  If Upfactor is greater than InterpolationFactor
%
%  /\ L -- H1(z) -- H2(z^LK)
%
%  /\ L -- H2(z^LK) -- H1(z)
%
%  H2(z^K) -- /\L -- H1(z)
%
%  H2(z^K) -- H1_interp_L(z)
%
%  K is often 1, in which case we dont need a multirate.

% Divide the interpolation factor into two parts.
L = gcd(upfactor, ifactor);
K = ifactor/L;  % Usually 1
final_upfactor = upfactor/L;

% Remove the extra ones from the numerator.  We will get these
% automatically when we use the interpolating structure.
b_h2 = Hd.Stage(2).Numerator(1:L:end);

% If there is "left over" interpolation put it into the first stage.
if K == 1
    Hm1 = Hd.Stage(2);
    set(Hm1, 'Numerator', b_h2);
else
    Hm1 = feval(filtstruct, K, K*b_h2);
end

% Create the 2nd stage.
if L == 1
    Hm2 = Hd.Stage(1);
else
    Hm2 = feval(filtstruct, L, L*Hd.Stage(1).Numerator);
end

Hm = mfilt.cascade(Hm1, Hm2);


% [EOF]
