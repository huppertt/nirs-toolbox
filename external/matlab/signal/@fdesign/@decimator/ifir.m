function varargout = ifir(this, varargin)
%IFIR   Design an interpolated FIR.

%   Author(s): J. Schickler
%   Copyright 2005-2014 The MathWorks, Inc.

% Make sure that IFIR is valid for this set of specifications.
if ~isdesignmethod(this, 'ifir')
    error(message('signal:fdesign:decimator:ifir:invalidMethod', 'IFIR', this.Specification));
end

% Parse the inputs for the multirate filter structure.
[filtstruct, varargin] = parsestruct(this,'firdecim', 'ifir', varargin{:});

filtstruct = ['mfilt.', filtstruct];

dfactor = get(this, 'DecimationFactor');

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

[Hd, upfactor] = ifir(this.CurrentFDesign, varargin{:}, 'rcf', dfactor);

fm = getfmethod(Hd);

if isa(Hd, 'dfilt.cascade'),
    % Integrate noble identities in Cascase
    [Hm,final_upfactor] = cascadenoble(Hd,upfactor,dfactor,filtstruct);
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
    [Hm1,final_upfactor] = cascadenoble(Hd.stage(s1),upfactor,dfactor,filtstruct);
    Hvec = decimdelaynoble(this,Hd.stage(s2),dfactor);
    Hvec = decimmergedelays(this,Hvec(end:-1:1));
    Hm2 = cascade(Hvec);
    Hm = parallel(Hm1,Hm2);
else
    % If we are returned a single section, the design was simple enough
    % that we can just use a single interpolator.
    Hm = feval(filtstruct, dfactor, Hd.Numerator);
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
function [Hm,final_upfactor] = cascadenoble(Hd,upfactor,dfactor,filtstruct)
% Integrate noble identities in Cascase

%  H1(z) -- H2(z^M) -- \/ MK
%
%  H1(z) -- \/ M -- H2(z) --- \/ K
%
%  H1_decim_M(z) -- H2_decim_K(z)
%
%  If MK is greater than the DecimationFactor requested.
%
%  H1(z) -- H2(z^MK) -- \/ M
%
%  H1(z) -- \/ M -- H2(z^K)
%
%  H1_decim_M(z) -- H2(z^k)
%
%
%  K is often 1, in which case we dont need a multirate.

% Divide the interpolation factor into two parts.
M = gcd(upfactor, dfactor);
K = dfactor/M;
final_upfactor = upfactor/M;

% Remove the extra ones from the numerator.  We will get these
% automatically when we use the interpolating structure.
b_h2 = Hd.Stage(2).Numerator(1:M:end);

b_h1 = get(Hd.Stage(1), 'Numerator');

% Create the 1st stage.
if M == 1
    if strcmpi(filtstruct, 'mfilt.firtdecim')
        Hm1 = dfilt.dffirt(b_h1);
    else
        Hm1 = Hd.Stage(1);
    end
else
    Hm1 = feval(filtstruct, M, b_h1);
end

% If there is "left over" decimation put it into the second stage.
if K == 1 %|| dfactor < upfactor
    if strcmpi(filtstruct, 'mfilt.firtdecim')
        Hm2 = dfilt.dffirt(b_h2);
    else
        Hm2 = Hd.Stage(2);
        set(Hm2, 'Numerator', b_h2);
    end
else
    Hm2 = feval(filtstruct, K, b_h2);
end

Hm = mfilt.cascade(Hm1, Hm2);


% [EOF]
