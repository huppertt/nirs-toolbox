function varargout = ifir(this, varargin)
%IFIR   Design an interpolated FIR filter.

%   Author(s): J. Schickler
%   Copyright 2005-2014 The MathWorks, Inc.

% Make sure that IFIR is valid for this set of specifications.
if ~isdesignmethod(this, 'ifir')
    error(message('signal:fdesign:rsrc:ifir:invalidMethod', 'IFIR', this.SpecificationType));
end

validatercf(this);

% Parse the inputs for the multirate filter structure.
[filtstruct, varargin] = parsestruct(this, 'firsrc', 'ifir', varargin{:});

filtstruct = ['mfilt.', filtstruct];

dfactor = get(this, 'DecimationFactor');
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

[Hd, upfactor] = ifir(this.CurrentFDesign, varargin{:}, ...
    'rcf', max(ifactor, dfactor));

fm = getfmethod(Hd);

if isa(Hd, 'dfilt.cascade')
    [Hm,final_upfactor] = cascadenoble(Hd,upfactor,ifactor,dfactor,filtstruct);
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
    [Hm1,final_upfactor] = cascadenoble(Hd.stage(s1),upfactor,ifactor,dfactor,filtstruct);
    if ifactor > dfactor,
        Hvec = interpdelaynoble(this,Hd.stage(s2),ifactor);
        Hvec = interpmergedelays(this,Hvec);
        % Merge downsampler
        if strcmpi(class(Hvec(end)),'mfilt.firinterp'),
            Hvec(end) = mfilt.firsrc(Hvec(end).InterpolationFactor,dfactor,...
                Hvec(end).Numerator);
        else
            % Must be a delay
            Hvec(end) = mfilt.firdecim(dfactor,[zeros(1,Hvec(end).Latency),1]);
        end
    else
        Hvec = decimdelaynoble(this,Hd.stage(s2),dfactor);
        Hvec = decimmergedelays(this,Hvec(end:-1:1));
          % Merge upsampler
        if strcmpi(class(Hvec(1)),'mfilt.firdecim'),
            Hvec(1) = mfilt.firsrc(ifactor,Hvec(1).DecimationFactor,...
                ifactor*Hvec(1).Numerator);
        else
            % Must be a delay
            Hvec(1) = mfilt.firinterp(ifactor,[zeros(1,Hvec(1).Latency),1]);
        end
    end     
    Hm2 = cascade(Hvec);
    Hm = parallel(Hm1,Hm2);
else

    % If we are returned a single section, the design was simple enough
    % that we can just use a single interpolator.
    Hm = feval(filtstruct, ifactor, dfactor, ifactor*Hd.Numerator);
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
function [Hm,final_upfactor] = cascadenoble(Hd,upfactor,ifactor,dfactor,filtstruct)

% If InterpolationFactor (LK) is greater than DecimationFactor (M).
%
%  /\ LK -- H1(z) -- H2(z^L) -- \/ M
%
%  /\ LK -- H2(z^L) -- H1(z) -- \/ M
%
%  /\ K -- H2(z) -- /\ L -- H1(z) -- \/ M
%
%  H2_interp_K(z) -- H1_src_L_M(z)
%
% If DecimationFactor (MK) is greater than InterpolationFactor (L)
%
%  /\ L -- H1(z) -- H2(z^M) -- \/ MK
%
%  /\ L -- H1(z) -- \/M -- H2(z) -- \/ K
%
%  H1_src_L_M(z) -- H2_decim_K(z)
%
%  J & K are often 1, in which case we dont need a multirate.

% When the decimation factor and upsampling factor do not share a
% denominator, we use the interpolation first code.
if ifactor > dfactor || gcd(upfactor, dfactor) == 1

    % Use GCD here bcause UPFACTOR and IFACTOR may not be integer
    % factors of each other, but we are sure to have some divisors
    % be equal.
    L = gcd(upfactor, ifactor);
    K = ifactor/L;
    final_upfactor = upfactor/L;

    % Remove the extra ones from the numerator.  We will get these
    % automatically when we use the interpolating structure.
    b_h2 = Hd.Stage(2).Numerator(1:L:end);

    % If there is "left over" interpolation put it into the first stage.
    if K == 1
        Hm1 = Hd.Stage(2);
        set(Hm1, 'Numerator', b_h2);
    else
        Hm1 = mfilt.firinterp(K, K*b_h2);
    end

    % Create the 2nd stage.
    if L == 1
        if dfactor == 1
            Hm2 = Hd.Stage(1);
        else
            Hm2 = mfilt.firdecim(dfactor, Hd.Stage(1).Numerator);
        end
    else
        Hm2 = feval(filtstruct, L, dfactor, L*Hd.Stage(1).Numerator);
    end
else
    M = gcd(upfactor, dfactor);
    K = dfactor/M;
    final_upfactor = upfactor/K;

    % Remove the extra ones from the numerator.  We will get these
    % automatically when we use the interpolating structure.
    b_h2 = Hd.Stage(2).Numerator(1:M:end);

    % Create the 1st stage.
    if M == 1
        if ifactor == 1
            Hm1 = Hd.Stage(1);
        else
            Hm1 = mfilt.firinterp(ifactor, ifactor*Hd.Stage(1).Numerator);
        end
    else
        Hm1 = feval(filtstruct, ifactor, M, ifactor*Hd.Stage(1).Numerator);
    end

    % If there is "left over" interpolation put it into the 2nd stage.
    if K == 1
        Hm2 = Hd.Stage(2);
        set(Hm2, 'Numerator', b_h2);
    else
        Hm2 = mfilt.firdecim(K, b_h2);
    end
end
Hm = mfilt.cascade(Hm1, Hm2);

% [EOF]
