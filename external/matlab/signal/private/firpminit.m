function [nfilt,ff,grid,des,wt,ftype,sign_val,hilbert,neg] = firpminit(order, ff, aa, varargin)
%FIRPMINIT   

%   Copyright 2006-2013 The MathWorks, Inc.

narginchk(3,6);

if all( ff(2:2:end)-ff(1:2:end) == 0),
	error(message('signal:firpminit:InvalidFreqVecZeroWidth'))
end

if order < 3
    error(message('signal:firpminit:InvalidRangeOrder'));
end
%
% Define default values for input arguments:
%
ftype = 'f';
wtx = ones(fix((1+length(ff))/2),1);
lgrid = 16;   % Grid density (should be at least 16)
%
% parse inputs and alter defaults
%
%  First find cell array and remove it if present
for i=1:length(varargin)
    if iscell(varargin{i})
        lgrid = varargin{i}{:};
        if lgrid < 16,
            warning(message('signal:firpminit:InvalidParamGridDensity'));
        end
        if lgrid < 1,
            error(message('signal:firpminit:MustBePositive'));
        end
        varargin(i) = [];
        break
    end
end
if length(varargin) == 1
    if ischar(varargin{1})
        ftype = varargin{1};
    else
        wtx = varargin{1};
    end
elseif length(varargin)==2
    wtx = varargin{1};
    ftype = varargin{2};
end

if isempty(ftype), ftype = 'f'; end

%
% Error checking
%
if rem(length(ff),2)
    error(message('signal:firpminit:MustBeEven'));
end
if any((ff < 0) | (ff > 1))
    error(message('signal:firpminit:InvalidRangeFreqs'));
end
df = diff(ff);
if (any(df < 0))
    error(message('signal:firpminit:InvalidFreqVecNonDecreasingFreqs'));
end
if length(wtx) ~= fix((1+length(ff))/2)
    error(message('signal:firpminit:InvalidDimensions'));
end

if (any(sign(wtx) == 1) && any(sign(wtx) == -1)) || any(sign(wtx) == 0),
    error(message('signal:firpminit:InvalidRangeWeights'));
end

% Determine "Frequency Response Function" (frf)

% The following comment, MATLAB compiler pragma, is necessary to allow the
% Compiler to find the FIRPMFRF private function. Don't remove.
%#function firpmfrf
%#function multiband
%#function lowpass
%#function highpass
%#function bandpass
%#function bandstop
%#function invsinc
%#function hilbfilt
%#function differentiator

if ischar(aa)
    frf = str2func(aa);
    frf_params = {};
elseif isa(aa, 'function_handle');
    frf = aa;
    frf_params = {};
elseif iscell(aa)
    frf = aa{1};
    if ischar(frf)
        frf = str2func(frf);
    end

    frf_params = aa(2:end);
else
    % Check for valid filter length. Ideally we would check in all cases,
    % for now we only do it when aa is a vector
    
    % Cast to enforce precision rules
    aa = double(aa);
    
    exception = 0;
    if any(strncmpi(ftype, {'differentiator','hilbert'}, length(ftype))),
        exception = 1;
    end
    [order,msg1,msg2,msgobj] = firchk(order,ff(end),aa,exception);
    if ~isempty(msg1)
      error(msgobj); 
    end
    
    if ~isempty(msg2),
        warning(message('signal:firpminit:OrderIncreased', ...
          msg2,'h','firpm(N,F,A,W,''h'')'));
    end
    % Grid and weights will be cast to double in firpm to enforce precision
    % rules.
    frf = @firpmfrf;
    frf_params = { aa, strcmpi(ftype(1),'d') };
end

% We need the following code to generate fcn handles to private functions.
% Otherwise the syntax h=firpm(n,f,{@firpmfrf,m},w); will not work
fcns = functions(frf);
if isempty(fcns.file)
    frf = str2func(func2str(frf));
end

%
% Determine symmetry of filter
%
sign_val = 1.0;
nfilt = order + 1;        % filter length
nodd = rem(nfilt,2);      % nodd == 1 ==> filter length is odd
% nodd == 0 ==> filter length is even

hilbert = 0;
if ftype(1) == 'h' || ftype(1) == 'H'
    ftype = 3;  % Hilbert transformer
    hilbert = 1;
    if ~nodd
        ftype = 4;
    end
elseif ftype(1) == 'd' || ftype(1) == 'D'
    ftype = 4;  % Differentiator
    sign_val = -1;
    if nodd
        ftype = 3;
    end
else
    % If symmetry was not specified, call the fresp function
    % with 'defaults' string and a cell-array of the actual
    % function call arguments to query the default value.
    try
        h_sym = frf('defaults', {order, ff, [], wtx, frf_params{:}} );
    catch
        h_sym = 'even';
    end

    if ~any(strcmp(h_sym,{'even' 'odd'})),
        error(message('signal:firpminit:InvalidParamSym', h_sym, frf, '''even''', '''odd'''));
    end

    switch h_sym
        case 'even'
            ftype = 1;  % Regular filter
            if ~nodd
                ftype = 2;
            end
        case 'odd'
            ftype = 3;  % Odd (antisymmetric) filter
            if ~nodd
                ftype = 4;
            end
    end
end

if (ftype == 3 || ftype == 4)
    neg = 1;  % neg == 1 ==> antisymmetric imp resp,
else
    neg = 0;  % neg == 0 ==> symmetric imp resp
end


%
% Create grid of frequencies on which to perform firpm exchange iteration
%
grid = firpmgrid(nfilt,lgrid,ff,neg,nodd);
while length(grid) <= nfilt,
    lgrid = lgrid*4;  % need more grid points
    grid = firpmgrid(nfilt,lgrid,ff,neg,nodd);
end
%
% Get desired frequency characteristics at the frequency points
% in the specified frequency band intervals.
%
% NOTE! The frf needs to see normalized frequencies in the range
% [0,1].
[des,wt] = frf(order, ff, grid, wtx, frf_params{:}); 

%--------------------------------------------------------------------------
function grid = firpmgrid(nfilt,lgrid,ff,neg,nodd)
% firpmgrid
%    Generate frequency grid
nfcns = fix(nfilt/2);
if nodd == 1 && neg == 0
    nfcns = nfcns + 1;
end
grid(1) = ff(1);
delf = 1/(lgrid*nfcns);
% If value at frequency 0 is constrained, make sure first grid point
% is not too small:
if neg ~= 0 && grid(1) < delf
    % Handle narrow bands
    if ff(1) > sqrt(eps)
       grid(1) = ff(1);
    elseif delf < ff(2)
        grid(1) = delf;
    else
        grid(1) = 0.5*(ff(2) - ff(1));
    end
end
j = 1;
l = 1;
while (l+1) <= length(ff)
    fup = ff(l+1);
    newgrid = (grid(j)+delf):delf:(fup+delf);
    if length(newgrid) < 11
        delf1 = ((fup+delf) - (grid(j)+delf)) / 10;
        newgrid = (grid(j)+delf1):delf1:(fup+delf1);
    end
    grid = [grid newgrid];
    jend = length(grid);
    if jend > 1,
        grid(jend-1) = fup;
        j = jend;
    else
        j = jend + 1;
    end

    l = l + 2;
    if (l+1) <= length(ff)
        grid(j) = ff(l);
    end
end
ngrid = j - 1;
% If value at frequency 1 is constrained, remove that grid point:
if neg == nodd && (grid(ngrid) > 1-delf)
    if ff(end-1) < 1-delf
        ngrid = ngrid - 1;
    else
        grid(ngrid) = ff(end-1);
    end
end
grid = grid(1:ngrid);



% [EOF]
