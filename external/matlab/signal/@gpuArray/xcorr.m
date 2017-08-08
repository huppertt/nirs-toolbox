function [c,lags] = xcorr(x,varargin)
%XCORR Cross-correlation function estimates.
%   C = XCORR(A,B), where A and B are length M vectors (M>1), returns
%   the length 2*M-1 cross-correlation sequence C. If A and B are of
%   different length, the shortest one is zero-padded.  C will be a
%   row vector if A is a row vector, and a column vector if A is a
%   column vector.
%
%   XCORR produces an estimate of the correlation between two random
%   (jointly stationary) sequences:
%          C(m) = E[A(n+m)*conj(B(n))] = E[A(n)*conj(B(n-m))]
%   It is also the deterministic correlation between two deterministic
%   signals.
%
%   XCORR(A), when A is a vector, is the auto-correlation sequence.
%   XCORR(A), when A is an M-by-N matrix, is a large matrix with
%   2*M-1 rows whose N^2 columns contain the cross-correlation
%   sequences for all combinations of the columns of A.
%   The zeroth lag of the output correlation is in the middle of the
%   sequence, at element or row M.
%
%   XCORR(...,MAXLAG) computes the (auto/cross) correlation over the
%   range of lags:  -MAXLAG to MAXLAG, i.e., 2*MAXLAG+1 lags.
%   If missing, default is MAXLAG = M-1.
%
%   [C,LAGS] = XCORR(...)  returns a vector of lag indices (LAGS).
%
%   XCORR(...,SCALEOPT), normalizes the correlation according to SCALEOPT:
%     'biased'   - scales the raw cross-correlation by 1/M.
%     'unbiased' - scales the raw correlation by 1/(M-abs(lags)).
%     'coeff'    - normalizes the sequence so that the auto-correlations
%                  at zero lag are identically 1.0.
%     'none'     - no scaling (this is the default).
%
%   See also XCORR, GPUARRAY.

%   Copyright 2012 The MathWorks, Inc.

%   References:
%     S.J. Orfanidis, "Optimum Signal Processing. An Introduction"
%     2nd Ed. Macmillan, 1988.

narginchk(1,4);

[x,nshift] = shiftdim(x);
[xIsMatrix,autoFlag,maxlag,scaleType] = parseinput(x,varargin{:});

%Guard against the case where only maxlag is on the GPU. In that case
%dispatch to the CPU version since that is where the data resides.
xisgpu = isa(x, 'gpuArray');
if (autoFlag && ~xisgpu) || (~autoFlag && ~xisgpu && ~isa(varargin{1}, 'gpuArray')),
    %dispatch to CPU
    varargin = cellfun(@gather, varargin, 'UniformOutput', false);
    x = shiftdim(x, -nshift);
    [c, lags] = xcorr(x, varargin{:});
    return;
else
    % Run GPU version. Put data in the best places
    maxlag = gather(maxlag);
    x = gpuArray(x);
end


ifftStyle = getifftStyle(x,autoFlag,varargin{:});

if xIsMatrix,
    [c,M,N] = matrixCorr(x, ifftStyle);
else
    [c,M,N] = vectorXcorr(x,autoFlag, ifftStyle, varargin{:});
end

lags = -maxlag:maxlag;


if maxlag >= M,
    %	c = [zeros(maxlag-M+1,N^2);c(end-M+2:end,:);c(1:M,:);zeros(maxlag-M+1,N^2)];
    c_bottom = substruct('()', {size(c,1)-M+2:size(c,1), ':'});
    c_top = substruct('()', {1:M, ':'});
    c = [gpuArray.zeros(maxlag-M+1,N^2); ...
        subsref(c, c_bottom); ...
        subsref(c, c_top); ...
        gpuArray.zeros(maxlag-M+1,N^2)];
else
    %c = [c(end-maxlag+1:end,:);c(1:maxlag+1,:)];
    c_bottom = substruct('()', { ((size(c,1) - maxlag + 1):size(c,1)), ':' });
    c_top = substruct( '()', {1:(maxlag+1), ':'});
    c = [subsref(c, c_bottom); subsref(c, c_top)];
end

% Scale as specified
c = scaleXcorr(c,xIsMatrix,scaleType,autoFlag,M,maxlag,lags,x,varargin{:});

% If first vector is a row, return a row
c = shiftdim(c,-nshift);

%----------------------------------------------------------------
function [c,M,N] = matrixCorr(x, ifftStyle)
% Compute all possible auto- and cross-correlations for a matrix input

[M,N] = size(x);

X = fft(x,2^nextpow2(2*M-1));

Xc = conj(X);

[~,NX] = size(X);

X = reshape(X, [], 1, NX); %make X an MX x 1 x NX cube
C = bsxfun(@times, X, Xc); %multiply each page by XC w/expansion

C = reshape(C,[], NX*NX); %Turn back into a matrix



c = ifft(C, ifftStyle);

%----------------------------------------------------------------
function [c,M,N] = vectorXcorr(x,autoFlag,ifftStyle, varargin)
% Compute auto- or cross-correlation for vector inputs

x = reshape(x, [], 1); %make a column vector

[M,N] = size(x);


if autoFlag,
    % Autocorrelation
    % Compute correlation via FFT
    fftSize = 2^nextpow2(2*M - 1);
    X = fft(x,fftSize);
    fh = @(w) abs(w).^2;
    c = ifft(arrayfun(fh,X), ifftStyle);
    
else
    % xcorrelation
    y = gpuArray(varargin{1}); %make sure y is on the GPU
    y = reshape(y,[],1); %make a column vector
    L = length(y);
    
    % Cache the length(x)
    Mcached = M;
    
    % Recompute length(x) in case length(y) > length(x)
    M = max(Mcached,L);
    fftSize = 2^nextpow2(2*M - 1);
    
    if (L ~= Mcached) && any([L./Mcached, Mcached./L] > 10),
        
        % Vector sizes differ by a factor greater than 10,
        % fftfilt is faster
        neg_c = conj(fftfilt(conj(x),flipud(y))); % negative lags
        pos_c = flipud(fftfilt(conj(y),flipud(x))); % positive lags
        
        % Make them of almost equal length (remove zero-th lag from neg)
        lneg = length(neg_c);
        lpos = length(pos_c);
        
        %neg_c_idxed = neg_c(1:(numel(neg_c)-1));
        neg_c_subs = substruct('()', {1:numel(neg_c)-1} );
        neg_c_idxed = subsref(neg_c, neg_c_subs);
        
        if ~isempty(neg_c_idxed),
            neg_c = [gpuArray.zeros(lpos-lneg,1); ...
                neg_c_idxed];
        else
            neg_c = gpuArray.zeros(lpos-lneg,1);
        end
        
        pos_c = [pos_c;gpuArray.zeros(lneg-lpos,1)];
        
        c = [pos_c;neg_c];
        
    else
        if L ~= Mcached,
            % Force equal lengths
            if L > Mcached
                x = [x;gpuArray.zeros(L-Mcached,1)];
                
            else
                y = [y;gpuArray.zeros(Mcached-L,1)];
            end
        end
        
        % Transform both vectors
        X = fft(x,fftSize);
        Y = fft(y,fftSize);
        
        % Compute cross-correlation
        c = ifft(X.*conj(Y), ifftStyle);
    end
end

%----------------------------------------------------------------
function c = scaleXcorr(c,xIsMatrix,scaleType,autoFlag,...
    M,maxlag,lags,x,varargin)
% Scale correlation as specified

switch scaleType,
    case 'none',
        return
    case 'biased',
        % Scales the raw cross-correlation by 1/M.
        c = c./M;
    case 'unbiased',
        % Scales the raw correlation by 1/(M-abs(lags)).
        scale = (M-abs(lags)).';
        scale(scale<=0)=1; % avoid divide by zero, when correlation is zero
        
        c = bsxfun(@rdivide, c, scale);
    case 'coeff',
        % Normalizes the sequence so that the auto-correlations
        % at zero lag are identically 1.0.
        if ~autoFlag,
            % xcorr(x,y)
            % Compute autocorrelations at zero lag
            fh = @(w) abs(w).^2;
            cxx0 = sum(arrayfun(fh,x));
            cyy0 = sum(arrayfun(fh, varargin{1} ));
            scale = sqrt(cxx0*cyy0);
            c = c./scale;
        else
            if ~xIsMatrix,
                % Autocorrelation case, simply normalize by c[0]
                %c(maxlag+1)
                c_subs = substruct('()', {maxlag+1});
                c0 = subsref(c, c_subs);
                
                c = c./c0;
            else
                % Compute the indices corresponding to the columns for which
                % we have autocorrelations (e.g. if c = n by 9, the autocorrelations
                % are at columns [1,5,9] the other columns are cross-correlations).
                [~,nc] = size(c);
                jkl = reshape(1:nc,sqrt(nc),sqrt(nc))';
                
                %c(maxlag+1,diag(jkl))
                subsref_c = substruct('()', {(maxlag+1), diag(jkl)});
                tmp = subsref(c, subsref_c);
                
                fh = @(u,v)sqrt(u)*sqrt(v);
                tmp = bsxfun(fh, reshape(tmp,[],1), tmp);
                tmp = reshape(tmp, 1, []);
                c = bsxfun(@rdivide, c, tmp);
            end
        end
end

%----------------------------------------------------------------
function [xIsMatrix,autoFlag,maxlag,scaleType] = parseinput(x,varargin)
%    Parse the input and determine optional parameters:
%
%    Outputs:
%       xIsMatrix - flag indicating if x is a matrix
%       AUTOFLAG  - 1 if autocorrelation, 0 if xcorrelation
%       maxlag    - Number or lags to compute
%       scaleType - String with the type of scaling wanted

% Set some defaults
scaleType = '';
autoFlag = 1; % Assume autocorrelation until proven otherwise
maxlag = [];
xIsMatrix = false;

switch nargin,
    case 2,
        % Can be (x,y), (x,maxlag), or (x,scaleType)
        if ischar(varargin{1}),
            % Second arg is scaleType
            scaleType = varargin{1};
            
        elseif (isnumeric(varargin{1}) || isa(varargin{1}, 'gpuArray'))
            % Can be y or maxlag
            if length(varargin{1}) == 1,
                maxlag = varargin{1};
            else
                autoFlag = 0;
                y = varargin{1};
            end
        else
            % Not recognized
            error(message('signal:xcorr:UnknInput'));
        end
    case 3,
        % Can be (x,y,maxlag), (x,maxlag,scaleType) or (x,y,scaleType)
        maxlagflag = 0; % By default, assume 3rd arg is not maxlag
        if ischar(varargin{2}),
            % Must be scaletype
            scaleType = varargin{2};
            
        elseif isnumeric(varargin{2}) || isa(varargin{2}, 'gpuArray'),
            % Must be maxlag
            maxlagflag = 1;
            maxlag = varargin{2};
            
        else
            % Not recognized
            error(message('signal:xcorr:UnknInput'));
        end
        
        if (isnumeric(varargin{1}) || isa(varargin{1}, 'gpuArray')) ,
            if maxlagflag,
                autoFlag = 0;
                y = varargin{1};
            else
                % Can be y or maxlag
                if length(varargin{1}) == 1,
                    maxlag = varargin{1};
                else
                    autoFlag = 0;
                    y = varargin{1};
                end
            end
        else
            % Not recognized
            error(message('signal:xcorr:UnknInput'));
        end
        
    case 4,
        % Must be (x,y,maxlag,scaleType)
        autoFlag = 0;
        y = varargin{1};
        
        maxlag = varargin{2};
        
        scaleType = varargin{3};
end

% Determine if x is a matrix or a vector
[xIsMatrix,m] = parse_x(x);



if ~autoFlag,
    % Test y for correctness
    maxlag = parse_y(y,m,xIsMatrix,maxlag);
end

maxlag = parse_maxlag(maxlag,m);


% Test the scaleType validity
scaleType = parse_scaleType(scaleType,autoFlag,m,varargin{:});


%-------------------------------------------------------------------
function [xIsMatrix,m] = parse_x(x)


xIsMatrix = (size(x,2) > 1);

m = size(x,1);


%-------------------------------------------------------------------
function maxlag = parse_y(y,m,xIsMatrix,maxlag)
[my,ny] = size(y);
if ~any([my,ny] == 1),
    % Second arg is a matrix, error
    error(message('signal:xcorr:BMustBeVector', 'B'));
end

if xIsMatrix,
    % Can't do xcorr(matrix,vector)
    error(message('signal:xcorr:MismatchedAB', 'B', 'A'));
end
if (length(y) > m) && isempty(maxlag),
    % Compute the default maxlag based on the length of y
    maxlag = length(y) - 1;
end

%-------------------------------------------------------------------
function maxlag = parse_maxlag(maxlag,m)
if isempty(maxlag),
    % Default hasn't been assigned yet, do so
    maxlag = m-1;
else
    % Test maxlag for correctness
    if  length(maxlag)>1
        error(message('signal:xcorr:MaxLagMustBeScalar', 'MAXLAG'));
    end
    if maxlag < 0,
        maxlag = abs(maxlag);
    end
    if maxlag ~= round(maxlag),
        error(message('signal:xcorr:MaxLagMustBeInteger', 'MAXLAG'));
    end
end

%--------------------------------------------------------------------
function ifftStyle = getifftStyle(x,autoFlag,varargin)
% If outputs should be real, force a symmetric transform

ifftStyle = 'nonsymmetric'; % Flag to determine whether we should force real corr

if (isreal(x) && autoFlag) || (isreal(x) && isreal(varargin{1})),
    ifftStyle = 'symmetric';
end

%--------------------------------------------------------------------
function scaleType = parse_scaleType(scaleType,autoFlag,m,varargin)
if isempty(scaleType),
    scaleType = 'none';
else
    scaleOpts = {'biased','unbiased','coeff','none'};
    indx = find(strncmpi(scaleType, scaleOpts, length(scaleType)));
    
    if isempty(indx),
        error(message('signal:xcorr:UnknInput'));
    else
        scaleType = scaleOpts{indx};
    end
    
    if ~autoFlag && ~strcmpi(scaleType,'none') && (m ~= length(varargin{1})),
        error(message('signal:xcorr:NoScale', 'SCALEOPT', 'none', 'A', 'B'));
    end
end

% EOF

