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
%   % Example:
%   %   Measure the delay between two correlated signals.
%
%   load sensorData s1 s2
%   [acor,lag] = xcorr(s2,s1);
%   [~,I] = max(abs(acor));
%   timeDiff = lag(I)         % sensor 2 leads sensor 1 by 350 samples
%   subplot(311); plot(s1); title('s1');
%   subplot(312); plot(s2); title('s2');
%   subplot(313); plot(lag,acor);
%   title('Cross-correlation between s1 and s2')
%
%   See also XCOV, CORRCOEF, CONV, CCONV, COV and XCORR2.

%   Copyright 1988-2013 The MathWorks, Inc.

%   References:
%     S.J. Orfanidis, "Optimum Signal Processing. An Introduction"
%     2nd Ed. Macmillan, 1988.

narginchk(1,4);

[x,nshift] = shiftdim(x);
[xIsMatrix,autoFlag,maxlag,scaleType] = parseinput(x,varargin{:});

if xIsMatrix,
    [c,M,N] = matrixCorr(x);
else
    [c,M,N] = vectorXcorr(x,autoFlag,varargin{:});
end

% Force correlation to be real when inputs are real
c = forceRealCorr(c,x,autoFlag,varargin{:});


lags = -maxlag:maxlag;


% Keep only the lags we want and move negative lags before positive lags
if maxlag >= M,
    c = [zeros(maxlag-M+1,N^2);c(end-M+2:end,:);c(1:M,:);zeros(maxlag-M+1,N^2)];
else
    c = [c(end-maxlag+1:end,:);c(1:maxlag+1,:)];
end

% Scale as specified
c = scaleXcorr(c,xIsMatrix,scaleType,autoFlag,M,maxlag,lags,x,varargin{:});

% If first vector is a row, return a row
c = shiftdim(c,-nshift);

%----------------------------------------------------------------
function [c,M,N] = matrixCorr(x)
% Compute all possible auto- and cross-correlations for a matrix input

[M,N] = size(x);

X = fft(x,2^nextpow2(2*M-1));

Xc = conj(X);

[MX,NX] = size(X);
C = zeros(MX,NX*NX);
for n =1:N,
    C(:,(((n-1)*N)+1):(n*N)) = repmat(X(:,n),1,N).*Xc;
end

c = ifft(C);

%----------------------------------------------------------------
function [c,M,N] = vectorXcorr(x,autoFlag,varargin)
% Compute auto- or cross-correlation for vector inputs

x = x(:);

[M,N] = size(x);

if autoFlag,
    % Autocorrelation
    % Compute correlation via FFT
    X = fft(x,2^nextpow2(2*M-1));
    c = ifft(abs(X).^2);
    
else
    % xcorrelation
    y = varargin{1};
    y = y(:);
    L = length(y);
    
    % Cache the length(x)
    Mcached = M;
    
    % Recompute length(x) in case length(y) > length(x)
    M = max(Mcached,L);
    
    if (L ~= Mcached) && any([L./Mcached, Mcached./L] > 10),
        
        % Vector sizes differ by a factor greater than 10,
        % fftfilt is faster
        neg_c = conj(fftfilt(conj(x),flipud(y))); % negative lags
        pos_c = flipud(fftfilt(conj(y),flipud(x))); % positive lags
        
        % Make them of almost equal length (remove zero-th lag from neg)
        lneg = length(neg_c); lpos = length(pos_c);
        neg_c = [zeros(lpos-lneg,1);neg_c(1:end-1)];
        pos_c = [pos_c;zeros(lneg-lpos,1)];
        
        c = [pos_c;neg_c];
        
    else
        if L ~= Mcached,
            % Force equal lengths
            if L > Mcached
                x = [x;zeros(L-Mcached,1)];
                
            else
                y = [y;zeros(Mcached-L,1)];
            end
        end
        
        % Transform both vectors
        X = fft(x,2^nextpow2(2*M-1));
        Y = fft(y,2^nextpow2(2*M-1));
        
        % Compute cross-correlation
        c = ifft(X.*conj(Y));
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
        
        if xIsMatrix,
            scale = repmat(scale,1,size(c,2));
        end
        c = c./scale;
    case 'coeff',
        % Normalizes the sequence so that the auto-correlations
        % at zero lag are identically 1.0.
        if ~autoFlag,
            % xcorr(x,y)
            % Compute autocorrelations at zero lag
            cxx0 = sum(abs(x).^2);
            cyy0 = sum(abs(varargin{1}).^2);
            scale = sqrt(cxx0*cyy0);
            c = c./scale;
        else
            if ~xIsMatrix,
                % Autocorrelation case, simply normalize by c[0]
                c = c./c(maxlag+1);
            else
                % Compute the indices corresponding to the columns for which
                % we have autocorrelations (e.g. if c = n by 9, the autocorrelations
                % are at columns [1,5,9] the other columns are cross-correlations).
                [mc,nc] = size(c);
                jkl = reshape(1:nc,sqrt(nc),sqrt(nc))';
                tmp = sqrt(c(maxlag+1,diag(jkl)));
                tmp = tmp(:)*tmp;
                cdiv = repmat(tmp(:).',mc,1);
                c = c ./ cdiv; % The autocorrelations at zero-lag are normalized to
                % one
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
xIsMatrix = false; %#ok<NASGU>

switch nargin,
    case 2,
        % Can be (x,y), (x,maxlag), or (x,scaleType)
        if ischar(varargin{1}),
            % Second arg is scaleType
            scaleType = varargin{1};
            
        elseif isnumeric(varargin{1}),
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
            
        elseif isnumeric(varargin{2}),
            % Must be maxlag
            maxlagflag = 1;
            maxlag = varargin{2};
            
        else
            % Not recognized
            error(message('signal:xcorr:UnknInput'));
        end
        
        if isnumeric(varargin{1}),
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
% Cast to enforce precision rules
maxlag = double(maxlag);

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
function c = forceRealCorr(c,x,autoFlag,varargin)
% Force correlation to be real when inputs are real

forceReal = 0; % Flag to determine whether we should force real corr

if (isreal(x) && autoFlag) || (isreal(x) && isreal(varargin{1})),
    forceReal = 1;
end


if forceReal,
    c = real(c);
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

