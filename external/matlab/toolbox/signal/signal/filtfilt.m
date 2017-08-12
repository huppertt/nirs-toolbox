function y = filtfilt(b,a,x)
%FILTFILT Zero-phase forward and reverse digital IIR filtering.
%   Y = FILTFILT(B, A, X) filters the data in vector X with the filter
%   described by vectors A and B to create the filtered data Y.  The filter
%   is described by the difference equation:
%
%     a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                           - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
%
%   The length of the input X must be more than three times the filter
%   order, defined as max(length(B)-1,length(A)-1).
%
%   Y = FILTFILT(SOS, G, X) filters the data in vector X with the
%   second-order section (SOS) filter described by the matrix SOS and the
%   vector G.  The coefficients of the SOS matrix must be expressed using
%   an Lx6 matrix where L is the number of second-order sections. The scale
%   values of the filter must be expressed using the vector G. The length
%   of G must be between 1 and L+1, and the length of input X must be more
%   than three times the filter order (input length must be greater than
%   one when the order is zero). You can use filtord(SOS) to get the
%   order of the filter. The SOS matrix should have the following form:
%
%   SOS = [ b01 b11 b21 a01 a11 a21
%           b02 b12 b22 a02 a12 a22
%           ...
%           b0L b1L b2L a0L a1L a2L ]
%
%   Y = FILTFILT(D, X) filters the data in vector X with the digital filter
%   D. You design a digital filter, D, by calling the <a href="matlab:help designfilt">designfilt</a> function.
%   The length of the input X must be more than three times the filter
%   order. You can use filtord(D) to get the order of the digital filter D.
%
%   After filtering in the forward direction, the filtered sequence is then
%   reversed and run back through the filter; Y is the time reverse of the
%   output of the second filtering operation.  The result has precisely
%   zero phase distortion, and magnitude modified by the square of the
%   filter's magnitude response. Startup and ending transients are
%   minimized by matching initial conditions.
%
%   Note that FILTFILT should not be used when the intent of a filter is to
%   modify signal phase, such as differentiators and Hilbert filters.
%
%   % Example 1:
%   %   Zero-phase filter a noisy ECG waveform using an IIR filter.
%  
%   load noisysignals x;                    % noisy waveform
%   [b,a] = butter(12,0.2,'low');           % IIR filter design
%   y = filtfilt(b,a,x);                    % zero-phase filtering
%   y2 = filter(b,a,x);                     % conventional filtering
%   plot(x,'k-.'); grid on ; hold on
%   plot([y y2],'LineWidth',1.5);
%   legend('Noisy ECG','Zero-phase Filtering','Conventional Filtering');
%
%   % Example 2:
%   %   Use the designfilt function to design a highpass IIR digital filter 
%   %   with order 4, passband frequency of 75 KHz, and a passband ripple 
%   %   of 0.2 dB. Sample rate is 200 KHz. Apply zero-phase filtering to a
%   %   vector of data. 
%  
%   D = designfilt('highpassiir', 'FilterOrder', 4, ...
%            'PassbandFrequency', 75e3, 'PassbandRipple', 0.2,...
%            'SampleRate', 200e3);
%
%   x = rand(1000,1);
%   y = filtfilt(D,x);
%
%   See also FILTER, SOSFILT.

%   References:
%     [1] Sanjit K. Mitra, Digital Signal Processing, 2nd ed.,
%         McGraw-Hill, 2001
%     [2] Fredrik Gustafsson, Determining the initial states in forward-
%         backward filtering, IEEE Transactions on Signal Processing,
%         pp. 988-992, April 1996, Volume 44, Issue 4

%   Copyright 1988-2014 The MathWorks, Inc.

narginchk(3,3)

% Only double precision is supported
if ~isa(b,'double') || ~isa(a,'double') || ~isa(x,'double')
    error(message('signal:filtfilt:NotSupported'));
end
if isempty(b) || isempty(a) || isempty(x)
    y = [];
    return
end

% If input data is a row vector, convert it to a column
isRowVec = size(x,1)==1;
if isRowVec
    x = x(:);
end
[Npts,Nchans] = size(x);

% Parse SOS matrix or coefficients vectors and determine initial conditions
[b,a,zi,nfact,L] = getCoeffsAndInitialConditions(b,a,Npts);

% Filter the data
if Nchans==1
    if Npts<10000
        y = ffOneChanCat(b,a,x,zi,nfact,L);
    else
        y = ffOneChan(b,a,x,zi,nfact,L);
    end
else
    y = ffMultiChan(b,a,x,zi,nfact,L);
end

if isRowVec
    y = y.';   % convert back to row if necessary
end

%--------------------------------------------------------------------------
function [b,a,zi,nfact,L] = getCoeffsAndInitialConditions(b,a,Npts)

[L, ncols] = size(b);
na = numel(a);

% Rules for the first two inputs to represent an SOS filter:
% b is an Lx6 matrix with L>1 or,
% b is a 1x6 vector, its 4th element is equal to 1 and a has less than 2
% elements. 
if ncols==6 && L==1 && na<=2,
    if b(4)==1,
        warning(message('signal:filtfilt:ParseSOS', 'SOS', 'G'));
    else
        warning(message('signal:filtfilt:ParseB', 'a01', 'SOS'));
    end
end
issos = ncols==6 && (L>1 || (b(4)==1 && na<=2));
if issos,
    %----------------------------------------------------------------------
    % b is an SOS matrix, a is a vector of scale values
    %----------------------------------------------------------------------
    g = a(:);
    ng = na;
    if ng>L+1,
        error(message('signal:filtfilt:InvalidDimensionsScaleValues', L + 1));
    elseif ng==L+1,
        % Include last scale value in the numerator part of the SOS Matrix
        b(L,1:3) = g(L+1)*b(L,1:3);
        ng = ng-1;
    end
    for ii=1:ng,
        % Include scale values in the numerator part of the SOS Matrix
        b(ii,1:3) = g(ii)*b(ii,1:3);
    end
    
    ord = filtord(b);
    
    a = b(:,4:6).';
    b = b(:,1:3).';
         
    nfact = max(1,3*ord); % length of edge transients    
    if Npts <= nfact % input data too short
        error(message('signal:filtfilt:InvalidDimensionsDataShortForFiltOrder',num2str(nfact)))
    end
    
    % Compute initial conditions to remove DC offset at beginning and end of
    % filtered sequence.  Use sparse matrix to solve linear system for initial
    % conditions zi, which is the vector of states for the filter b(z)/a(z) in
    % the state-space formulation of the filter.
    zi = zeros(2,L);
    for ii=1:L,
        rhs  = (b(2:3,ii) - b(1,ii)*a(2:3,ii));
        zi(:,ii) = ( eye(2) - [-a(2:3,ii),[1;0]] ) \ rhs;
    end
    
else
    %----------------------------------------------------------------------
    % b and a are vectors that define the transfer function of the filter
    %----------------------------------------------------------------------
    L = 1;
    % Check coefficients
    b = b(:);
    a = a(:);
    nb = numel(b);
    nfilt = max(nb,na);   
    nfact = max(1,3*(nfilt-1));  % length of edge transients
    if Npts <= nfact      % input data too short
        error(message('signal:filtfilt:InvalidDimensionsDataShortForFiltOrder',num2str(nfact)));
    end
    % Zero pad shorter coefficient vector as needed
    if nb < nfilt
        b(nfilt,1)=0;
    elseif na < nfilt
        a(nfilt,1)=0;
    end
    
    % Compute initial conditions to remove DC offset at beginning and end of
    % filtered sequence.  Use sparse matrix to solve linear system for initial
    % conditions zi, which is the vector of states for the filter b(z)/a(z) in
    % the state-space formulation of the filter.
    if nfilt>1
        rows = [1:nfilt-1, 2:nfilt-1, 1:nfilt-2];
        cols = [ones(1,nfilt-1), 2:nfilt-1, 2:nfilt-1];
        vals = [1+a(2), a(3:nfilt).', ones(1,nfilt-2), -ones(1,nfilt-2)];
        rhs  = b(2:nfilt) - b(1)*a(2:nfilt);
        zi   = sparse(rows,cols,vals) \ rhs;
        % The non-sparse solution to zi may be computed using:
        %      zi = ( eye(nfilt-1) - [-a(2:nfilt), [eye(nfilt-2); ...
        %                                           zeros(1,nfilt-2)]] ) \ ...
        %          ( b(2:nfilt) - b(1)*a(2:nfilt) );
    else
        zi = zeros(0,1);
    end
end
%--------------------------------------------------------------------------
function y = ffOneChanCat(b,a,y,zi,nfact,L)

for ii=1:L
    % Single channel, data explicitly concatenated into one vector
    y = [2*y(1)-y(nfact+1:-1:2); y; 2*y(end)-y(end-1:-1:end-nfact)]; %#ok<AGROW>
    
    % filter, reverse data, filter again, and reverse data again
    y = filter(b(:,ii),a(:,ii),y,zi(:,ii)*y(1));
    y = y(end:-1:1);
    y = filter(b(:,ii),a(:,ii),y,zi(:,ii)*y(1));
    
    % retain reversed central section of y
    y = y(end-nfact:-1:nfact+1);
end

%--------------------------------------------------------------------------
function y = ffOneChan(b,a,xc,zi,nfact,L)
% Perform filtering of input data with no phase distortion
%
% xc: one column of input data
% yc: one column of output, same dimensions as xc
% a,b: IIR coefficients, both of same order/length
% zi: initial states
% nfact: scalar
% L: scalar

% Extrapolate beginning and end of data sequence using reflection.
% Slopes of original and extrapolated sequences match at end points,
% reducing transients.
%
% yc is length numel(x)+2*nfact
%
% We wish to filter the appended sequence:
% yc = [2*xc(1)-xc(nfact+1:-1:2); xc; 2*xc(end)-xc(end-1:-1:end-nfact)];
%
% We would use the following code:
% Filter the input forwards then again in the backwards direction
%   yc = filter(b,a,yc,zi*yc(1));
%   yc = yc(length(yc):-1:1); % reverse the sequence
%
% Instead of reallocating and copying, just do filtering in pieces
% Filter the first part, then xc, then the last part
% Retain each piece, feeding final states as next initial states
% Output is yc = [yc1 yc2 yc3]
%
% yc2 can be long (matching user input length)
% yc3 is short (nfilt samples)
%
for ii=1:L
    
    xt = -xc(nfact+1:-1:2) + 2*xc(1);
    [~,zo] = filter(b(:,ii),a(:,ii), xt, zi(:,ii)*xt(1)); % yc1 not needed
    [yc2,zo] = filter(b(:,ii),a(:,ii), xc, zo);
    xt = -xc(end-1:-1:end-nfact) + 2*xc(end);
    yc3 = filter(b(:,ii),a(:,ii), xt, zo);
    
    % Reverse the sequence
    %   yc = [flipud(yc3) flipud(yc2) flipud(yc1)]
    % which we can do as
    %   yc = [yc3(end:-1:1); yc2(end:-1:1); yc1(end:-1:1)];
    %
    % But we don't do that explicitly.  Instead, we wish to filter this
    % reversed result as in:
    %   yc = filter(b,a,yc,zi*yc(1));
    % where yc(1) is now the last sample of the (unreversed) yc3
    %
    % Once again, we do this in pieces:
    [~,zo] = filter(b(:,ii),a(:,ii), yc3(end:-1:1), zi(:,ii)*yc3(end));
    yc5 = filter(b(:,ii),a(:,ii), yc2(end:-1:1), zo);
    
    % Conceptually restore the sequence by reversing the last pieces
    %    yc = yc(length(yc):-1:1); % restore the sequence
    % which is done by
    %    yc = [yc6(end:-1:1); yc5(end:-1:1); yc4(end:-1:1)];
    %
    % However, we only need to retain the reversed central samples of filtered
    % output for the final result:
    %    y = yc(nfact+1:end-nfact);
    %
    % which is the reversed yc5 piece only.
    %
    % This means we don't need yc4 or yc6.  We need to compute yc4 only because
    % the final states are needed for yc5 computation.  However, we can omit
    % yc6 entirely in the above filtering step.
    xc = yc5(end:-1:1);
end
y = xc;

%--------------------------------------------------------------------------
function y = ffMultiChan(b,a,xc,zi,nfact,L)
% Perform filtering of input data with no phase distortion
%
% xc: matrix of input data
% yc: matrix of output data, same dimensions as xc
% a,b: IIR coefficients, both of same order/length
% zi: initial states
% nfact: scalar
% L: scalar

% Same comments as in ffOneChan, except here we need to use bsxfun.
% Instead of doing scalar subtraction with a vector, we are doing
% vector addition with a matrix.  bsxfun replicates the vector
% for us.
%
% We also take care to preserve column dimensions
%
sz = size(xc);
xc = reshape(xc,sz(1),[]);%converting N-D matrix to 2-D.

for ii=1:L
    xt = bsxfun(@minus, 2*xc(1,:), xc(nfact+1:-1:2,:));
    [~,zo] = filter(b(:,ii),a(:,ii),xt,zi(:,ii)*xt(1,:)); % outer product
    [yc2,zo] = filter(b(:,ii),a(:,ii),xc,zo);
    xt = bsxfun(@minus, 2*xc(end,:), xc(end-1:-1:end-nfact,:));
    yc3 = filter(b(:,ii),a(:,ii),xt,zo);
    
    [~,zo] = filter(b(:,ii),a(:,ii),yc3(end:-1:1,:),zi(:,ii)*yc3(end,:)); % outer product
    yc5 = filter(b(:,ii),a(:,ii),yc2(end:-1:1,:), zo);
    
    xc = yc5(end:-1:1,:);
end
y = xc;
y = reshape(y,sz);% To match the input size.
