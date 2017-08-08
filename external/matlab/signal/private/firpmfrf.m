function [DH,DW] = firpmfrf(N, F, GF, W, A, diff_flag)
%FIRPMFRF Frequency Response Function for FIRPM.
%  FIRPM(N,F,A, ...) or
%  FIRPM(N,F,{@firpmfrf,A}, ...) designs a linear-phase FIR filter
%  using FIRPM.
%
%  The symmetry option SYM defaults to 'even' if unspecified in the
%  call to FIRPM.
%
%  See also FIRPM.

%   Copyright 1988-2004 The MathWorks, Inc.

%  [DH,DW]=FIRPMFRF(N,F,GF,W,A,diff_flag)
%      N: filter order (length minus one)
%      F: vector of band edges
%     GF: vector of interpolated grid frequencies
%      W: vector of weights, one per band
%      A: vector of amplitudes of desired frequency response at band edges F
% diff_flag: ==1 for differentiator (1/f) weights, ==0 otherwise
%
%     DH: vector of desired filter response (mag & phase)
%     DW: vector of weights (positive)
%
% NOTE: DH(GF) and DW(GF) are specified as functions of frequency

% Support query by FIRPM for the default symmetry option:
if nargin==2,
  % Return symmetry default:
  if strcmp(N,'defaults'),
    DH = 'even';   % can be 'even' or 'odd'
    return
  end
end

if nargin < 6
    diff_flag = 0;
end

% Prevent discontinuities in desired function
for k=2:2:length(F)-2
    if F(k) == F(k+1)
        error(message('signal:firpmfrf:InvalidFreqVec'))
    end
end
if length(F) ~= length(A)
    error(message('signal:firpmfrf:InvalidDimensions'))
end

nbands = length(A)/2;
l = 1;
while (l+1)/2 <= nbands
    sel = find( GF>=F(l) & GF<=F(l+1) );
    % desired magnitude is line connecting A(l) to A(l+1)
    if F(l+1)~=F(l)   % 
        slope=(A(l+1)-A(l))/(F(l+1)-F(l));
        DH(sel) = polyval([slope A(l)-slope*F(l)],GF(sel));
    else   % zero bandwidth band 
        DH(sel) = (A(l)+A(l+1))/2;  
    end
    DW(sel) = W((l+1)/2) ./ (1 +(diff_flag & A(l+1) >= .0001)*(GF(sel)/2 - 1));
    l = l + 2;
end

