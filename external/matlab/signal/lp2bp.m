function [at,bt,ct,dt] = lp2bp(a,b,c,d,wo,bw)
%LP2BP Lowpass to bandpass analog filter transformation.
%   [NUMT,DENT] = LP2BP(NUM,DEN,Wo,Bw) transforms the lowpass filter
%   prototype NUM(s)/DEN(s) with unity cutoff frequency to a
%   bandpass filter with center frequency Wo and bandwidth Bw.
%   [AT,BT,CT,DT] = LP2BP(A,B,C,D,Wo,Bw) does the same when the
%   filter is described in state-space form.
%
%   % Example:
%   %   Design a Elliptic analog lowpass filter with 5dB of ripple in the
%   %   passband and a stopband 90 decibels down. Transform this filter to
%   %   a bandpass filter with center frequency 24Hz and Bandwidth of 10Hz.
%
%   Wo= 24; Bw=10                   % Define Center frequency and Bandwidth
%   [z,p,k]=ellipap(6,5,90);        % Lowpass filter prototype
%   [b,a]=zp2tf(z,p,k);             % Specify filter in polynomial form
%   [num,den]=lp2bp(b,a,Wo,Bw);     % Convert LPF to BPF
%   freqs(num,den)                  % Frequency response of analog filter
%
%   See also BILINEAR, IMPINVAR, LP2LP, LP2BS and LP2HP

%   Author(s): J.N. Little and G.F. Franklin, 8-4-87
%   Copyright 1988-2013 The MathWorks, Inc.


if nargin == 4		% Transfer function case
    % handle column vector inputs: convert to rows
    if size(a,2) == 1
        a = a(:).';
    end
    if size(b,2) == 1
        b = b(:).';
    end
    
    % Cast to enforce precision rules
    wo = signal.internal.sigcasttofloat(c,'double','lp2bp','Wo',...
      'allownumeric');
    bw = signal.internal.sigcasttofloat(d,'double','lp2bp','Bw',...
      'allownumeric');
    if any([signal.internal.sigcheckfloattype(b,'single','lp2bp','DEN')...
        signal.internal.sigcheckfloattype(a,'single','lp2bp','NUM')])
      b = single(b);
      a = single(a);
    end
    % Transform to state-space
    [a,b,c,d] = tf2ss(a,b);

elseif (nargin == 6)
  % Cast to enforce precision rules
  if any([signal.internal.sigcheckfloattype(a,'single','lp2bp','A')...
      signal.internal.sigcheckfloattype(b,'single','lp2bp','B')...
      signal.internal.sigcheckfloattype(c,'single','lp2bp','C')...
      signal.internal.sigcheckfloattype(d,'single','lp2bp','D')])
    a = single(a);
    b = single(b);   
    c = single(c);
    d = single(d);
  end
  wo = signal.internal.sigcasttofloat(wo,'double','lp2bp','Wo',...
    'allownumeric');
  bw = signal.internal.sigcasttofloat(bw,'double','lp2bp','Bw',...
    'allownumeric');
else
  error(message('signal:lp2bp:MustHaveNInputs'));
end

error(abcdchk(a,b,c,d));

nb = size(b, 2);
[mc,ma] = size(c);

% Transform lowpass to bandpass
q = wo/bw;
at = wo*[a/q eye(ma); -eye(ma) zeros(ma)];
bt = wo*[b/q; zeros(ma,nb)];
ct = [c zeros(mc,ma)];
dt = d;

if nargin == 4		% Transfer function case
    % Transform back to transfer function
    zinf = ltipack.getTolerance('infzero',true);
    [z,k] = ltipack.sszero(at,bt,ct,dt,[],zinf);
    num = k * poly(z);
    den = poly(at);
    at = num;
    bt = den;
end
