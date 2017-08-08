function varargout=ss2sos(A,B,C,D,varargin)
%SS2SOS State-space to second-order sections model conversion.
%   [SOS,G]=SS2SOS(A,B,C,D) finds a matrix SOS in second-order sections
%   form and a gain G which represent the same system as the one with
%   single-input, single-output state space matrices A, B, C, and D. 
%   The zeros and poles of the system A, B, C, D must be in complex
%   conjugate pairs. 
%
%   [SOS,G] = SS2SOS(A,B,C,D,IU) uses the IUth input of the multi-input,
%   single-output state space matrices A, B, C and D.
% 
%   SOS is an L by 6 matrix with the following structure:
%       SOS = [ b01 b11 b21  1 a11 a21 
%               b02 b12 b22  1 a12 a22
%               ...
%               b0L b1L b2L  1 a1L a2L ]
%
%   Each row of the SOS matrix describes a 2nd order transfer function:
%                 b0k +  b1k z^-1 +  b2k  z^-2
%       Hk(z) =  ----------------------------
%                  1 +  a1k z^-1 +  a2k  z^-2
%   where k is the row index.
%
%   G is a scalar which accounts for the overall gain of the system. If
%   G is not specified, the gain is embedded in the first section. 
%   The second order structure thus describes the system H(z) as:
%       H(z) = G*H1(z)*H2(z)*...*HL(z)
%
%   NOTE: Embedding the gain in the first section when scaling a
%   direct-form II structure is not recommended and may result in erratic
%   scaling. To avoid embedding the gain, use ss2sos with two outputs. 
%
%   SS2SOS(...,DIR_FLAG) specifies the ordering of the 2nd order
%   sections. If DIR_FLAG is equal to 'UP', the first row will contain
%   the poles closest to the origin, and the last row will contain the
%   poles closest to the unit circle. If DIR_FLAG is equal to 'DOWN', the
%   sections are ordered in the opposite direction. The zeros are always
%   paired with the poles closest to them. DIR_FLAG defaults to 'UP'.
%
%   SS2SOS(...DIR_FLAG,SCALE) specifies the desired scaling of the
%   gain and the numerator coefficients of all 2nd order sections. SCALE
%   can be either 'NONE', Inf or 2 which correspond to no scaling,
%   infinity norm scaling and 2-norm scaling respectively. SCALE defaults
%   to 'NONE'. The filter must be stable in order to scale in the 2-norm
%   or inf-norm sense.  Using infinity-norm scaling in conjunction with 'UP'
%   ordering will minimize the probability of overflow in the realization.
%   On the other hand, using 2-norm scaling in conjunction with 'DOWN' 
%   ordering will minimize the peak roundoff noise.
%
%   NOTE: Infinity-norm and 2-norm scaling are appropriate only for direct
%   form II structures. 
%
%   % Example:
%   %   Find a second-order section form of a Butterworth lowpass filter.
%
%   [A,B,C,D] = butter(5,0.2);  % Butterworth filter design
%   sos = ss2sos(A,B,C,D)       % Second-order sections form
%   fvtool(sos)                 % Visualize filter
%
%   See also ZP2SOS, SOS2ZP, SOS2TF, SOS2SS, tf2SOS, CPLXPAIR.

%   NOTE: restricted to real coefficient systems (poles  and zeros 
%             must be in conjugate pairs)

%   References:
%     [1] L. B. Jackson, DIGITAL FILTERS AND SIGNAL PROCESSING, 3rd Ed.
%              Kluwer Academic Publishers, 1996, Chapter 11.
%     [2] S.K. Mitra, DIGITAL SIGNAL PROCESSING. A Computer Based Approach.
%              McGraw-Hill, 1998, Chapter 9.
%     [3] P.P. Vaidyanathan. ROBUST DIGITAL FILTER STRUCTURES. Ch 7 in
%              HANDBOOK FOR DIGITAL SIGNAL PROCESSING. S.K. Mitra and J.F.
%              Kaiser Eds. Wiley-Interscience, N.Y.

%   Author(s): R. Losada 
%   Copyright 1988-2002 The MathWorks, Inc.

error(nargchk(4,7,nargin,'struct'))

if nargin == 4,
   varargin = {};
end

[IU,dir_flag,scale] = parse_options(varargin{:});

 
if IU > size(D,2), % A, B and C can be empty, if the system is simply a gain given by D
   error(message('signal:ss2sos:InvalidParam', size( D, 2 )));
end

% Find Poles and Zeros
[z,p,k] = ss2zp(A,B,C,D,IU);

[varargout{1:max(1,nargout)}] = zp2sos(z,p,k,dir_flag,scale);

%------------------------------------------------------------------------------------
function [IU,dir_flag,scale] = parse_options(varargin)
%PARSE_OPTIONS   Parse optional args of SS2SOS.

switch nargin,
case 0,
   IU = 1;
   dir_flag = 'up';
   scale = 'none';
case 1,
   if isnumeric(varargin{1}),
      IU = varargin{1};
      dir_flag = 'up';
      scale = 'none';
   else
      IU = 1;
      dir_flag = varargin{1};
      scale = 'none';
   end
case 2,
   if isnumeric(varargin{1}),
      IU = varargin{1};
      dir_flag = varargin{2};
      scale = 'none';
   else
      IU = 1;
      dir_flag = varargin{1};
      scale = varargin{2};
   end
case 3,
   IU = varargin{1};
   dir_flag = varargin{2};
   scale = varargin{3};
end
