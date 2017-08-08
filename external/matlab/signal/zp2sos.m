function [sos,g] = zp2sos(varargin)
%ZP2SOS Zero-pole-gain to second-order sections model conversion.
%   [SOS,G] = ZP2SOS(Z,P,K) finds a matrix SOS in second-order sections 
%   form and a gain G which represent the same system H(z) as the one
%   with zeros in vector Z, poles in vector P and gain in scalar K.
%   The poles and zeros must be in complex conjugate pairs. 
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
%   scaling. To avoid embedding the gain, use zp2sos with two outputs. 
%
%   ZP2SOS(Z,P,K,DIR_FLAG) specifies the ordering of the 2nd order
%   sections. If DIR_FLAG is equal to 'UP', the first row will contain
%   the poles closest to the origin, and the last row will contain the
%   poles closest to the unit circle. If DIR_FLAG is equal to 'DOWN', the
%   sections are ordered in the opposite direction. The zeros are always
%   paired with the poles closest to them. DIR_FLAG defaults to 'UP'.
%
%   ZP2SOS(Z,P,K,DIR_FLAG,SCALE) specifies the desired scaling of the gain
%   and the numerator coefficients of all 2nd order sections. SCALE can be
%   either 'NONE', Inf or 2 which correspond to no scaling, infinity
%   norm scaling and 2-norm scaling respectively. SCALE defaults to 'NONE'.
%   The filter must be stable in order to scale in the 2-norm or inf-norm sense.
%   Using infinity-norm scaling in conjunction with 'UP' ordering will
%   minimize the probability of overflow in the realization. On the other
%   hand, using 2-norm scaling in conjunction with 'DOWN' ordering will
%   minimize the peak roundoff noise.
%
%   ZP2SOS(Z,P,K,DIR_FLAG,SCALE,KRZFLAG) specifies whether or not to keep
%   real zeros that are the negative of each other together rather than
%   ordering according to their proximity to poles. If KRZFLAG is true,
%   this is done and the result is a numerator with a middle coefficient
%   equal to zero. The default is false.
%
%   NOTE: Infinity-norm and 2-norm scaling are appropriate only for direct
%   form II structures. 
%
%   % Example:
%   %    Find a second-order section form of a Butterworth lowpass filter.
%
%   [z,p,k] = butter(5,0.2);    % 5th order Butterworth filter (Wn =0.2)
%   sos = zp2sos(z,p,k)         % Obtain second-order sections form
%   fvtool(sos)                 % Visualize the filter
%
%   See also TF2SOS, SOS2ZP, SOS2TF, SOS2SS, SS2SOS, CPLXPAIR.

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

%   Copyright 1988-2013 The MathWorks, Inc.

narginchk(2,6)

[z,p,k,direction_flag,scale,krzflag] = parseinput(varargin{:});

   
% Order the poles and zeros in complex conj. pairs and return the number of
% poles and zeros
[z,p,lz,lp,L] = orderzp(z,p);

% Break up conjugate pairs and real poles and zeros
[z_conj,z_real,p_conj,p_real] = breakconjrealzp(z,lz,p,lp);


% Order poles according to proximity to unit circle
[new_p,p_conj,p_real] = orderp(p_conj,p_real);


% Order zeros according to proximity to pole pairs
new_z = orderz(z_conj,p_conj,z_real,p_real,krzflag);

% Form SOS matrix
sos = formsos(lz,lp,new_z,new_p,L);

% Change direction if requested
if strcmp(direction_flag,'down'),
   sos = flipud(sos);
end

% At this point no scaling has been performed
% The leading coefficients of both num and den are one.

% Perform appropriate scaling
[sos,k] = scaling(sos,k,L,scale);


% Embed the gain if only one output argument was specified
if nargout == 1,
   sos(1,1:3) = k*sos(1,1:3);   
else
   g = k;
end

%------------------------------------------------------------------------------------
function [z,p,k,direction_flag,scale,krzflag] = parseinput(varargin)

z = varargin{1};
p = varargin{2};

if all(size(z) > 1) || all(size(p) > 1)
  error(message('signal:zp2sos:InvalidDimiensions'))
end

% Make sure z,p inputs are converted to columns
z = z(:);
p = p(:);

% Setup default values
k = 1; 
direction_flag = 'up';
scale = 'none';
krzflag = false;

% Replace with given values
if nargin > 2,
   if ~isempty(varargin{3}),
      k = varargin{3};
   end
   if nargin > 3,
      if ~isempty(varargin{4}),
         direction_flag = varargin{4};
      end
      if nargin > 4,
         if ~isempty(varargin{5}),
            scale = varargin{5};
         end
         if nargin > 5,
             if ~isempty(varargin{6}),
                 krzflag = varargin{6};
             end
         end
      end
   end
end

% Input check
diropts = {'up','down'};
indx1 = find(strncmpi(direction_flag, diropts, length(direction_flag)));
if isempty(indx1),
   error(message('signal:zp2sos:DirFlag', 'DIR_FLAG', 'UP', 'DOWN'));
end
direction_flag = diropts{indx1};

% Scale must be one of 'none', 2 or inf. For backwards compatibility we allow 'two'
% and 'inf' as well.
if ischar(scale),
    scaleopts = {'none','inf','two'};
    indx2 = find(strncmpi(scale, scaleopts, length(scale)));
    if isempty(indx2),
        error(message('signal:zp2sos:BadScale', 'SCALE', 'NONE', 'inf', 2));
    end
    scale = scaleopts{indx2};
    
    if strcmpi(scale,'two'),
        scale = 2;
    end
    
    if strcmpi(scale,'inf'),
        scale = inf;
    end
    
else
    if ~(isinf(scale) || scale == 2),
        error(message('signal:zp2sos:BadScale', 'SCALE', 'NONE', 'inf', 2));
    end
end

if ~islogical(krzflag),
    error(message('signal:zp2sos:BadLogical', 'KRZFLAG'));
end
    

%------------------------------------------------------------------------------------
function [z,p,lz,lp,L] = orderzp(z,p)
% Order the poles and zeros in complex conj. pairs and return the number of poles and
% zeros

% Order the poles and zeros in complex conj. pairs
z = cplxpair(z);
p = cplxpair(p);

% Get the number of poles and zeros
lz = length(z);
lp = length(p);
L = ceil(lp/2);

if lz > lp,
   error(message('signal:zp2sos:TooManyZeros'));
end



%------------------------------------------------------------------------------------
function [z_conj,z_real,p_conj,p_real] = breakconjrealzp(z,lz,p,lp)
% break up conjugate pairs and real poles and zeros

ind = find(abs(imag(p))>0);
p_conj = p(ind);   % the poles that have conjugate pairs
ind_complement = 1:lp;
if ~isempty(ind)
    ind_complement(ind) = [];
end
p_real = p(ind_complement);    % the poles that are real

% break up conjugate pairs and real zeros
ind = find(abs(imag(z))>0);
z_conj = z(ind);   % the zeros that have conjugate pairs
ind_complement = 1:lz;
if ~isempty(ind)
    ind_complement(ind) = [];
end
z_real = z(ind_complement);    % the zeros that are real

%------------------------------------------------------------------------------------
function [new_p,p_conj,p_real] = orderp(p_conj,p_real)
% Order poles according to proximity to unit circle

% order the conjugate pole pairs according to proximity to unit circle
[temp,ind] = sort(abs(p_conj - exp(1i*angle(p_conj)))); %#ok
p_conj = p_conj(ind);

% order the real poles according to proximity to unit circle too
[temp,ind] = sort(abs(p_real - sign(p_real))); %#ok
p_real = p_real(ind);

% Save the ordered poles
new_p = [p_conj;p_real];

%------------------------------------------------------------------------------------
function new_z = orderz(z_conj,p_conj,z_real,p_real,krzflag)
% Order zeros according to proximity to pole pairs

% order the conjugate zero pairs according to proximity to pole pairs
new_z = [];
for i = 1:length(z_conj)/2,
    if ~isempty(p_conj),
        [temp,ind1] = min(abs(z_conj-p_conj(1))); %#ok
        new_z = [new_z; z_conj(ind1)]; %#ok<*AGROW>
        z_conj(ind1) = [];
        [temp,ind2] = min(abs(z_conj-p_conj(2))); %#ok
        new_z = [new_z; z_conj(ind2)];
        z_conj(ind2) = [];
        p_conj([1 2]) = [];
    elseif ~isempty(p_real),
        [temp,ind] = min(abs(z_conj-p_real(1))); %#ok
        new_z = [new_z; z_conj(ind); z_conj(ind+1)];
        z_conj([ind ind+1]) = [];
        p_real(1) = [];
        if ~isempty(p_real),
            p_real(1) = [];
        end
    else
        new_z = [new_z; z_conj];
        break
    end
end

% If KRZFLAG is true, keep real zeros along with their negatives if present
if krzflag && ~isempty(z_real)
  remaining_z = [];
  while ~isempty(z_real)    
      realindx = find(z_real == -z_real(1));
      if ~isempty(realindx) && z_real(1) ~= 0 
        new_z = [new_z; z_real(1); z_real(realindx(1))];
        z_real(realindx(1)) = [];
        z_real(1) = [];
      else
        remaining_z = [remaining_z z_real(1)];        
        z_real(1) = [];        
      end
  end
else
  remaining_z = z_real;
end

% order remaining real zeros according to proximity to pole pairs too
for i = 1:length(remaining_z),
    if ~isempty(p_conj),
        [temp,ind] = min(abs(remaining_z-p_conj(1))); %#ok
        new_z = [new_z; remaining_z(ind)];
        remaining_z(ind) = [];
        p_conj(1) = [];
    elseif ~isempty(p_real),
        [temp,ind] = min(abs(remaining_z-p_real(1))); %#ok
        new_z = [new_z; remaining_z(ind)];
        remaining_z(ind) = [];
        p_real(1) = [];
    else
        new_z = [new_z; remaining_z];
        break
    end
end

%------------------------------------------------------------------------------------
function sos = formsos(lz,lp,new_z,new_p,L)
% Form SOS matrix

% Initialize SOS matrix
sos = [];

if lz == 0,
   if lp == 0,
      sos = [1 0 0 1 0 0];
   elseif ~rem(lp,2),
      sos = sosfun(1,2*L-1,new_p,sos); %even number of poles
   else
      sos = sosfun(1,2*(L-1)-1,new_p,sos); %odd number of poles
      sos = last_pole([],new_p(end),sos);%handle the last pole separately
   end   
else
   if ~rem(lz,2),
      sos = sosfun2(1,lz-1,new_z,new_p,sos); %even number of zeros
      % Now continue for the excess poles if any
      if ~rem(lp,2),
         sos = sosfun(lz+1,lp,new_p,sos); %even number of poles
      else
         sos = sosfun(lz+1,lp-1,new_p,sos); %odd number of poles
         sos = last_pole([],new_p(end),sos);%handle the last pole separately
      end
   else
      sos = sosfun2(1,lz-1,new_z,new_p,sos); %odd number of zeros
      %handle the last zero separately
      if lz == lp, %if number of poles = number of zeros
         sos = last_pole(new_z(lz),new_p(lz),sos);%handle the last pole separately
      else % more poles than zeros      
         [num,den] = zp2tf(new_z(lz),new_p(lz:lz+1),1);
         sos = [num den;sos];
         % Now continue for the excess poles if any
         if ~rem(lp,2),
            sos = sosfun(lz+2,lp,new_p,sos); %even number of poles
         else
            sos = sosfun(lz+2,lp-1,new_p,sos); %odd number of poles
            sos = last_pole([],new_p(end),sos);%handle the last pole separately
         end
      end
   end
end

%------------------------------------------------------------------------------------
function sos = last_pole(z,p,sos)
% Handle the last pole separately
[num,den] = zp2tf(z,p,1); 
num = [num 0];
den = [den 0];
sos = [num den;sos];

%------------------------------------------------------------------------------------
function sos = sosfun(start,stop,p,sos)
% This function was made to not repeat code in several places
for m = start:2:stop, 
   [num,den] = zp2tf([],p(m:m+1),1);
   sos = [num den;sos];
end

%------------------------------------------------------------------------------------
function sos = sosfun2(start,stop,z,p,sos)
% This function was made to not repeat code in several places
for m = start:2:stop, 
   [num,den] = zp2tf(z(m:m+1),p(m:m+1),1);
   sos = [num den;sos];
end

%------------------------------------------------------------------------------------
function [sos,g] = scaling(sos,k,L,scale) 
% SCALING, scale the cascaded sos filters using the infinity norm
% or the 2 norm as specified.

if strcmpi(scale,'none'),
    g = k;
    return
end

% Find the L scaling factors s(m) and perform the scaling
Fnum = 1;
Fden = 1;
den = sos(1,4:6);

% Compute the first scaling constant separetely
s(1) = filternorm(1,den,scale);

% Now compute the other scaling constants
for m = 2:L
   den = sos(m,4:6);
   Fnum = conv(Fnum,sos(m-1,1:3));
   Fden = conv(Fden,sos(m-1,4:6));
   Fden2 = conv(Fden,den);
   
   s(m) = filternorm(Fnum,Fden2,scale);

   % Now perform the scaling for the first L-1 sos sections
   sos(m-1,1:3) = s(m-1)./s(m).*sos(m-1,1:3);
end

% And now scale the last section
sos(end,1:3) = k*s(end)*sos(end,1:3);
g = 1/s(1);

% [EOF]
