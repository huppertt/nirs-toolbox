function [y,z,tapidx] = firinterpfilter(q,L,p,x,z,tapidx,nx,nchans,ny) %#ok<INUSL>
%FIRINTERPFILTER   Filtering method for fir interpolator.

%   Copyright 1999-2013 The MathWorks, Inc.

% Quantize input
x = quantizeinput(q,x);
if isreal(p)
    if isreal(x) && isreal(z),
        [y,z,tapidx] = firinterpfilter(p,L,x,uint32(ny),z,tapidx);
    else
         % Complex output
        [y,z,tapidx] = getcmplxoutput(p,L,x,z,tapidx,ny);
    end
else
    preal = real(p);
    pimag = imag(p);
    
    lcltapidx = tapidx+0;       % force new memory allocation
    lclz = z+0;                 % force new memory allocation
    
    % Filter input with real part of coefficients (output1)
    [y1,zf1,lcltapidx] = getcmplxoutput(preal,L,x,lclz,lcltapidx,ny); %#ok<NASGU>
       
    % Filter input with  imaginary part of coefficients (output2)
    [y2,~,tapidx] = getcmplxoutput(pimag,L,x,z,tapidx,ny);
    
    % Combine both outputs. Output2 should be multiplied by i. Cannot use
    % COMPLEX function since both outputs could be complex numbers by themselves    
    y = y1+1i*y2;
    z = zf1;   
    
end

%--------------------------------------------------------------------------
function  [y,z,tapidx] = getcmplxoutput(p,L,x,z,tapidx,ny)
% GETCMPLXOUTPUT Get complex output

% Filter real part of input with real initial conditions
% Copy tap index
lcltapidx = tapidx+0; % Force new memory allocation
zr = real(z);
[yr,zr,lcltapidx] = firinterpfilter(p,L,real(x),uint32(ny),zr,lcltapidx); %#ok<NASGU>

% Now filter imag part of input with imaginary initial conditions
zi = imag(z);
[yi,zi,tapidx] = firinterpfilter(p,L,imag(x),uint32(ny),zi,tapidx);
y = complex(yr,yi);
z = complex(zr,zi);

% [EOF]
