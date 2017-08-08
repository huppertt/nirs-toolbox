function [y,z,tapidx] = firsrcfilter(q,L,M,p,x,z,tapidx,im,inOffset,Mx,Nx,My) %#ok<INUSL>
%FIRINTERPFILTER   Filtering method for fir interpolator.

%   Copyright 1999-2013 The MathWorks, Inc.

% Quantize input
x = quantizeinput(q,x);

if isreal(p)
    if isreal(x) && isreal(z),        
        [y,z,tapidx] = firsrcfilter(p,L,M,x,uint32(My),z,tapidx,im,inOffset);
    else        
        [y,z,tapidx] = getcmplxoutput(p,L,M,x,My,z,tapidx,im,inOffset);
    end
else
    preal = real(p);
    pimag = imag(p);
    
    lcltapidx = tapidx+0;       % force new memory allocation
    lclz = z+0;                 % force new memory allocation
    
    % Filter input with real part of coefficients (output1)
    [y1,z1,lcltapidx] = getcmplxoutput(preal,L,M,x,My,lclz,lcltapidx,im,inOffset);  %#ok<NASGU>
    
    % Filter input with imaginary part of coefficients (output1)
    [y2,~,tapidx]    = getcmplxoutput(pimag,L,M,x,My,z,tapidx,im,inOffset);
    
    % Combine both outputs. Output2 should be multiplied by i. Cannot use
    % COMPLEX function since both outputs may be complex numbers by themselves    
    y = y1+1i*y2;
    z = z1;      
end

% ------------------------------------------------------------------------
function  [y,z,tapidx] = getcmplxoutput(p,L,M,x,My,z,tapidx,im,inOffset)
% GETCMPLXOUTPUT Get complex output

% Filter real part of input with real initial conditions
% Copy tap index
tapidxc = tapidx+0; % Force new memory allocation
zr = real(z);
[yr,zr,tapidx] = firsrcfilter(p,L,M,real(x),uint32(My),zr,tapidx,im,inOffset);

% Now filter imag part of input with imaginary initial conditions
zi = imag(z);
[yi,zi,tapidxc] = firsrcfilter(p,L,M,imag(x),uint32(My),zi,tapidxc,im,inOffset); %#ok<NASGU>
y = complex(yr,yi);
z = complex(zr,zi);
    
% [EOF]
