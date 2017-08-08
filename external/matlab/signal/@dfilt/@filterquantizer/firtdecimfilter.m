function [y, zf, acc, phaseidx] = firtdecimfilter(q,M,p,x,zi,acc,phaseidx,nx,nchans,ny) %#ok<INUSL>
%FIRTDECIMFILTER   

%   Copyright 1988-2013 The MathWorks, Inc.

% Quantize input
x = quantizeinput(q,x);

if isreal(p)    
    if isreal(x) && isreal(zi) && isreal(acc),
        % Real output
        % Call C-mex
        [y,zf,acc,phaseidx] = firtdecimfilter(p,M,x,zi,acc,phaseidx,uint32(ny));        
    else
        % Complex output
        [y,zf,acc,phaseidx] = getcmplxoutput(p,M,x,zi,acc,phaseidx,ny);         
    end
else      
    % Separate real and imaginary parts of the coefficients
    preal = real(p);
    pimag = imag(p);   
    
    lclphaseidx = phaseidx+0;   % force new memory allocation
    lclzi = zi+0;               % force new memory allocation    
    
    
    % Filter input with real part of coefficients (output1)
    [y1,zf1,acc1,lclphaseidx] = getcmplxoutput(preal,M,x,lclzi,acc,lclphaseidx,ny); %#ok<NASGU>
       
    % Filter input with  imaginary part of coefficients (output2)
    [y2,zf2,acc2,phaseidx] = getcmplxoutput(pimag,M,x,zi,acc,phaseidx,ny);
    
    % Combine both outputs. Output2 should be multiplied by i. Cannot use
    % COMPLEX function since both outputs could be complex numbers by themselves
    
    y = y1+1i*y2;
    zf = zf1+1i*zf2;
    acc = acc1+1i*acc2;           
end

%--------------------------------------------------------------------------
function [y,zf,acc,phaseidx] = getcmplxoutput(p,M,x,zi,acc,phaseidx,ny)
% GETCMPLXOUTPUT Get complex output

% Complex output
lclphaseidx = phaseidx+0; % Force new memory allocation
% Real part
xr = real(x);
zir = real(zi);
accr = real(acc);
[yr,zir,accr,lclphaseidx] = firtdecimfilter(p,M,xr,zir,accr,lclphaseidx,uint32(ny)); %#ok<NASGU>

% Imaginary part
xi = imag(x);
zii = imag(zi);
acci = imag(acc);
[yi,zii,acci,phaseidx] = firtdecimfilter(p,M,xi,zii,acci,phaseidx,uint32(ny));

% y = complex(yr,yi);
% zf = complex(zir,zii);
% acc = complex(accr,acci);
y = yr + 1i*yi;
zf = zir + 1i*zii;
acc = accr + 1i*acci;


