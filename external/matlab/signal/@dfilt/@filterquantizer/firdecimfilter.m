function [y, zf, acc, phaseidx, tapidx] = firdecimfilter(q,M,p,x,zi,acc,phaseidx,tapidx,nx,nchans,ny) %#ok<INUSL>
%FIRDECIMFILTER   

%   Copyright 1988-2013 The MathWorks, Inc.

% Quantize input
x = quantizeinput(q,x);

if isreal(p)   
    if isreal(x) && isreal(zi) && isreal(acc),
        % Real output                
        % Call C-mex
        [y,zf,acc,phaseidx,tapidx] = firdecimfilter(p,M,x,zi,acc,phaseidx,tapidx,uint32(ny));        
    else   
        % Complex output
        [y,zf,acc,phaseidx,tapidx] = getcmplxoutput(p,M,x,zi,acc,phaseidx,tapidx,ny);       
    end
else
    % Separate real and imaginary parts of the coefficients
    preal = real(p);
    pimag = imag(p);      
    
    lclphaseidx = phaseidx+0;   % force new memory allocation
    lcltapidx   = tapidx+0;     % force new memory allocation
    lclzi = zi+0;               % force new memory allocation
    
    % Filter input with real part of coefficients (output1)
    [y1,zf1,acc1,lclphaseidx,lcltapidx] = getcmplxoutput(preal,M,x,zi,acc,lclphaseidx,lcltapidx,ny);   %#ok<NASGU,ASGLU>
   
    % Filter input with  imaginary part of coefficients (output2)
    [y2,~,acc2,phaseidx,tapidx] = getcmplxoutput(pimag,M,x,lclzi,acc,phaseidx,tapidx,ny);    
    
    % Combine both outputs. Output2 should be multiplied by i. Cannot use
    % COMPLEX function since both outputs could be complex numbers by themselves
    
    y    = y1+1i*y2;
    zf   = zf1;
    acc  = acc1+1i*acc2;
    
end
          
%--------------------------------------------------------------------------
function [y,zf,acc,phaseidx,tapidx] = getcmplxoutput(p,M,x,zi,acc,phaseidx,tapidx,ny)
% GETCMPLXOUTPUT Get complex output

lclphaseidx = phaseidx+0; % Force new memory allocation
lcltapidx = tapidx+0; % Force new memory allocation
% Real part
xr = real(x);
zir = real(zi);
accr = real(acc);
[yr,zir,accr,lclphaseidx,lcltapidx] = firdecimfilter(p,M,xr,zir,accr,lclphaseidx,lcltapidx,uint32(ny)); %#ok<NASGU,ASGLU>

% Imaginary part
xi = imag(x);
zii = imag(zi);
acci = imag(acc);
[yi,zii,acci,phaseidx,tapidx] = firdecimfilter(p,M,xi,zii,acci,phaseidx,tapidx,uint32(ny));

y = complex(yr,yi);
zf = complex(zir,zii);
acc = complex(accr,acci);

