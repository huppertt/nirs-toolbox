function [h,err,res] = firpmmex(order, ff, aa, varargin)
%FIRPMMEX   

%   Author(s): R. Losada
%   Copyright 2006 The MathWorks, Inc.

[nfilt,ff,grid,des,wt,ftype,sign_val,hilbert,neg] = firpminit(order, ff, aa, varargin{:});

[h,err,iext,ret_code,iters] = remezmex(nfilt,ff,grid,des,wt,ftype);
err = abs(err);
h = h * sign_val;

if (ret_code ~= 0) 
    handleErrorCode (ret_code, iters(1));
elseif any(err < 100*max(abs(des))*eps)
    warning(message('signal:firpmmex:machineAccuracy'));
end

%
% arrange 'results' structure
%
if nargout > 2
    res.fgrid = grid(:);
    res.H = freqz(h,1,res.fgrid*pi);
    if neg  % asymmetric impulse response
        linphase = exp(sqrt(-1)*(res.fgrid*pi*(order/2) - pi/2));
    else
        linphase = exp(sqrt(-1)*res.fgrid*pi*(order/2));
    end
    if hilbert == 1  % hilbert
        res.error = real(des(:) + res.H.*linphase);
    else
        res.error = real(des(:) - res.H.*linphase);
    end
    res.des = des(:);
    res.wt = wt(:);
    res.iextr = iext(1:end);
    res.fextr = grid(res.iextr);  % extremal frequencies
    res.fextr = res.fextr(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Handle error codes 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handleErrorCode (ret_code, iters)
if ret_code == -1
    error(message('signal:firpmmex:CantApproximate'));
elseif ret_code == -2     
    warning(message('signal:firpmmex:notConverging', iters, 'err'));
elseif ret_code == -5
    error(message('signal:firpmmex:SignalErr'));
end

% [EOF]
