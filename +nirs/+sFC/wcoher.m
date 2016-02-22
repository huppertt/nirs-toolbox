function [r,p,dfe]=wcoher(data,freq,wname,robust_flag)

if(~isempty(strfind(class(data),'.core.Data')))
    Fs=data.Fs;
    data=data.data;
else
    Fs=1;
end

if(nargin<4)
    robust_flag=false;
end

if(nargin<3 || isempty(wname))
   wname='morl';
end

if(nargin<2 || isempty(freq))
   freq=[0 1]*Fs/2;
end

if(robust_flag)
    [r,p] = nirs.math.robust_wavelet_coher(data,Fs,freq,wname);
else
    [r,p] = nirs.math.wavelet_coher(data,Fs,freq,wname);
end

dfe=size(data,1);
