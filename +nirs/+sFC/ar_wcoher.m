function [r,p,dfe]=ar_wcoher(data,modelorder,freq,wname,robust_flag)

if(~isempty(strfind(class(data),'.core.Data')))
    Fs=data.Fs;
    data=data.data;
else
    Fs=1;
end

if(nargin<5)
    robust_flag=false;
end

if(nargin<4 || isempty(wname))
   wname='morl';
end

if(nargin<3 || isempty(freq))
   freq=[0 1]*Fs/2;
end

if(nargin<2 || isempty(modelorder))
    modelorder=4;
end


if(isa(data,'nirs.core.Data'))
    Fs=data.Fs;
    data=data.data;
else
   Fs=1;
end


if(~isreal(data))
    warning('code does not support masked data yet');
    data=real(data);
end


if(isstr(modelorder))
    p = Fs*str2num(modelorder(1:strfind(modelorder,'x')-1));
else
    p=modelorder;
end


[data,f] = nirs.math.innovations(data,p);
maxf=0;
for i=1:length(f)
    maxf=max(maxf,length(f{i}));
end
% Mask out boundary values
data(1:maxf,:) = [];

if(robust_flag)
    [r,p] = nirs.math.robust_wavelet_coher(data,Fs,freq,wname);
else
    [r,p] = nirs.math.wavelet_coher(data,Fs,freq,wname);
end

dfe=size(data,1)-2;
