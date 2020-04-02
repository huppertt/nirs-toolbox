function [plv,p] = phaselocking(Y,FS,wname,freq)
% This function is a wrapper for the wcoher function and computes the
% phaselocking between two signals (a,b) over the frequency bin defined by freq
%
% inputs:
%   a- signal 1
%   b- signal 2
%   Fs-  sample rate
%   freq - range of frequencies to use in calc.  If not defined uses whole spectrum
%   waven-  wavelet name.  Default = real valued morlet


if(nargin<3 || isempty(wname))
    %wname='db45';
    wname='hilbert';
end
if(nargin<2 || isempty(FS))
    FS=1;
end
if(nargin<4 || isempty(freq))
    freq = [FS*2/length(Y) FS/2];
end

Y=Y-ones(size(Y,1),1)*mean(Y,1);
Y=Y-ones(size(Y,1),1)*mean(Y,1);

ang=compute_angle(Y,wname,FS,freq);

plv=zeros(size(Y,2),size(Y,2),size(ang,3));
for i=1:size(Y,2)
    plv(i,i,:)=1;
    for j=i+1:size(Y,2)
        plv(i,j,:)=squeeze(abs(mean(exp(1i*(ang(:,i,:)-ang(:,j,:))),1)));
        plv(j,i,:)=plv(i,j,:);
    end
end

% run Monte Carlo
niter=5;

Znull=zeros(niter,1,size(ang,3));
for iter=1:niter
    Yrnd=randn(size(Y,1),2);
    ang=compute_angle(Yrnd,wname,FS,freq);
    Znull(iter,1,:)=abs(mean(exp(1i*(ang(:,1,:)-ang(:,2,:))),1));
end
Znull=5*log((1+Znull)./(1-Znull));
Znull=max(min(Znull,8),-8); 
Znull=mean(Znull,3);

dist=fitdist(Znull,'normal');

Z=.5*log((1+plv)./(1-plv));
Z=max(min(Z,8),-8); 
Z=mean(Z,3);
plv=tanh(Z);

p=1-dist.cdf(abs(Z));
p=min(p+eye(size(p)),1);

end

function ang=compute_angle(Y,wname,FS,freq)

if(~strcmp(wname,'hilbert'))
    
    % TODO: add code to deal with other wavelets, but for now lets use the
    % continuious ones that matlab supports
    scales=[2:256];
    f=scal2frq(scales,wname,1/FS);
    scales = find(f<.5*FS & f>min(freq) & f<max(freq) & f>3/(FS*size(Y,1)));
    f=scal2frq(scales,wname,1/FS);
    
    for i=1:size(Y,2)
        cfs_s1(:,:,i)    = cwt(Y(:,i),scales,wname);     
    end
    cfs_s1=permute(cfs_s1,[2 3 1]);
else
    for i=1:size(Y,2)
        cfs_s1(:,i,1)=hilbert(Y(:,i));
    end
end

ang=unwrap(angle(cfs_s1));

end



