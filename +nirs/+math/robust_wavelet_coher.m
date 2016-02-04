function [r,p] = wavelet_coher(Y,FS,freq,wname,nscales)
% This function is a wrapper for the wcoher function and computes the
% complex coherence between two signals (a,b) over the frequency bin defined by freq
%
% inputs:
%   a- signal 1
%   b- signal 2
%   Fs-  sample rate
%   freq - range of frequencies to use in calc.  If not defined uses whole spectrum
%   waven-  wavelet name.  Default = real valued morlet


if(nargin<6)
    nscales=128;
end
if(nargin<5)
    %wname='db45';
    wname='morl';
end
if(nargin<3 || isempty(FS))
    FS=1;
end
if(nargin<4 || isempty(freq))
    freq = [Fs*2/length(Y) Fs/2];
end

Y=Y-ones(size(Y,1),1)*mean(Y,1);
Y=Y-ones(size(Y,1),1)*mean(Y,1);



% TODO: add code to deal with other wavelets, but for now lets use the
% continuious ones that matlab supports
f=scal2frq([1:256],wname,1/FS);

scales = find(f<.5*FS & f>min(freq) & f<max(freq) & f>3/(FS*size(Y,1)));
f=scal2frq(scales,wname,1/FS);

F=1./(scales');

wav=zeros(length(scales),size(Y,1),size(Y,2));
wav_sm=zeros(length(scales),size(Y,1),size(Y,2));

for i=1:size(Y,2)
   % -------- Wavelet transforms
    wav(:,:,i) = nirs.math.robust_cwt(Y(:,i),scales,wname);
    %-------- Wavelet smoothing
    wav_sm(:,:,i) = conv2(abs(squeeze(wav(:,:,i))).^2,F,'same');
end


r=eye(size(Y,2));
for i=1:size(Y,2)
    for j=i+1:size(Y,2)
        % -------- Cross wavelet 
        wavAB=squeeze(wav(:,:,i)).*conj(squeeze(wav(:,:,j)));
        wavAB_sm = conv2(abs(wavAB).^2,F,'same');

        %--------- Wavelet coherence
        Rsq=wavAB_sm./(squeeze(wav_sm(:,:,i)).*squeeze(wav_sm(:,:,j)));
        
        
        lst=find(f>=min(freq) & f<max(freq)); % get the band of interest
        Rsq=Rsq(lst,:);
        
        R=sqrt(Rsq).*sign(wavAB(lst,:));
        Z = .5*log((1+R)./(1-R));
        
        r(i,j)=tanh(mean(Z(:)));
        r(j,i)=r(i,j);
    end
end

% ????? how many DFE are in this model.  Its not the length(a)
dfe=size(Y,1);
t=r./sqrt((1-r.^2)/(dfe-2));
p=2*tcdf(-abs(t),max(dfe));
p=min(p+eye(size(p)),1);
