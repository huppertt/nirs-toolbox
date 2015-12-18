function [r,p] = wavelet_coher(a,b,FS,freq,wname)
% This function is a wrapper for the wcoher function and computes the
% complex coherence between two signals (a,b) over the frequency bin defined by freq
%
% inputs:
%   a- signal 1
%   b- signal 2
%   Fs-  sample rate
%   freq - range of frequencies to use in calc.  If not defined uses whole spectrum
%   waven-  wavelet name.  Default = real valued morlet

if(nargin<5)
    %wname='db45';
    wname='morl';
end
if(nargin<3 || isempty(FS))
    FS=1;
end
if(nargin<4 || isempty(freq))
    freq = [Fs*2/length(a) Fs/2];
end

a=a-mean(a);
b=b-mean(b);


% TODO: add code to deal with other wavelets, but for now lets use the
% continuious ones that matlab supports
scales=[1:1024];


f=scal2frq(scales,wname,1/FS);

% -------- Wavelet transforms
wavA = cwt(a,scales,wname);
wavB = cwt(b,scales,wname);

%-------- Wavelet smoothing
F=1./(scales');
wavA_sm = conv2(abs(wavA).^2,F,'same');
wavB_sm = conv2(abs(wavB).^2,F,'same');

% -------- Cross wavelet 
wavAB=wavA.*conj(wavB);
wavAB_sm = conv2(abs(wavAB).^2,F,'same');

%--------- Wavelet coherence
Rsq=wavAB_sm./(wavA_sm.*wavB_sm);

lst=find(f>=min(freq) & f<max(freq)); % get the band of interest
f=f(lst);
Rsq=Rsq(lst,:);

R=sqrt(Rsq).*sign(wavAB(lst,:));
Z = .5*log((1+R)./(1-R));

rtime=tanh(mean(Z,1));
r=tanh(mean(Z(:)));

% ????? how many DFE are in this model.  Its not the length(a)
dfe=size(a,1);
t=r./sqrt((1-r.^2)/(dfe-2));
p=2*tcdf(-abs(t),max(dfe));



