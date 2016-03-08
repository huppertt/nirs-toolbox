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


if(nargin<4)
    %wname='db45';
    wname='morl';
end
if(nargin<2 || isempty(FS))
    FS=1;
end
if(nargin<3 || isempty(freq))
    freq = [Fs*2/length(Y) Fs/2];
end

Y=Y-ones(size(Y,1),1)*mean(Y,1);
Y=Y-ones(size(Y,1),1)*mean(Y,1);

% TODO: add code to deal with other wavelets, but for now lets use the
% continuious ones that matlab supports
scales=[2:256];
f=scal2frq(scales,wname,1/FS);
scales = find(f<.5*FS & f>min(freq) & f<max(freq) & f>3/(FS*size(Y,1)));
f=scal2frq(scales,wname,1/FS);

flag_SMOOTH = true;
NTW = fix(FS);
NSW = fix(length(find(f>=min(freq) & f<max(freq)))/5);

for i=1:size(Y,2)
    cfs_s1(:,:,i)    = nirs.math.robust_cwt(Y(:,i),scales,wname);
    cfs_s10(:,:,i)   = cfs_s1(:,:,i);
    cfs_s1(:,:,i)    = smoothCFS(abs(squeeze(cfs_s1(:,:,i))).^2,flag_SMOOTH,NSW,NTW);
    cfs_s1(:,:,i)    = sqrt(cfs_s1(:,:,i));
end

disp('Progress');
n=size(Y,2)*(size(Y,2)-1)/2;
cnt=1;
str='  0';
fprintf('%s %',str(end-2:end));

r=eye(size(Y,2));
lst=find(f>=min(freq) & f<max(freq)); % get the band of interest
for i=1:size(Y,2)
    for j=i+1:size(Y,2)
        str=['   ' num2str(round(100*cnt/n))];
        fprintf('\b\b\b\b%s %',str(end-2:end));
        % -------- Cross wavelet
        cfs_cross = conj(squeeze(cfs_s10(lst,:,i))).*squeeze(cfs_s10(lst,:,j));
        cfs_cross = smoothCFS(cfs_cross,flag_SMOOTH,NSW,NTW);
        WCOH      = cfs_cross./(squeeze(cfs_s1(lst,:,i)).*squeeze(cfs_s1(lst,:,j)));
        
        Z = .5*log((1+WCOH)./(1-WCOH));
        Z=max(min(Z,8),-8);  %avoid having Inf since one value will screw up everything
        r(i,j)=tanh(median(Z(:)));
        r(j,i)=r(i,j);
        cnt=cnt+1;
    end
end
disp('completed')
      

dist=montecarlo(50,length(Y),scales,wname,flag_SMOOTH,NSW,NTW);
Z=.5*log((1+r)./(1-r));
Z=max(min(Z,8),-8); 
p=1-dist.cdf(abs(Z));
p=min(p+eye(size(p)),1);


% 
% % ????? how many DFE are in this model.  Its not the length(a)
% dfe=size(Y,1);
% t=r./sqrt((1-r.^2)/(dfe-2));
% p=2*tcdf(-abs(t),max(dfe));
% p=min(p+eye(size(p)),1);


end


function dist=montecarlo(niter,sizeY,scales,wname,flag_SMOOTH,NSW,NTW);

cnt=1;  Znull=[];
for iter=1:niter
    Y=randn(sizeY,2);
    for i=1:size(Y,2)
        cfs_s1(:,:,i)    = cwt(Y(:,i),scales,wname);
        cfs_s10(:,:,i)   = cfs_s1(:,:,i);
        cfs_s1(:,:,i)    = smoothCFS(abs(squeeze(cfs_s1(:,:,i))).^2,flag_SMOOTH,NSW,NTW);
        cfs_s1(:,:,i)    = sqrt(cfs_s1(:,:,i));
    end
    for i=1:size(Y,2)
        for j=i+1:size(Y,2)
            % -------- Cross wavelet
            cfs_cross = conj(squeeze(cfs_s10(:,:,i))).*squeeze(cfs_s10(:,:,j));
            cfs_cross = smoothCFS(cfs_cross,flag_SMOOTH,NSW,NTW);
            WCOH      = cfs_cross./(squeeze(cfs_s1(:,:,i)).*squeeze(cfs_s1(:,:,j)));
            
            Z = .5*log((1+WCOH)./(1-WCOH));
            Z=max(min(Z,8),-8);  %avoid having Inf since one value will screw up everything
            Znull(cnt)=median(Z(:));
            cnt=cnt+1;
        end
    end
    
end
dist=fitdist(Znull','normal');

end



%----------------------------------------------------------------------
function CFS = smoothCFS(CFS,flag_SMOOTH,NSW,NTW)

if ~flag_SMOOTH , return; end
if ~isempty(NTW)
    len = NTW;
    F   = ones(1,len)/len;
    CFS = conv2(CFS,F,'same');
end
if ~isempty(NSW)
    len = NSW;
    F   = ones(1,len)/len;    
    CFS = conv2(CFS,F','same');
end
%----------------------------------------------------------------------
end


