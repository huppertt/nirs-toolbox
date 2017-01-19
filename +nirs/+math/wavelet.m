function cfs = wavelet(d,FS,scales,wname,flag_normalize,flag_SMOOTH,NSW,NTW)
if(nargin<4 || isempty(wname))
    wname='morl';
end
if(nargin<3 || isempty(scales))
    scales=[2:128];
     f=scal2frq(scales,wname,1/FS);
     scales = find(f<.5*FS & f>3/(FS*size(d,1)));
end


if(nargin<5 || isempty(flag_normalize))
    flag_normalize=false;
end
if(nargin<6 || isempty(flag_SMOOTH))
    flag_SMOOTH=true;
end

if(nargin<7 || isempty(NSW))
    f=scal2frq(scales,wname,1/FS);
   NSW = fix(length(f)/10);
end
if(nargin<8 || isempty(NTW))
     NTW = fix(FS);
end





if(flag_normalize)
    maxcomp=30;
    n=round(size(d,1)/30);
    [X]=convmtx(d,n);
    [U,S,V]=nirs.math.mysvd(X(n+1:end-n,:));
    
    nn=min(find(cumsum(diag(S))/sum(diag(S))>.99));
    nn=min(nn,maxcomp);
    [icasig, A, W] = fastica(X(n+1:end-n,:)', 'numOfIC',nn);

    cfs=nirs.math.wavelet(icasig(:,1),FS,scales,wname,false,flag_SMOOTH,NSW,NTW);
    cfs=cfs/sqrt(sum(cfs(:).*conj(cfs(:))));
    for i=2:size(icasig,1);
        tmp=nirs.math.wavelet(icasig(:,1),FS,scales,wname,false,flag_SMOOTH,NSW,NTW);
        cfs=cfs+tmp/sqrt(sum(tmp(:).*conj(tmp(:))));
    end
    c=nirs.math.wavelet(d,FS,scales,wname,false,flag_SMOOTH,NSW,NTW);
    cfs=cfs/sqrt(sum(cfs(:).*conj(cfs(:))))*sqrt(sum(c(:).*conj(c(:))));
    return
else
    
    cfs(:,:) = cwt(d,scales,wname);
    if(flag_SMOOTH)
        cfs(:,:)    = smoothCFS(abs(squeeze(cfs(:,:))).^2,flag_SMOOTH,NSW,NTW);
        cfs(:,:)    = sqrt(cfs(:,:));
    end
end
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