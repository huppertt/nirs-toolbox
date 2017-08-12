numsubjects = 5;
numchannels = 12;
numtimepts = 3000;
Fs = 200;

% let's simulate some data
corefreq =[4 8 12 30];
BW=.1;  % bandwidth of the jitter.  0=pure sine wave
spatial=randn(numchannels,length(corefreq));
for i=1:numsubjects
    Ymeg{i}=zeros(numchannels,numtimepts);
    for j=1:length(corefreq)
        a = randn(numtimepts,1);
        [fa,fb]=butter(4,.01);
        a=filtfilt(fa,fb,a);
        a=a-mean(a);
        a=a./std(a);
        a=sin((corefreq(j)+BW*a).*[1:numtimepts]'/Fs*2*pi);
        Ymeg{i}=Ymeg{i} + spatial(:,j)*a';
    end
    Ymeg{i}=Ymeg{i}+randn(numchannels,numtimepts);
end

% Now do the wavelet
a0 = 2^(1/32);
scales = a0.^(1:7*32);
f=scal2frq(scales,'cmor1-1.5',1/Fs);
scales(find(f>Fs/2))=[];


MEG = zeros(length(scales),numtimepts,numchannels,numsubjects);
for i=1:numsubjects
    for j=1:numchannels
        [cfs,sc,frequencies] = cwt(Ymeg{i}(j,:),scales,'cmor1-1.5',1/Fs,'scal');
        MEG(:,:,j,i)=sc;
        fprintf(1,'.');
    end
    disp('');
    disp(['Subject ' num2str(i) ' completed']);
end
MEG=permute(MEG,[1 3 2 4]);
MEG = reshape(MEG,size(MEG,1),size(MEG,2),size(MEG,3)*size(MEG,4));

[ssX,Corco] = pftest(3,MEG,5);

model = parafac(MEG,3);