
[US,V]=nirs.math.hogSVD({JacobEEG JacobNIRS});

job=nirs.modules.OpticalDensity;
dOD=job.run(nir_raw);

L=[US{2}.hbo US{2}.hbr];
hbX=pinv(L)*dOD.data';
hbo=hbX(1:end/2,:);
hbr=hbX(end/2+1:end,:);

eeg=pinv(US{1}.eeg)*eeg_raw.data';

hbo_wht=nirs.math.innovations(hbo',10);
hbr_wht=nirs.math.innovations(hbr',10);
eeg_wht=nirs.math.innovations(eeg',10);

mintime=max(min(nir_raw.time),min(eeg_raw.time));
maxtime=min(max(nir_raw.time),max(eeg_raw.time));
ntime=sort(unique([nir_raw.time(find(nir_raw.time>=mintime & nir_raw.time<=maxtime));...
    eeg_raw.time(find(eeg_raw.time>=mintime & eeg_raw.time<=maxtime))]));

eeg_conv=[]; hbo_conv=[]; hbr_conv=[];
twin=2;
for i=1:size(eeg_wht,2)
    a=conv(eeg_wht(:,i),ones(fix(eeg_raw.Fs*twin),1),'same');
    eeg_conv(:,i)=interp1(eeg_raw.time,a,ntime);
end
for i=1:size(hbo_wht,2)
    a=conv(hbo_wht(:,i),ones(fix(nir_raw.Fs*twin),1),'same');
    hbo_conv(:,i)=interp1(nir_raw.time,a,ntime);
end
for i=1:size(hbr_wht,2)
    a=conv(hbr_wht(:,i),ones(fix(nir_raw.Fs*twin),1),'same');
    hbr_conv(:,i)=interp1(nir_raw.time,a,ntime);
end

d=[];
d(:,:,1)=resample(eeg_conv,1,40);
d(:,:,2)=resample(hbo_conv,1,40);
d(:,:,3)=resample(hbr_conv,1,40);

Fac=25;

Options(1)=1E-12;  % Convergence
Options(2)=0;  % Initialization
Options(3)=0;  %ploting
Options(4)=0;  %scaling
Options(5)=10;  % display update
Options(6)=200; %max iters

const = [0 1 0];  % L1/non-negative

OldLoad=[];  %initial guess of loading
FixMode=[];  
 
[ssX,Corco] = pftest(1,d,Fac,Options,const,OldLoad,FixMode);
model = parafac(d,Fac,Options,const,OldLoad,FixMode);
[A,B,C] = fac2let(model);

aa=zeros(size(d));
for i=1:25; for j=1:3; aa(:,:,j)=aa(:,:,j)+(A(:,i)*B(:,i)')*C(j,i); end; end;
