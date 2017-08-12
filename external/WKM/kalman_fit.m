function [CMRO2,CBF]=kalman_fit(nlgr,data)
   
t=get(data,'SamplingInstants');
dt=get(data,'Ts');
ntps=length(t);

fileName=get(nlgr,'FileName');

x=ones(length(nlgr.InitialStates),ntps);
for i=1:length(nlgr.InitialStates)
    x(i,1)=nlgr.InitialStates(i).Value;
end

x_aposteriori =ones(length(nlgr.InitialStates),ntps);
x_apriori     =ones(length(nlgr.InitialStates),ntps);
x_apriori(:,1)=x(:,1);
x_aposteriori(:,1)=x(:,1);

P_aposteriori =zeros(length(nlgr.InitialStates),length(nlgr.InitialStates),ntps+1);
P_apriori     =zeros(length(nlgr.InitialStates),length(nlgr.InitialStates),ntps);

y=get(data,'OutputData');
I=eye(length(nlgr.InitialStates));

u=zeros(nlgr.Order.nu,ntps);

dy=diff(y,1);
Q=zeros(length(nlgr.InitialStates));
P=zeros(length(nlgr.InitialStates));
for i=1:nlgr.Order.nu
    Q(i,i)=1.2536*mad(dy(:),1).^2;
    P(i,i)=1000;
end

[yf]=nirs.math.innovations(y,20);
R=(1.2536*eye(2).*diag(max(mad(yf,1),1E-6))).^2;

% Take the param values from the model instead of the object in case I want
% to add my code to fit data with the model
for i=1:length(nlgr.Parameters)
    param{i}=nlgr.Parameters(i).Value;
end


for i=2:ntps
    dx = feval(fileName,t, x(:,i-1),u(:,i), param{:});
    x(:,i)=x(:,i-1)+dx*dt;
    
    [~, yhat,F,H] = feval(fileName,t, x(:,i),u(:,i), param{:});
    P=F*P*F'+Q;
   
    x_apriori(:,i)=x(:,i);
    P_apriori(:,:,i)=P;
   
    innov=y(i,:)-yhat';
    K=P*H'*pinv(H*P*H'+R);
    
    x(:,i)=x(:,i)+K*innov';
    P=(I-K*H)*P;
    
    x_aposteriori(:,i)=x(:,i);
    P_aposteriori(:,:,i)=P;
    
end


% Run the RTS smoother
x=zeros(length(nlgr.InitialStates),ntps);
x(:,end)=x_aposteriori(:,end);

for i=ntps-1:-1:1
    [~,~,F] = feval(fileName,t, x(:,i+1),u(:,i+1), param{:});
    C=squeeze(P_aposteriori(:,:,i))*F'*pinv(squeeze(P_apriori(:,:,i+1)));
    x(:,i)=x_aposteriori(:,i)+C*(x(:,i+1)-x_apriori(:,i+1));
    P=squeeze(P_aposteriori(:,:,i))+C*(P-P_apriori(:,:,i+1))*C';
end

for idx=1:length(nlgr.InitialStates)
    StateNames{idx}=nlgr.InitialStates(idx).Name;
end

i=find(ismember(StateNames,'CBF'));
CBF=x(i,:);

i=find(ismember(StateNames,'OEF'));
if(~isempty(i))
    CMRO2=x(i,:).*CBF;
else
    i=find(ismember(StateNames,'q'));
    i2=find(ismember(StateNames,'CBV'));
    CMRO2=x(i,:)./x(i2,:).*CBF;
end

CBF=CBF-1;
CMRO2=CMRO2-1;


i=find(ismember(StateNames,'q'));
i2=find(ismember(StateNames,'CBV'));

E0=nlgr.Parameters(end-1).Value;
HbT0=nlgr.Parameters(end).Value;

yhat=[];
yhat(1,:)=(x(i2,:)-1)*HbT0-(x(i,:)-1)*HbT0*E0;  %HbO2
yhat(2,:)=(x(i,:)-1)*HbT0*E0;  % Hb

yhat=yhat';

return;

