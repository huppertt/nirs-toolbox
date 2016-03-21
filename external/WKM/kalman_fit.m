function [CMRO2,CBF]=kalman_fit(nlgr,data)

t=get(data,'SamplingInstants');
dt=get(data,'Ts');
ntps=length(t);

x_aposteriori=ones(length(nlgr.InitialStates),ntps+1);
x_apriori=ones(length(nlgr.InitialStates),ntps);
for i=1:length(nlgr.InitialStates)
    x_aposteriori(i,1)=nlgr.InitialStates(i).Value;
    x_x_apriori(i,1)=nlgr.InitialStates(i).Value;
end

P_aposteriori=zeros(length(nlgr.InitialStates),length(nlgr.InitialStates),ntps+1);
P_apriori=zeros(length(nlgr.InitialStates),length(nlgr.InitialStates),ntps);

P_aposteriori(:,:,1)=1000*eye(length(nlgr.InitialStates));

y=get(data,'OutputData');
I=eye(length(nlgr.InitialStates));

dy=diff(y,1);
Q=(1.2536*eye(length(nlgr.InitialStates)).*diag(mad(dy(:),1))).^2;
Q(1,1)=0;
[yf]=nirs.math.innovations(y,20);
R=0*(1.2536*eye(2).*diag(max(mad(yf,1),1E-6))).^2;

% Take the param values from the model instead of the object in case I want
% to add my code to fit data with the model
for i=1:length(nlgr.Parameters)
    param{i}=nlgr.Parameters(i).Value;
end


for i=1:ntps
    dx = WKM_m(t, x_aposteriori(:,i),0, param{:});
    
    x_apriori(:,i+1)=x_aposteriori(:,i)+dx*dt;
    
    [dx, yhat,F,H] = WKM_m(t, x_apriori(:,i+1),0, param{:});
    yy(:,i)=yhat;
    P_apriori(:,:,i+1)=F*squeeze(P_aposteriori(:,:,i))*F'+Q;
    
    innov=y(i,:)-yhat';
    S=H*squeeze(P_apriori(:,:,i+1))*H'+R;
    K=squeeze(P_apriori(:,:,i+1))*H'*inv(S);
    
    x_aposteriori(:,i+1)=x_apriori(:,i+1)+K*innov';
    P_aposteriori(:,:,i+1)=(I-K*H)*P_apriori(:,:,i+1);
end


% Run the RTS smoother
x=ones(length(nlgr.InitialStates),ntps+1);
x(:,end)=x_aposteriori(:,end);
P=eye(length(nlgr.InitialStates));


for i=ntps:-1:1
    [~,~,F] = WKM_m(t, x(:,i),0, param{:});
    C=squeeze(P_aposteriori(:,:,i))*F'*pinv(squeeze(P_apriori(:,:,i+1)));
    x(:,i)=x_aposteriori(:,i)+C*(x(:,i+1)-x_apriori(:,i+1));
    P=squeeze(P_aposteriori(:,:,i))+C*(P-P_apriori(:,:,i+1))*C';
end

CBF=x(2,:);
OEF=x(4,:)./x(3,:);
CMRO2=OEF.*CBF;

CBF=CBF-1;
CMRO2=CMRO2-1;

