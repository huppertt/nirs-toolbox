function [CMRO2,CBF]=nlfit_WKM(nlgr,data)

t=get(data,'SamplingInstants');
dt=get(data,'Ts');
ntps=length(t);


y=get(data,'OutputData');

% Take the param values from the model instead of the object in case I want
% to add my code to fit data with the model
for i=1:4
    param{i}=nlgr.Parameters(i).Value;
end

x=ones(4,ntps+1);

for i=2:ntps
    dx = WKM_m(t, x(:,i-1),zeros(2,1), param{:});
    x(:,i)=x(:,i-1)+dx*dt;
    iter=1;
    while(1)
        [~, yhat,F,H] = WKM_m(t, x(:,i),zeros(2,1), param{:});
        yy(i,:)=yhat';
        innov=(y(i,:)'-yhat);
        dx=pinv(H'*H)*H'*innov;
        x(:,i)=x(:,i)+dx;
        if(norm(abs(dx))<1E-6 | iter>50)
            disp(i);
            disp(x(:,i))
            break
        end
        iter=iter+1;
    end
end

CBF=(x(3,1:ntps)-1);
CMRO2=(x(4,1:ntps)-1);