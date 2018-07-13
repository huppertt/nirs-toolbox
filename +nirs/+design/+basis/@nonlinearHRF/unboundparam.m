function pIn=unboundparam(obj,pOut,LB,UB);


if(size(pOut,1)~=length(LB))
    pOut=pOut';
end

for idx=1:size(pOut,1)
    if(~isnan(LB(idx)) & ~isnan(UB(idx)))
        %Sigmoid shape
        md=(LB(idx)+UB(idx))/2;
        range=(UB(idx)-LB(idx));
        rho=.05;
         for idx2=1:size(pOut,2)
             if(pOut(idx,idx2)<md)
                pIn(idx,idx2)=1/rho*log((pOut(idx,idx2)-LB(idx))/(range/2))+md;
               %pOut(idx,idx2)=LB+range/2*exp(rho*(pIn(idx,idx2)-md));
             else
                 pIn(idx,idx2)=-1/rho*log(1-(pOut(idx,idx2)-md)/(range/2))+md;
                %pOut(idx,idx2)=md+range/2*(1-exp(-rho*(pIn(idx,idx2)-md)));
             end
         end
    elseif(~isnan(LB(idx)) & isnan(UB(idx))) 
        %Half-sigmoid, above midpt slope=1;
         md=LB(idx)+abs(LB(idx)*2);
        range=max(abs(LB(idx)),1);
        rho=.05;
        a = 1/rho*log(2/rho/range);
         for idx2=1:size(pOut,2)
             if(pOut(idx,idx2)<LB(idx)+range/2*exp(rho*a))
                 pIn(idx,idx2)=1/rho*log((pOut(idx,idx2)-LB(idx))/(range/2))+md-a;
%                 pOut(idx,idx2)=LB+range/2*exp(rho*(pIn(idx,idx2)-md+a));
             else
                 pIn(idx,idx2)=pOut(idx,idx2)+LB(idx)-range/2*exp(rho*a);
%                 pOut(idx,idx2)=pIn(idx,idx2)-LB+range/2*exp(rho*a);
             end
         end
    elseif(isnan(LB(idx)) & ~isnan(UB(idx)))
        %Half-sigmoid, below midpt slope=1;
        md=UB(idx)-abs(UB(idx)/2);
        range=abs(UB(idx));
        rho=.05;
        a = -1/rho*log(2/rho/range);
         for idx2=1:size(pOut,2)
             if(pOut(idx,idx2)<md+range/2*(1-exp(-rho*a)))
                 pIn(idx,idx2)=pOut(idx,idx2)-range/2*(1-exp(-rho*a));
%                 pOut(idx,idx2)=pIn(idx,idx2)+range/2*(1-exp(-rho*a));
             else
               pIn(idx,idx2)=-1/rho*log(1-(pOut(idx,idx2)-md)/(range/2))+md-a;
                 
%                 pOut(idx,idx2)=md+range/2*(1-exp(-rho*(pIn(idx,idx2)-md+a)));
             end
         end
    else
        pIn(idx,:)=pOut(idx,:);
    end
end

