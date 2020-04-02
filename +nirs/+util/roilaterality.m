function tbl = roilaterality(Stats,ROI1,ROI2,flag)
% computes a laterality index defined as
% LI = (ROI1-ROI2)/(ROI1+ROI2)

if(nargin<4)
    flag=false;
end

Diff = nirs.util.roi_math(ROI1,'-',ROI2);
Sum = nirs.util.roi_math(ROI1,'+',ROI2);

if(flag)
    for(i=1:length(Stats))
        Stats(i).beta=abs(Stats(i).beta);
    end
end

tblD = nirs.util.roiAverage(Stats,Diff);
tblS = nirs.util.roiAverage(Stats,Sum);

tbl=tblD;
tbl.ROI=repmat(cellstr('LatIndex'),height(tbl),1); 
tbl.Beta=tblD.Beta./tblS.Beta;
tbl.SE = (abs(tblS.Beta).*tblD.SE+abs(tblD.Beta).*tblS.SE)./(tblS.Beta.^2);
tbl.T = tbl.Beta./tbl.SE;
tbl.DF=tblS.DF+tblD.DF;
tbl.p = 2*tcdf(-abs(tbl.T),tbl.DF);


tbl.q = nirs.math.fdr( tbl.p );
[~,tbl.power] = nirs.math.MDC(tbl,.8,.05);
return


