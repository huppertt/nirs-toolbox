% This demo is a work in progress.  I haven't moved the ROC methods into
% the same testing framework as the other code.  The draw code is still in
% rough for the connectivity stats objects.  I also want to implement the
% joint testing (HbO/Hb) as part of this since currently the models always
% include the HbO-Hb cross terms in the multi-variate model.  The model is
% still a bit sensitive to motion-artifacts (which should be solved once I
% include a robust regression iteration within the Grangers model)
% Finally, this doesn't handle the 0th lag (correlation) case under the Granger's model
% yet (to do so is not typical of Granger's anyway, but I think it might
% needed for fNIRS).  

truth=[];
for i=1:2;
    [data(i),truth] = nirs.testing.simData_connectivity([],truth);
end

j = nirs.modules.OpticalDensity();
dOD=j.run(data);

%This runs the mult-variate version of Grangers
j = nirs.modules.Connectivity();
ConnStats=j.run(dOD);

j=nirs.modules.MixedEffectsConnectivity();
GroupStats = j.run(ConnStats); 



a=(sum(truth(:,:,2:end),3)>0)*1;

lst=find(a==1);
lstN=find(a~=1);

% Run the pairwise version for comparision
g3=[];

for idx=1:length(data)
    d=dOD(idx).data;
    n=size(d,2);
    gr=zeros(n,n);
    for i=1:n;
        for j=1:n
            if(i~=j)
                gr(i,j)=nirs.math.grangers(d(:,i),d(:,j),5);
            end
        end
    end
    g3=[g3; gr(lst); gr(lstN)];
    disp(idx);
end



g=[GroupStats.Grangers(lst); GroupStats.Grangers(lstN)];
g2=[GroupStats.Pearsons(lst); GroupStats.Pearsons(lstN)];
t=[ones(length(lst),1); zeros(length(lstN),1)];

plotroc(t',g','MV-Grangers',t',g2','Pearsons',repmat(t,length(data),1)',g3','Grangers');