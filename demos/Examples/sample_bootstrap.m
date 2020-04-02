data = nirs.testing.simDataSet;

j=nirs.modules.Resample;
j=nirs.modules.OpticalDensity(j);
j=nirs.modules.BeerLambertLaw(j);
hb=j.run(data);

j=nirs.modules.Hyperscanning;
ScanA = [1:2:length(hb)]';  % The list of all the "A" files
ScanB = [2:2:length(hb)]';  % The list of all the "B" files

OffsetA = zeros(size(ScanA));  % The time shift of the "A" files (in sec)
OffsetB = zeros(size(ScanB));  % The time shift of the "B" files (in sec)
hyperlink = table(ScanA,ScanB,OffsetA,OffsetB);
j.link=hyperlink;

% do this fast for testing purposes
j.corrfcn=@(data)nirs.sFC.corr(data,false);

ConnStats =  j.run(hb);
j=nirs.modules.MixedEffectsConnectivity;
GroupConnStats = j.run(ConnStats);


%now setup the null model
j=nirs.modules.Hyperscanning;
ScanA = reshape(repmat(1:length(hb),length(hb),1),[],1);
ScanB = reshape(repmat(1:length(hb),length(hb),1)',[],1);

OffsetA = zeros(size(ScanA));  % The time shift of the "A" files (in sec)
OffsetB = zeros(size(ScanB));  % The time shift of the "B" files (in sec)
link = table(ScanA,ScanB,OffsetA,OffsetB);

link(link.ScanA>=link.ScanB,:)=[];
link(ismember(link,hyperlink),:)=[];

link=link(1:60,:);  % normally run the whole list, but I am being lazy
j.link=link;

% do this fast for testing purposes
j.corrfcn=@(data)nirs.sFC.corr(data,false);
NullData = j.run(hb);


% the Bootstrap class works like the ROC class
BootStrap=nirs.testing.permutation.ChannelStatsBootStrap;
BootStrap.pipeline=nirs.modules.MixedEffectsConnectivity;
BootStrap.output='Z';  % what field to use for output
BootStrap.samplefunction=@()NullData(randi(length(NullData),30,1));  %sampling function


BootStrap=BootStrap.run(10);
% run 10 more
BootStrap=BootStrap.run(10);


% now compute the stats
BootStrap.pval(GroupConnStats.Z)



