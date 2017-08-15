function probe1020=registerfromfiff(fiff,probe)

% This function requires that the order of the optodes in the
% probe.optodes field matches the order in which they were digitized.
% note- you can reorder the optode list prior to running this code
% E,g, probe.optodes-
%     Name           X      Y     Z       Type       Units
%     _______________    ___    ___    _    __________    _____
% 
%     'Source-0001'       70      0    0    'Source'      'mm' 
%     'Source-0002'       50     20    0    'Source'      'mm' 
%     'Source-0003'       75     65    0    'Source'      'mm' 
%     'Source-0004'       50     90    0    'Source'      'mm' 
%     'Source-0005'      -50     20    0    'Source'      'mm' 
%     'Source-0006'      -70      0    0    'Source'      'mm' 
%     'Source-0007'      -50     90    0    'Source'      'mm' 
%     'Source-0008'      -75     65    0    'Source'      'mm' 
%     'Detector-0001'     70     35    0    'Detector'    'mm' 
%     'Detector-0002'     50     50    0    'Detector'    'mm' 
%     'Detector-0003'     70    110    0    'Detector'    'mm' 
%     'Detector-0004'     50    115    0    'Detector'    'mm' 
%     'Detector-0005'    -50     50    0    'Detector'    'mm' 
%     'Detector-0006'    -70     35    0    'Detector'    'mm' 
%     'Detector-0007'    -50    115    0    'Detector'    'mm' 
%     'Detector-0008'    -70    110    0    'Detector'    'mm' 
%     'Detector-0009'     90      0    0    'Detector'    'mm' 
%     'Detector-0010'     30     20    0    'Detector'    'mm' 
%     'Detector-0011'     96     65    0    'Detector'    'mm' 
%     'Detector-0012'     30     90    0    'Detector'    'mm' 
%     'Detector-0013'    -30     20    0    'Detector'    'mm' 
%     'Detector-0014'    -90      0    0    'Detector'    'mm' 
%     'Detector-0015'    -30     90    0    'Detector'    'mm' 
%     'Detector-0016'    -95     65    0    'Detector'    'mm' 

digpts=eeg.io.readFIFFdigpts(fiff);

pts1020=nirs.util.list_1020pts('?');

[loca,locb]=ismember(pts1020.Name,digpts.kind);
lst=find(loca);
xyz=[pts1020.X(loca(lst)) pts1020.Y(loca) pts1020.Z(loca(lst))];
xyz2=[digpts.X(locb(lst)) digpts.Y(locb(lst)) digpts.Z(locb(lst))];
xyz(:,4)=1;
xyz2(:,4)=1;
T=xyz2\xyz;
T(3,3)=1;

p = [digpts.X digpts.Y digpts.Z ones(size(digpts.X))]*T;

BEM=nirs.registration.Colin27.BEM;

[r,t]=icp(BEM.mesh(1).nodes',p(:,[1:3])');
p2=(r * p(:,1:3)' + t*ones(1,size(p,1)))';

lst=[8:8+height(probe.optodes)-1];

% remove the short distances 
sdist=probe.distances<25;

lst(size(probe.srcPos,1)+probe.link.detector(sdist))=lst(probe.link.source(sdist));


optodes_reg=probe.optodes;
optodes_reg.X=p2(lst,1);
optodes_reg.Y=p2(lst,2);
optodes_reg.Z=p2(lst,3);

fid=table(digpts.kind(1:7),p2(1:7,1),p2(1:7,2),p2(1:7,3),repmat({'FID-anchor'},7,1),repmat({'mm'},7,1),...
    'VariableNames',{'Name','X','Y','Z','Type','Units'});



probe1020=nirs.core.Probe1020;
probe1020.optodes=probe.optodes;
probe1020.link=probe.link;
probe1020.optodes_registered=[optodes_reg; fid];

BEM.mesh(1).fiducials.Draw(:)=false;
probe1020=probe1020.regsister_mesh2probe(BEM.mesh);
probe1020.defaultdrawfcn='3D mesh';



