function probe = simProbe()
% This function returns a default regstered NIRS probe
        
    srcPos(:,1) = (-80:20:80)';
    srcPos(:,2:3) = 0;
    
    detPos(:,1) = (-70:20:70)';
    detPos(:,2) = 25;
    detPos(:,3) = 0;
    
    probe = nirs.core.Probe(srcPos,detPos);
    
    link = [1	1	690
        2	1	690
        2	2	690
        3	2	690
        3	3	690
        4	3	690
        4	4	690
        5	4	690
        5	5	690
        6	5	690
        6	6	690
        7	6	690
        7	7	690
        8	7	690
        8	8	690
        9	8	690
        1	1	830
        2	1	830
        2	2	830
        3	2	830
        3	3	830
        4	3	830
        4	4	830
        5	4	830
        5	5	830
        6	5	830
        6	6	830
        7	6	830
        7	7	830
        8	7	830
        8	8	830
        9	8	830];
    
    link = sortrows(link);
    
    probe.link = table(link(:,1), link(:,2), link(:,3), ...
        'VariableNames', {'source', 'detector', 'type'});
    
    
Name{1}='FpZ';
xyz(1,:)=[0 0 0];
Type{1}='FID-anchor';  % This is an anchor point
Units{1}='mm';

%Now let's add a few more
Name{2}='Cz';
xyz(2,:)=[0 100 0];
Type{2}='FID-attractor';  % This is an attractor
Units{2}='mm';

Name{3}='T7';
xyz(3,:)=[-200 0 0];
Type{3}='FID-attractor';  % This is an attractor
Units{3}='mm';

Name{4}='T8';
xyz(4,:)=[200 0 0];
Type{4}='FID-attractor';  % This is an attractor
Units{4}='mm';
    
fid=table(Name',xyz(:,1),xyz(:,2),xyz(:,3),Type',Units',...
    'VariableNames',{'Name','X','Y','Z','Type','Units'});
% and concatinate it to the probe
probe.optodes=[probe.optodes; fid];
probe=nirs.util.registerprobe1020(probe);

lambda=unique(probe.link.type);
fwdBEM=nirs.registration.Colin27.BEM(lambda);
probe=probe.regsister_mesh2probe(fwdBEM.mesh);

end