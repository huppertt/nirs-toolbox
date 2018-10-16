function data = loadNIRSLabDOTMAT(filename,tplfile)
% this function loads a NIRS lab output file into the toolbox format
% filename should be one of "cnt_??" (data file), "mnt_???" (probe file) or "mrk_???" (events file) and all 3
% files will be loaded based on this naming scheme.

if(~iscellstr(filename))
    filename=cellstr(filename);
end

% use the first mnt file to define the probe
[p,f,e]=fileparts(filename{1});
rootfile=[f(min(strfind(f,'_')):end) e];
mntfile=['mnt' strtok(rootfile,'.') '.mat'];
disp('creating registered probe');
probe1020=nirs.io.TPL2probe_Info(tplfile,mntfile);

n=1;
for i=1:length(filename)
    [p,f,e]=fileparts(filename{i});
    rootfile=[f(min(strfind(f,'_')):end) e];
    
    if(~exist(fullfile(p,['cnt' rootfile]),'file') | ~exist(fullfile(p,['mnt' rootfile]),'file') | ~exist(fullfile(p,['mrk' rootfile]),'file')) 
        warning(['unable to find all files for ' filename{i}]);
        continue;
    end
    try
        data(n,1)=nirs.core.Data;
    
        % get the data
        tmp=load(fullfile(p,['cnt' rootfile]));
        cnt=tmp.(['cnt' strtok(rootfile,'.')]);
        tmp=load(fullfile(p,['mrk' rootfile]));
        mrk=tmp.(['mrk' strtok(rootfile,'.')]);
        
        for id=1:length(mrk)
            stim=nirs.design.StimulusEvents;
            stim.name=mrk(id).className{1};
            stim.onset=mrk(id).time'/cnt.oxy.fs;   % THESE TIMES DO NOT MAKE SENSE
            stim.amp=mrk(id).y';
            stim.dur=mrk(id).event.desc;
            data(n).stimulus(stim.name)=stim;
        end
        
        
        % the center of the measurements are marked (not source/det pairs)
        data(n).data=[cnt.oxy.x cnt.deoxy.x]; 
        data(n).time=[0:size(data(n).data,1)-1]/cnt.oxy.fs;    
        data(n).probe=probe1020;
        link=data(n).probe.link;
        link=link(ismember(link.type,link.type{1}),:);
        types=[repmat({'hbo'},height(link),1); repmat({'hbr'},height(link),1)];
        link=[link; link];
        link.type=types;
        data(n).probe.link=link;
        
        
    n=n+1;
    end
end