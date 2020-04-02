function export_nirstorm(data,folder)
% this function will export nirs data from my toolbox into the
% BrainStorm/NIRStorm format

dbFile = bst_get('BrainstormDbFile');
bstOptions = load(dbFile);


bsfolder=bstOptions.BrainStormDbDir;
system(['mkdir -p ' fullfile(bsfolder,folder,'data')]);
system(['mkdir -p ' fullfile(bsfolder,folder,'anat')]);

i=1;
demo=nirs.createDemographicsTable(data);
if(ismember(demo.Properties.VariableNames,'subject'))
    subj=demo.subject;
else
    subj=repmat(cellstr('UnknownSubject'),length(data),1);
end

pp= which('nirs.io.export_nirstorm');
pp=fileparts(pp);

us=unique(subj);
for i=1:length(us)
   lst=find(ismember(subj,us{i}));
   save_BS(data(lst),fullfile(bsfolder,folder),us{i});
   save_BS_anat( data(lst),fullfile(bsfolder,folder),us{i});
   system(['cp -r ' fullfile(pp,'private','nirstorm','data','subj','@*') ' ' fullfile(bsfolder,folder,'data',us{i})]);
end

system(['cp -r ' fullfile(pp,'private','nirstorm','data','@*') ' ' fullfile(bsfolder,folder,'data')]);

system(['cp -r ' fullfile(pp,'private','nirstorm','anat','@*') ' ' fullfile(bsfolder,folder,'anat')]);


brainstorm('startjava');
org.brainstorm.file.Pack.zip(fullfile(bsfolder,folder),fullfile(bsfolder,[folder '.zip']));
system(['rm -rf ' fullfile(bsfolder,folder)]); 
 
ProtocolName = 'NIRSToolbox';
if ~brainstorm('status')
    brainstorm('start');
end
gui_brainstorm('DeleteProtocol', ProtocolName);
gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);

brainstorm('start');
import_subject(fullfile(bsfolder,[folder '.zip']));
system(['rm ' fullfile(bsfolder,[folder '.zip'])]);

end

function save_BS_anat(data,folder,subjID)

system(['mkdir -p ' fullfile(folder,'anat',subjID)]);

folder=fullfile(folder,'anat',subjID);

m=data(1).probe.getmesh;
scalp=makesurf(m(1),fullfile(folder,'tess_scalp.mat'));
cortex=makesurf(m(end),fullfile(folder,'tess_cortex.mat'));
innerskull='';
outerskull='';

if(~exist(fullfile(folder,'brainstormsubject.mat'),'file'))
    bsstudy=struct;
    bsstudy.Comments='';
    bsstudy.Scalp=scalp;
    bsstudy.InnerSkull=innerskull;
    bsstudy.OuterSkull=outerskull;
    bsstudy.UseDefaultAnat=0;
    bsstudy.UseDefaultChannel=0;
    bsstudy.Anatomy='';
    bsstudy.Cortex=cortex;
    save(fullfile(folder,'brainstormsubject.mat'),'-STRUCT','bsstudy');
end



end


function save_BS(data,folder,subjID)

system(['mkdir -p ' fullfile(folder,'data',subjID)]);     


for i=1:length(data)
    
    if(1)
        fold=fullfile(folder,'data',subjID,'@dataNIRS');
        name=['data_' num2str(i-1) 'NIRS'];
        savelink(data(i),fold,name);
    end
end
 
 
 
% 
% data_hproj_180403_1058.mat
% headmodel_nirs_mcx_fluence.mat
% results_Sensitivities_-_WL1_180501_1644.mat
% results_Sensitivities_-_WL2_180501_1644.mat
% results_Summed_sensitivities_-_WL1_180501_1644.mat
% results_Summed_sensitivities_-_WL2_180501_1644.mat
% brainstormstudy.mat
% channel_nirsbrs.mat
% data_180403_1413.mat
% data_Visual_average_180403_1425.mat



end


function savelink(data,folder,name)

system(['mkdir -p ' folder]);
i=1;
if(~exist(fullfile(folder,'brainstormstudy.mat'),'file'))
    bsstudy=struct;
    bsstudy.DataOfStudy = datestr(now);
    bsstudy.Name=name;
    bsstudy.BadTrials=[];
    save(fullfile(folder,'brainstormstudy.mat'),'-STRUCT','bsstudy');
end
% 
% if(~exist(fullfile(folder,'brainstormsubject.mat'),'file'))
%     bsstudy=struct;
%     bsstudy.DataOfStudy = datestr(now);
%     bsstudy.Name=name;
%     bsstudy.BadTrials=[];
%     save(fullfile(folder,'brainstormsubject.mat'),'-STRUCT','bsstudy');
% end

%  brainstormstudy.mat
%      DateOfStudy: '03-Apr-2018'
%            Name: '@rawVisualTask_BlockDesign'
%       BadTrials: []

if(~exist(fullfile(folder,'channel_nirsbrs.mat'),'file'))
    nirsbrs=struct;
    nirsbrs.Comment = data(i).description;
    nirsbrs.MegRefCoef=[];
    nirsbrs.Projector=[];
    nirsbrs.TransfMeg{1}=eye(4);
    nirsbrs.TransfMegLabels={'Native=>Brainstorm/CTF'};
    nirsbrs.TransfEeg{1}=eye(4);
    nirsbrs.TransfEegLabels={'Native=>Brainstorm/CTF'};
    
    m=data.probe.getmesh;
    lstC=find(ismember(m(1).fiducials.Name,{'lpa','rpa','nas'}));
    lstNC=find(~ismember(m(1).fiducials.Name,{'lpa','rpa','nas'}));
    
    nirsbrs.HeadPoints=struct;
    nirsbrs.HeadPoints.Loc=[m(1).fiducials.X([lstC; lstNC])...
        m(1).fiducials.Y([lstC; lstNC])...
        m(1).fiducials.Z([lstC; lstNC])]/1000;
    nirsbrs.HeadPoints.Label=vertcat(m(1).fiducials.Name(lstC),m(1).fiducials.Name(lstNC));
    nirsbrs.HeadPoints.Type=vertcat(repmat({'CARDINAL'},length(lstC),1),...
        repmat({'FID'},length(lstNC),1));
    
    %            Loc: [3x3 double]
    %            Label: {'Nasion'  'LeftEar'  'RightEar'}
    %            Type: {'CARDINAL'  'CARDINAL'  'CARDINAL'}
    
    nirsbrs.Channel=struct;
    p=data.probe.swap_reg;
    
    for i=1:height(data.probe.link)
        type=['WL' num2str(p.link.type(i))];
        nirsbrs.Channel(i).Name=['S' num2str(data.probe.link.source(i)) ...
            'D' num2str(data.probe.link.detector(i)) type];
        nirsbrs.Channel(i).Type='NIRS';
        nirsbrs.Channel(i).Loc=[p.srcPos(data.probe.link.source(i),:)'...
            p.detPos(data.probe.link.detector(i),:)']/1000;
        nirsbrs.Channel(i).Orient=[];
        nirsbrs.Channel(i).Weight=1;
        nirsbrs.Channel(i).Comment=[];
        nirsbrs.Channel(i).Group=type;
    end
    nirsbrs.IntraElectrodes=[];
    nirsbrs.History={'export from NIRS-toolbox'};
    
    nirsbrs.SCS=struct;
    ii=find(ismember(m(1).fiducials.Name,{'nas'}));
    nirsbrs.SCS.NAS=[m(1).fiducials.X(ii) m(1).fiducials.Y(ii) m(1).fiducials.Z(ii)]/1000;
    ii=find(ismember(m(1).fiducials.Name,{'lpa'}));
    nirsbrs.SCS.LPA=[m(1).fiducials.X(ii) m(1).fiducials.Y(ii) m(1).fiducials.Z(ii)]/1000;
    ii=find(ismember(m(1).fiducials.Name,{'rpa'}));
    nirsbrs.SCS.RPA=[m(1).fiducials.X(ii) m(1).fiducials.Y(ii) m(1).fiducials.Z(ii)]/1000;
    nirsbrs.SCS.R=eye(3);
    nirsbrs.SCS.T=[0 0 0]';
    nirsbrs.SCS.Origin=[0 0 0]';
    nirsbrs.Nirs=struct;
    nirsbrs.Nirs.Wavelengths=unique(p.link.type);
    save(fullfile(folder,'channel_nirsbrs.mat'),'-STRUCT','nirsbrs');
end


rawdata=struct;
rawdata.F=data.data';
rawdata.Std=[];
[~,rawdata.Comment]=fileparts(data.description);
rawdata.ChannelFlag=ones(height(data.probe.link),1);
rawdata.Time=[data.time(1) data.time(end)];
rawdata.Datatype='recordings';
rawdata.Device='unknown';
rawdata.nAvg=1;
rawdata.Events=struct;
rawdata.ColormapType=[];
rawdata.DisplayUnits=[];
rawdata.History={'import from nirs toolbox'};

% F=struct;
% F.filename = data.description;
% F.format='NIRS-BRS';
% F.device='unknown';
% F.condition='';
% F.comment='';
% F.byteorder='l';
% F.header=[];
% F.channelflag=ones(height(data.probe.link),1);
% F.acq_data=[];
% F.fid=[];

F.prop=struct;
F.prop.times=[data.time(1) data.time(end)];
F.prop.samples=[0 length(data.time)-1];
F.prop.sfreq=data.Fs;
F.prop.nAvg=1;
F.prop.currCtfComp=[];
F.prop.destCtfComp=[];

F.epochs=struct('label',{},'samples',[],'times',[],'nAvg',[],...
    'select',[],'bad',[],'channelflag',[]);

F.events=struct;
st=data.stimulus;
for i=1:st.count
    s=st(st.keys{i});
    F.events(i).label=s.name;
    F.events(i).color=[1 0 0];
    F.events(i).epochs=s.amp;
    F.events(i).times=s.onset;
    F.events(i).samples=find(ismember(data.time,s.onset));
    F.events(i).reactTimes=[];
    F.events(i).select=0;
end

rawdata.Events=F.events;

system(['mkdir -p ' folder]);
save(fullfile(folder,name),'-Struct','rawdata');

end

function file=makesurf(m,file);

[p,f,e]=fileparts(file);
[~,p]=fileparts(p);
s=struct;
s.Vertices=m.nodes/1000;
s.Faces=m.faces;
s.VertConn=[];
s.VertNormals=[];
s.Curvature=[];
s.SulciMap=[];
s.Atlas=struct('Name',[],'Scouts',[]);
s.iAtlas=1;
s.tess2mri_interp=[];
s.Reg=struct;
s.Reg.Sphere=struct('Vertices',[]);
s.History={'exported from nirs toolbox'};
s.Comment=f(strfind(f,'_')+1:end);

save(file,'-STRUCT','s');

file=[p filesep f e];
end