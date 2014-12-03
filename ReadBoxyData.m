function ISSdata = ReadBoxyData(fileName,WF)
%  ISSdata = ReadBoxyData(fileName)
%  
%  This function reads the BOXY text file into a structured data class
%  
%  Written by T. Huppert
%  April 7th, 2007
global wavelengthlst;

wavelengthlst{1}=[1 16 17 32];  %676
wavelengthlst{2}=[2 15 18 31];  %690
wavelengthlst{3}=[3 14 19 30]; %750
wavelengthlst{4}=[4 13 20 29 8 9 24 25]; %788
wavelengthlst{5}=[5 12 21 28]; %800
wavelengthlst{6}=[6 11 22 27]; %808
wavelengthlst{7}=[7 10 23 26]; %830




%Open the file and verify its contents
fid = fopen(fileName,'r');
if fid==0
    error('Data file could not be opened');
end

%Initialize the data structure
ISSdata=[];

ISSdata.SD.Lambda=[676 690 750 788 800 808 830];

%BOXY files are delimited by '#' at the start of each block 
%Read in each block
[BlockStrACQinfo errorflag] = GetBlock(fid,'#ACQ INFORMATION');
if(errorflag)
    error('Failed to read aquistion info settings');
end
[ISSdata] = ParseInfo(ISSdata,BlockStrACQinfo,'ACQinfo');


[BlockStrFileInfo errorflag] = GetBlock(fid,'#FILE INFORMATION');
if(errorflag)
    error('Failed to read file info settings');
end
[ISSdata] = ParseInfo(ISSdata,BlockStrFileInfo,'FILEinfo');


[BlockStrCalibration errorflag] = GetBlock(fid,'#CALIBRATION STATE');
if(errorflag)
    error('Failed to read calibration settings');
end
[ISSdata] = ParseInfo(ISSdata,BlockStrCalibration,'CALIBRATIONinfo');


%***********
%Calibration info
[BlockStrAuxCalibration errorflag] = GetBlock(fid,'#AUX CALIBRATION VALUES');
if(errorflag)
    error('Failed to read auxillary calibration settings');
end
ISSdata = ParseCalibrationInfo(ISSdata,BlockStrAuxCalibration,'Aux');

[BlockStrWFCalibration errorflag] = GetBlock(fid,'#WF CALIBRATION VALUES');
if(errorflag)
    error('Failed to read WF calibration settings');
end
ISSdata = ParseCalibrationInfo(ISSdata,BlockStrWFCalibration,'WF');



%***********
%Distance info
[BlockStrDistance errorflag] = GetBlock(fid,'#DISTANCE SETTINGS');
if(errorflag)
    warning('Failed to read distance settings');
else
    ISSdata = ParseDistance(ISSdata,BlockStrDistance);
end

%************
%Finally, the Data
[BlockStrDataA errorflag] = GetBlock(fid,'#DATA BEGINS');
[BlockStrData errorflag] = GetBlock(fid,'#DATA BEGINS',true);
if(errorflag)
    error('Failed to read data');
end
BlockStrData{1}=BlockStrDataA{1};
ISSdata = ParseData(ISSdata,BlockStrData);

if(exist('WF'))
    %remove calibration (if was applied)
    ntp=length(ISSdata.Data.time);
    
    if(~isempty(ISSdata.CALIBRATIONinfo))
    
    factor=reshape(ISSdata.CalibrationValues.WF.Factor(:,1:end/3)',[],1)*ones(1,ntp);
    term=reshape(ISSdata.CalibrationValues.WF.Factor(:,1:end/3)',[],1)*ones(1,ntp);
    ISSdata.Data.AC=ISSdata.Data.AC./factor-term;
    
    factor=reshape(ISSdata.CalibrationValues.WF.Factor(:,1+end/3:2*end/3)',[],1)*ones(1,ntp);
    term=reshape(ISSdata.CalibrationValues.WF.Factor(:,1+end/3:2*end/3)',[],1)*ones(1,ntp);
    ISSdata.Data.DC=ISSdata.Data.DC./factor-term;
    
    term=reshape(ISSdata.CalibrationValues.WF.Factor(:,1+2*end/3:end)',[],1)*ones(1,ntp);
    ISSdata.Data.Phase=ISSdata.Data.Phase-term;
    end
  
    if(~isempty(WF))
        ISSdata.CalibrationValues.WF=WF;
        
        factor=reshape(ISSdata.CalibrationValues.WF.Factor(:,1:end/3)',[],1)*ones(1,ntp);
        term=reshape(ISSdata.CalibrationValues.WF.Term(:,1:end/3)',[],1)*ones(1,ntp);
        ISSdata.Data.AC=(ISSdata.Data.AC+term).*factor;
        
        factor=reshape(ISSdata.CalibrationValues.WF.Factor(:,1+end/3:2*end/3)',[],1)*ones(1,ntp);
        term=reshape(ISSdata.CalibrationValues.WF.Term(:,1+end/3:2*end/3)',[],1)*ones(1,ntp);
        ISSdata.Data.DC=(ISSdata.Data.DC+term).*factor;
        
        term=reshape(ISSdata.CalibrationValues.WF.Term(:,1+2*end/3:end)',[],1)*ones(1,ntp);
        ISSdata.Data.Phase=ISSdata.Data.Phase+term;
    end
    
    
end


%Setup SD info
ISSdata.SD.MeasurementLst=ISSdata.Data.MeasurementList;
ISSdata.SD.NumDet=length(unique(ISSdata.SD.MeasurementLst(:,2)));
ISSdata.SD.MeasurementLst(:,3)=1;

for idx=1:length(ISSdata.SD.Lambda)
    lst=find(ismember(ISSdata.SD.MeasurementLst(:,1),wavelengthlst{idx}));
    ISSdata.SD.MeasurementLst(lst,4)=idx;
end
ISSdata.SD.ModFreq=110;

ISSdata.Data.Phase=ISSdata.Data.Phase*2*pi/360;

return




%**************************************************************************
%**************************************************************************
function [BlockStr errorflag] = GetBlock(fid,str,con2num)
% This subfunction reads in a block of the BOXY text file begining with the string "#str"
% until the next # header

errorflag=true;

try
    %Make sure we start at the begining of the file
    frewind(fid);

    %skip the initial header info
    while(1)
        nextLine=fgetl(fid);
        if ~isempty(strfind(nextLine,'#'))
            break;
        end
    end

    %Now look for the block of interest
    while(isempty(strfind(nextLine,str)))
        while(1)
            nextLine=fgetl(fid);
            if ~isempty(strfind(nextLine,'#'))
                break;
            end
            if(feof(fid))
                warning(['Failed to find ' str]);
                %Hit the end of the file
                BlockStr=[];
                return
            end
        end
        
        
        
    end

    %Save the block of strings
    BlockStr={};
    cnt=1;
    while(1)
        nextLine=fgetl(fid);
        if ~isempty(strfind(nextLine,'#'))
            break;
        else
            if(exist('con2num') && con2num)
                BlockStr{cnt}=sscanf(nextLine,'%f');
            else
                BlockStr{cnt}=nextLine;
            end
            cnt=cnt+1;
        end
    end
end

errorflag=false;  %This will only reset if no errors occured
return



%***************************************************************************
%***************************************************************************
function [ISSdata] = ParseInfo(ISSdata,BlockStr,FieldStr)
%This subfunction parses the string infomation

try
    terms{1} = 'Detector Channels';
    terms{2} = 'External MUX Channels';
    terms{3} = 'Auxiliary Channels';
    terms{4} = 'Waveform (CCF) Frequency (Hz)';
    terms{5} = 'Waveforms Skipped';
    terms{6} = 'Waveforms Averaged';
    terms{7} = 'Cycles Averaged';
    terms{8} = 'Acquisitions per Waveform';
    terms{9} = 'Updata Rate (Hz)';
    terms{10} = 'External MUX Channel results are PARSED';
    terms{11} = 'AC DC and Phase are NOT grouped';
    terms{12} = 'Companion Program Settings File CREATED';
    terms{13} = 'AC data not excluded';
    terms{14} = 'DC data not excluded';
    terms{15} = 'Phase data not excluded';
    terms{16} = 'Auxillary Chn. data not excluded';
    terms{17}='Waveform Calibration Values APPLIED';


    FieldValues =[];
    for termIdx=1:length(terms)
        for strIdx=1:length(BlockStr)
            if(~isempty(strfind(BlockStr{strIdx},terms{termIdx})))
                AddField=terms{termIdx};
                AddField(strfind(AddField,' '))='_';
                AddField(strfind(AddField,'('))='_';
                AddField(strfind(AddField,')'))='_';
                AddField(strfind(AddField,'.'))='_';

                [ValueField]=strtok(BlockStr{strIdx},' ');
                if strcmp(ValueField,'TRUE')
                    ValueField=true;
                elseif strcmp(ValueField,'FALSE')
                    ValueField=false;
                else
                    ValueField = str2num(ValueField);
                end
                FieldValues=setfield(FieldValues,AddField,ValueField);
            end
        end
    end
    ISSdata=setfield(ISSdata,FieldStr,FieldValues);
catch
    error(['Failed to read ' FieldStr ' info']);
end

return


%***************************************************************************
%***************************************************************************
function ISSdata = ParseCalibrationInfo(ISSdata,BlockStr,FieldName)
%This subfunction parses the calibration info

try
    if ~isfield(ISSdata,'CalibrationValues')
        ISSdata=setfield(ISSdata,'CalibrationValues',[]);
    end
    FieldValues=[];

    nDet = ISSdata.ACQinfo.Detector_Channels;
    nSrc = ISSdata.ACQinfo.External_MUX_Channels;
    if(isfield(ISSdata.ACQinfo,'Auxiliary_Channels'))
    nAux = ISSdata.ACQinfo.Auxiliary_Channels;
    else
        ISSdata.ACQinfo.Auxiliary_Channels=0;
        nAux=0;
    end

    if(strcmp(FieldName,'Aux'))
        if nAux==0
            %no aux data
            return
        end
        for StrIdx=1:length(BlockStr)
            if ~isempty(strfind(BlockStr{StrIdx},'Term'))
                Term = str2num(BlockStr{StrIdx}(length('Term')+1:end));
                FieldValues = setfield(FieldValues,'Term',Term);
            elseif ~isempty(strfind(BlockStr{StrIdx},'Factor'))
                strTmp=BlockStr{StrIdx}(length('Factor')+1:end);
                lst=strfind(strTmp,'X');
                strTmp(lst)='1';
                tmp=str2num(strTmp);
                strTmp=BlockStr{StrIdx}(length('Factor')+1:end);
                strTmp(lst)='2';
                tmp2=str2num(strTmp);
                lst=find(tmp-tmp2~=0);
                tmp(lst)=NaN;
                Factor = tmp;
                FieldValues = setfield(FieldValues,'Factor',Factor);
            end
        end
    else
        for idx=1:nDet
            DetectorNames{idx} = char(idx+double('A')-1);
        end
        
        Factor = zeros(nDet,nSrc*3);
        Term   = zeros(nDet,nSrc*3);
        Starts=zeros(nDet,1);
        %Figure out the block edges
        for DetIdx=1:nDet
            for StrIdx=1:length(BlockStr)
                if(~isempty(strfind(BlockStr{StrIdx},[DetectorNames{DetIdx} ' = Detector Channel'])))
                    Starts(DetIdx)=StrIdx;     
                end
            end
        end
        %Now read the data
        [Starts,DetOrder]=sort(Starts);
        Starts(end+1)=length(BlockStr);
        for idx=1:length(DetOrder);
            DetIdx=DetOrder(idx);
            for StrIdx=Starts(idx):Starts(idx+1)
                if ~isempty(strfind(BlockStr{StrIdx},'Term'))
                    Term(DetIdx,:) = str2num(BlockStr{StrIdx}(length('Term')+1:end));
                elseif ~isempty(strfind(BlockStr{StrIdx},'Factor'))
                    strTmp=BlockStr{StrIdx}(length('Factor')+1:end);
                    lst=strfind(strTmp,'X');
                    strTmp(lst)='1';
                    tmp=str2num(strTmp);
                    strTmp=BlockStr{StrIdx}(length('Factor')+1:end);
                    strTmp(lst)='2';
                    tmp2=str2num(strTmp);
                    lst=find(tmp-tmp2~=0);
                    tmp(lst)=NaN;
                    Factor(DetIdx,:) = tmp;
                end
            end
        end
        FieldValues = setfield(FieldValues,'Term',Term);
        FieldValues = setfield(FieldValues,'Factor',Factor);
    end
    ISSdata.CalibrationValues=setfield(ISSdata.CalibrationValues,FieldName,FieldValues);

catch
    error(['Failed to read calibration info for ' FieldName ]);
end
return


function ISSdata = ParseDistance(ISSdata,BlockStr)
try
    nDet = ISSdata.ACQinfo.Detector_Channels;
    nSrc = ISSdata.ACQinfo.External_MUX_Channels;
    
    for idx=1:nDet
        DetectorNames{idx} = char(idx+double('A')-1);
    end
    Distances = zeros(nDet,nSrc);
    Starts=zeros(nDet,1);
    
    %Figure out the block edges
    for DetIdx=1:nDet
        for StrIdx=1:length(BlockStr)
            if(~isempty(strfind(BlockStr{StrIdx},[DetectorNames{DetIdx} ' = Detector Channel'])))
                Starts(DetIdx)=StrIdx;
            end
        end
    end
    %Now read the data
    [Starts,DetOrder]=sort(Starts);
    Starts(end+1)=length(BlockStr);
    for idx=1:length(DetOrder);
        DetIdx=DetOrder(idx);
        for StrIdx=Starts(idx):Starts(idx+1)
            [tmp, ok]=str2num(BlockStr{StrIdx});
            if ok
                Distances(DetIdx,:)=tmp;
            end
        end
    end
    ISSdata=setfield(ISSdata,'Distances',Distances);
    
catch
    error('Failed to parse distance files');
end
return

%**************************************************************************
%**************************************************************************
function ISSdata = ParseData(ISSdata,BlockStr)
%This function reads in the data

global wavelengthlst;

try
    terms{1}='time';
    terms{2}='group';
    terms{3}='step';
    terms{4}='mark';
    terms{5}='flag';
    
    ISSdata = setfield(ISSdata,'Data',[]);
    
%     Data =[];
%     for StrIdx=1:length(BlockStr)
%         if(ischar(BlockStr{StrIdx}))
%         [tmp, ok]=str2num(BlockStr{StrIdx});
% %         else
% %             tmp=BlockStr{StrIdx};
% %             ok=true;
%         end
%         if ok
%             Data =[Data; tmp];
%         else
%             if ~isempty(BlockStr{StrIdx})
%                 ColumnNames=BlockStr{StrIdx};
%             end
%         end
%     end
    
    ColumnNames=BlockStr{1};
    Data=horzcat(BlockStr{3:end})';
    
    %Parse the columnnames
    cnt=1; 
    Cnames={};
    [Cnames{cnt},rest]=strtok(ColumnNames); cnt=cnt+1;
    while(~isempty(rest))
        [Cnames{cnt},rest]=strtok(rest); cnt=cnt+1;
    end
    
    for idxTerm=1:length(terms)
        for idxColumn=1:length(Cnames)
            if(~isempty(strfind(terms{idxTerm},Cnames{idxColumn})))
                ISSdata.Data = setfield(ISSdata.Data,terms{idxTerm},Data(:,idxColumn));
            end
        end
    end
    
    %Do the detectors
    nDet = ISSdata.ACQinfo.Detector_Channels;
    nSrc = ISSdata.ACQinfo.External_MUX_Channels;
    nAux = ISSdata.ACQinfo.Auxiliary_Channels;

    DataDC=zeros(nDet*nSrc,size(Data,1));
    DataAC=zeros(nDet*nSrc,size(Data,1));
    DataPH=zeros(nDet*nSrc,size(Data,1));
    MeasurementList = [repmat([1: nSrc]',nDet,1) reshape(repmat([1: nDet],nSrc,1),[],1)];
    
    for widx=1:length(wavelengthlst); 
        lst=find(ismember(MeasurementList(:,1),wavelengthlst{widx}));
        MeasurementList(lst,4)=widx;
    end
    
    
    for idx=1:nDet
        DetectorNames{idx} = char(idx+double('A')-1);
    end
    
    for MLidx=1:size(MeasurementList,1)
        DCname = [num2str(DetectorNames{MeasurementList(MLidx,2)}) '-DC' num2str(MeasurementList(MLidx,1))];
        ACname = [num2str(DetectorNames{MeasurementList(MLidx,2)}) '-AC' num2str(MeasurementList(MLidx,1))];
        PHname = [num2str(DetectorNames{MeasurementList(MLidx,2)}) '-Ph' num2str(MeasurementList(MLidx,1))];
        for idx=1:length(Cnames)
            if(~isempty(strfind(DCname,Cnames{idx})))
                DataDC(MLidx,:)=Data(:,idx)';
            elseif (~isempty(strfind(ACname,Cnames{idx})))
                DataAC(MLidx,:)=Data(:,idx)';
            elseif (~isempty(strfind(PHname,Cnames{idx})))
                DataPH(MLidx,:)=Data(:,idx)';
            end
        end
    end

    ISSdata.Data = setfield(ISSdata.Data,'Phase',DataPH);
    ISSdata.Data = setfield(ISSdata.Data,'AC',DataAC);
    ISSdata.Data = setfield(ISSdata.Data,'DC',DataDC);
    ISSdata.Data = setfield(ISSdata.Data,'MeasurementList',MeasurementList);
    
    
    %Finally, do the auxillary
    if nAux~=0
        aux = zeros(nAux,size(Data,1));
        for AuxIdx=1:nAux
            AuxName = ['aux-' num2str(AuxIdx)];
            for idx=1:length(Cnames)
                if (~isempty(strfind(AuxName,Cnames{idx})))
                    aux(AuxIdx,:)=Data(:,idx)';
                end
            end
        end
    else
        aux=[];
    end
    ISSdata.Data = setfield(ISSdata.Data,'Aux',aux);
catch
    error('Failed to parse data structure');
end
return