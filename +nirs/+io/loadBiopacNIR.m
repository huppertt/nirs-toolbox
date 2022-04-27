function raw = loadBiopacNIR(nir_filename)
% Reads .nir files from FNIR Devices model 1200 to 2000 devices
% 	Merges with associated marker files
%   auto loads .log files for Subject information
%   attempts to load .mrk and _C.mrk files if available or _Mark1.mrk
%       Up to marker 5 (for files recorded in COBI studio modern
%   Markers are marked as Mark[SourceNumber]_[MarkerValue]
%       Legacy manual markers are recorded as Souce #100

% Developed by Adrian Curtin (adrian.b.curtin@drexel.edu)

if nargin < 1 % No Arguments - Open fNIRS and mrk file
	[nir_filename, pathname] = uigetfile({'*.nir';'*.*'},'Select fNIR Devices (BioPac) .nir filename');
	nir_filename=[pathname nir_filename];
    
    if(isempty(nir_filename)||~isstr(nir_filename))
		raw=[];
       return; 
    end
end



if ~isstr(nir_filename)
	error('Input must be a string representing a filename');
end

if(~isempty(nir_filename))
	fid = fopen(nir_filename);
end

if fid==-1
	error('Data %s not found or permission denied',nir_filename);
end

fname=nir_filename;

fileroot=fname(1:strfind(fname,'.nir')-1);

%Default names for markers, manual markers, and log file

if(contains(fileroot,'_Dev'))% new marker type
    fileroot=fname(1:strfind(fname,'_Dev')-1);
    mrk_filename{1}=sprintf('%s_Mark1.mrk',fileroot); %device markers
    mrk_filename{2}=sprintf('%s_Mark2.mrk',fileroot); %device markers
    mrk_filename{3}=sprintf('%s_Mark3.mrk',fileroot); %device markers
    mrk_filename{4}=sprintf('%s_Mark4.mrk',fileroot); %device markers
    mrk_filename{5}=sprintf('%s_Mark5.mrk',fileroot); %device markers
    log_filename=sprintf('%s.log',fileroot);
    
    [~,~,~,~,subDeviceNum]=regexpi(fname,'^.*_dev(\d{1,2})_p[0-9]{1,2}.*.nir','match');
    [~,~,~,~,subProbeNum]=regexpi(fname,'^.*_dev\d{1,2}_p(\d{1,2}).*.nir','match');
    
    if(~isempty(subDeviceNum))
       subDeviceNum=str2double(subDeviceNum{1}); 
    elseif(length(subDeviceNum)>1)
        warning('Appears to have multiple device identifiers in filename. Assuming Device 1');
        subDeviceNum=1;
    else
        warning('Unable to determine device from filename. Assuming Device 1');
        subDeviceNum=1;
    end
    
    if(~isempty(subProbeNum))
       subProbeNum=str2double(subProbeNum{1}); 
    elseif(length(subProbeNum)>1)
        warning('Appears to have multiple probe identifiers in filename. Assuming Device 1');
        subProbeNum=1;
    else
        warning('Unable to determine probe from filename. Assuming Device 1');
        subProbeNum=1;
    end
   
else % Old Marker names
    mrk_filename{1}=sprintf('%s.mrk',fileroot); %manual markers
    mrk_filename{2}=sprintf('%s_C.mrk',fileroot); %device markers
    mrk_filename{2}=sprintf('%s_B.mrk',fileroot); %calculated markers
    log_filename=sprintf('%s.log',fileroot);
    
    subProbeNum=1; % Unused by old device files
    subDeviceNum=1;
end



log_info=importCOBIlog(log_filename,subDeviceNum,subProbeNum);


lineF=fgetl(fid);
while(ischar(lineF))
    if(contains(lineF, 'Start Code:'))
        sC=strsplit(lineF,'\t');
        if(length(sC)>1&&~isempty(sC{2}))
            header.startCodeAlt=str2double(sC{2});
        end
        if(length(sC)>2&&~isempty(sC{3}))
            header.startCode=str2double(sC{3})/1000;
            if((header.startCodeAlt-header.startCode)>1)
                %%warning(sprintf('StartCode Diff %.2f\nMay affect manual marker integrity',(header.startCodeAlt-header.startCode)));
            end
        else
            header.startCode=header.startCodeAlt;
        end
        
    end
    if(~isempty(strfind(lineF, 'Freq Code:')))
        header.freqCode=sscanf(lineF,'Freq Code:\t%f\n',1);  %record frequency code
    end
    if(~isempty(strfind(lineF, 'Current:')))
        header.current=sscanf(lineF,'Current:\t%f\n',1);  %record LED current used
    end
    if(~isempty(strfind(lineF, 'Gains:')))
        header.gains=sscanf(lineF,'Gains:\t%f\n',1);  %record Detector Gain used
    end
    if(~isempty(strfind(lineF, 'Other:')))
        header.other=sscanf(lineF,'Other:\t%s\n',1);  %record Detector Gain used
    end
    
    if(~isempty(strfind(lineF,'Start Time:')))
       temp=strsplit(lineF,'\t');
       temp=temp{2};
       temp=strsplit(temp,' ');
       for i=1:length(temp)
          if(~isempty(strfind(temp,':')))
             header.Time=temp;
          end
       end
    end
    
    lineF=fgetl(fid); %Get Next Line
    
    if(contains(lineF, 'Baseline Started'))
        countCheckFlag=true;
        lineF=fgetl(fid); %Get Next Line

        break;
    end
end


if(iscell(mrk_filename))
	
    tempData=[];
    for i=1:length(mrk_filename)
		markerCell{i}=importMrk(mrk_filename{i},header.startCode,i);
        if(isfield(markerCell{i},'data'))
            tempData=[tempData;markerCell{i}.data];
        end
    end
    
    
    markers=markerCell{1};
    if(~isempty(tempData)) % Merge and sort marker values (.info sections are only preserved for first marker file)
        [~,idx]=sort(tempData(:,1));
       markers.data=tempData(idx,:); 
    end
    
else
	markers=importMrk(mrk_filename,header.startCode);
end


spaceParsingMode=false;

while(ischar(lineF))
    
    if(countCheckFlag)
       countCheckFlag=false;
       numVar=sum(lineF(:)=='	');    
       
       if(numVar==0)
            numVar=sum(lineF(:)==' ')+1;
            spaceParsingMode=true;
       end
       baseline=nan(1000,numVar);
       blLineCount=0;
    end
    if(contains(lineF, 'Baseline end'))
        baseline(isnan(baseline(:,1)),:)=[]; %trim NaN rows;
        baseline(baseline(:,1)<=0,:)=[]; %remove markerrows and zero rows
        lineF=fgetl(fid); %Get Next Line
       break; 
    end
    
    if(~countCheckFlag)
        blLineCount=blLineCount+1;
        [baseline(blLineCount,:), numVar]=sscanf(lineF,'%f',[1 numVar]);
    end
    lineF=fgetl(fid); %Get Next Line
end




if(lineF==-1)
   return;
else
    if(~spaceParsingMode)
        line1=str2double(strsplit(lineF,'\t'));
    else
        line1=str2double(strsplit(lineF,' '));
    end
end

while (ischar(lineF)&&length(line1)>1)&&(line1(1)<0) 
            % Keeps searching for first line with values if it can't find
            % it for some reason
    if(~spaceParsingMode)
        line1=str2double(strsplit(lineF,'\t'));
    else
        line1=str2double(strsplit(lineF,' '));
    end
    lineF=fgetl(fid); %Get Next Line
end

numVar=length(line1(~isnan(line1))); %Count Elements


data=nan(5e5,numVar); %overinitialize array



lineCount=1;
while(lineF~=-1)
    if(~spaceParsingMode)
        numTabs=sum(lineF(:)=='	');
        if(lineF(end)~='\t') %.nir files are terminared with \t\n, but not always true
            numTabs=numTabs+1;
        end
    else
        numTabs=sum(lineF(:)==' ');
        if(lineF(end)~=' ') %.nir files are terminared with \t\n, but not always true
            numTabs=numTabs+1;
        end
    end
    if(numTabs<numVar)
        data(lineCount,:)=zeros(1,numVar);
        count=numTabs;
    else
        [data(lineCount,:), count]=sscanf(lineF,'%f',[1 numVar]);
    end
    if(count==numVar)
        lineCount=lineCount+1;
    end
    lineF=fgetl(fid);
end
data(isnan(data(:,1)),:)=[]; %trim nan rows
data(data(:,1)<=0,:)=[]; %trim zero or negative rows

data=[[baseline,zeros(size(baseline,1),size(data,2)-size(baseline,2))];data];

fclose(fid);
clear fid line1 lineCount line



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Building into nirs-toolbox format


% Now put all the data together
raw=nirs.core.Data;

raw.time=data(:,1);

newColTest=diff(data(:,end));
newColTest(newColTest<-63000)=0; % counts up to about 65k and then resets sample count

if(sum(newColTest<0)==0) % newer COBI, last column is sample time
    raw.data=data(:,2:end-1);
else
    raw.data=data(:,2:end);
end





if(~isempty(markers)&&size(markers.data,2)>1)
    [ua,~,uc]=unique(markers.data(:,2:3),'rows');
    stim=cell(length(ua),1);

    for i=1:length(ua)
        idx=(uc==i);
        stimCount=nnz(idx);
        stim{i}=nirs.design.StimulusEvents(...
            sprintf('Mark%i_%i',ua(i,2),ua(i,1)),... %name
            markers.data(idx,1)',... %onset
            zeros(stimCount,1),... %duration
            ones(stimCount,1)); %amplitude
        if(isfield(log_info,'MarkerDict'))
            dictIdx=strcmp(log_info.MarkerDict.Source,sprintf('Mark%i',ua(i,2)))...
                &strcmp(log_info.MarkerDict.Value,sprintf('%i',ua(i,1)));
            stim{i}.metadata=log_info.MarkerDict(dictIdx,:);
        end
        raw.stimulus(stim{i}.name)=stim{i};
    end
end

%probe=nirs.core.Probe(srcpos,detpos,link);


if(~isempty(log_info)&&isfield(log_info,'SubjectID'))
    raw.demographics('ID')=log_info.SubjectID;
end
if(~isempty(log_info)&&isfield(log_info,'Experimenter'))
    raw.demographics('Experimenter')=log_info.Experimenter;
end
    
if(~isempty(log_info)&&isfield(log_info,'ExperimentID'))
    raw.demographics('Session')=log_info.ExperimentID;
end
    
if(~isempty(log_info)&&isfield(log_info,'Sex'))
    raw.demographics('Gender')=log_info.Sex;
end

if(~isempty(log_info)&&isfield(log_info,'Age')&&~isnan(log_info.Age))
    raw.demographics('Age')=log_info.Age;
end

if(~isempty(log_info)&&isfield(log_info,'Description'))
    raw.demographics('Description')=log_info.Description;
end

if(~isempty(log_info)&&isfield(log_info,'Comments'))
    raw.demographics('Comment')=log_info.Comments;
end
    
%raw.demographics('Patient Info')=info{1}.Patient_Information';

raw.description=fullfile(pwd,fname);

%%%%%%

% Now deal with the probe
SrcPos=[]; DetPos=[]; link=table;

probe=guessprobefromchannelcount(size(raw.data,2));
l=probe.link;
l.source=l.source+size(SrcPos,1);
l.detector=l.detector+size(DetPos,1);
SrcPos=[SrcPos; probe.srcPos];
DetPos=[DetPos; probe.detPos];
link=[link; l];


raw.probe=nirs.core.Probe(SrcPos,DetPos,link);


end


function probe=guessprobefromchannelcount(numch)

% figure out what probe this is
switch(numch)
    
    case(48)  % 3wv  2 x 8 (16ch x3wv)  fNIR Devices 1100-1200 
       m=2;
       n=8;
       numOpt=16;
       
       DetPosX	=	[0,0,3.53553390593274,3.53553390593274,7.07106781186548,7.07106781186548,10.6066017177982,10.6066017177982,14.1421356237310,14.1421356237310];
       DetPosY	=	[3.53553390593274,0,3.53553390593274,0,3.53553390593274,0,3.53553390593274,0,3.53553390593274,0];
       DetPosZ  =   zeros(size(DetPosY));
       dI=[1,2,3,4,3,4,5,6,5,6,7,8,7,8,9,10];
       SrcPosX =	[1.76776695296637,5.30330085889911,8.83883476483184,12.3743686707646]; % Source positions per channel in cm
       SrcPosY	=	[1.76776695296637,1.76776695296637,1.76776695296637,1.76776695296637];
       SrcPosZ	=	[zeros(size(SrcPosY))];
       sI=[1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4];
       
       SrcPos=[SrcPosX',SrcPosY',SrcPosZ']*10; %convert to mm
       DetPos=[DetPosX',DetPosY',DetPosZ']*10;
       WL=[730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850];	% wavelength in nm, 0 indicates dark or ambient 
        
       link=table([repelem(sI,3)'],[repelem(dI,3)'],WL(1:numOpt*3)','VariableNames',{'source','detector','type'});
       

      probe=nirs.core.Probe(SrcPos,DetPos,link);
      %probe.link=probe.link(probe.distances==25,:);
    case(54)    %3wv 2 x 8  (16ch) + 2 ss x3wc fNIR Devices 2000
       m=2;
       n=8;
       ss=2;
       numOpt=18;
       
       DetPosX	=	[0,0,3.53553390593274,3.53553390593274,7.07106781186548,7.07106781186548,10.6066017177982,10.6066017177982,14.1421356237310,14.1421356237310,0.76776695296637,13.37437];
       DetPosY	=	[3.53553390593274,0,3.53553390593274,0,3.53553390593274,0,3.53553390593274,0,3.53553390593274,0,1.76776695296637,1.76776695296637];
       DetPosZ  =   zeros(size(DetPosY));
       dI=[1,2,3,4,3,4,5,6,5,6,7,8,7,8,9,10,11,12];
       SrcPosX =	[1.76776695296637,5.30330085889911,8.83883476483184,12.3743686707646]; % Source positions per channel in cm
       SrcPosY	=	[1.76776695296637,1.76776695296637,1.76776695296637,1.76776695296637];
       SrcPosZ	=	[zeros(size(SrcPosY))];
       sI=[1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,1,4];
       
       SrcPos=[SrcPosX',SrcPosY',SrcPosZ']*10; %convert to mm
       DetPos=[DetPosX',DetPosY',DetPosZ']*10;
       WL=[730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850];	% wavelength in nm, 0 indicates dark or ambient 
        
       link=table([repelem(sI,3)'],[repelem(dI,3)'],WL(1:numOpt*3)','VariableNames',{'source','detector','type'});
       

      probe=nirs.core.Probe(SrcPos,DetPos,link);
  case(24)    %3wv 4 x 1 (2 probe)   (4ch x 3 wv x 2 probes)  fNIR Devices Split Probe x2
       m=2;
       n=4;
       numOpt=8;
       
       DetPosX	=	[0,0,3.53553390593274,3.53553390593274,7.07106781186548,7.07106781186548,10.6066017177982,10.6066017177982];
       DetPosY	=	[3.53553390593274,0,3.53553390593274,0,3.53553390593274,0,3.53553390593274,0];
       DetPosZ  =   zeros(size(DetPosY));
       dI=[1,2,3,4];
       SrcPosX =	[1.76776695296637,8.83883476483184]; % Source positions per channel in cm
       SrcPosY	=	[1.76776695296637,1.76776695296637];
       SrcPosZ	=	[zeros(size(SrcPosY))];
       sI=[1,1,1,1];
       
       SrcPos=[SrcPosX',SrcPosY',SrcPosZ']*10; %convert to mm
       DetPos=[DetPosX',DetPosY',DetPosZ']*10;
       WL=[730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850];	% wavelength in nm, 0 indicates dark or ambient 
        
       link=table([repelem(sI,3)'],[repelem(dI,3)'],WL(1:numOpt*3)','VariableNames',{'source','detector','type'});
       

      probe=nirs.core.Probe(SrcPos,DetPos,link);
  case(24)    %3wv 4 x 1 (2 probe)   (4ch x 3 wv x 2 probes)  fNIR Devices Split Probe x2
       m=2;
       n=4;
       numOpt=8;
       
       DetPosX	=	[0,0,3.53553390593274,3.53553390593274,7.07106781186548,7.07106781186548,10.6066017177982,10.6066017177982];
       DetPosY	=	[3.53553390593274,0,3.53553390593274,0,3.53553390593274,0,3.53553390593274,0];
       DetPosZ  =   zeros(size(DetPosY));
       dI=[1,2,3,4];
       SrcPosX =	[1.76776695296637,8.83883476483184]; % Source positions per channel in cm
       SrcPosY	=	[1.76776695296637,1.76776695296637];
       SrcPosZ	=	[zeros(size(SrcPosY))];
       sI=[1,1,1,1];
       
       SrcPos=[SrcPosX',SrcPosY',SrcPosZ']*10; %convert to mm
       DetPos=[DetPosX',DetPosY',DetPosZ']*10;
       WL=[730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850];	% wavelength in nm, 0 indicates dark or ambient 
        
       link=table([repelem(sI,3)'],[repelem(dI,3)'],WL(1:numOpt*3)','VariableNames',{'source','detector','type'});
       

      probe=nirs.core.Probe(SrcPos,DetPos,link);
    case(15)    %3wv 2 x 2 +1   (4ch x 3 wv)  fNIR Devices Split Probe
       m=1;
       n=4;
       numOpt=5;
       
       DetPosX	=	[0,0,3.53553390593274,3.53553390593274,2.76776695296637];
       DetPosY	=	[3.53553390593274,0,3.53553390593274,0,1.76776695296637];
       DetPosZ  =   zeros(size(DetPosY));
       dI=[1,2,3,4,5];
       SrcPosX =	[1.76776695296637]; % Source positions per channel in cm
       SrcPosY	=	[1.76776695296637];
       SrcPosZ	=	[zeros(size(SrcPosY))];
       sI=[1,1,1,1,1];
       
       SrcPos=[SrcPosX',SrcPosY',SrcPosZ']*10; %convert to mm
       DetPos=[DetPosX',DetPosY',DetPosZ']*10;
       WL=[730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850];	% wavelength in nm, 0 indicates dark or ambient 
        
       link=table([repelem(sI,3)'],[repelem(dI,3)'],WL(1:numOpt*3)','VariableNames',{'source','detector','type'});
       

      probe=nirs.core.Probe(SrcPos,DetPos,link);
        case(15)    %3wv 2 x 2 +1   (4ch x 3 wv)  fNIR Devices Split Probe
       m=1;
       n=4;
       numOpt=5;
       
       DetPosX	=	[0,0,3.53553390593274,3.53553390593274,2.76776695296637];
       DetPosY	=	[3.53553390593274,0,3.53553390593274,0,1.76776695296637];
       DetPosZ  =   zeros(size(DetPosY));
       dI=[1,2,3,4,5];
       SrcPosX =	[1.76776695296637]; % Source positions per channel in cm
       SrcPosY	=	[1.76776695296637];
       SrcPosZ	=	[zeros(size(SrcPosY))];
       sI=[1,1,1,1,1];
       
       SrcPos=[SrcPosX',SrcPosY',SrcPosZ']*10; %convert to mm
       DetPos=[DetPosX',DetPosY',DetPosZ']*10;
       WL=[730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850];	% wavelength in nm, 0 indicates dark or ambient 
        
       link=table([repelem(sI,3)'],[repelem(dI,3)'],WL(1:numOpt*3)','VariableNames',{'source','detector','type'});
       

      probe=nirs.core.Probe(SrcPos,DetPos,link);
    case(18)    %3wv 1 x 4 +2   (6ch x 3 wv)  fNIR Devices Linear probe
       m=4;
       n=1;
       numOpt=6;
       
       DetPosX	=	[0,5,7,12,3.5,8.5];
       DetPosY	=	[0,0,0,0,0,0];
       DetPosZ  =   zeros(size(DetPosY));
       dI=[1,2,3,4,5,6];
       SrcPosX =	[2.5,9.5]; % Source positions per channel in cm
       SrcPosY	=	[0,0];
       SrcPosZ	=	[zeros(size(SrcPosY))];
       sI=[1,1,2,2,1,2];
       
       SrcPos=[SrcPosX',SrcPosY',SrcPosZ']*10; %convert to mm
       DetPos=[DetPosX',DetPosY',DetPosZ']*10;
       WL=[730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850,730,0,850];	% wavelength in nm, 0 indicates dark or ambient 
        
       link=table([repelem(sI,3)'],[repelem(dI,3)'],WL(1:numOpt*3)','VariableNames',{'source','detector','type'});
       

      probe=nirs.core.Probe(SrcPos,DetPos,link);
    otherwise
        % I don't want to just assume I can do this based on the mode
        error('This is a different probe design');
end




end



function markers=importMrkFile(mrkid,mrkSrcID)
    if(nargin<2)
        mrkSrcID=0;
    end
	frewind(mrkid);
	for i=1:5
		   lineF=fgetl(mrkid);
		   markers.info.headers{i,1}=lineF;
			  
		   if(contains(lineF,'Start Time:'))
			   temp=strsplit(lineF,'\t');
			   temp=temp{2};
			   temp=strsplit(temp,' ');
			   for i=1:length(temp)
				  if(contains(temp,':'))
					 markers.info.Time=temp;
				  end
			   end
		   end
	end

		
	linecount=0;
	lineF=fgetl(mrkid);
	markers.data=[];
	while(lineF~=-1)
	   linecount=linecount+1;
	   
	   temp=sscanf(lineF,'%f\t%f\t%f')';
	   if(length(temp)>2)
			markers.data(linecount,:)=[temp(1:2),mrkSrcID];
	   elseif(length(temp)>1)
		   markers.data(linecount,:)=[temp(1:2),mrkSrcID];
	   else
		   markers.data(linecount,:)=[zeros(1,2),0];
	   end
	   lineF=fgetl(mrkid);
	   
    end
    
    markers.data(markers.data(:,2)==0&markers.data(:,1)==0,:)=[];
	fclose(mrkid);
end


function markers=importManualMrkFile(mrkid,startCode,useStrSplit)
    if(nargin<3)
        useStrSplit=true;
    end

	frewind(mrkid);
    
    markers.info.startcode=startCode; % time noting start of fNIRS file 
        % manual marker values are subtracted from this
	
    if(~useStrSplit) %faster but more error prone
	
        [times, count]= fscanf(mrkid,'%f\t%d\t%*d\t%*s\t%*s %*s %*s %*s %*s\n', [1 inf]);
        markers.data=nan(ceil(count/2),1);
        markers.data(:,3)=times(1:2:end)'; % marker time
        markers.data(:,2)=times(2:2:end)'; % marker code
        markers.data(:,1)=markers.data(:,3)-startCode; %adjusted marker time
        markers.data(:,3)=100; %We don't use marker time here, so just use 100 as the ID

    else
        linecount=0;
        lineF=fgetl(mrkid);
        markers.data=nan(1,3);
        while(lineF~=-1)
           temp=strsplit(lineF,'\t');
           if(length(temp)>=3)
               linecount=linecount+1;
                markers.data(linecount,3)=str2num(temp{1}); % marker time
                markers.data(linecount,2)=str2num(temp{2}); % marker code
                markers.data(linecount,1)=markers.data(linecount,3)-startCode; %adjusted marker time
                if(length(temp)>=5)
                    markers.info.markerDateTime{linecount}=temp{5};
                end
           else
               continue;
           end
           lineF=fgetl(mrkid);

        end
        
        markers.data(:,3)=100;
    end

   clear count;
   clear times;
   fclose(mrkid);
end


function markers=importMrk(mrk_filename,startCode,mrkSourceID)
    if(nargin<3)
        mrkSourceID=0;
    end
	if(isstr(mrk_filename)&&~isempty(mrk_filename))
		mrkid = fopen(mrk_filename);
	else
		mrkid=-1;
	end
	
	markers=[];

	if (mrkid==-1&&~isempty(mrk_filename))
		%disp('Marker nir_filename not found or permission denied: Loading without markers');
		return;
	elseif(mrkid==-1)
		%no file provided so just return
		return
	end
	
	
	lineF=fgetl(mrkid);
	if(lineF~=-1)
		[times, count]= sscanf(lineF,'%f\t%*d\t%*d\t%*s\t%*s %*s %*s %*s %*s\n', [1 inf]);
	else
		count=0;
		lineF='';
		fclose(mrkid);
		mrkid=-1;
	end
	
   if(count~=0)  % Then its a manual marker file
		markers=importManualMrkFile(mrkid,startCode,true);
	elseif contains(lower(lineF),'marker') % then its an automatic marker file
		markers=importMrkFile(mrkid,mrkSourceID);
	else
	   markers=[];
	end

	if(~isempty(markers))
		markers.info.fname=mrk_filename;
	end

end

function loginfo=importCOBIlog(log_filename,subDevNum,subProbeNum)
    if(nargin<2)
        subProbeNum=1;
        warning('Unable to Determine Probe Number. Assuming Probe 1');
    end
    if(nargin<3)
        subDevNum=1;
        warning('Unable to Determine Device Number. Assuming Device 1');
    end
    if(isstr(log_filename)&&~isempty(log_filename))
            logfid = fopen(log_filename);
        else
            logfid=-1;
    end
        
    loginfo=[];

    if (logfid==-1&&~isempty(log_filename))
        warning('COBI log file not found, loading without log file');
        return;
    elseif(logfid==-1)
        %no file provided so just return
        return
    end
    
    linecount=1;
    
     lineF=fgetl(logfid);

    while(ischar(lineF))
       
       if(contains(lineF,'COBI log file'))
           loginfo.version=sscanf(lineF,'COBI log file %f');
       elseif(contains(lineF,'Experimenter:'))
           lineF=fgetl(logfid);
           if(lineF~=-1)
               loginfo.Experimenter=lineF;
           end
       elseif(contains(lineF,'SubjectID:'))
           lineF=fgetl(logfid);
           linecount=linecount+1;
           if(lineF~=-1)
               loginfo.SubjectID=lineF;
           end
       elseif(contains(lineF,'ExperimentID:'))
           lineF=fgetl(logfid);
           linecount=linecount+1;
           if(lineF~=-1)
               loginfo.ExperimentID=lineF;
           end
       elseif(contains(lineF,'Description:'))
           lineF=fgetl(logfid);
           linecount=linecount+1;
           loginfo.Description='';
           while(lineF~=-1)
               if(isempty(loginfo.Description))
                   loginfo.Description=lineF;
               else
                   loginfo.Description=sprintf('%s\n%s',loginfo.Description,lineF);
               end
               lineF=fgetl(logfid);
               linecount=linecount+1;
           end
      elseif(contains(lineF,'Data Sources:'))
           lineF=fgetl(logfid);
           loginfo.DataSources='';
           linecount=linecount+1;
           while(lineF~=-1)
               if(isempty(loginfo.DataSources))
                   loginfo.DataSources=lineF;
               else
                   loginfo.DataSources=sprintf('%s\n%s',loginfo.DataSources,lineF);
               end
               lineF=fgetl(logfid);
               linecount=linecount+1;
           end
     elseif(contains(lineF,'Marker Sources:'))
           lineF=fgetl(logfid);
           loginfo.MarkerSources='';
           linecount=linecount+1;
           while(lineF~=-1)
               if(isempty(loginfo.MarkerSources))
                   loginfo.MarkerSources=lineF;
               else
                   loginfo.MarkerSources=sprintf('%s\n%s',loginfo.MarkerSources,lineF);
               end
               lineF=fgetl(logfid);
               linecount=linecount+1;
           end
      elseif(contains(lineF,'Broadcasts:'))
           lineF=fgetl(logfid);
           loginfo.Broadcasts='';
           linecount=linecount+1;
           while(lineF~=-1)
               if(isempty(loginfo.Broadcasts))
                   loginfo.Broadcasts=lineF;
               else
                   loginfo.Broadcasts=sprintf('%s\n%s',loginfo.Broadcasts,lineF);
               end
               lineF=fgetl(logfid);
               linecount=linecount+1;
           end
           
       elseif(contains(lineF,'Subject Info:'))
           lineF=fgetl(logfid);
           linecount=linecount+1;
           if(lineF~=-1)
               loginfo.SubjectInfo=lineF;
           end
           if(~isempty(loginfo.SubjectInfo))
               subInfoPart=strsplit(loginfo.SubjectInfo,'-');
               if(length(subInfoPart)>1)
                   if(contains(subInfoPart{2},'Male'))
                       loginfo.Sex='M';
                   elseif(contains(subInfoPart{2},'Female'))
                       loginfo.Sex='F';
                   else
                      loginfo.Sex=''; 
                   end
               else
                   loginfo.Sex='';
               end
               
               if(~isempty(subInfoPart))
                   loginfo.Age=str2double(subInfoPart{1});
               end
           
           end
       elseif(contains(lineF,'Participants:'))
           participantString=sprintf('Dev%i_P%i',subDevNum,subProbeNum);
           
           lineF=fgetl(logfid);
           linecount=linecount+1;
           if(lineF~=-1)
               tableHeaders=strsplit(lineF,'\t');
               lastItem=tableHeaders{end};
               if(lastItem(end)==':')
                   lastItem=lastItem(1:end-1);
                   tableHeaders{end}=lastItem;
               end
               
           else
              continue; 
           end
           lineF=fgetl(logfid);
           linecount=linecount+1;
           
           participantDict=cell(1,length(tableHeaders));
           numPartDictItems=1;
           
           while(lineF~=-1)
               items=strsplit(lineF,'\t','CollapseDelimiters',false);
               if(~isempty(items)&&iscell(items))
                   participantDict(numPartDictItems,1:length(items))=items;
                   numPartDictItems=numPartDictItems+1;
               end
               lineF=fgetl(logfid);
               linecount=linecount+1;
           end
           loginfo.ParticipantDict=cell2table(participantDict(:,1:length(tableHeaders)),'VariableNames',tableHeaders);
           
           particpantIndex=contains(loginfo.ParticipantDict.Source,participantString);
           
           if(sum(particpantIndex)==0)
              warning('Unable to find participant with source %s',participantString); 
           end
           
           if(sum(particpantIndex)>0&&size(loginfo.ParticipantDict,1)>=find(particpantIndex==1))
               if(iscell(loginfo.ParticipantDict.Gender(particpantIndex))&&~isempty(loginfo.ParticipantDict.Gender{particpantIndex}))
                   if(contains(loginfo.ParticipantDict.Gender{particpantIndex},'Male'))
                       loginfo.Sex='M';
                   elseif(contains(loginfo.ParticipantDict.Gender{particpantIndex},'Female'))
                       loginfo.Sex='F';
                   else
                      loginfo.Sex=''; 
                   end
               else
                   loginfo.Sex='';
               end
               
               if(~isempty(loginfo.ParticipantDict.Age(particpantIndex)))
                   loginfo.Age=str2double(loginfo.ParticipantDict.Age(particpantIndex));
               end 
               
               if(~isempty(loginfo.ParticipantDict.SubjectID(particpantIndex)))
                   loginfo.SubjectID=loginfo.ParticipantDict.SubjectID(particpantIndex);
               end
               
               if(~isempty(loginfo.ParticipantDict.Source(particpantIndex)))
                   loginfo.Source=loginfo.ParticipantDict.Source(particpantIndex);
               end
           end
          
       elseif(contains(lineF,'Comments:'))
           lineF=fgetl(logfid);
           linecount=linecount+1;
           loginfo.Comments='';
           while(lineF~=-1)
               if(isempty(loginfo.Comments))
                   loginfo.Comments=lineF;
               else
                   loginfo.Comments=sprintf('%s\n%s',loginfo.Comments,lineF);
               end
               lineF=fgetl(logfid);
               linecount=linecount+1;
           end
       elseif(contains(lineF,'Marker Dictionary:'))
           lineF=fgetl(logfid);
           linecount=linecount+1;
           if(lineF~=-1)
               tableHeaders=strsplit(lineF,'\t');
               lastItem=tableHeaders{end};
               if(lastItem(end)==':')
                   lastItem=lastItem(1:end-1);
                   tableHeaders{end}=lastItem;
               end
               
           else
              continue; 
           end
           lineF=fgetl(logfid);
           linecount=linecount+1;
           
           cellMarkerDict=cell(1,length(tableHeaders));
           numMrkDictItems=1;
           
           while(lineF~=-1)
               items=strsplit(lineF,'\t');
               if(~isempty(items)&&iscell(items))
                   cellMarkerDict(numMrkDictItems,1:length(items))=items;
                   numMrkDictItems=numMrkDictItems+1;
               end
               lineF=fgetl(logfid);
               linecount=linecount+1;
           end
           loginfo.MarkerDict=cell2table(cellMarkerDict(:,1:length(tableHeaders)),'VariableNames',tableHeaders);
           
       end
       lineF=fgetl(logfid);
       linecount=linecount+1;
       %disp(lineF)

    end
    fclose(logfid);
       

end



 

    
