function validate_session(folder,saveflag)
% This function makes a report (PDF) about the data quality of files in a
% folder.

if(nargin==0)
    folder=uigetdir([], 'Pick a Subject Data Directory');
end

if(nargin<2)
    saveflag=false;
end

if(isa(folder,'nirs.core.Data'))
    raw=folder;
    folder=fileparts(raw(1).description);
else
    raw=nirs.io.loadDirectory(folder,{});
end

f=figure;
set(f,'MenuBar','none')
set(f,'DockControls','off')
set(f,'Units','normalized');
set(f,'Position',[0 0 1 1]);
set(f,'Color','w')

% make the summary table
for i=1:length(raw)
    filenames{i,1}=raw(i).description;
    [~,filenames{i,1}]=fileparts(filenames{i,1});
    
    data_length{i,1} = raw(i).time(end);
    s=nirs.getStimNames(raw(i));
    if(~isempty(s))
        ss=s{1};
        for j=2:length(s)
            ss=sprintf('%s : %s',ss,s{j});
        end
        stimulus_events{i,1}=ss;
    else
         stimulus_events{i,1}=[];
    end
    
    j=nirs.modules.FixNaNs;
    raw=j.run(raw);
    
    SNI{i}=nirs.math.structnoiseindex(raw(i).data);
    n=length(find(SNI{i}>3));
    types=unique(raw(i).probe.link.type);
    lst1=find(raw(i).probe.link.type==types(1));
    
    lst2=find(raw(i).probe.link.type==types(2));
    
    Number_of_Good_Channels{i,1}=sprintf('%d%s (%d of %d)',100*n/length(SNI{i}),'%',n,length(SNI{i}));
    DataQuality_830{i,1}=sprintf('%0.1f [%0.1f-%0.1f]',median(SNI{i}(lst1)),min(SNI{i}(lst1)),max(SNI{i}(lst1)));
    DataQuality_690{i,1}=sprintf('%d1.f [%0.1f-%0.1f]',median(SNI{i}(lst2)),min(SNI{i}(lst2)),max(SNI{i}(lst2)));
    
end

    


T=table(filenames,data_length,Number_of_Good_Channels,DataQuality_830,DataQuality_690,stimulus_events);

a=subplot(length(raw)+1,4,[1 2 3 4]);
axis off
s = uipanel('Title','Summary Info','FontSize',12,...
                'BackgroundColor','white',...
                'Position',get(a,'Position'));

% uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
%     'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
TString = evalc('disp(T)');
% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
annotation(s,'Textbox','String',TString,'Interpreter','Tex',...
    'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);


c=hot(64);
mS=max(horzcat(SNI{:}));
sp=linspace(0,mS,64);

c(find(sp<=2),1)=1;
c(find(sp<=2),2)=0;
c(find(sp<=2),3)=0;


c(find(sp>2 & sp<=5),1)=1;
c(find(sp>2 & sp<=5),2)=1;
c(find(sp>2 & sp<=5),3)=0;

c(find(sp>5),1)=0;
c(find(sp>5),2)=1;
c(find(sp>5),3)=0;


for i=1:length(raw)
    s=subplot(length(raw)+1,4,[4+4*(i-1)+1 4+4*(i-1)+2]);
    raw(i).draw;
    title([filenames{i}]);
     
     
    lst1=find(raw(i).probe.link.type==types(1));
    lst=dsearchn(sp',SNI{i}(lst1)');
    s=subplot(length(raw)+1,4,[4+4*(i-1)+3]);
    raw(i).probe.draw(c(lst,:),[],s)
    colormap(c);
    set(s,'clim',[0 mS])
    colorbar
    title(['SNI ' num2str(types(1))]);
    
    
    lst1=find(raw(i).probe.link.type==types(2));
    lst=dsearchn(sp',SNI{i}(lst1)');
    s=subplot(length(raw)+1,4,[4+4*(i-1)+4]);
    raw(i).probe.draw(c(lst,:),[],s)
    colormap(c);
    set(s,'clim',[0 mS])
    colorbar
       title(['SNI ' num2str(types(2))]);
end

if(saveflag)
    fileout=fullfile(folder,'summary');
    export_fig(fileout,'-jpeg');
end
    