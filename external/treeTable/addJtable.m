function jtable = addJtable(data,parent)

if(~isempty(data))

demo=nirs.createDemographicsTable(data.raw);

dd=get(parent,'userdata');
if(isempty(dd))
    dd.data=[];
end

flds=fields(data);
flds={flds{~ismember(flds,'raw')}};

DataTable={}; cnt=1;
groups=unique(demo.group);
for gI=1:length(groups)
    subjs=unique(demo(ismember(demo.group,groups{gI}),:).subject);
    for sI=1:length(subjs)
        files=unique(demo(ismember(demo.group,groups{gI}) &...
            ismember(demo.subject,subjs{sI}),:).file);
        for fI=1:length(files)
            DataTable{cnt,1}=groups{gI};
            DataTable{cnt,2}=subjs{sI};
            DataTable{cnt,3}=files{fI};
            DataTable{cnt,4}='raw';
            cnt=cnt+1;
            
            for ii=1:length(flds)
                demo2=nirs.createDemographicsTable(data.(flds{ii}));
                if(ismember(groups{gI},demo2.group) &...
                        ismember(subjs{sI},demo2.subject) &...
                        ismember(files{fI},demo2.file))
                    DataTable{cnt,1}=groups{gI};
                    DataTable{cnt,2}=subjs{sI};
                    DataTable{cnt,3}=files{fI};
                    DataTable{cnt,4}=flds{fI};
                    cnt=cnt+1;
                end
                
            end
        end
    end
end
jtable = treeTable(parent,{'Group','Subject','File','Type'},DataTable,...
    'ColumnTypes',{'char','char','char','char'},'groups',[true true true false]);
else
    jtable = javaObjectEDT('javax.swing.JTextArea', ...
        sprintf('Load data from pulldown menu or drag files here.\n\n'));
    
    % Create Java Swing JScrollPane
    jScrollPane = javaObjectEDT('javax.swing.JScrollPane', jtable);
    jScrollPane.setVerticalScrollBarPolicy(jScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
    
    % Add Scrollpane to figure
    [~,hContainer] = javacomponent(jScrollPane,[],parent);
    set(hContainer,'Units','normalized','Position',[0 0 1 1]);
    DataTable={};
end

dndcontrol.initJava();
dnd=dndcontrol(jtable);
dnd.DropFileFcn=@dndcallback;
dnd.DropStringFcn=@dndcallback;

dd.table=DataTable;

set(parent,'userdata',dd);


function dndcallback(varargin)

return
 