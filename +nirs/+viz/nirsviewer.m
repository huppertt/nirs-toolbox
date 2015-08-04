function nirsviewer(data)

h=nirs_viewer;
handles=guihandles(h);

%% Populate the filename/subject/group info

table=nirs.createDemographicsTable(data);
delete(findobj('tag','uitree_cont'));
delete(findobj('tag','uitree_obj'));

import javax.swing.*
import javax.swing.tree.*;
root = uitreenode('v0','Subjects', 'Subjects', [], false);

groups=unique(table.group);
for gIdx=1:length(groups)
    g = uitreenode('v0', groups{gIdx},  groups{gIdx}, [], false);
    lst=find(ismember(table.group, groups{gIdx}));
    subj=table.subject;
    subj=unique({subj{lst}});
    for sIdx=1:length(subj)

        s = uitreenode('v0', subj{sIdx},  subj{sIdx}, [], false);
        lstfiles=find(ismember(table.group, groups{gIdx}) & ismember(table.subject, subj{sIdx}));
        for cf=1:length(lstfiles)
            [~,file,ext]=fileparts(data(lstfiles(cf)).description);
            f = uitreenode('v0', num2str(lstfiles(cf)),[file ext],[],true);
            s.add(f);
        end
        g.add(s);
    end
    root.add(g);
end

[mtree,container] = uitree('v0', 'Root', root);
set(container,'units',get(handles.uipanel1,'units'));
set(container,'tag','uitree_cont');

set(container,'position',get(handles.uipanel1,'position'));
set(container,'visible','on');
set(container,'userdata',mtree);

mtree.expand(mtree.getRoot);
mtree.expand(s);
mtree.setSelectedNode(f);
set(mtree,'NodeSelectedCallback',@updatewindow);


set(handles.figure_nirsview,'Userdata',data);
updatewindow;
return

function updatewindow(varargin)

handles=guihandles(findobj('tag','figure_nirsview'));

a=findobj('tag','uitree_cont');
a=get(a,'Userdata');
node=get(a.Tree,'LastSelectedPathComponent');

filename=node.getName;
val=str2num(node.getValue);
data=get(handles.figure_nirsview,'Userdata');

cla(handles.axes_SDG);
axes(handles.axes_SDG);
hold on;
SDGhandlesBase=data(val).probe.draw;
set(SDGhandlesBase,'Color',[.8 .8 .8]);
SDGhandles=data(val).probe.draw;

cla(handles.axes_maindata);
axes(handles.axes_maindata);
linehandles=data(val).draw;

typesAll=arrayfun(@(x){num2str(x)},data(val).probe.link.type);
types=unique(typesAll);
lstdisp=[1];

SDcolors=makeSDcolors(data(val).probe.link);

set(linehandles,'visible','off');
for idx=1:length(lstdisp)
    lst=find(ismember(typesAll,types{lstdisp(idx)}));
    set(linehandles(lst),'visible','on');
    for idx2=1:length(lst)
        set(linehandles(lst(idx2)),'color',SDcolors(idx2,:));
    end
    LstAll(:,idx)=lst;
end
global linelinks; 
clear linelinks;
for idx=1:length(SDGhandles)
    set(SDGhandles(idx),'tag',['SDG' num2str(idx)]);
    set(SDGhandlesBase(idx),'tag',['SDG' num2str(idx)]);
    set(SDGhandles(idx),'color',SDcolors(idx,:));
    linelinks(idx)=linkprop([linehandles(LstAll(idx,:)),SDGhandles(idx)],{'Visible','Color'});
end
for idx=1:size(LstAll,2)
    linelinks(end+1)=linkprop(linehandles(LstAll(:,idx)),'LineStyle');
end

set(SDGhandles,'ButtonDownFcn','set(gcbo,''visible'',''off'')');
set(SDGhandlesBase,'ButtonDownFcn','set(findobj(''tag'',get(gcbo,''tag'')),''visible'',''on'')');
setappdata(handles.figure_nirsview,'linelinks',linelinks)
return
