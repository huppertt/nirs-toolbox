function f = gui(obj);

if(length(obj)>1)
    warning('GUI only supports single stats entries: using only first value in array');
    obj=obj(1);
end

obj = obj.sorted();


f=figure;
set(f,'menu','none');

handles.axes = axes('parent',f,'units','normalized');
set(handles.axes,'Position',[.05 .3 .9 .65],'box','on','XTick',[],'YTick',[]);


try
    types=arrayfun(@(x){num2str(x)},obj(1).probe.link.type);
catch
    types=arrayfun(@(x){x{1}},obj(1).probe.link.type);
end
handles.selecttype = uicontrol('Style','popupmenu','String',unique(types),...
    'Units','normalized','Position',[.57 .18 .3 .1]);
set(handles.selecttype,'callback',@changetype);
uicontrol('Style','text','String','Data Type','Units','normalized',...
    'position',[.48 .18 .1 .1]);

cond=nirs.getStimNames(obj);
handles.selectcond = uicontrol('Style','popupmenu','String',cond,...
    'Units','normalized','Position',[.57 .1 .3 .1]);
a=uicontrol('Style','text','String','Contrast','Units','normalized',...
    'position',[.48 .1 .1 .1]);
set(handles.selectcond,'callback',@changetype);

handles.editcond = uicontrol('Style','pushbutton','Units','normalized','String','Edit','Position',[.87 .15 .1 .05]);
set(handles.editcond,'callback',@editcond);


handles.selectbeta = uicontrol('Style','popupmenu','String',{'tstat','beta'},...
    'Units','normalized','Position',[.57 .02 .3 .1]);
a=uicontrol('Style','text','String','Variable','Units','normalized',...
    'position',[.48 .02 .1 .1]);
set(handles.selectbeta,'callback',@changetype);


handles.selectthres = uicontrol('Style','popupmenu','String',{'p<','q<'},...
    'Units','normalized','Position',[.15 .1 .12 .1]);
set(handles.selectthres,'callback',@changetype);
a=uicontrol('Style','text','String','Threshold','Units','normalized',...
    'position',[.05 .1 .1 .1]);

handles.editthres = uicontrol('Style','edit','String','0.05',...
    'Units','normalized','Position',[.28 .15 .1 .05]);
set(handles.editthres,'callback',@changetype);

handles.checkFWE = uicontrol('Style','radiobutton','String','FWE-correct all data',...
            'units','normalized','Position',[.05 .05 .4 .1]);
set(handles.checkFWE,'callback',@changetype);

handles.data=obj;
set(f,'userdata',handles);

changetype(f);




return


function changetype(varargin)

handles=get(gcf,'userdata');
type=get(handles.selecttype,'String');
type=type{get(handles.selecttype,'Value')};

cond=get(handles.selectcond,'String');
cond=cond{get(handles.selectcond,'Value')};

stat=get(handles.selectbeta,'String');
stat=stat{get(handles.selectbeta,'Value')};

values = handles.data.(stat);

vmax    = max(abs(values(:)));
vrange  = vmax*[-1 1];

threstype= get(handles.selectthres,'String');
threstype=threstype{get(handles.selectthres,'value')};
threstype=threstype(1);
threstype = handles.data.(threstype);

if(~get(handles.checkFWE,'Value'))
    threstype(~ismember(handles.data.table.cond,cond))=NaN;
end

mask = (threstype<str2num(get(handles.editthres,'String')));


[~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',128) )');
z = linspace(vrange(1), vrange(2), size(cmap,1))';

lst = find(ismember(handles.data.variables.type,type) &...
      ismember(handles.data.variables.cond,cond)); 
  
  vals = values(lst);
  
  % this mask
  m = mask(lst);
  
  % map to colors
  idx = bsxfun(@minus, vals', z);
  [~, idx] = min(abs(idx), [], 1);
  
  colors = cmap(idx, :);
  
  % line styles
  lineStyles = {};
  for i = 1:length(idx)
      if m(i)
          lineStyles(i,:) = {'LineStyle', '-', 'LineWidth', 8};
      else
          lineStyles(i,:) = {'LineStyle', '--', 'LineWidth', 4};
      end
  end
  cla(handles.axes);
  handles.data.probe.draw(colors, lineStyles,handles.axes);
  
  c = colorbar('Eastoutside'); colormap(cmap); caxis(vrange);

  
  
return

function editcond(varargin)

handles=get(gcf,'userdata');
cont=nirs.getStimNames(handles.data);

prompt={'Create new contrast:'};
name='Edit Contrast';
numlines=1;
defaultanswer={cont{1}};
str=inputdlg(prompt,name,numlines,defaultanswer);
 
   
cont={cont{:} str{1}};
cont=unique(cont);

lastwarn('');
handles.data=handles.data.ttest(cont);
[a,b]=lastwarn;
if(isempty(a))
    cont=nirs.getStimNames(handles.data);
    set(handles.selectcond,'String',cont);
else
    warning('invalid contrast');
end

set(gcf,'userdata',handles);
changetype(gcf);


return