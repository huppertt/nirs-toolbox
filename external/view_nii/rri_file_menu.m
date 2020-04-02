%  Imbed a file menu to any figure. If file menu exist, it will append
%  to the existing file menu. This file menu includes: Copy to clipboard,
%  print, save, close etc.
%
%  Usage: rri_file_menu(fig);
%
%         rri_file_menu(fig,0) means no 'Close' menu.
%
%  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)
%
%--------------------------------------------------------------------

function rri_file_menu(action, varargin)

   if isnumeric(action)
      fig = action;
      action = 'init';
   end

   %  clear the message line,
   %
   h = findobj(gcf,'Tag','MessageLine');
   set(h,'String','');

   if ~strcmp(action, 'init')
      set(gcbf, 'InvertHardcopy','off');
%      set(gcbf, 'PaperPositionMode','auto');
   end

   switch action
      case {'init'}
         if nargin > 1
            init(fig, 1);		% no 'close' menu
         else
            init(fig, 0);
         end
      case {'print_fig'}
         printdlg(gcbf);
      case {'copy_fig'}
         copy_fig;
      case {'export_fig'}
         export_fig;
   end

   return					% rri_file_menu


%------------------------------------------------
%
%  Create (or append) File menu
%
function init(fig, no_close)

   %  search for file menu
   %
   h_file = [];
   menuitems = findobj(fig, 'type', 'uimenu');

   for i=1:length(menuitems)
      filelabel = get(menuitems(i),'label');

      if strcmpi(strrep(filelabel, '&', ''), 'file')
         h_file = menuitems(i);
         break;
      end
   end

   set(fig, 'menubar', 'none');

   if isempty(h_file)
      if isempty(menuitems)
         h_file = uimenu('parent', fig, 'label', 'File');
      else
         h_file = uimenu('parent', fig, 'label', 'Copy Figure');
      end

      h1 = uimenu('parent', h_file, ...
         'callback','rri_file_menu(''copy_fig'');', ...
         'label','Copy to Clipboard');
   else
      h1 = uimenu('parent', h_file, ...
         'callback','rri_file_menu(''copy_fig'');', ...
         'separator','on', ...
         'label','Copy to Clipboard');
   end

   h2 = uimenu(h_file, ...
      'callback','pagesetupdlg(gcbf);', ...
      'label','Page Setup...');

   h2 = uimenu(h_file, ...
      'callback','printpreview(gcbf);', ...
      'label','Print Preview...');

   h2 = uimenu('parent', h_file, ...
      'callback','printdlg(gcbf);', ...
      'label','Print Figure ...');

   h2 = uimenu('parent', h_file, ...
      'callback','rri_file_menu(''export_fig'');', ...
      'label','Save Figure ...');
  
  
  
%Yi Sui file browser
    h2 = uimenu('parent', h_file, ...
      'callback',@open_nii, ...
      'label','Open nii file...');
  
     h2 = uimenu('parent', h_file, ...
      'callback',@file_browser, ...
      'label','File Browser');
   arch = computer;
   if ~strcmpi(arch(1:2),'PC')
      set(h1, 'enable', 'off');
   end

   if ~no_close
      h1 = uimenu('parent', h_file, ...
         'callback','close(gcbf);', ...
         'separator','on', ...
         'label','Close');
   end

   return;					% init

%------------------------------------------------
function file_browser(varargin)

Obj = varargin{1};
h_nii_view=gcbf;

if strcmp( get(Obj,'Selected') , 'on')
    set(Obj,'Label','File Browser',...
        'Selected', 'off');
    file_browser_objs = findobj(h_nii_view,'Tag','file_browser');
    delete(file_browser_objs);
    pos = [0.05 0.05 0.95 0.9];
    option.setarea=pos;
    option.usestretch=0;
    view_nii(h_nii_view,option);
    
else % 'off'
    set(Obj,'Label','Close File Browser',...
        'Selected', 'on');
    
    
    pos = [0.3 0.05 0.7 0.9];
    option.setarea=pos;
    option.usestretch=0;
    
    view_nii(h_nii_view,option);
    root = pwd;
    [t container] = uitree('v0','Root', [root filesep],'SelectionChangeFcn',@uitree_selection,'parent',h_nii_view);
    set(container,'Tag','file_browser','parent',h_nii_view);
    set(t,'Units','normalized','Position',[0 0.1 pos(1)*0.98 0.89 ]);
    t.expand(t.getRoot);
    btn = uicontrol('style','pushbutton','String','Change Directory',...
        'Callback',@(varargin) pushbutton1_callback(varargin{:},t),...
        'Units','Normalized',...
        'Position',[0 0 pos(1)*0.98 0.09],...
        'Tag','file_browser');
end
            
function pushbutton1_callback(varargin)
tree = varargin{end};
folder_name=uigetdir;
folder_name = [folder_name filesep];
rootnode = uitreenode('v0',folder_name,folder_name,[matlabroot, '/toolbox/matlab/icons/foldericon.gif'], 0);
tree.setRoot(rootnode);
tree.expand(rootnode);
        
            
            


function uitree_selection(hObj,data)
%  h_nii_view = getappdata(gcf,'h_nii_view');
fpath = hObj.Tree.getSelectionPath.getPath;
fname = fpath(end).getValue; 
fname = fname(1:end-1); 
% 
if exist(fname,'file')==2 

%     nii = load_untouch_nii(fname);
  try
    nii = load_nii(fname);
    opt.command = 'updateimg';
    opt.command = 'updatenii';
    pos = [0.4 0.05 0.55 0.9];
    opt.setarea=pos;
    opt.usestretch=0;
    view_nii(gcf,nii,opt);
  catch err
  end
end




%------------------------------------------------
function listbox1_callback(varargin)

hObject = varargin{1};
    index_selected = get(hObject,'Value');
    file_list = get(hObject,'String');
    filename = file_list{index_selected};
    curdir = getappdata(hObject,'curdir');
    
    fullname = fullfile(curdir,filename);
    
    if  isdir(fullname)
        
        if strcmp(filename,'..')
            curdir = fileparts(curdir);
            setappdata(hObject,'curdir',curdir);
%         cd (filename)
            load_listbox(curdir,hObject)
        elseif strcmp(filename,'.')
        else
            setappdata(hObject,'curdir',fullname);
%         cd (filename)
            load_listbox(fullname,hObject)
        end
    end
% ------------------------------------------------------------
function load_listbox(dir_path,handles)
% cd (dir_path)
dir_struct = dir(dir_path);
% [sorted_names,sorted_index] = sortrows({dir_struct.name}');
% dir_struct(1).desc='';
% dir_struct(2).desc='';
% for k=3:numel(dir_struct)
%    dir_struct(k).desc='';
%    if dir_struct(k).isdir
%        
%        dcmfiles = dir(fullfile(dir_path,dir_struct(k).name));
%        for j =3:length(dcmfiles);
%            if dcmfiles(j).isdir
%                continue;
%            end
%            if strcmp(dcmfiles(j).name,'DICOMDIR')
%                continue;
%            end
%            fname=fullfile(dir_path,dir_struct(k).name,dcmfiles(j).name);
%            if isdicom(fname)
%                dinfo=dicominfo(fname);
%                try
%                dir_struct(k).desc=['    ''' dinfo.SeriesDescription ''''];
%                catch
%                    continue;
%                end
%                break;
%            end
%        end
%    end
% end


[sorted_names,sorted_index] = sort_nat({dir_struct.name}');
% handles.file_names = sorted_names;
% handles.file_desc = {dir_struct(sorted_index).desc}';
% handles.is_dir = [dir_struct.isdir];
% handles.sorted_index = sorted_index;
% guidata(handles.figure1,handles);
% set(handles.listbox1,'String', cellfun(@(x,y)[x y],handles.file_names,handles.file_desc,'UniformOutput',false),...
% 	'Value',1)
set(handles,'String', sorted_names,...
	'Value',1)
% set(handles.NofFiles,'String',handles.dfolder)




%------------------------------------------------
function button1_callback(varargin)
hdl = getappdata(gcf,'hdl');
load_listbox(pwd,hdl.listbox1);
setappdata(hdl.listbox1,'curdir',pwd);






%------------------------------------------------
function open_nii(varargin)
[filename, pathname] = uigetfile({'*.nii' '*.nii.gz'}, 'Select NII file');
nii = load_nii(fullfile(pathname,filename));
% opt.command = 'updateimg';
% view_nii(gcf,nii.img,opt);
opt.command = 'updatenii';
view_nii(gcf,nii,opt);

%------------------------------------------------
%
%  Copy to clipboard
%
function copy_fig

   arch = computer;
   if(~strcmpi(arch(1:2),'PC'))
      error('copy to clipboard can only be used under MS Windows');
      return;
   end

   print -noui -dbitmap;

   return					% copy_fig


%------------------------------------------------
%
%  Save as an image file
%
function export_fig

   curr = pwd;
   if isempty(curr)
      curr = filesep;
   end

   [selected_file, selected_path] = rri_select_file(curr,'Save As');

   if isempty(selected_file) | isempty(selected_path)
      return;
   end

   filename = [selected_path selected_file];

   if(exist(filename,'file')==2)		% file exist

      dlg_title = 'Confirm File Overwrite';
      msg = ['File ',filename,' exist. Are you sure you want to overwrite it?'];
      response = questdlg(msg,dlg_title,'Yes','No','Yes');

      if(strcmp(response,'No'))
         return;
      end

   end

   old_pointer = get(gcbf,'pointer');
   set(gcbf,'pointer','watch');

   try
      saveas(gcbf,filename);
   catch
      msg = 'ERROR: Cannot save file';
      set(findobj(gcf,'Tag','MessageLine'),'String',msg);
   end

   set(gcbf,'pointer',old_pointer);

   return;					% export_fig

