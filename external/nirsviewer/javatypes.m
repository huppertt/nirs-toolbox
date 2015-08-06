function prop=javatypes(types,varargin)

%This is needed because all the editors need unique ID names
uniq=datestr(now,'yyyymmddHHMMSSFFF');

warning('off','MATLAB:hg:PossibleDeprecatedJavaSetHGProperty');
prop = com.jidesoft.grid.DefaultProperty();

switch(types)
    case('double')
         javatype=javaclass('double');
        set(prop,'Type',javatype);
    case('logical')
        javatype=javaclass('logical');
        editor=com.jidesoft.grid.BooleanCheckBoxCellEditor;
        context = com.jidesoft.grid.EditorContext(['comboboxeditor' uniq]);
        com.jidesoft.grid.CellEditorManager.initDefaultEditor();
        com.jidesoft.grid.CellEditorManager.registerEditor(javatype, editor,context);
        set(prop,'Type',javatype,'EditorContext',context);
    case('int')
        javatype=javaclass('int32');
        set(prop,'Type',javatype);
    case('float')
        javatype=javaclass('double');
        set(prop,'Type',javatype);
    case('char')
        javatype=javaclass('char',1);
        set(prop,'Type',javatype);
    case('enum')
        javatype = javaclass('char', 1);
       model = javax.swing.SpinnerListModel(varargin{1});
        
        editor = com.jidesoft.grid.SpinnerCellEditor(model);
        context = com.jidesoft.grid.EditorContext(['spinnereditor' uniq]);
        com.jidesoft.grid.CellEditorManager.registerEditor(javatype, editor, context);
        set(prop,'Type',javatype,'EditorContext',context);
    case('spinner')
        javatype=javaclass('int32');
        val=varargin{1};
        spinner = javax.swing.SpinnerNumberModel(val(1),val(2),val(3),val(4));
        editor = com.jidesoft.grid.SpinnerCellEditor(spinner);
        context = com.jidesoft.grid.EditorContext(['spinnereditor' uniq]);
        com.jidesoft.grid.CellEditorManager.registerEditor(javatype, editor, context);
        set(prop,'Type',javatype,'EditorContext',context);

end

%      case 'logical'  % logical array of true and false values
%             jclassname = 'java.lang.Boolean';
%         case 'char'  % character array
%             jclassname = 'java.lang.Character';
%         case {'int8','uint8'}  % 8-bit signed and unsigned integer array
%             jclassname = 'java.lang.Byte';
%         case {'int16','uint16'}  % 16-bit signed and unsigned integer array
%             jclassname = 'java.lang.Short';
%         case {'int32','uint32'}  % 32-bit signed and unsigned integer array
%             jclassname = 'java.lang.Integer';
%         case {'int64','uint64'}  % 64-bit signed and unsigned integer array
%             jclassname = 'java.lang.Long';
%         case 'single'  % single-precision floating-point number array
%             jclassname = 'java.lang.Single';
%         case 'double'  % double-precision floating-point number array
%             jclassname = 'java.lang.Double';
%         case 'complexsingle'
%             jclassname = 'hu.bme.aut.matlab.ComplexF';
%         case 'complexdouble'
%             jclassname = 'hu.bme.aut.matlab.Complex';

% 
% warning('off','MATLAB:hg:PossibleDeprecatedJavaSetHGProperty');
% switch(optval.type)
%     case('enum')
%         javatype = javaclass('char', 1);
%         options = optval.typesub;
%         
%         model = javax.swing.SpinnerListModel(options);
%         
%         editor = com.jidesoft.grid.SpinnerCellEditor(model);
%         context = com.jidesoft.grid.EditorContext(['spinnereditor' num2str(cnt)]);
%         com.jidesoft.grid.CellEditorManager.registerEditor(javatype, editor, context);
%         set(prop,'Type',javatype,'EditorContext',context);
%         
%     case('old-enum')
%         %This broke on Matlab 2014b on Mac OS Yosimite.  Fix this later.
%         javatype = javaclass('char', 1);
%         options = optval.typesub;
%         editor = com.jidesoft.grid.ListComboBoxCellEditor(options);
%         context = com.jidesoft.grid.EditorContext(['comboboxeditor' num2str(cnt)]);
%         com.jidesoft.grid.CellEditorManager.initDefaultEditor();
%         com.jidesoft.grid.CellEditorManager.registerEditor(javatype, editor,context);
%         % set(prop,'Type',javatype, 'EditorContext',context);
%         set(prop,'Type',javatype);
%     case('logical')
%         javatype=javaclass('logical');
%         editor=com.jidesoft.grid.BooleanCheckBoxCellEditor;
%         context = com.jidesoft.grid.EditorContext(['comboboxeditor' num2str(cnt)]);
%         com.jidesoft.grid.CellEditorManager.initDefaultEditor();
%         com.jidesoft.grid.CellEditorManager.registerEditor(javatype, editor,context);
%         set(prop,'Type',javatype,'EditorContext',context);
%         %set(prop,'Type',javatype);
%     case('int')
%         javatype=javaclass('int32');
%         set(prop,'Type',javatype);
%     case('float')
%         javatype=javaclass('double');
%         set(prop,'Type',javatype);
%     case('spinner')
%         javatype=javaclass('int32');
%         spinner = javax.swing.SpinnerNumberModel(optval.default,optval.typesub(1),optval.typesub(end),mean(diff(optval.typesub)));
%         editor = com.jidesoft.grid.SpinnerCellEditor(spinner);
%         context = com.jidesoft.grid.EditorContext(['spinnereditor' num2str(cnt)]);
%         com.jidesoft.grid.CellEditorManager.registerEditor(javatype, editor, context);
%         set(prop,'Type',javatype,'EditorContext',context);
%     case('char')
%         javatype=javaclass('char',1);
%         set(prop,'Type',javatype);
% end
% 
% set(prop,'Name',optval.name,'Value',optval.default);
% set(prop,'Category',optval.category);
% set(prop,'Description',optval.help);
% opts=setfield(opts,optval.name,optval.default);
% end