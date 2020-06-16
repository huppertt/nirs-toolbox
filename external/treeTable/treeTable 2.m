function jtable = treeTable(varargin)
% treeTable - create Java-based tree-table based on javax.swing.JTable and javax.swing.JTree
%
% Syntax:
%    jtable = treeTable (pnContainer, headers, data, 'PropName',PropValue, ...)
%
% Input Parameters:
%    pnContainer - optional handle to container uipanel or figure. If empty/unsupplied then current figure will be used
%    headers     - optional cell array of column header strings. If unsupplied then = {'A','B','C'}
%    data        - optional vector/matrix (either scalar or cell array) of data values
%    'PropName',PropValue - 
%                  optional list of property pairs (e.g., 'iconFilenames',{'a.gif','b.jpg'},'columnTypes',{'char','logical'})
%                  Note: All optional parameters of treeTable may be specified using PropName/PropValue pairs,
%                        case-insensitive, in whichever order (see the bottom example below):
%                        - 'Container'      => HG handle             (default=gcf)
%                        - 'Headers'        => cell array of labels  (default={'A','B',...} based on data size)
%                        - 'Data'           => 2D cell/numeric array (default=[])
%                        - 'IconFilenames'  => filepath strings      (default={leafIcon,folderClosedIcon,folderOpenIcon})
%                        - 'ColumnTypes'    => 'char'/'logical'/{}   (default={'char','char',...})
%                              Note 1: {'a','b','c'} indicates a non-editable combo-box with the specified data items
%                              Note 2: {'a','b','c', ''} indicates *editable* combo-box with the specified data items
%                        - 'ColumnEditable' => array of true/false   (default=[true,true,...])
%                        - 'Groupable'           => logical flag     (default=true; if false, display as non-groupable table)
%                        - 'InteractiveGrouping' => logical flag     (default=false; if true, enables interactive grouping as in Outlook)
%
% Output parameters:
%    jtable      - handle to Java tree-table object
%
% Examples:
%    jtable = treeTable;  % show demo with default parameter values
%
%    jtable = treeTable(gcf, 'column name');
%
%    data = {1,'M11',true,true,false; 1,'M12',true,false,true; 1,'M13',false,true,true; 2,'M21',true,true,false; 2,'M22',false,true,false;};
%    jtable = treeTable(figure,{'Group','Panel','Mask','Object','ID'},data);
%    jtable = treeTable(figure,{'Group','Panel','Mask','Object','ID'},data,'ColumnTypes',{'','char','logical','logical','logical'});
%
%    data = {1,'M11',true,true,'Yes'; 1,'M12',true,false,'No'; 1,'M13',false,true,'No'; 2,'M21',true,true,'Yes'; 2,'M22',false,true,'Maybe';};
%    jtable = treeTable(figure,{'Group','Panel','Mask','Object','ID'},data,'ColumnTypes',{'','char','logical','logical',{'No','Yes',''}});
%
% Usage:
%    The table is sortable.
%    The table automatically resizes to fill the pnContainer (you may modify this via the 'Position' property).
%    The table automatically sets the columns' cell editor and renderer based on the supplied data. Logical values are
%       given a checkbox, strings are left-aligned (numbers are right-aligned). You can always override the defaults.
%    You can change column widths by dragging the column borders.
%    You can exchange columns by simply dragging the column header right or left.
%    For additional tips about how to set multiple aspects of the table, refer to:
%       <a href="http://java.sun.com/docs/books/tutorial/uiswing/components/table.html">http://java.sun.com/docs/books/tutorial/uiswing/components/table.html</a>
%
% Bugs and suggestions:
%    Please send to Yair Altman (altmany at gmail dot com)
%
% See also:
%    uitable, uitree, java, javaclasspath
%
% Release history:
%    1.0 2011-01-02: Initial version
%    1.1 2011-01-04: Added leaf/node icons
%    1.2 2013-06-21: Adaptations for R2013b & HG2, supported uiextras.Panel parent, enabled multi-column sorting, added groupable flag, supported numeric data
%    1.3 2013-08-04: Initial version posted on File Exchange: http://www.mathworks.com/matlabcentral/fileexchange/index?term=authorid%3A27420
%    1.4 2013-08-06: Added InteractiveGrouping option

% License to use and modify this code is granted freely to all interested, as long as the original author is
% referenced and attributed as such. The original author maintains the right to be solely associated with this work.

% Programmed and Copyright by Yair M. Altman: altmany(at)gmail.com
% $Revision: 1.4 $  $Date: 2013/08/06 14:31:16 $

  %try
      % Ensure that java swing is enabled...
      if ~usejava('swing')
          error('treeTable:NeedSwing','Java tables require Java Swing.');
      end
      import javax.swing.*

      % Process optional arguments
      paramsStruct = processArgs(varargin{:});

      if isa(handle(paramsStruct.container), 'figure')
          pnContainerPos = getpixelposition(paramsStruct.container,0);  % Fix for Matlab 7.0.4 as per Sebastian Hölz
          pnContainerPos(1:2) = 0;
      else
          try
              pnContainerPos = getpixelposition(paramsStruct.container,1);  % Fix for Matlab 7.0.4 as per Sebastian Hölz
          catch
              pnContainerPos = getpixelposition(paramsStruct.container);  % Fix for uiextras.Panel
          end
      end

      % Get handle to parent figure
      hFig = ancestor(paramsStruct.container,'figure');

      % Get the uitable's required position within the container
      margins = [1,1,0,0];
      tablePosition = pnContainerPos + margins;    % Relative to the figure

      % Create a sortable uitable within the container
      try
          % Use JideTable if available on this system
          com.mathworks.mwswing.MJUtilities.initJIDE;

          % Prepare the tree-table with the requested data & headers
          %model = javax.swing.table.DefaultTableModel(paramsStruct.data, paramsStruct.headers);
          try
              model = MultiClassTableModel(paramsStruct.data, paramsStruct.headers);  %(model)
          catch
              try
                  javaaddpath(fileparts(mfilename('fullpath')));
                  model = MultiClassTableModel(paramsStruct.data, paramsStruct.headers);  %(model)
              catch
                  % Revert to the default table model
                  % (which has problematic numeric sorting since it does not recognize numeric data columns)
                  err = lasterror;
                  model = javax.swing.table.DefaultTableModel(paramsStruct.data, paramsStruct.headers);
              end
          end
          jtable = eval('com.jidesoft.grid.GroupTable(model);');  % prevent JIDE alert by run-time (not load-time) evaluation
          jtable = handle(javaObjectEDT(jtable), 'CallbackProperties');
          jtable.setRowAutoResizes(true);
          jtable.setColumnAutoResizable(true);
          jtable.setColumnResizable(true);
          jtable.setShowGrid(false);

          % Wrap the standard model in a JIDE GroupTableModel
          %model = jtable.getModel;
          model = com.jidesoft.grid.DefaultGroupTableModel(model);
          %model = StyledGroupTableModel(jtable.getModel);

          % Automatically group by the first column (only if it has multiple value)
          shouldGroup = paramsStruct.groupable;
          %{
          col0Data = paramsStruct.data(:,1);
          try
              shouldGroup = length(unique(col0Data)) > 1;
          catch
              col0Data = cellfun(@num2str,paramsStruct.data(:,1),'Uniform',false);
              shouldGroup = length(unique(col0Data)) > 1;
          end
          %}
          
          
          if(isfield(paramsStruct,'groups') && ~isempty(paramsStruct.groups))
              for i=1:length(paramsStruct.groups)
                  if(paramsStruct.groups(i))
                      model.addGroupColumn(i-1);
                     
                      %paramsStruct.columntypes(1) = [];
                  end
                   model.groupAndRefresh;
              end
          else
              if shouldGroup
                  model.addGroupColumn(0);
                  model.groupAndRefresh;
                  paramsStruct.columntypes(1) = [];
              end
          end
          
          jtable.setModel(model);

          % Enable multi-column sorting
          jtable.setSortable(true);

          % Automatically resize all columns - this can be extremely SLOW for >1K rows!
          %jideTableUtils = eval('com.jidesoft.grid.TableUtils;');  % prevent JIDE alert by run-time (not load-time) evaluation
          %jideTableUtils.autoResizeAllColumns(jtable);

          % Set default cell renderers/editors based on the requested column types
          for colIdx = 0 : length(paramsStruct.columntypes)-1
              setColumnRenderersEditors(jtable,colIdx,paramsStruct);
          end
          jtable.setSelectionBackground(java.awt.Color(0.9608*.8,0.9608*.8, 0.9608));  % light-blue

          % Modify the group style (doesn't work on new Matlab releases)
          %{
          try
              iconPath = paramsStruct.iconfilenames{2};
              groupStyle = model.getGroupStyle;
              groupStyle.setBackground(java.awt.Color(.7,.7,.7));  % light-gray
              if ~isempty(iconPath)
                  icon = javax.swing.ImageIcon(iconPath);
                  groupStyle.setIcon(icon);
              end
          catch
              %fprintf(2, 'Invalid group icon: %s (%s)\n', char(iconPath), lasterr);
              a=1;   % never mind - probably an invalid icon
          end
          %}
          try
              jtable.setExpandedIcon (javax.swing.ImageIcon(paramsStruct.iconfilenames{2}));
              jtable.setCollapsedIcon(javax.swing.ImageIcon(paramsStruct.iconfilenames{3}));
          catch
              %fprintf(2, 'Invalid group icon: %s (%s)\n', char(iconPath), lasterr);
              a=1;   % never mind - probably an invalid icon
          end

          % Attach a GroupTableHeader so that we can use Outlook-style interactive grouping
          try
              jTableHeader = com.jidesoft.grid.GroupTableHeader(jtable);
              jtable.setTableHeader(jTableHeader);
              if paramsStruct.interactivegrouping
                  jTableHeader.setGroupHeaderEnabled(true);
              end
          catch
              warning('YMA:treeTable:InteractiveGrouping','InteractiveGrouping is not supported - try using a newer Matlab release');
          end

          % Present the tree-table within a scrollable viewport on-screen
          scroll = javaObjectEDT(JScrollPane(jtable));
          hParent = paramsStruct.container;
          try
              % HG2 sometimes needs double(), sometimes not, so try both of them...
              [jhscroll,hcontainer] = javacomponent(scroll, tablePosition, double(hParent));
          catch
              [jhscroll,hcontainer] = javacomponent(scroll, tablePosition, hParent);
          end
          set(hcontainer,'units','normalized','pos',[0,0,1,1]);  % this will resize the table whenever its container is resized
          pause(0.05);
      catch
          err = lasterror;
          hcontainer = [];
      end

      % Fix for JTable focus bug : see http://bugs.sun.com/bugdatabase/view_bug.do;:WuuT?bug_id=4709394
      % Taken from: http://xtargets.com/snippets/posts/show/37
      jtable.putClientProperty('terminateEditOnFocusLost', java.lang.Boolean.TRUE);

      % Store the uitable's handle within the container's userdata, for later use
      try
          % add to parent userdata, so we have a handle for deletion
          set(paramsStruct.container,'userdata',[get(paramsStruct.container,'userdata'), jtable]);
      catch
          try set(paramsStruct.container,'userdata',jtable); catch, end
      end

      % Enable multiple row selection, auto-column resize, and auto-scrollbars
      %scroll = mtable.TableScrollPane;
      scroll.setVerticalScrollBarPolicy(scroll.VERTICAL_SCROLLBAR_AS_NEEDED);
      scroll.setHorizontalScrollBarPolicy(scroll.HORIZONTAL_SCROLLBAR_AS_NEEDED);
      jtable.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);

      % Comment the following line in order to prevent column resize-to-fit
      jtable.setAutoResizeMode(jtable.java.AUTO_RESIZE_SUBSEQUENT_COLUMNS)
      %jtable.setAutoResizeMode(jtable.java.AUTO_RESIZE_OFF)

      % Set the jtable name based on the containing panel's tag
      try
          basicTagName = get(paramsStruct.container,'tag');
          jtable.setName([basicTagName 'Table']);
      catch
          % never mind...
      end

      % Move the selection to first table cell (if any data available)
      if (jtable.getRowCount > 0)
          jtable.changeSelection(0,0,false,false);
      end

      % Set default editing & selection callbacks
      % {
      try
          oldWarnState = warning('off','MATLAB:hg:JavaSetHGProperty');
          %set(handle(getOriginalModel(jtable),'CallbackProperties'), 'TableChangedCallback', {@tableChangedCallback, jtable});
          set(handle(jtable.getSelectionModel,'CallbackProperties'), 'ValueChangedCallback', {@selectionCallback,    jtable});
          warning(oldWarnState);
      catch
          a=1;  % never mind...
      end
      %}
      
      % Process optional args
      processParams(paramsStruct,hcontainer,jtable);

      % Fix for deployed app - show docking control
      set(hFig,'DockControls','on');
  %catch
      % Insert your code here
      %handleError;
  %end
end  % treeTable

%% Set the column renderers/editors
function setColumnRenderersEditors(jtable,colIdx,paramsStruct)

    import javax.swing.*

    % Set the column's cellRenderer and editor based on the declared ColumnTypes
    try
        colType = paramsStruct.columntypes{colIdx+1};
        cellRenderer = getDefaultCellRenderer();  % default cellRenderer = label

        % Cell array => combo-box
        if iscellstr(colType)
            % Combo-box editor (no need for a special renderer - use default label)
            emptyIdx = cellfun('isempty',colType);
            colType(emptyIdx) = [];
            editableFlag = any(emptyIdx);
            cb = JComboBox(colType);
            cb.setEditable(editableFlag);
            cbe = DefaultCellEditor(cb);
            try
                % If non-editable, disable combo-box
                if ~paramsStruct.columneditable{colIdx+1}
                    cbe.setClickCountToStart(1e6);
                end
            catch
            end
            jtable.getColumnModel.getColumn(colIdx).setCellEditor(cbe);
            jtable.getColumnModel.getColumn(colIdx).setCellRenderer(cellRenderer);

        % 'logical' => checkbox
        elseif strcmpi(colType,'logical')
            % Checkbox editor
            cbe = DefaultCellEditor(JCheckBox);
            cbe.getComponent.setHorizontalAlignment(SwingConstants.CENTER);
            try
                % If non-editable, disable checkbox
                if ~paramsStruct.columneditable{colIdx+1}
                    cbe.setClickCountToStart(1e6);
                end
            catch
            end
            jtable.getColumnModel.getColumn(colIdx).setCellEditor(cbe);
            
            % Checkbox renderer
            cellRenderer = javaObject('javax.swing.JTable$BooleanRenderer');
            cellRenderer.setHorizontalAlignment(SwingConstants.CENTER);
            jtable.getColumnModel.getColumn(colIdx).setCellRenderer(cellRenderer);

        % 'label' or 'char' => label
        elseif strcmpi(colType,'label') || strcmpi(colType,'char')
            try
                % If non-editable, disable checkbox
                if ~paramsStruct.columneditable{colIdx+1}
                    jtf = JTextField;
                    jtf.setEditable(false);
                    jte = DefaultCellEditor(jtf);
                    jte.setClickCountToStart(intmax);
                    jtable.getColumnModel.getColumn(colIdx).setCellEditor(jte);
                end
            catch
            end
            jtable.getColumnModel.getColumn(colIdx).setCellRenderer(cellRenderer);
            
        %else
            % never mind - leave as-is (label)
        end
    catch
        err = lasterror;  % never mind
        a=1;
    end
    %return;
    
    % The first column should have an icon
    if colIdx == 0 && paramsStruct.groupable
        try
            iconPath = paramsStruct.iconfilenames{1};
            if ~isempty(iconPath)
                icon = javax.swing.ImageIcon(iconPath);
                cellRenderer.setIcon(icon);
            end
        catch
            fprintf(2, 'Invalid leaf icon: %s (%s)\n', char(iconPath), lasterr);
            a=1;   % never mind - probably an invalid icon
        end
    end
end  % setColumnRenderersEditors

%% Get the basic JTable data model
function originalModel = getOriginalModel(jtable)
    originalModel = jtable.getModel;
    try
        while(true)
            originalModel = originalModel.getActualModel;
        end;
    catch
        a=1;  % never mind - bail out...
    end
end  % getOriginalModel

%% Process optional arguments
function paramsStruct = processArgs(varargin)

    % Fix edge-case
    if nargin>=2 && ischar(varargin{2})
        varargin{2} = {varargin{2}};
    end

    % Get the properties in either direct or P-V format
    [regParams, pvPairs] = parseparams(varargin);

    % Now process the optional P-V params
    try
        % Initialize
        paramName = [];
        paramsStruct = [];
        paramsStruct.container = [];
        paramsStruct.headers = {' '};  % 5 columns by default
        paramsStruct.data = {};
        paramsStruct.iconfilenames = {fullfile(matlabroot,'/toolbox/matlab/icons/greenarrowicon.gif'), ...
                                      fullfile(matlabroot,'/toolbox/matlab/icons/file_open.png'), ...
                                      fullfile(matlabroot,'/toolbox/matlab/icons/foldericon.gif'), ...
                                      };
        paramsStruct.columntypes = {}; %{'label','logical',{'True','False',''},{'Yes','No'}};
        paramsStruct.columneditable = {true, true, true, true, true};
        paramsStruct.dockinggroup = 'Figures';
        paramsStruct.extra = {};
        paramsStruct.groupable = true;
        paramsStruct.interactivegrouping = false;

        paramsStruct.groups=[];
        
        for i=1:length(varargin)-1
            if(strcmp(varargin{i},'groups'))
                paramsStruct.groups=varargin{i+1};
            end
        end
        
        % Parse the regular (non-named) params in recption order
        if length(regParams)>0,  paramsStruct.container = regParams{1};  end  %#ok
        if length(regParams)>1,  paramsStruct.headers   = regParams{2};  end
        if length(regParams)>2,  paramsStruct.data      = regParams{3};  end

        % Parse the optional param PV pairs
        supportedArgs = fieldnames(paramsStruct);  % ={'container', 'headers', 'data', 'iconfilenames', 'columntypes', 'columneditable', 'dockinggroup'};
        while ~isempty(pvPairs)

            % Ensure basic format is valid
            paramName = '';
            if ~ischar(pvPairs{1})
                error('YMA:treeTable:invalidProperty','Invalid property passed to treeTable');
            elseif length(pvPairs) == 1
                error('YMA:treeTable:noPropertyValue',['No value specified for property ''' pvPairs{1} '''']);
            end

            % Process parameter values
            paramName  = pvPairs{1};
            paramValue = pvPairs{2};
            pvPairs(1:2) = [];
            if any(strncmpi(paramName,supportedArgs,length(paramName)))
                paramsStruct.(lower(paramName)) = paramValue;
            else
                paramsStruct.extra = {paramsStruct.extra{:} paramName paramValue};
            end
        end  % loop pvPairs

        % Create a panel spanning entire figure area, if container handle was not supplied
        if isempty(paramsStruct.container) || (~ishandle(paramsStruct.container) && ~isa(paramsStruct.container,'uiextras.Panel'))
            paramsStruct.container = uipanel('parent',gcf,'tag','TablePanel');
        end

        % Set default header names, if not supplied
        if isempty(paramsStruct.headers)
            if isempty(paramsStruct.data)
                paramsStruct.headers = {' '};
            else
                paramsStruct.headers = cellstr(char('A'-1+(1:size(paramsStruct.data,2))'))';
            end
        elseif ischar(paramsStruct.headers)
            paramsStruct.headers = {paramsStruct.headers};
        end

        % Convert data to cell-format (if not so already)
        paramsStruct.data           = cellizeData(paramsStruct.data);
        paramsStruct.headers        = cellizeData(paramsStruct.headers);
        paramsStruct.columntypes    = cellizeData(paramsStruct.columntypes);
        paramsStruct.columneditable = cellizeData(paramsStruct.columneditable);

        % Start with dummy data, just so that uitable can be initialized (or use supplied data, if available)
        if isempty(paramsStruct.data)
            selector = {'One','Two','Many'};
            paramsStruct.columntypes = {'label','label','char','logical',selector,'double'};
            paramsStruct.headers = {'ID','Label','Logical1','Logical2','Selector','Numeric'};  % 5 columns by default
            paramsStruct.data = {1,'M11',true, false,'One', 1011;  ...
                                 1,'M12',true, true, 'Two', 12;   ...
                                 1,'M13',false,false,'Many',13.4; ...
                                 2,'M21',true, false,'One', 21;  ...
                                 2,'M22',true, true, 'Two', 22;   ...
                                 3,'M31',true, true, 'Many',31;   ...
                                 3,'M32',false,true, 'One', -32;  ...
                                 3,'M33',false,false,'Two', 33; ...
                                 3,'M34',true, true, 'Many',34;  ...
                                 };
        end

        % Ensure we have valid column types & editable flags for all columns
        numCols = size(paramsStruct.data,2);
        [paramsStruct.headers{end+1:numCols}]        = deal(' ');
        [paramsStruct.columntypes{end+1:numCols}]    = deal('char');
        [paramsStruct.columneditable{end+1:numCols}] = deal(true);

        % TBD - Ensure icon filenames are readable and in the correct format

        a=1;
    catch
        if ~isempty(paramName),  paramName = [' ''' paramName ''''];  end
        err = lasterror;
        error('YMA:treeTable:invalidProperty',['Error setting treeTable property' paramName ':' char(10) lasterr]);
    end
end  % processArgs

%% Convert a numeric matrix to a cell array (if not so already)
function data = cellizeData(data)
    if ~iscell(data)
        %data = mat2cell(data,ones(1,size(data,1)),ones(1,size(data,2)));
        data = num2cell(data);
    end
end  % cellizeData

%% Process optional arguments on the newly-created table object
function processParams(paramsStruct,container,jtable)
    try
        % Process regular extra parameters
        paramName = '';
        th = jtable.getTableHeader;
        %container = get(mtable,'uicontainer');
        drawnow; pause(0.05);
        for argIdx = 1 : 2 : length(paramsStruct.extra)
            if argIdx<2
                % We need this pause to let java complete all table rendering
                % TODO: We should really use calls to awtinvoke() instead, though...
                pause(0.05);
            end
            if (length(paramsStruct.extra) > argIdx)   % ensure the arg value is there...
                paramsStruct.extra{argIdx}(1) = upper(paramsStruct.extra{argIdx}(1));  % property names always start with capital letters...
                paramName  = paramsStruct.extra{argIdx};
                paramValue = paramsStruct.extra{argIdx+1};
                propMethodName = ['set' paramName];
                
                % First try to modify the container
                try
                    set(container, paramName, paramValue);
                catch
                    try % if ismethod(mtable,propMethodName)
                        % No good, so try the mtable...
                        set(mtable, paramName, paramValue);
                    catch %elseif ismethod(jtable,propMethodName)
                        try
                            % sometimes set(t,x,y) failes but t.setX(y) is ok...
                            javaMethod(propMethodName, mtable, paramValue);
                        catch
                            try
                                % Try to modify the underlying JTable itself
                                set(jtable, paramName, paramValue);
                            catch
                                try
                                    javaMethod(propMethodName, jtable, paramValue);
                                catch
                                    try
                                        % Try to modify the table header...
                                        set(th, paramName, paramValue);
                                    catch
                                        javaMethod(propMethodName, th, paramValue);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end  % for argIdx
        drawnow; pause(0.05);

        % Process requested image columns
        try
            if ~isempty(which('ImageCellRenderer')) && ~isempty(paramsStruct.iconfilenames)
                %TODO TBD
                cr = ImageCellRenderer(paramsStruct.imagetooltipheight, paramsStruct.imagerescandelay*1000);
                if ischar(paramsStruct.iconfilenames)
                    % Maybe a header name?
                    jtable.getColumn(paramsStruct.imagecolumns).setCellRenderer(cr);
                elseif iscellstr(paramsStruct.imagecolumns)
                    % Cell array of header names
                    for argIdx = 1 : length(paramsStruct.imagecolumns)
                        jtable.getColumn(paramsStruct.imagecolumns{argIdx}).setCellRenderer(cr);
                        drawnow;
                    end
                else
                    % Try to treat as a numeric index array
                    for argIdx = 1 : length(paramsStruct.imagecolumns)
                        colIdx = paramsStruct.imagecolumns(argIdx) - 1;  % assume 1-based indexing
                        %jtable.setEditable(colIdx,0);  % images are editable!!!
                        jtable.getColumnModel.getColumn(colIdx).setCellRenderer(cr);
                        drawnow;
                    end
                end
                drawnow;
            elseif ~isempty(paramsStruct.imagecolumns)  % i.e., missing Renderer
                warning('YMA:treeTable:missingJavaClass','Cannot set image columns: ImageCellRenderer.class is missing from the Java class path');
            end
            jtable.repaint;
        catch
        end

        % Process UIContextMenu
        cm = get(container,'uicontextmenu');
        if ~isempty(cm)
            popupMenu = jtable.getRowHeaderPopupMenu;
            %popupMenu.list;
            popupMenu.removeAll; drawnow; pause(0.1);
            cmChildren = get(cm,'child');
            itemNum = 0;
            for cmChildIdx = length(cmChildren) : -1 : 1
                %{
                if itemNum == 6
                    % add 2 hidden separators which will be removed by the Matlab mouse listener...
                    popupMenu.addSeparator;
                    popupMenu.addSeparator;
                    popupMenu.getComponent(5).setVisible(0);
                    popupMenu.getComponent(6).setVisible(0);
                    itemNum = 8;
                end
                % Add a possible separator
                if strcmpi(get(cmChildren(cmChildIdx),'Separator'),'on')
                    popupMenu.addSeparator;
                    itemNum = itemNum + 1;
                end
                if itemNum == 6
                    % add 2 hidden separators which will be removed by the Matlab mouse listener...
                    popupMenu.addSeparator;
                    popupMenu.addSeparator;
                    popupMenu.getComponent(5).setVisible(0);
                    popupMenu.getComponent(6).setVisible(0);
                    itemNum = 8;
                end
                %}
                % Ramiro's fix:
                % "Though your code supports it, right now it has a little bug that if the user
                %  has more than 6 entries in the context menu, only 1 separator is shown at a
                %  fixed position. I was able to go around this problem using this code:"
                if itemNum==1 || itemNum==2 || itemNum==8 || itemNum==9
                    popupMenu.addSeparator;
                end
                % End Ramiro's fix

                % Add the main menu item
                jMenuItem = javax.swing.JMenuItem(get(cmChildren(cmChildIdx),'Label'));
                set(jMenuItem,'ActionPerformedCallback',get(cmChildren(cmChildIdx),'Callback'));
                popupMenu.add(jMenuItem);
                itemNum = itemNum + 1;
            end
            for extraIdx = itemNum+1 : 7
                popupMenu.addSeparator;
                popupMenu.getComponent(extraIdx-1).setVisible(0);
            end
            drawnow;
        end

        % Process docking group
        drawnow; pause(0.05);
        group = paramsStruct.dockinggroup;
        if ~strcmp(group,'Figures')
            try
                jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
                currentGroupNames = cell(jDesktop.getGroupTitles);
                if ~any(strcmp(group,currentGroupNames))
                    jDesktop.addGroup(group);
                end

                % Temporarily dock first figure into the group, to ensure container creation
                % Note side effect: group container becomes visible
                hFig = ancestor(paramsStruct.container,'figure');
                try
                    jFrame = get(mtable,'ParentFigureValidator');
                catch
                    % old Matlab versions used a different property name
                    jFrame = get(mtable,'UserParentFigure');
                end
                set(jFrame,'GroupName',group);
                %oldStyle = get(hFig,'WindowStyle');
                %set(hFig,'WindowStyle','docked');  drawnow
                %set(hFig,'WindowStyle',oldStyle);  drawnow
                %commandwindow;
                drawnow; pause(0.02);
                try
                    jDesktop.setGroupDocked(group,false);
                    jDesktop.showGroup(group,true);
                catch
                    % never mind...
                end
                figure(hFig);
            catch
                warning('YMA:treeTable:Docking',['Cannot dock figure: ' lasterr]);
            end
        end
    catch
        if ~isempty(paramName),  paramName = [' ''' paramName ''''];  end
        err = lasterror;
        error('YMA:treeTable:invalidProperty',['Error setting treeTable property' paramName ' (line #' num2str(err.stack(1).line) '):' char(10) lasterr]);
    end
end  % processParams

%% Get the default cell renderer object
function cr = getDefaultCellRenderer()
    try
        % Custom cell renderer (striping, cell FG/BG color, cell tooltip)
        cr = CustomizableCellRenderer;
        cr.setRowStriping(false);
    catch
        % Use the standard JTable cell renderer
        %cr = [];
        cr = javax.swing.table.DefaultTableCellRenderer;
    end
end  % getDefaultCellRenderer

%% Sample row insertion function
function addRow(jtable,newData)
    [data,headers] = getTableData(jtable);
    data = [data; newData];
    setTableData(jtable,data,headers);
end  % addRow

%% Sample row deletion function
function deleteRow(jtable,rowIdx)
    jtable.getModel.removeRow(rowIdx);
    jtable.repaint;
end  % addRow

%% Sample table-editing callback
function tableChangedCallback(hModel,hEvent,jtable)
    % Get the modification data
    modifiedRow = get(hEvent,'FirstRow');
    modifiedCol = get(hEvent,'Column');
    label   = hModel.getValueAt(modifiedRow,1);
    newData = hModel.getValueAt(modifiedRow,modifiedCol);

    % Now do something useful with this info
    fprintf('Modified cell %d,%d (%s) to: %s\n', modifiedRow+1, modifiedCol+1, char(label), num2str(newData));
end  % tableChangedCallback

%% Sample table-selection callback
function selectionCallback(hSelectionModel,hEvent,jtable)
    % Get the selection data
    MinSelectionIndex  = get(hSelectionModel,'MinSelectionIndex');
    MaxSelectionIndex  = get(hSelectionModel,'MaxSelectionIndex');
    LeadSelectionIndex = get(hSelectionModel,'LeadSelectionIndex');

    % Now do something useful with this info
    fprintf('Selected rows #%d-%d\n', MinSelectionIndex+1, MaxSelectionIndex+1);
end  % selectionCallback
