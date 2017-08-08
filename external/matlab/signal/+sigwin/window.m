classdef (CaseInsensitiveProperties=true,Abstract) window < hgsetget & matlab.mixin.Copyable & matlab.mixin.SetGet
  %sigwin.window class
  %    sigwin.window properties:
  %       Name - Property is of type 'string' (read only)
  %
  %    sigwin.window methods:
  %       coefficientvariables - Window coefficient variables.
  %       disp - Display a window object
  %       dispstr - Returns the strings to display the window
  %       exportdata - Extract data to export.
  %       exportinfo - Export information for SIGWIN objects.
  %       getinfoheader - Return the header to the info.
  %       getparamnames - Default method of the base class
  %       info - Information about a window
  %       isminordersupported - Default method of the base class
  %       loadobj -   Load this object.
  %       propstoadd - Method that lists the properties associated with this object.
  %       propstoaddtospectrum - Return properties to be added to the SPECTRUM object.
  %       saveobj -   Save this object.
  %       sigwinname - SIGWIN object name.
  %       sigwinvariable - SIGWIN object variable name.
  %       windowname -  Window name.
  %       winwrite - Write a window file.
  
  
  properties (Access=protected, AbortSet, SetObservable, GetObservable)
    %VERSION Property is of type 'double'
    Version = 1.0;
  end
  
  properties (SetAccess=protected, AbortSet, SetObservable, GetObservable)
    %NAME Property is of type 'string' (read only)
    Name = '';
  end
  
  
  methods
    function set.Version(obj,value)
      % DataType = 'double'
      validateattributes(value,{'double'}, {'scalar'},'','Version')
      obj.Version = value;
    end
    
    function set.Name(obj,value)
      % DataType = 'string'
      validateattributes(value,{'char'}, {'row'},'','Name')
      obj.Name = set_name(obj,value);
    end
    
  end   % set and get functions
  
  methods  %% public methods
    function c = coefficientvariables(this) %#ok
      %COEFFICIENTVARIABLES Window coefficient variables.
      
      %   This should be a private method.
      
      c = {'Win'};
    end
    
    function disp(hWIN)
      %DISP Display a window object
      
      
      disp(reorderstructure(get(hWIN), 'Name'));
      
    end
    
    function strs = dispstr(hWin)
      %DISPSTR Returns the strings to display the window
      
      
      strs = sprintf('%g\n', generate(hWin));
    end
    
    function data = exportdata(Hwin)
      %EXPORTDATA Extract data to export.
      
      % This should be a private method.
      
      
      data = {generate(Hwin)};
      
    end
    
    function s = exportinfo(Hwin)
      %EXPORTINFO Export information for SIGWIN objects.
      
      %   This should be a private method.
      
      
      % Both coefficientnames & coefficientvariables return cell arrays.
      s.variablelabel = sigwinname(Hwin);
      s.variablename =  coefficientvariables(Hwin);
      
      % SIGWINs can be exported as both objects and arrays. Filters can also be
      % exported as System objects so we need to add this option to the popup
      % strings.
      s.exportas.tags = {'Coefficients','Objects','System Objects'};
      
      % SIGWIN object specific labels and names
      s.exportas.objectvariablelabel = sigwinname(Hwin);
      s.exportas.objectvariablename  = sigwinvariable(Hwin);
      
      % Optional fields (destinations & constructors) if exporting to destinations other
      % than the built-in 'Workspace','Text-file', or, 'MAT-file';
      s.destinations  = {'Workspace','Text-File','MAT-File'};
      s.constructors  = {'','sigio.xp2winfile',''};
    end
    
    function str = getinfoheader(h)
      %GETINFOHEADER Return the header to the info.
      
      % This should be a private method.
      
      
      str = [h.Name, ' Window'];
    end
    
    function ParamNames = getparamnames(this) %#ok
      %GETPARAMNAMES Default method of the base class
      
      %   Author(s): V.Pellissier
      %   Copyright 1988-2004 The MathWorks, Inc.
      
      ParamNames = '';
      
    end
    
    function s = info(h)
      %INFO Information about a window
      %   S = INFO(Hwin) returns a string matrix with information about the window.
      %
      %   See also SIGWIN.
      
      [p, v] = thisinfo(h);
      titlestr = getinfoheader(h);
      infostrs = { ...
        getinfoheader(h), ...
        repmat('-', 1, size(titlestr, 2)), ...
        [strvcat(p{:}) repmat('  : ', length(p), 1), strvcat(v{:})], ...
        }; %#ok
      s = strvcat(infostrs{:}); %#ok
      
    end
    
    function bool = isminordersupported(this) %#ok
      %ISMINORDERSUPPORTED Default method of the base class
      
      bool = 0;
      
    end
    
    function props = propstoadd(this)
      %PROPSTOADD Method that lists the properties associated with this object.
      %
      % This method can be used in conjunction with the utility function
      % ADDPROPS, which dynamically adds properties from one object to another.
      % This is useful in the case of one object containing another object.
      
      props = fieldnames(this);
      
    end
    
    function props2add = propstoaddtospectrum(this)
      %PROPSTOADDTOSPECTRUM Return properties to be added to the SPECTRUM object.
      
      %   Author(s): P. Pacheco
      %   Copyright 1988-2003 The MathWorks, Inc.
      
      allProps = propstoadd(this);
      
      % Exclude the properties listed in the cell array.
      props2add = setdiff(allProps,{'Name','Length'});
      
    end
    
    function s = saveobj(this)
      %SAVEOBJ   Save this object.
      
      % Get all the public properties.
      s = rmfield(get(this), 'Name');
      
      % Get the class.
      s.class = class(this);
      
    end
    
    function c = sigwinname(this) %#ok
      %SIGWINNAME SIGWIN object name.
      
      % This should be a private method.
      
      c = {'Window'};
      
    end
    
    
    function c = sigwinvariable(this) %#ok
      %SIGWINVARIABLE SIGWIN object variable name.
      
      % This should be a private method.
      
      
      c = {'Hwin'};
      
    end
    
    function c = windowname(this) %#ok
      %WINDOWNAME  Window name.
      
      c = {'Window'};
      
    end
    
    function winwrite(h,filename,hCD)
      %WINWRITE Write a window file.
      %   WINWRITE(Hwin) writes an ASCII-file with window weights.  The window
      %   values are extracted from the SIGWIN window object, Hwin. Hwin maybe a
      %   single SIGWIN object or a vector of SIGWINs. A dialog box is displayed
      %   to fill in a file name. The default file name is 'untitled.wf'.
      %
      %   WINWRITE(Hwin,FILENAME) writes the file to a disk file called
      %   FILENAME in the present working directory.
      %
      %   The extension '.wf' will be added to FILENAME if it doesn't already
      %   have an extension.
      %
      %   See also INFO, SIGWIN.
      
      narginchk(1,3);
      if nargin < 2 
        filename = []; 
      end
      
      hetArr = false;
      if nargin == 3
        hetArr = true;
      end
      
      if isempty(filename),
        filename = 'untitled.wf';
        dlgStr = 'Export Window to .WF file';
        [filename,pathname] = uiputfile('*.wf', dlgStr, filename);
      else
        % File will be created in present directory
        s = pwd;
        pathname = [s filesep];
      end
      
      if ~isempty(filename),
        if isempty(findstr(filename,'.')), filename=[filename '.wf']; end
        filename = [pathname filename];
      end
      
      if ~any(filename == 0),
        if hetArr
          save2txtfile(hCD,filename,hetArr);
        else
          save2txtfile(h,filename,hetArr);
        end
      end
      
    end
    
    
    
  end  %% public methods
  
  methods (Hidden) %% possibly private or hidden
    function hout = copyobj(hWIN)
      %COPYOBJ
      
      hout = feval(str2func(class(hWIN)));
      s = rmfield(get(hWIN), 'Name');
      set(hout, s);
      
    end
    
    function thisloadobj(this, s) %#ok
      %THISLOADOBJ
      
      % NO OP, meant to be overloaded.
    end
    
    function c = windowvariable(this) %#ok
      %WINDOWVARIABLE
      %
      %   Output:
      %       C  - cell array
      
      %   This should be a private method.
      
      c = {'Win'};
      
    end
    
  end  %% possibly private or hidden
  
  
  methods (Static) %% static methods
    function this = loadobj(s)
      %LOADOBJ   Load this object.
      
      this = feval(s.class);
      
      % All windows have a length
      set(this,'Length', s.Length);
      
      thisloadobj(this, s);
      
    end
  end  %% static methods
  
end  % classdef

function str = set_name(~,str)

if ~license('checkout','signal_toolbox')
  error(message('signal:sigwin:window:schema:LicenseRequired'));
end
end  % set_name

%--------------------------------------------------------------
function save2txtfile(h,file,hetArr)

% Write the coefficients out to a file.
fid = fopen(file,'w');

% Display header information
fprintf(fid,'%s\n',sptfileheader);

txt = getfiletxt(h,hetArr);

sz = size(txt);
for j = 1:sz(1), % Rows
  fprintf(fid, '%s\n', num2str(txt(j,:),10));
end
fprintf(fid, '\n');

fclose(fid);

% Launch the MATLAB editor (to display the generated file)
edit(file);

end

% -------------------------------------------------------------------------
function txt = getfiletxt(Hb,hetArr)
% txt is a character array

if ~hetArr && ~iscell(Hb)
  Hb = {Hb};
end

strs = cell(length(Hb)*4, 1);
for idx = 1:length(Hb)  
  strs{idx*4-3} = info(Hb{idx});
  strs{idx*4-2} = sprintf('\n');
  strs{idx*4-1} = dispstr(Hb{idx});
  strs{idx*4}   = sprintf('\n');
end
txt = strvcat(strs{:}); %#ok

end
% [EOF]








