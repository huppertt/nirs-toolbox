function generatehdl(filterobj, varargin)
%GENERATEHDL Generate HDL.

%   Copyright 2003-2012 The MathWorks, Inc.

  % check for Filter Design HDL Coder
  fdhdlcInstallCheck;

  [cando, ~, errObj] = ishdlable(filterobj);
  if ~cando
    error(errObj);
  end

  % Add filter object variable name to varargin 
  if ~any(strcmpi(varargin,'name'))
    varargin(end+1) = {'name'};
    if ~isempty(inputname(1))
        dutname = inputname(1);
        varargin(end+1) = {dutname};
    else
        error(message('signal:dfilt:basefilter:generatehdl:genhdlcalledwithconst'));
    end
  end
    
  % Add testbench name to varargin when testbench is requested but name not
  % provided -- make up a name based on filter variable name
  indices = strcmpi(varargin, 'generatehdltestbench');
  pos = 1:length(indices);
  pos = pos(indices);
  
  indices_name = strcmpi(varargin, 'name');
  posname = 1:length(indices_name);
  posname = posname(indices_name);
  
  if (~isempty(pos) && ~strcmpi(varargin{pos+1},'off')) && ... % tb is requested
          ~any(strcmpi(varargin,'testbenchname')) % but name not provided
    varargin(end+1) = {'testbenchname'};
    if ~isempty(pos) % if dut name is provided then derive the tb name
        varargin(end+1) = {[varargin{posname+1}, '_tb']}; 
    else
        varargin(end+1) = {[inputname(1) '_tb']}; 
    end
  end
  
  if any(strcmpi(fieldnames(get(filterobj)),'Arithmetic')) %all but cascades
      if (strcmpi(get(filterobj, 'Arithmetic'), 'fixed') && filterobj.InputWordLength == 1)
          error(message('signal:dfilt:basefilter:generatehdl:OneBitInputNotSupported'));
      end
  elseif strcmpi(get(filterobj.Stage(1), 'Arithmetic'), 'fixed') && filterobj.Stage(1).InputWordLength == 1
      error(message('signal:dfilt:basefilter:generatehdl:OneBitInputNotSupported'));
  end
  
  privgeneratehdl(filterobj,varargin{:});

% [EOF]

