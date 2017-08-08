classdef (CaseInsensitiveProperties=true,Abstract) samplingflagwin < sigwin.variablelength & hgsetget
  %sigwin.samplingflagwin class
  %   sigwin.samplingflagwin extends sigwin.variablelength.
  %
  %    sigwin.samplingflagwin properties:
  %       Name - Property is of type 'string' (read only)
  %       Length - Property is of type 'spt_uint32 user-defined'
  %       SamplingFlag - Property is of type 'signalSignalwindowsWindowSampling_flag enumeration: {'symmetric','periodic'}'
  %
  %    sigwin.samplingflagwin methods:
  %       thisinfo - Information for this class.
  %       thisloadobj -   Load this object.
  
  
  properties (AbortSet, SetObservable, GetObservable)
    %SAMPLINGFLAG Property is of type 'signalSignalwindowsWindowSampling_flag enumeration: {'symmetric','periodic'}'
    SamplingFlag = 'symmetric';
  end
  
  properties (Transient, SetObservable, GetObservable, Hidden)
    %SAMPLING_FLAG Property is of type 'signalSignalwindowsWindowSampling_flag enumeration: {'symmetric','periodic'}' (hidden)
    sampling_flag = 'symmetric';
  end
  
  
  methods
    function set.SamplingFlag(obj,value)
      % Enumerated DataType = 'signalSignalwindowsWindowSampling_flag enumeration: {'symmetric','periodic'}'
      value = validatestring(value,{'symmetric','periodic'},'','SamplingFlag');
      obj.SamplingFlag = value;
    end
    
    function value = get.sampling_flag(obj)
      value = getsampling_flag(obj,obj.sampling_flag);
    end
    
    function set.sampling_flag(obj,value)
      % Enumerated DataType = 'signalSignalwindowsWindowSampling_flag enumeration: {'symmetric','periodic'}'
      value = validatestring(value,{'symmetric','periodic'},'','sampling_flag');
      obj.sampling_flag = setsampling_flag(obj,value);
    end
    
    function [p, v] = thisinfo(h)
      %THISINFO Information for this class.
      
      [pvl, vvl] = varlen_thisinfo(h);
      p = {pvl{:}, getString(message('signal:dfilt:info:SamplingFlag'))};
      v = {vvl{:},get(h, 'SamplingFlag')};
    end
    
    function thisloadobj(this, s)
      %THISLOADOBJ   Load this object.
      
      set(this,'SamplingFlag',s.SamplingFlag);
    end
    
    function varargout = set(obj,varargin)
      
      varargout = signal.internal.signalset(obj,varargin{:});
      varargout = {varargout};
      
    end
    
    function values = getAllowedStringValues(~,prop)
      % This function gives the the valid string values for object properties.
      
      switch prop
        case 'SamplingFlag'
          values = {...
            'symmetric'
            'periodic'};
          
        case 'sampling_flag'
          values = {...
            'none'
            'symmetric'
            'periodic'};
          
        otherwise
          values = {};
      end
      
    end
    
  end  %% public methods
  
end  % classdef

function sf = setsampling_flag(this, sf)

warning(message('signal:sigwin:samplingflagwin:schema:DeprecatedProperty'));

set(this, 'SamplingFlag', sf);
end  % setsampling_flag


% -------------------------------------------------------------------------
function sf = getsampling_flag(this, sf) %#ok

warning(message('signal:sigwin:samplingflagwin:schema:DeprecatedProperty'));

sf = get(this, 'SamplingFlag');
end  % getsampling_flag


% [EOF]
