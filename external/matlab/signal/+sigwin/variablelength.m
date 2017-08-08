classdef (CaseInsensitiveProperties=true, Abstract) variablelength < sigwin.window 
  %sigwin.variablelength class
  %   sigwin.variablelength extends sigwin.window.
  %
  %    sigwin.variablelength properties:
  %       Name - Property is of type 'string' (read only)
  %       Length - Property is of type 'spt_uint32 user-defined'
  %
  %    sigwin.variablelength methods:
  %       thisinfo - Information for this class.
  %       varlen_thisinfo - Information for variablelength class.
  
  
  properties (AbortSet, SetObservable, GetObservable)
    %LENGTH Property is of type 'spt_uint32 user-defined'
    Length = 64;
  end
  
  
  methods
    function set.Length(obj,value)
      % User-defined DataType = 'spt_uint32 user-defined'
      obj.Length = value;
    end
    
    function [p, v] = thisinfo(h)
      %THISINFO Information for this class.
      
      % This should be a private method.
      
      [p, v] = varlen_thisinfo(h);
      
    end
    
    function [p, v] = varlen_thisinfo(h)
      %VARLEN_THISINFO Information for variablelength class.
      
      % This should be a private method.
      
      p = {getString(message('signal:dfilt:info:Length'))};
      v = {sprintf('%g', h.Length)};
    end
    
  end  %% public methods
  
end  % classdef

