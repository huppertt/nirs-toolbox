classdef (Sealed) SPTCustomSettings < handle
%SPTCustomSettings define SPT settings
% This is a Singleton persistent class that, after being instantiated, will
% exist throughout the MATLAB session.
  
  % Copyright 2013 The MathWorks, Inc.
  properties
    % This property indicates whether MATLAB is running on MOTW. This flag
    % is set to false at startup of a MOTW session.
    DDGSupport = true;
  end
  methods (Access = private)
    function obj = SPTCustomSettings
      % Private constructor, this class can only be instantiated by calling
      % one of its methods.
    end
  end
  methods (Static)
    function singleObj = getInstance
      mlock
      persistent localObj
      if isempty(localObj) || ~isvalid(localObj)
        localObj = signal.internal.SPTCustomSettings;
      end
      singleObj = localObj;
    end
    
    % If an instant of this class does not exist, it will be
    % automatically created when calling one of these methods.
    function setDDGSupportFlag(flag)      
      obj = signal.internal.SPTCustomSettings.getInstance;
      obj.DDGSupport = flag;
    end
    function flag = isDDGSupported
      obj = signal.internal.SPTCustomSettings.getInstance;
      flag = obj.DDGSupport;
    end
  end
end

 