classdef nnfcnFunctionType < nnfcnType
%NNFUNCTIONTYPE Modular function type info.

% Copyright 2010 The MathWorks, Inc.

  properties
    supportedVersions = [];
    toolboxFolder = '';
  end
  
  methods
    
    function x = nnfcnFunctionType(mname,name,version,versions,folder)
      x = x@nnfcnType(mname,name,version);
      x.supportedVersions = versions;
      x.toolboxFolder = folder;
    end
    
    function disp(x)
      sv = '[';
      for i=1:length(x.supportedVersions)
        if (i>1), sv = [sv ', ']; end
        sv = [sv num2str(x.supportedVersions(i))];
      end
      sv = [sv ']'];
      
      isLoose = strcmp(get(0,'FormatSpacing'),'loose');
      disp@nnfcnType(x)
      if (isLoose), fprintf('\n'), end
      disp([nnlink.prop2link('supportedVersions') sv]);
      disp([nnlink.prop2link('toolboxFolder') 'nnet' filesep x.toolboxFolder]);
      if (isLoose), fprintf('\n'), end
    end
    
  end
  
end
