classdef (CaseInsensitiveProperties=true,Abstract) parameterizewin < dynamicprops & sigwin.variablelength 
  %sigwin.parameterizewin class
  %   sigwin.parameterizewin extends sigwin.variablelength.
  %
  %    sigwin.parameterizewin properties:
  %       Name - Property is of type 'string' (read only)
  %       Length - Property is of type 'spt_uint32 user-defined'
  %
  %    sigwin.parameterizewin methods:
  %       createdynamicprops - CREATEDYNAMICSPROPS Create dynamic properties
  %       getparamnames - Get the name of the dynamic properties
  %       thisinfo - Information for this class (without the name)
  %       thisloadobj -   Load this object.
  
  
  properties (Access=protected, AbortSet, SetObservable, GetObservable)
    %DYNAMICPROP Property is of type 'handle vector'
    DynamicProp = [];
  end
  
  
  methods
    function set.DynamicProp(obj,value)
      % DataType = 'handle vector'
      validateattributes(value,{'handle'}, {'vector'},'','DynamicProp')
      obj.DynamicProp = value;
    end
    
    function createdynamicprops(hWIN, propName, propType, propDes)
      %CREATEDYNAMICSPROPS Create dynamic properties
      
      if iscell(propName),
        for i=1:length(propName)
          % p(i) = schema.prop(hWIN, propName{i}, propType{i}); %#ok
          % set(p(i),'Description',propDes);
          name = propName{idx};
          p(i) = hWIN.addprop(name);
          p(i).NonCopyable = false;
        end
      else
        % p = schema.prop(hWIN, propName, propType);
        % set(p,'Description',propDes);
        p = hWIN.addprop(propName);
        p.Description = propDes;
        p.NonCopyable = false;
      end
      
      if isempty(hWIN.DynamicProp)
        hWIN.DynamicProp = p;
      else
        hWIN.DynamicProp = [hWIN.DynamicProp p];
      end
    end
    
    function [ParamNames, des] = getparamnames(hWIN)
      %GETPARAMNAMES Get the name of the dynamic properties
      
      p = hWIN.DynamicProp;
      np = length(p); % number of dynamic properties.
      if np==1 % to be compatible with previous outputs.
        ParamNames = p.Name;
        des = p.Description;
      else
        ParamNames = cell(np, 1);
        des = cell(np, 1);
        for pI = 1:np
          ParamNames{pI} = p(pI).Name;
          des{pI} = p(pI).Description;
        end
      end
      
    end
    
    function [p, v] = thisinfo(h)
      %THISINFO Information for this class (without the name)
      
      [pvl, vvl] = varlen_thisinfo(h);
      [param, des] = getparamnames(h);
      if ~iscell(param)
        param = {param};
      end
      if ~iscell(des)
        des = {des};
      end
      p = {pvl{:}, des{:}};
      ndp = length(h.DynamicProp);    % number of dynamic properties
      v = cell(1, ndp+1);
      v(1) = vvl;
      for dpI = 1:ndp
        v(dpI+1) = {sprintf('%g', get(h, param{dpI}))};
      end
      
    end
    
    function thisloadobj(this, s)
      %THISLOADOBJ   Load this object.
      
      ParamNames = getparamnames(this);
      
      if iscell(ParamNames)
        for paramI = 1:length(ParamNames)
          oneName = ParamNames{paramI};
          set(this, oneName, s.(oneName));
        end
      else
        set(this,ParamNames,s.(ParamNames));
      end
      
    end
    
  end  %% public methods
  
end  % classdef

