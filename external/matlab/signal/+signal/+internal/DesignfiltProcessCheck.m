classdef DesignfiltProcessCheck < handle
%DesignfiltProcessCheck Parse and design filters for designfilt function

  %   Copyright 2013 The MathWorks, Inc.

%#ok<*AGROW>
  
properties (Hidden)
  % Filter Object Database variables
  InputResponseTypes
  ValidDesignOpts
  SymbolDefs
  FilterResponse = '';
  FilterType = '';
  PropertyNames = {};
  PropertyValues = {};  
end

methods (Hidden)
  function obj = DesignfiltProcessCheck
    % Constructor
    % Initialize the Filter Database
    initFiltObj(obj);
  end
end

methods (Hidden)
  function [err,filterFcn,parseParams] = checkConstraints(obj,filterFcn,userProps,userVals,userValsNames)
    % Checks the Constraints of the specified Filter
    
    parseParams.ResponseOnly = false;
    
    % Check if the response specified by user is Valid
    [err,errid,filterFcn] = verifyInputResponse(obj,filterFcn);
    
    if ~isempty(err)
      parseParams.ErrorId = errid;
      return
    end
    
    if isempty(err) && isempty(userProps) && isempty(userVals)
      % No Input PV-Pairs Specified.
      parseParams.ResponseOnly = true;
      return
    end
    
    % Look for valid filter specifications, design options, design methods,
    % and sample rate.
    
    % Get the entire list of properties and design options
    % available for the designfilt function regardless of the input
    % constraints. Make sure all input names are un-ambiguous.
    designOptsList = obj.ValidDesignOpts(:,1);    
    propList = obj.SymbolDefs(:,2);
    
    fullList = [designOptsList; propList; {'SampleRate'}; {'DesignMethod'}];        
    parseParams.FullInputNameList = fullList;
    
    % validateStringList outputs a cell array with correct set of strings
    % (it completes incomplete strings, and capitalizes strings according
    % to valid list in fullList). 
    [userProps,~,err,errType,errid] = obj.validateStringList(fullList,userProps);
    if errType == 2
      parseParams.ResponseOnly = true;
      parseParams.ErrorId = errid;
      return
    else
      err = '';
    end
    
    % Put different inputs in different buckets
    
    % Design options
    [cachedDesOptInputs,cachedOptVal,cachedOptValNames,...
      remInputs,remVals,remValsNames] = ...
      cacheInputList(designOptsList,userProps,userVals,userValsNames);
    
    % Sample rate
    [cachedSampleRate,cachedSampleRateVal,...
      cachedSampleRateValName,remInputs,remVals,remValsNames] = ...
      cacheInputList({'SampleRate'},remInputs,remVals,remValsNames);
    
    % Design method
    [cachedDesignMethod, cachedDesignMethodVal,...
      cachedDesignMethodValName,remInputs,remVals,remValsNames] = ...
      cacheInputList({'DesignMethod'},remInputs,remVals,remValsNames);
    
    if ~isempty(cachedDesignMethod) && ~ischar(cachedDesignMethodVal{1})
      error(message('signal:designfilt:DesignMethodMustBeString'));
    end
    
    % Make sure design method names are valid and not ambiguous
    designMethodIsAmbiguousFlag = false;
    if ~isempty(cachedDesignMethod)
      validDesMethodList = obj.getDesignMethodsList;
      [cachedDesignMethodVal,~,err,errType,~] = obj.validateStringList(validDesMethodList,cachedDesignMethodVal);
      if errType == 2
        designMethodIsAmbiguousFlag = true;              
      end
    end
    
    % Filter specs
    [userProps, userVals,userValsNames,remInputs,~,~] = ...
      cacheInputList(propList,remInputs,remVals,remValsNames);
    
    unrecognizedInputs = remInputs;
        
    if isempty(userProps) && isempty(userVals)
      % No Input PV-Pairs Specified.
      parseParams.ResponseOnly = true;
      if isempty(unrecognizedInputs)
        parseParams.ErrorId = 'signal:designfilt:MustSpecifyValidSpecifications';
        parseParams.ErrorMessage = getString(message(parseParams.ErrorId));
      elseif numel(unrecognizedInputs) == 1
        parseParams.ErrorId = 'signal:designfilt:MustSpecifyValidSpecifications_ParameterNotRecognized';
        str1 = getString(message('signal:designfilt:MustSpecifyValidSpecifications'));
        str1 = [str1(1:end-1) ', '];
        str2 = getString(message('signal:designfilt:ParameterNotRecognized',unrecognizedInputs{:}));
        parseParams.ErrorMessage = [str1 str2];
      else
        parseParams.ErrorId = 'signal:designfilt:MustSpecifyValidSpecifications_ParametersNotRecognized';
        str1 = getString(message('signal:designfilt:MustSpecifyValidSpecifications'));
        str1 = [str1(1:end-1) ', '];
        str2 = getString(message('signal:designfilt:ParametersNotRecognized',cell2str(unrecognizedInputs,', ',true)));
        parseParams.ErrorMessage = [str1 str2];
      end
        err = parseParams.ErrorMessage;
      return
    end
    
    % Check valid filter specifications. A hard error occurs if a property
    % has an invalid value regardless of the fact that the user passed in a
    % valid set of specs or not. 
    [parseParams,~] = identifyConstraintSetFromSpecifiedProperties( ...
      obj,filterFcn,userProps,userVals,cachedDesignMethodVal);
    
    % Check valid design option names. Do it for all possible design
    % options, later on we will check that these design options are valid
    % for the specified set of specifications and design method. Do this
    % check after we have checked the filter spec values inside
    % identifyConstraintSetFromSpecifiedProperties. String properties are
    % replaced with corrected string value. 
    for idx = 1:numel(cachedDesOptInputs)
      cachedOptVal{idx} = validateDesingOptValue(obj,cachedDesOptInputs{idx},cachedOptVal{idx});
    end
        
    parseParams.ResponseOnly = false;
    parseParams.FullInputNameList = fullList;
    
    % Add all user inputs to the parseParams structure
    parseParams.PropertyNames = userProps;
    parseParams.PropertyValNames = userValsNames;
    parseParams.DesignOptsProps = cachedDesOptInputs;
    parseParams.DesignOptsVals = cachedOptVal;
    parseParams.DesignOptsValNames = cachedOptValNames;
    parseParams.DesignMethodProps = cachedDesignMethod;
    parseParams.DesignMethodVals = cachedDesignMethodVal;
    parseParams.DesignMethodValNames = cachedDesignMethodValName;
    parseParams.SampleRateProps = cachedSampleRate;
    parseParams.SampleRateVals = cachedSampleRateVal;
    parseParams.SampleRateValNames = cachedSampleRateValName;
    
    isfir = ~isempty(strfind(filterFcn,'fir'));
    
    if ~isempty(parseParams.FdesCode)
      % Choose a design method from input or the default - design method
      % validity is checked later if there are no design specification parsing
      % errors.
      if ~isempty(cachedDesignMethodVal) && ~designMethodIsAmbiguousFlag
        convertedName = obj.convertdesignmethodnames(cachedDesignMethodVal{:},isfir,'designfiltToFdesign');
        if any(strcmpi(parseParams.FdesMethods,convertedName))
          parseParams.FdesMethods = convertedName;
        else
          % Use default method if requested method is not valid
          parseParams.FdesMethods = getDefaultDesignMethod(parseParams.FdesMethods);          
        end
      else
        % Use default method
        parseParams.FdesMethods = getDefaultDesignMethod(parseParams.FdesMethods);
      end
      
      if parseParams.FinalizeFdesCode
        % Add sample rate if any then close parenthesis on fdesign code
        % string. If parseParams.FinalizeFdesCode is false, the code has
        % already been completed.
        if isempty(cachedSampleRate)
          parseParams.FdesCode = [parseParams.FdesCode ',''FromDesignfilt'')'];
        else
          parseParams.FdesCode = [parseParams.FdesCode ',' ...
            num2str(cachedSampleRateVal{:}) ',''FromDesignfilt'')'];
        end
      end
    end
    
    % If design specification parsing errors were found then add them to
    % the error output. Concatenate parsing errors with unrecognized input
    % errors.
    if ~isempty(parseParams.ErrorMessage)
      if ~isempty(unrecognizedInputs)
        if numel(unrecognizedInputs) > 1
          parseParams.ErrorMessage = [parseParams.ErrorMessage(1:end-1) ', ' ...
            getString(message('signal:designfilt:ParametersNotRecognized',cell2str(unrecognizedInputs,', ',true)))];          
          parseParams.ErrorId = [parseParams.ErrorId '_ParametersNotRecognized'];                    
        else
          parseParams.ErrorMessage = [parseParams.ErrorMessage(1:end-1) ', ' ...
            getString(message('signal:designfilt:ParameterNotRecognized',unrecognizedInputs{1}))];
          parseParams.ErrorId = [parseParams.ErrorId '_ParameterNotRecognized'];
        end
        err = parseParams.ErrorMessage;
        if ~isempty(parseParams.InputRecommendation)          
          parseParams.InputRecommendation(1) = lower(parseParams.InputRecommendation(1));
          if numel(unrecognizedInputs) > 1
            parseParams.InputRecommendation(1) = lower(parseParams.InputRecommendation(1));
            parseParams.InputRecommendation = [getString(message('signal:designfilt:ParametersNotRecognized1',cell2str(unrecognizedInputs,', ',true))) ' ' parseParams.InputRecommendation];
          else
            parseParams.InputRecommendation = [getString(message('signal:designfilt:ParameterNotRecognized1',unrecognizedInputs{1})) ' ' parseParams.InputRecommendation];
          end
        end
      else
        err = parseParams.ErrorMessage;
      end
      return;
    end
    
    % If we make it this far, then it means no specification errors were found.
    % Now we need to check design method and design option inputs.
    
    % Verify correct design method, add error if design method input is invalid
    cstr = parseParams.Constraint{:};
    designMethods = getDesignMethodsForConstraintStr(obj,filterFcn,cstr); % fdesign names
    designMethods = obj.convertdesignmethodnames(designMethods,isfir,'fdesignToDesignfilt'); % designfilt names
    
    if isempty(cachedDesignMethodVal)
      % Use default if not specified
      cachedDesignMethodVal = {getDefaultDesignMethod(designMethods)};
    elseif designMethodIsAmbiguousFlag
      % Ambiguous design method
      parseParams.ErrorId = 'signal:designfilt:AmbiguousDesignMethodName';
      parseParams.ErrorMessage = getString(message(parseParams.ErrorId,cachedDesignMethodVal{:}));         
      err = parseParams.ErrorMessage;      
    else
      if ~any(strcmpi(designMethods,cachedDesignMethodVal{:}))
        propStr = cell2str(userProps,', ');
        parseParams.ErrorId = 'signal:designfilt:InvalidDesignMethod';
        err = getString(message(parseParams.ErrorId,cachedDesignMethodVal{:},filterFcn,propStr));        
      end
    end
    
    % Verify correct design options based on valid design method. Do not check
    % design options if design method was invalid. If any design option is
    % valid for the design method that will be used to launch the design
    % assistant, its value will be added by the filterbuilder code.
    parseParams.DesignOptsStruct = [];
    if isempty(err)
      fdesDesOpts = struct;
      if ~isempty(cachedDesOptInputs)
        designOpts = getValidDesignOptions(obj,filterFcn,cstr,cachedDesignMethodVal{1});
        for idx = 1:length(cachedDesOptInputs)
          if any(strcmpi(designOpts,cachedDesOptInputs{idx}))
            % Get valid design options and create a strucutre that we can pass
            % to the fdesign code that will design the filter.
            
            % Make sure values are valid - if value is a string, then
            % correct string is placed in cachedOptVal output, otherwise,
            % output returns same value as cachedOptVal input.
            desMethName = obj.convertdesignmethodnames(parseParams.FdesMethods,isfir,'fdesignToDesignfilt'); 
            cachedOptVal{idx} = validateDesingOptValue(obj,cachedDesOptInputs{idx},cachedOptVal{idx},desMethName);
                        
            name = obj.convertdesignoptnames(cachedDesOptInputs{idx},'designfiltToFdesign');
            fdesDesOpts.(name) = cachedOptVal{idx};
          else
            parseParams.ErrorId = 'signal:designfilt:InvalidDesignOption';
            err = getString(message(parseParams.ErrorId,cachedDesOptInputs{idx},cachedDesignMethodVal{:}));            
          end
        end
        if isempty(err)
          parseParams.DesignOptsStruct = fdesDesOpts;
        end
      end
      
    end
    
    % Add errors or design the filter
    if ~isempty(unrecognizedInputs)
      if isempty(err)
        if numel(unrecognizedInputs) > 1
          parseParams.ErrorId = 'signal:designfilt:ParametersNotRecognized2';
          parseParams.ErrorMessage = getString(message(parseParams.ErrorId,cell2str(unrecognizedInputs,', ',true)));         
        else
          parseParams.ErrorId = 'signal:designfilt:ParameterNotRecognized2';
          parseParams.ErrorMessage = getString(message(parseParams.ErrorId,unrecognizedInputs{1}));          
        end
        err = parseParams.ErrorMessage;
      else
        % Concatenate errors
        if numel(unrecognizedInputs) > 1
          parseParams.ErrorMessage = sprintf('%s %s',...
            getString(message('signal:designfilt:ParametersNotRecognized1',cell2str(unrecognizedInputs,', ',true))),err);
          parseParams.ErrorId = [parseParams.ErrorId '_ParametersNotRecognized1'];
        else
          parseParams.ErrorMessage = sprintf('%s %s',...
            getString(message('signal:designfilt:ParameterNotRecognized1',unrecognizedInputs{1})),err);
          parseParams.ErrorId = [parseParams.ErrorId '_ParameterNotRecognized1'];
        end
        err = parseParams.ErrorMessage;
      end
    elseif ~isempty(err)
      parseParams.ErrorMessage = err;
    else
      % No parsing errors, so design the filter. If FDESIGN has errors
      % (like invalid frequency values for example, throw error as if the
      % parser was the caller).
      try
        fdes = eval(parseParams.FdesCode);
      catch ME
        throwAsCaller(ME)        
      end
      try
        parseParams.DfiltObj = design(fdes,parseParams.FdesMethods,fdesDesOpts);
      catch ME
        throwAsCaller(ME)        
      end        
      parseParams.DesignedFilter = digitalFilter(parseParams.DfiltObj,parseParams);
    end
  end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
methods (Access = protected)
  %------------------------------------------------------------------------
  function cstrs = getConstraintStrsForFilterFcn(obj,filterFcn)
    % Return cell string with all constraint strings for the given filter
    % response.
    
    filtType = filterFcn(1:end-3);
    respType = filterFcn(end-2:end);
    
    h = getFdesignObject(obj,filtType);
    specs = h.getspeclist;
    
    cstrs = {};
    
    for idx = 1:numel(specs)
      h.Specification = specs{idx};
      m = designmethods(h,respType,'signalonly');
      if ~isempty(m)
        cstrs = [cstrs ; specs{idx}];
      end
    end    
  end
  %--------------------------------------------------------------------------
  function pCell = getAllReqdPropertyNamesFromValidConstraintStr(obj,resp,constraintStr,fullTextFlag,varargin)
    % Return cell vector of property name strings based on
    % comma-separated list of constraint symbols. 
    %
    % This function is called by the filter design assistant GUI to get
    % property names. 
    
    % Inputs:
    %  resp          - filter response (lowpass, ...) without fir or iir 
    %                  sufix
    %  constraintStr - constraint string
    %  fullTexrFlag  - if true, output is a description of each property
    %                  instead of the property name. If set to 'fdesign'
    %                  then output is the set of fdesign property names.
    
    assert(~iscellstr(constraintStr));
    
    if nargin < 4
      fullTextFlag = false;
    end   
    
    % Put each comma-delimited constraint into its own cell
    cCell = csl2cell(constraintStr);
    
    % Translate each into full property strings
    pCell = cell(size(cCell));
    defs = obj.SymbolDefs;
    allConstraintStrs = defs(:,1);
    
    N = numel(cCell);
    for i = 1:N
      idx = strcmpi(cCell{i},allConstraintStrs);
      assert(~isempty(idx),...
        getString(message('signal:designfilt:UnrecognizedConstraintSymbol',cCell{i})));
      if ~ischar(fullTextFlag) && fullTextFlag
        pCell{i} = defs{idx,3};
      elseif ischar(fullTextFlag) && strcmp(fullTextFlag,'fdesign')
        pCell{i} = defs{idx,4};
      else
        pCell{i} = defs{idx,2};
      end
    end
    
    % arbag response with multiple bands
    if strcmp(resp,'arbmag') && strcmpi(constraintStr,'n,b,f,a')
      pCell = setdiff(pCell,{'Frequencies','Amplitudes'});
      idx = strcmp(obj.PropertyNames,'NumBands');
      numBands = 10;
      if any(idx)
        numBands = obj.PropertyValues{idx};        
      end      
      for idx = 1:numBands
        bandStr = num2str(idx);
        pCell = [pCell; {sprintf('BandFrequencies%s',bandStr); ...
          sprintf('BandAmplitudes%s',bandStr)}];
      end
    end            
  end  
  %--------------------------------------------------------------------------
  function designMethods = getDesignMethodsForConstraintStr(obj,filterFcn,cstr)
    %getDesignMethodsForConstraintStr Get design methods for filter function.
    % If input cstr is passed, then we find design methods for that
    % particular constraint set. designMethods output has fdesign format
    
    filtType = filterFcn(1:end-3);
    respType = filterFcn(end-2:end);
    
    h = getFdesignObject(obj,filtType);
    
    if nargin == 3
      h.Specification = cstr;
      specs = {cstr};
    else
      specs = h.getspeclist;
    end
    
    designMethods = [];
    
    for idx = 1:length(specs)
      h.Specification = specs{idx};
      m = designmethods(h,respType,'signalonly');
      designMethods = [designMethods ; m];
    end
  end
  %--------------------------------------------------------------------------
  function [parseParams,earlyRet] = identifyConstraintSetFromSpecifiedProperties(obj,filterFcn,propertyNames,propertyValues,desMethod)
    % Determine if specified properties match property constraint sets for
    % the given filter response. Return match condition and advice as
    % needed.
    %
    % parseParams is a struct:
    %   .FilterFcn: the filter response string as passed by user
    %
    %   .PropertyNames: property list as specified by user
    %
    %   .Match: string indicating whether specified properties
    %       constitute a supported constraint set, or if they are
    %       close enough to provide useful feedback to user.
    %       'singleComplete'
    %       'singlePartial'
    %       'singleOver'
    %       'multiplePartial'
    %       'multipleOver'
    %       'between'
    %       'none'
    %   .Constraint:
    %      The target constraint string for 'singleComplete',
    %      'singlePartial', and 'singleOver'. A cell string of
    %      multiple constraint strings for 'multiplePartial' or
    %      'multipleOver'. An empty string for 'none'.
    %   .MatchingProperties
    %      A cell vector with the subset of specified properties
    %      that match the target contraint string if
    %      'singleComplete', 'singlePartial', 'singleOver'.  Also a
    %      cell vector for a limited case of 'multiplePartial'.
    %      Empty string if 'none', 'multipleOver', and most cases
    %      of 'multiplePartial'.
    %   .OtherProperties
    %      A cell vector of properties that must be specified in
    %      order to fulfill the constraint set if 'singlePartial';
    %      the properties to remove if 'singleOver'; empty for all
    %      other match strings.
    %   .ErrorMessage
    %      An error message for the caller to provide to the user
    %      if an error is warranted; this will be non-empty when
    %      .match is anything but 'singleComplete'.
    
    parseParams.FilterFcn = filterFcn;
    parseParams.PropertyNames = propertyNames;
    parseParams.PropertyValues = propertyValues;
    parseParams.Match = '';
    parseParams.Constraint = '';
    parseParams.MatchingProperties = '';
    parseParams.OtherProperties = '';
    parseParams.ErrorMessage = '';
    parseParams.ErrorId = '';
    parseParams.FdesCode = '';
    parseParams.FdesMethods = '';
    parseParams.InputRecommendation = '';
    parseParams.FinalizeFdesCode = true;    
    
    earlyRet = false;
    cstrs = getConstraintStrsForFilterFcn(obj,filterFcn);
    
    obj.FilterResponse = filterFcn(1:end-3);
    obj.FilterType = filterFcn(end-2:end);
    obj.PropertyNames = propertyNames;
    obj.PropertyValues = propertyValues;
    
    [invalidPropMsg,invalidPropMsgId,errIdx,propertyNames] = checkSpecifiedPropertiesAreValid(obj,parseParams,cstrs);
    % If we find valid property names that do not apply to the current
    % response type, then remove them and try to get fdesign code with as
    % many matches to the user inputs as possible. Cache the error and
    % throw this error after the parser is done. If all properties do not
    % apply, then return with an error.
    propertyNames(errIdx) = [];
    propertyValues(errIdx) = [];
    propertyNames = unique(propertyNames,'stable');    
    parseParams.PropertyNames = propertyNames;
    obj.PropertyNames = propertyNames;
    obj.PropertyValues = propertyValues;

    if isempty(propertyNames)
      parseParams.ErrorMessage = invalidPropMsg;
      parseParams.ErrorId = invalidPropMsgId;
      return;
    end
    
    % Get results from comparing specified set of property names to each
    % constraint set for this filter. Comparison results for each test will
    % be one of: 'match','subset','superset','none','partial'
    [compCell,compIdx] = identifyPropertyComparisons(obj,propertyNames,cstrs);
    [parseParams.Match,cmatch,csub,csuper,cpart] = identifyDetermineMatch(obj,compCell);
    [parseParams,idx] = identifyDetermineConstraint(obj,parseParams,cstrs,cmatch,csub ,csuper,cpart);
    parseParams = identifyDetermineMatchProps(obj,parseParams,propertyNames,compIdx,csub,idx);
    parseParams = identifyDetermineOther(obj,parseParams,propertyNames);
    parseParams = setupFdesignCode(obj,parseParams,desMethod);
    
    if ~isempty(invalidPropMsg)
      parseParams.InputRecommendation = '';
      parseParams.ErrorMessage = invalidPropMsg;
      parseParams.ErrorId = invalidPropMsgId;
    end
    
  end
  
  %--------------------------------------------------------------------------
  function [compCell,compIdx] = identifyPropertyComparisons(obj,propertyNames,cstrs)
    % Compares the property names to the constraint set symbols
    
    Ncs = numel(cstrs);     % constraint sets for this filter
    compCell = cell(Ncs,1); % hold all comparison results
    compIdx = cell(Ncs,1);  % hold specified indices of matching
    % props when 'single*'
    
    for i = 1:Ncs
      % Convert i'th constraint set symbols to property names
      c_i = cstrs{i};
      p_i = getAllReqdPropertyNamesFromValidConstraintStr(obj,obj.FilterResponse,c_i); 
      [compCell{i},compIdx{i}] = comparePropertyNames(obj,propertyNames,p_i);
    end
  end
  %------------------------------------------------------------------------
  function [comp,sidx,ridx] = comparePropertyNames(~,spec,req)
    % Compares specified and required property names.
    % Allows spec to be a string or a cellstring.
    % Returns one of the following conditions:
    %   'match'
    %      all specified match required, AND
    %      all required match specified
    %   'subset'
    %      all specified match required, AND
    %      at least one required does not match specified
    %   'superset'
    %      all required match specified, AND
    %      at least one specified does not match required
    %   'empty'
    %      no specified match required
    %   'partial'
    %      otherwise!  which is:
    %      at least one specified matches required, AND
    %      at least one specified does not match required, AND
    %      at least one required does not match specified
    %
    % Optionally returns indices of specified properties that match
    % required properties when comparison is 'match', 'subset' or
    % 'superset'.  Otherwise returns an empty vector.  Also returns indices
    % of required properties that matchd specified properties under same
    % conditions.
    
    spec = lower(spec);
    req = lower(req);
    sr = ismember(spec,req);
    rs = ismember(req,spec);
    all_sr = all(sr); % careful! goes true for no user specs
    all_rs = all(rs);
    if all_sr && all_rs
      comp = 'match';
      sidx = find(sr);
      ridx = find(rs);
    elseif all_sr && ~all_rs && ~isempty(sr)
      comp = 'subset';
      sidx = find(sr);
      ridx = find(rs);
    elseif all_rs && ~all_sr
      comp = 'superset';
      sidx = find(sr);
      ridx = find(rs);
    elseif ~any(sr)
      comp = 'empty';
      sidx = [];
      ridx = [];
    else
      comp = 'partial';
      sidx = [];
      ridx = [];
    end
  end
  
  %------------------------------------------------------------------------
  function [m,cmatch,csub,csuper,cpart] = identifyDetermineMatch(~,compCell)
    % Determine .Match:
    %
    % The vector of comparisons leads us to choose one of the
    % following overall match results, testing in this order:
    %   'singleComplete'
    %        if any single comparison is 'match'
    %   'multiplePartial'
    %        if two or more comparisons are 'subset'
    %        and no match is 'superset'
    %   'multipleOver'
    %        if two or more comparisons are 'superset'
    %        and no match is 'subset'
    %   'singlePartial'
    %        if any single comparison is 'subset'
    %        and no match is 'superset'
    %   'singleOver'
    %        if any single comparison is 'superset'
    %        and no match is 'subset'
    %   'between'
    %        if at least one 'superset' and
    %        at least one 'subset' is present, or
    %        any 'partial'
    %   'none'
    %        otherwise
    cmatch = strcmpi(compCell,'match');
    csub = [];
    csuper = [];
    cpart = [];
    
    if any(cmatch)
      m = 'singleComplete';
    else
      csub   = strcmpi(compCell,'subset');
      csuper = strcmpi(compCell,'superset');
      if any(csub) || any(csuper)
        if ~any(csuper)
          if sum(csub)>1
            m = 'multiplePartial';
          else
            m = 'singlePartial';
          end
        elseif ~any(csub)
          if sum(csuper)>1
            m = 'multipleOver';
          else
            m = 'singleOver';
          end
        else % if any(csub) && any(csuper)
          % If we have subsets and supersets,give priority to partial
          % completions as they can have default ripple and attenuation
          % inputs missing.
          if sum(csub)>1
            m = 'multiplePartial';
          else
            m = 'singlePartial';
          end                                              
        end
      else
        cpart = strcmpi(compCell,'partial');
        if any(cpart)
          m = 'between';
        else
          m = 'none';
        end
      end
    end
  end
  
  %------------------------------------------------------------------------
  function [parseParams,idx] = identifyDetermineConstraint(~,parseParams,cstrs,cmatch,csub,csuper,cpart)
    % Determine .Constraint:
    % The target constraint string for 'singleComplete', 'singlePartial',
    % and 'singleOver'. A cell string for 'multiplePartial' or
    % 'multipleOver'. An empty string for 'none'.
    
    idx = [];
    if strncmpi(parseParams.Match,'single',6)
      switch parseParams.Match
        case 'singleComplete'
          idx = find(cmatch);
        case 'singlePartial'
          idx = find(csub);
        case 'singleOver'
          idx = find(csuper);
        otherwise
          assert(false); % no other string should match
      end
      c = cstrs(idx);
    elseif strncmpi(parseParams.Match,'multiple',8)
      % Multiple - logical indexing, return cell with 2 or more
      % entries
      if sum(csub)>1
        c = cstrs(csub);
      else
        c = cstrs(csuper);
      end
    elseif strcmpi(parseParams.Match,'between')
      % Possibly multiple
      c = cstrs(cpart);
    else
      % No constraint info is relevant for no match ('none')
      c = {};
    end
    
    parseParams.Constraint = c;
  end
  %------------------------------------------------------------------------
  function parseParams = identifyDetermineMatchProps(~,parseParams,propertyNames,compIdx,csub,idx)
    % Determine .matchingProperties:
    %   .matchingProperties
    %      A cell vector with the subset of specified properties
    %      that match the target contraint string if
    %      'singleComplete', 'singlePartial', 'singleOver'.
    %      Empty string if 'none', 'multiplePartial',
    %      'multipleOver', or 'between'.
    if strncmpi(parseParams.Match,'single',6)
      parseParams.MatchingProperties = propertyNames(compIdx{idx});
    elseif strcmpi(parseParams.Match,'multiplePartial')
      % It's tempting to want to give advice in this case!
      % Only in one subset of this case can that occur.
      %
      % When the SAME subset of specified properties are a subset of
      % multiple required properties, we have some nice guidance to give.
      % When if A,B,C,D,E are specified, and B,C,D are a subset of multiple
      % constraints, we say so. The user is close.
      %
      % Otherwise it's too confusing: if A,B,C,D,E are specified, it could
      % be that A,B are a subset of properties for one constraint, while
      % B,C forms a subset of another, and D,E a subset of yet another.
      % What advice do we give? None.  The user is close to many options --
      % none are more likely than others.
      %
      % See if the same subset of specified props is matching a subset of
      % all constraint sets that are a partial match. The indices within
      % the required constraint set partial matches may vary, but the same
      % indices into the specified properties must remain the same.  Whew.
      %
      % We can rely on c being cstrs(csub) here. It's the set of constraint
      % strings that we partially match.
      %
      % Check compIdx entries for 1 and 5. These vectors of indices must
      % all be the same. If so, we have the special case we seek. any one
      % of those (identical) vectors indicates the specified properties
      % that are in common across all constraint sets.  We report those as
      % .matchingProperties.
      allToCheck = compIdx(csub);
      inCommon = isequal(allToCheck{:}); % works for multiple input args (>2!)
      if inCommon
        specIdx = allToCheck{1}; % or any of them - they are all equal
        parseParams.MatchingProperties = propertyNames(specIdx);
      end
    end
  end
  %------------------------------------------------------------------------
  function parseParams = identifyDetermineOther(obj,parseParams,propertyNames)
    % Determine .otherProperties: A cel of cell vectors of properties that
    % must be specified in order to fulfill the constraint if
    % 'singlePartial'; the properties to remove if 'singleOver'; empty for
    % all other match strings.
    switch parseParams.Match
      case 'singlePartial'
        % Get "canonical" capitalization, so use index "ia" and not the
        % first LHS from setdiff
        p_i = getAllReqdPropertyNamesFromValidConstraintStr(obj,obj.FilterResponse,parseParams.Constraint{1});
        [~,ia] = setdiff(lower(p_i),lower(parseParams.MatchingProperties));
        o = {p_i(ia)};
        
      case 'singleOver'
        o = {setxor(propertyNames,parseParams.MatchingProperties).'};
        
      case 'multiplePartial'
        mp = parseParams.MatchingProperties;
        if ~isempty(mp)
          % We have common props identified. Construct one cell per
          % constraint, indicating the props yet to be specified for
          % each.
          
          Nc = numel(parseParams.Constraint);
          o = cell(Nc,1); % one entry per possible matching constraint
          for i = 1:Nc
            p_i = getAllReqdPropertyNamesFromValidConstraintStr(obj,obj.FilterResponse,parseParams.Constraint{i});
            o{i} = setdiff(p_i,parseParams.MatchingProperties);
          end
        else
          o = {{''}};
        end
      otherwise
        o = {{''}};
    end
    
    % o is a cell array of cells
    parseParams.OtherProperties = o;
  end
  %--------------------------------------------------------------------------
  function parseParams  = setupFdesignCode(obj,parseParams,desMethod)
    %setupFdesignCode Builds the fdesign code and stores into the structure
    
    invalidDesMethodForDefaultValues = false;
    isfir = ~isempty(strfind(parseParams.FilterFcn,'fir'));
    
    % Setup fdesign code using the perfect match or the best match
    if strcmpi(parseParams.Match,'none')
      parseParams.ErrorId = 'signal:designfilt:MustSpecifyValidSpecifications';
      parseParams.ErrorMessage = getString(message(parseParams.ErrorId));
    else
      
      hf = getFdesignObject(obj,parseParams.FilterFcn(1:end-3));
      
      findBestMatch = true;
      if any(strcmpi(parseParams.Match,{'singlePartial','multiplepartial'}))
        % Find spec where only attenuation values are missing
        
        bestMatchCandidates = {};
        bestMatchCandidatesAll = {};
        
        otherprops = parseParams.OtherProperties;        
        for idx = 1:numel(parseParams.Constraint)
          matchIdx = false(numel(otherprops{idx}),1);
          matchIdx = (matchIdx | strncmpi(otherprops{idx},'PassbandRipple',length('PassbandRipple')));
          matchIdx = (matchIdx | strncmpi(otherprops{idx},'StopbandAttenuation',length('StopbandAttenuation')));
          otherprops{idx}(matchIdx) = [];
          
          parseParams.ConstraintNames{idx} = getAllReqdPropertyNamesFromValidConstraintStr(obj,obj.FilterResponse,parseParams.Constraint{idx});
          
          if isempty(otherprops{idx})
            if isempty(desMethod)
              bestMatchCandidates = [bestMatchCandidates parseParams.Constraint{idx}];
            else
              bestMatchCandidatesAll = [bestMatchCandidates parseParams.Constraint{idx}];
              hf.Specification =  parseParams.Constraint{idx};
              dm = designmethods(hf,'signalonly',parseParams.FilterFcn(end-2:end));
              specDesMethod = obj.convertdesignmethodnames(desMethod,isfir,'designfiltToFdesign');
              if any(ismember(specDesMethod,dm))
                bestMatchCandidates = [bestMatchCandidates parseParams.Constraint{idx}];
              end
            end
          end
        end
        
        if ~isempty(desMethod) && ~isempty(bestMatchCandidatesAll) && isempty(bestMatchCandidates)
          % There were matches but design method is not correct for any of the
          % matches
          bestMatch = bestMatchCandidatesAll{1};
          invalidDesMethodForDefaultValues = true;
          p = getAllReqdPropertyNamesFromValidConstraintStr(obj,obj.FilterResponse,bestMatch);
          hf.Specification = bestMatch;
          findBestMatch = false;
          
        elseif ~isempty(bestMatchCandidates)
          bestMatch = bestMatchCandidates{1};
          p = getAllReqdPropertyNamesFromValidConstraintStr(obj,obj.FilterResponse,bestMatch);
          hf.Specification = bestMatch;
          findBestMatch = false;
        end
        
      end
      
      if findBestMatch
        
        % Chose the match that contains the first specified property first.
        % If design method was passed in, choose a constraint set that
        % supports the design method if possible.
        
        % If design method was specified, reduce the constraint list to
        % specs that support the design method.
        if ~isempty(desMethod)
          constraintList = parseParams.Constraint;
          newConstraintIdx =false(1,length(constraintList));
          for idx = 1:length(constraintList)
            hf.Specification = constraintList{idx};
            % fdesign names
            dm = designmethods(hf,'signalonly',parseParams.FilterFcn(end-2:end));
            % Convert user input to fdesign name
            specDesMethod = obj.convertdesignmethodnames(desMethod,isfir,'designfiltToFdesign');
            if any(strcmp(dm,specDesMethod))
              newConstraintIdx(idx) = true;
            end
          end
          
          if any(newConstraintIdx)
            % Reduce the list if specified design method matches the design
            % method of one of the specification sets            
            if numel(parseParams.OtherProperties) == numel(parseParams.Constraint)
              % Make sure OtherProperties is a cell of a cell array
              parseParams.OtherProperties = parseParams.OtherProperties(newConstraintIdx);
            end
            parseParams.Constraint = parseParams.Constraint(newConstraintIdx);
          end
          
        end
        
        chk = [];
        for idx = 1:length(parseParams.Constraint)
          p{idx} = getAllReqdPropertyNamesFromValidConstraintStr(obj,obj.FilterResponse,parseParams.Constraint{idx});
          parseParams.ConstraintNames{idx} = p{idx};
          chk(idx) = find(ismember(parseParams.PropertyNames,p{idx})==1,1);
        end
        
        if ~isempty(chk)
          [~,matchIdx] = min(chk);
          bestMatch = parseParams.Constraint{matchIdx(1)};
          p = p{matchIdx};
        else
          bestMatch = parseParams.Constraint{1};
          p = getAllReqdPropertyNamesFromValidConstraintStr(obj,obj.FilterResponse,bestMatch);
        end
        hf.Specification = bestMatch;
      end
      
      parseParams.FdesMethods = designmethods(hf,parseParams.FilterFcn(end-2:end),'signalonly');
      
      if invalidDesMethodForDefaultValues
        % There were matches but design method is not correct for any of
        % the matches
        
        parseParams.InputRecommendation = '';
        parseParams.ErrorId = 'signal:designfilt:SpecifiedDesignMethodNotValidForInput';
        wText = getString(message(parseParams.ErrorId));        
        parseParams.ErrorMessage = wText;
        
      elseif ~strcmp(parseParams.Match,'singleComplete') && findBestMatch
        msg = getInputRecommendation(obj,parseParams);
        parseParams.InputRecommendation = msg;
        parseParams.ErrorId = 'signal:designfilt:MustSpecifyValidSpecifications';
        wText = getString(message(parseParams.ErrorId));
        parseParams.ErrorMessage = wText;
      else
        parseParams.FullPropertyNames = getAllReqdPropertyNamesFromValidConstraintStr(obj,obj.FilterResponse,parseParams.Constraint{1},true).';
      end
      
      % Generate fdesign code with complete set of specs if no errors have
      % been caught so far. Otherwise the fdesign code will only contain
      % the specstring.
      
      % Add input values if not in error mode
      if isempty(parseParams.ErrorMessage)
        fdesCode = sprintf('fdesign.%s(''%s'',',parseParams.FilterFcn(1:end-3),bestMatch);
        pfdes = getAllReqdPropertyNamesFromValidConstraintStr(obj,obj.FilterResponse,bestMatch,'fdesign');
        
        for idx = 1:length(p)
          propIdx = (strcmp(parseParams.PropertyNames,p{idx})== true);
          if ~any(propIdx)
            %Use default value if property was not specified
            v = hf.(pfdes{idx});
          else
            v = parseParams.PropertyValues{propIdx};
          end
          if idx == length(p)
            fdesCode =[fdesCode mat2str(v)];
          else
            fdesCode =[fdesCode mat2str(v) ','];
          end
        end
      else
        fdesCode = sprintf('fdesign.%s(''%s'')',parseParams.FilterFcn(1:end-3),bestMatch);
        % Tell the caller that the fdesign code string has been completed
        % and that there is no need to add anything else to it.
        parseParams.FinalizeFdesCode = false;
      end
      parseParams.FdesCode = fdesCode;
    end
  end
  %--------------------------------------------------------------------------
  function msg = getInputRecommendation(~,parseParams)
    % Gives the input recommendation to the user to specify the correct
    % parameter set as inputs to the designfilt function.   
    msg = '';
    
    switch lower(parseParams.Match)
      case 'singlepartial'        
        str = cell2str(parseParams.OtherProperties{1},', ');        
        msgStr1 = getString(message('signal:designfilt:HaveSpecifiedFewParameters',parseParams.FilterFcn));
        msgStr2 = getString(message('signal:designfilt:AddFollowingParameterToInputs'));        
        msg = sprintf('%s %s \n  - %s\n',msgStr1,msgStr2,str);
      case 'singleover'
        str = cell2str(parseParams.OtherProperties{1},', ');        
        msgStr1 = getString(message('signal:designfilt:TooManyParameters',parseParams.FilterFcn));
        msgStr2 = getString(message('signal:designfilt:RemoveFollowingParameterFromInputs'));        
        msg = sprintf('%s %s \n  - %s\n',msgStr1,msgStr2,str);
      case 'multiplepartial'         
        for idx = 1:length(parseParams.ConstraintNames)
          str{idx} = cell2str(parseParams.ConstraintNames{idx},', ');
        end        
        msgStr1 = getString(message('signal:designfilt:HaveSpecifiedFewParameters',parseParams.FilterFcn));
        msgStr2 = getString(message('signal:designfilt:ValidParameterSetsCloseToInputs'));                        
        msg = sprintf('%s\n%s\n',msgStr1,msgStr2);
        for idx = 1:length(str)
          msg = [msg sprintf('  - %s\n',str{idx})];
        end
      case 'multipleover'
        for idx = 1:length(parseParams.ConstraintNames)
          str{idx} = cell2str(parseParams.ConstraintNames{idx},', ');
        end        
        msgStr1 = getString(message('signal:designfilt:TooManyParameters',parseParams.FilterFcn));
        msgStr2 = getString(message('signal:designfilt:ValidParameterSetsCloseToInputs'));                        
        msg = sprintf('%s\n%s\n',msgStr1,msgStr2);
        for idx = 1:length(str)
          msg = [msg sprintf('  - %s\n',str{idx})];
        end
    end
  end
  %--------------------------------------------------------------------------
  function [err,errid,errIdx,propertyNames] = checkSpecifiedPropertiesAreValid(obj,parseParams, cstrs)
    % Confirms specified properties are in symbol list.
    
    err = '';
    errid = '';
    errIdx = [];
    invalidPropList = {};
    
    propertyNames = parseParams.PropertyNames;
    propertyValues = parseParams.PropertyValues;
    
    % Check valid values
    for idx = 1:length(propertyNames)
      propIdx = strcmp(obj.SymbolDefs(:,2),propertyNames{idx});
      propAttributes = obj.SymbolDefs{propIdx,5};
      validateattributes(propertyValues{idx},propAttributes{1},propAttributes{2},propertyNames{idx});
    end
    
    % Check valid properties for the current filter function
    p = {};
    for idx = 1:length(cstrs)
      plist = getAllReqdPropertyNamesFromValidConstraintStr(obj,obj.FilterResponse,cstrs{idx})';
      p = [p plist{:}];
    end
    p = unique(p);
    for idx = 1:length(propertyNames)
      d = setdiff(propertyNames{idx},p);
      if ~isempty(d)
        invalidPropList = [invalidPropList propertyNames(idx)];
        errIdx = [errIdx idx];
      end
    end
    
    if ~isempty(invalidPropList)
      if numel(invalidPropList) == 1
        errid = 'signal:designfilt:InvalidPropertyForResponse';
        msgStr = getString(message(errid,invalidPropList{1},parseParams.FilterFcn));
      else
        errid = 'signal:designfilt:InvalidPropertiesForResponse';
        proplistStr = cell2str(invalidPropList,' ,',true);
        msgStr = getString(message(errid,proplistStr,parseParams.FilterFcn));        
      end        
      err = (sprintf('%s',msgStr));
    end
    
  end
  %--------------------------------------------------------------------------
  function [err,errid,requestedResponse] = verifyInputResponse(obj,filterFcn)
    % Confirm filterFcn is valid. Optionally return filter index into
    % constraint set database. CheckInputResponseType - check
    
    requestedResponse = '';    
    errid = '';
    
    if isempty(filterFcn)
      errid = 'signal:designfilt:SpecifyFilterReponseAndPairs';
      err = sprintf(getString(message(errid)));
    else
      %validateStringList returns a cell array
      [outString, ~, err, errType] = obj.validateStringList(obj.InputResponseTypes,filterFcn);
      
      if isempty(err)
        requestedResponse = outString{:};
      elseif errType == 1
        errid = 'signal:designfilt:FilterResponseIsNotValid';
        err = getString(message(errid));
      else
        errid = 'signal:designfilt:FilterResponseInputAmbiguous';
        err = getString(message(errid));
      end
    end
  end
  
  %--------------------------------------------------------------------------
  function [dopts,doptsValues] = getValidDesignOptions(obj,filterFcn,cstr,designMethod)
    % This function gives the valid design options for the particular Filter
    % type and Response type. The designMethod input must be a name in the
    % designfilt format.
    
    filtType = filterFcn(1:end-3);
    respType = filterFcn(end-2:end);
    
    h = getFdesignObject(obj,filtType);
    
    % Get options for specified set of specs or for all the possible specs
    if nargin > 2
      specs = {cstr};
    else
      specs = h.getspeclist;
    end
    
    dopts = {};
    doptsValues = {};
    
    for idx = 1:length(specs)
      h.Specification = specs{idx};
      % Get options for specified design method or for all design methods
      % available for the spec at hand.
      if nargin == 4
        convertedName = signal.internal.DesignfiltProcessCheck.convertdesignmethodnames(designMethod,strcmpi(respType,'fir'),'designfiltToFdesign');
        m = {convertedName};
      else
        m = designmethods(h,respType,'signalonly');
      end
      
      for i = 1:length(m)
        
        opts = designopts(h,m{i},'signalonly');
        p = struct2cell(opts);
        
        doptsValues = [doptsValues; p];
        dopts = [dopts; fieldnames(opts)];
        
      end
      
      % arbmag multiband response
      if strcmp(filtType,'arbmag') && strcmpi(specs{idx},'n,b,f,a')        
        numBandsIdx = strcmp(obj.PropertyNames,'NumBands');
        numBands = 10;
        if any(numBandsIdx)
          numBands = obj.PropertyValues{numBandsIdx};
        end
        for kk = 2:numBands
          bandStr = num2str(kk);
          dopts = [dopts; {sprintf('B%sWeights',bandStr)}];
          doptsValues = [doptsValues; {1}];
        end
      end      
    end
  
  [dopts,IA] = unique(dopts);
  doptsValues = doptsValues(IA);
      
  % Convert names to designfilt format
  dopts = signal.internal.DesignfiltProcessCheck.convertdesignoptnames(dopts,'fdesignToDesignfilt');
  end
  %------------------------------------------------------------------------
  function outVal = validateDesingOptValue(obj,desOptName,desOptValue,desMethod)
    %validateDesingOptValue Validate design option value
    
    % If valid value is a string, we return the corrected value, otherwise,
    % we just return the same value as the input value
         
    outVal = desOptValue;
    
    idx = strcmp(obj.ValidDesignOpts(:,1),desOptName);
    if ~any(idx)
      % Invalid design option, this will be detected somewhere else in the
      % class.
      return;      
    end
            
    propAttributes = obj.ValidDesignOpts{idx,3};
    
    % Correct options that depend on desMethod if a desMethod has been
    % input
    if strcmpi(desOptName,'Window') 
      if ischar(desOptValue)
      try
        validWins = obj.getWindowsList;
        outVal = validatestring(desOptValue, validWins,'designfilt',desOptName); 
      catch        
        error(message('signal:designfilt:InvalidWindowString',desOptValue));
      end
      elseif iscell(desOptValue)
        if numel(desOptValue) < 2 || (~isa(desOptValue{1},'function_handle') && ~ischar(desOptValue{1}))
          error(message('signal:designfilt:InvalidWindowCell'));
        end        
        if ischar(desOptValue{1})
          outVal{1} = validatestring(desOptValue{1},obj.getWindowsList,'Window name');
        end
      end
    end
          
    if nargin == 4
      if strcmpi(desOptName,'MatchExactly') && ~isempty(desMethod) && any(strcmp({'butter','cheby1','cheby2'},desMethod))
        propAttributes = {{'passband','stopband'}};
      end
    end
    
    if strcmp(obj.ValidDesignOpts(idx,2),'validateattributes')
      validateattributes(desOptValue, propAttributes{1},propAttributes{2},'designfilt',desOptName);
    else      
      outVal = validatestring(desOptValue, propAttributes{1},'designfilt',desOptName);      
    end  
    
  end          
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
methods (Access = private)
  
  function initFiltObj(obj)
    %initFiltObj Defines Filter Database
    
    obj.InputResponseTypes = {'arbmagfir',...
      'bandpassfir','bandpassiir',...
      'bandstopfir','bandstopiir',...
      'differentiatorfir',...
      'highpassfir','highpassiir',...
      'hilbertfir',...
      'lowpassfir', 'lowpassiir'};
    
    % SymbolDefs Columns
    % 1-fdesign specstring symbol
    % 2-designfilt property names corresponding to specstring symbol
    % 3-description in plain english
    % 4-fdesign property name
    % 5-validate attribute inputs
    obj.SymbolDefs = { ...
      'N', 'FilterOrder','filter order','FilterOrder',{{'double'},{'scalar','integer','positive'}}; ...
      'Nb', 'NumeratorOrder','numerator order','NumOrder',{{'double'},{'scalar','integer','positive'}}; ...
      'Na', 'DenominatorOrder','denominator order','DenOrder',{{'double'},{'scalar','integer','positive'}}; ...
      'Fp', 'PassbandFrequency','passband frequency','Fpass',{{'double'},{'scalar','positive'}}; ...
      'Fp1', 'PassbandFrequency1','first passband frequency','Fpass1',{{'double'},{'scalar','positive'}}; ...
      'Fp2', 'PassbandFrequency2','second passband frequency','Fpass2',{{'double'},{'scalar','positive'}}; ...
      'Fst', 'StopbandFrequency','stopband frequency','Fstop',{{'double'},{'scalar','positive'}}; ...
      'Fst1', 'StopbandFrequency1','first stopband frequency','Fstop1',{{'double'},{'scalar','positive'}}; ...
      'Fst2', 'StopbandFrequency2','second stopband frequency','Fstop2',{{'double'},{'scalar','positive'}}; ...
      'Fc', 'CutoffFrequency','cutoff frequency','Fcutoff',{{'double'},{'scalar','positive'}}; ...
      'Fc1', 'CutoffFrequency1','first cutoff frequency','Fcutoff1',{{'double'},{'scalar','positive'}}; ...
      'Fc2', 'CutoffFrequency2','second cutoff frequency','Fcutoff2',{{'double'},{'scalar','positive'}}; ...
      'F3db', 'HalfPowerFrequency','3dB frequency','F3dB',{{'double'},{'scalar','positive'}}; ...
      'F3db1', 'HalfPowerFrequency1','first 3dB frequency','F3dB1',{{'double'},{'scalar','positive'}}; ...
      'F3db2', 'HalfPowerFrequency2','second 3dB frequency','F3dB2',{{'double'},{'scalar','positive'}}; ...
      'Ap', 'PassbandRipple','passband ripple','Apass',{{'double'},{'scalar','positive'}}; ...
      'Ap1', 'PassbandRipple1','first passband ripple','Apass1',{{'double'},{'scalar','positive'}}; ...
      'Ap2', 'PassbandRipple2','second passband ripple','Apass2',{{'double'},{'scalar','positive'}}; ...
      'Ast', 'StopbandAttenuation','stopband attenuation','Astop',{{'double'},{'scalar','positive'}}; ...
      'Ast1', 'StopbandAttenuation1','first stopband attenuation','Astop1',{{'double'},{'scalar','positive'}}; ...
      'Ast2', 'StopbandAttenuation2','second stopband attenuation','Astop2',{{'double'},{'scalar','positive'}}; ...
      'TW', 'TransitionWidth','transition width','TransitionWidth',{{'double'},{'scalar','positive'}}; ...
      'B', 'NumBands','number of bands','NBands',{{'double'},{'scalar','integer','positive'}}; ...
      'F', 'Frequencies','vector of frequencies','Frequencies',{{'double'},{'vector'}}; ...
      'A', 'Amplitudes','vector of amplitudes','Amplitudes',{{'double'},{'vector'}}};
      
    for idx = 1:10
      band = num2str(idx);
      obj.SymbolDefs = [obj.SymbolDefs ; {'',sprintf('BandFrequencies%s',band),sprintf('vector of frequencies for band %s',band),sprintf('B%sFrequencies',band),{{'double'},{'vector'}}}];   
    end
    
   for idx = 1:10
      band = num2str(idx);
      obj.SymbolDefs = [obj.SymbolDefs ; {'',sprintf('BandAmplitudes%s',band),sprintf('vector of amplitudes for band %s',band),sprintf('B%sAmplitudes',band),{{'double'},{'vector'}}}];   
   end    
    
   % Valid design options colums:
   % 1 - Design opt property name
   % 2 - Design opt type - can be 'validatestrings' or 'validateattributes'
   % 3 - If type is validatestrings, then cell array with inputs to validatestrings
   %     If type is validateattributes, then cell array with inputs to validateattributes
    obj.ValidDesignOpts = {...
      'Weights','validateattributes', {{'double'},{'vector','positive','finite','real'}}; ...
      'Window', 'validateattributes', {{'double','function_handle','cell','char'},{}}; ...
      'ScalePassband', 'validateattributes', {{'logical'},{'scalar'}}; ...
      'Zerophase', 'validateattributes', {{'logical'},{'scalar'}}; ...
      'PassbandOffset', 'validateattributes', {{'double'},{'vector','nonnegative','finite','real'}}; ...
      'PassbandWeight', 'validateattributes', {{'double'},{'scalar','positive','finite','real'}}; ...
      'PassbandWeight1', 'validateattributes', {{'double'},{'scalar','positive','finite','real'}}; ...
      'PassbandWeight2', 'validateattributes', {{'double'},{'scalar','positive','finite','real'}}; ...      
      'StopbandWeight', 'validateattributes', {{'double'},{'scalar','positive','finite','real'}}; ...
      'StopbandWeight1', 'validateattributes', {{'double'},{'scalar','positive','finite','real'}}; ...
      'StopbandWeight2', 'validateattributes', {{'double'},{'scalar','positive','finite','real'}}; ...
      'MatchExactly', 'validatestrings', {{'passband','stopband','both'}}};
   
   for idx = 1:10
      band = num2str(idx);
      obj.ValidDesignOpts = [obj.ValidDesignOpts ; {sprintf('BandWeights%s',band),'validateattributes',{{'double'},{'vector','positive','finite','real'}}}];   
   end    
   
  end
  
  %------------------------------------------------------------------------
  function h = getFdesignObject(~,filtType)
    % Returns fdesign object for specified filter type. Required at
    % compile time.
    switch filtType
      case 'arbmag'
        h = fdesign.arbmag;
      case 'bandpass'
        h = fdesign.bandpass;
      case 'bandstop'
        h = fdesign.bandstop;
      case 'differentiator'
        h = fdesign.differentiator;
      case 'highpass'
        h = fdesign.highpass;
      case 'hilbert'
        h = fdesign.hilbert;
      case 'lowpass'
        h = fdesign.lowpass;
    end
  end  
end
methods(Static)
  % Methods are also used in filterbuilder code
  %------------------------------------------------------------------------
  function desmethods = convertdesignmethodnames(desmethods,isfir,type,outFormat)
    %CONVERTDESIGNMETHODNAME Convert design method name from fdesign to
    %designfilt and vice versa
    
    if nargin < 4
      outFormat = 'short';
    end
    
    if nargin < 3
      type = 'fdesignToDesignfilt';
    end
    
    if nargin < 2
      isfir = true;
    end
    
    notCell = false;
    if ~iscell(desmethods)
      desmethods = {desmethods};
      notCell = true;
    end
    
    for idx = 1:length(desmethods)
      method = desmethods{idx};
      if strcmp(type,'designfiltToFdesign')        
        % These are names used in designfilt which have to be converted to
        % fdesign names.
        switch method
          case 'ls'
            if isfir
              method = 'firls';
            else
              method = 'iirls';
            end
          case 'cls'
            if isfir
              method = 'fircls';
            end
          otherwise
            method = lower(method);
        end
      elseif strcmp(type,'fdesignToDesignfilt')
        if strcmp(outFormat,'short')
          switch method
            case {'firls','iirls'}
              method = 'ls';
            case 'fircls'
              method = 'cls';
            otherwise
              method = lower(method);          
          end
        else
          switch method
            case 'equiripple'
              method = 'Equiripple';
            case 'ellip'
              method = 'Elliptic';
            case 'butter'
              method = 'Butterworth';
            case 'cheby1'
              method = 'Chebyshev type I';
            case 'cheby2'
              method = 'Chebyshev type II';
            case 'kaiserwin'
              method = 'Kaiser window';
            case {'firls','iirls'}
              method = 'Least-squares';
            case 'fircls'
              method = 'Constrained least-squares';
            case 'freqsamp'
              method = 'Frequency-sampling';
            case 'maxflat'
              method = 'Maximally flat';
            otherwise
              method = lower(method);          
          end                                                              
        end
      else        
        error(message('signal:designfilt:InvalidInputArgument','TYPE'))
      end
      desmethods{idx} = method;
    end
    
    if notCell
      desmethods = desmethods{:};
    end
  end
  %------------------------------------------------------------------------
  function names = convertdesignoptnames(names,type)
    %CONVERTDESIGNOPTNAMES Convert design option names from those in
    %fmethod to those in the filter designer object and vice versa
    
    if nargin == 1
      type = 'fdesignToDesignfilt';
    end
    
    notCell = false;
    if ~iscell(names)
      names = {names};
      notCell = true;
    end
    
    for idx = 1:length(names)
      name = names{idx};
      if strcmp(type,'fdesignToDesignfilt')
        switch name
          case 'Zerophase'
            name = 'ZeroPhase';
          case 'Wstop'
            name = 'StopbandWeight';
          case 'Wstop1'
            name = 'StopbandWeight1';
          case 'Wstop2'
            name = 'StopbandWeight2';
          case 'Wpass'
            name = 'PassbandWeight';
          case 'Wpass1'
            name = 'PassbandWeight1';
          case 'Wpass2'
            name = 'PassbandWeight2';
          case {'B1Weights','B2Weights','B3Weights','B4Weights','B5Weights',...
              'B6Weights','B7Weights','B8Weights','B9Weights','B10Weights'}
            if strcmp(name(3),'0')
              numStr = '10';
            else
              numStr = name(2);
            end
            name = ['BandWeights' numStr];
        end
      elseif strcmp(type,'designfiltToFdesign')
        switch name
          case 'ZeroPhase'
            name = 'Zerophase';
          case 'StopbandWeight'
            name = 'Wstop';
          case 'StopbandWeight1'
            name = 'Wstop1';
          case 'StopbandWeight2'
            name = 'Wstop2';
          case 'PassbandWeight'
            name = 'Wpass';
          case 'PassbandWeight1'
            name = 'Wpass1';
          case 'PassbandWeight2'
            name = 'Wpass2';
          case {'BandWeights1','BandWeights2','BandWeights3','BandWeights4',...
              'BandWeights5','BandWeights6','BandWeights7','BandWeights8',...
              'BandWeights9','BandWeights10'}
            numStr = name(end);
            if strcmp(numStr,'0')
              numStr = '10';
            end
            name = ['B' numStr 'Weights'];
        end
      else
        error(message('signal:designfilt:InvalidInputArgument','TYPE'))
      end
      names{idx} = name;
    end
    
    if notCell
      names = names{:};
    end
    
  end
  %----------------------------------------------------------------------
  function props = convertpropnames(props,type,opt)
    %CONVERTPROPNAMES Convert property names from filter designer to
    %fspecs or vice versa
    
    if nargin == 1
      type = 'fdesignToDesignfilt';
      opt = 'short';
    elseif nargin == 2
      opt = 'short';
    end
    
    notCell = false;
    if ~iscell(props)
      props = {props};
      notCell = true;
    end
    
    isShortName = strcmpi(opt,'short');
    
    % From fspec to filter designer names
    if strcmp(type,'fdesignToDesignfilt')
      for idx = 1:length(props)
        switch props{idx}
          case 'Fs'
            if isShortName
              props{idx} = 'SampleRate';
            else
              props{idx} = 'Sample rate';
            end
          case 'Fpass'
            if isShortName
              props{idx} = 'PassbandFrequency';
            else
              props{idx} = 'Passband frequency';
            end
          case 'Fpass1'
            if isShortName
              props{idx} = 'PassbandFrequency1';
            else
              props{idx} = 'First passband frequency';
            end
          case 'Fpass2'
            if isShortName
              props{idx} = 'PassbandFrequency2';
            else
              props{idx} = 'Second passband frequency';
            end
          case 'Fstop'
            if isShortName
              props{idx} = 'StopbandFrequency';
            else
              props{idx} = 'Stopband frequency';
            end
          case 'Fstop1'
            if isShortName
              props{idx} = 'StopbandFrequency1';
            else
              props{idx} = 'First stopband frequency';
            end
          case 'Fstop2'
            if isShortName
              props{idx} = 'StopbandFrequency2';
            else
              props{idx} = 'Second stopband frequency';
            end
          case 'Fcutoff'
            if isShortName
              props{idx} = 'CutoffFrequency';
            else
              props{idx} = 'Cutoff frequency';
            end
          case 'Fcutoff1'
            if isShortName
              props{idx} = 'CutoffFrequency1';
            else
              props{idx} = 'First cutoff frequency';
            end
          case 'Fcutoff2'
            if isShortName
              props{idx} = 'CutoffFrequency2';
            else
              props{idx} = 'Second cutoff frequency';
            end
          case 'F3dB'
            if isShortName
              props{idx} = 'HalfPowerFrequency';
            else
              props{idx} = 'Half power frequency';
            end
          case 'F3dB1'
            if isShortName
              props{idx} = 'HalfPowerFrequency1';
            else
              props{idx} = 'First half power frequency';
            end
          case 'F3dB2'
            if isShortName
              props{idx} = 'HalfPowerFrequency2';
            else
              props{idx} = 'Second half power frequency';
            end
          case 'Apass'
            if isShortName
              props{idx} = 'PassbandRipple';
            else
              props{idx} = 'Passband ripple';
            end
          case 'Apass1'
            if isShortName
              props{idx} = 'PassbandRipple1';
            else
              props{idx} = 'First passband ripple';
            end
          case 'Apass2'
            if isShortName
              props{idx} = 'PassbandRipple2';
            else
              props{idx} = 'Second passband ripple';
            end
          case 'Astop'
            if isShortName
              props{idx} = 'StopbandAttenuation';
            else
              props{idx} = 'Stopband attenuation';
            end
          case 'Astop1'
            if isShortName
              props{idx} = 'StopbandAttenuation1';
            else
              props{idx} = 'First stopband attenuation';
            end
          case 'Astop2'
            if isShortName
              props{idx} = 'StopbandAttenuation2';
            else
              props{idx} = 'Second stopband attenuation';
            end
          case 'NumOrder'
            if isShortName
              props{idx} = 'NumeratorOrder';
            else
              props{idx} = 'Numerator order';
            end
          case 'DenOrder'
            if isShortName
              props{idx} = 'DenominatorOrder';
            else
              props{idx} = 'Denominator order';
            end
          case {'NBands'}
            if isShortName
              props{idx} = 'NumBands';
            else
              props{idx} = 'Numer of bands';
            end
          case {'B1Frequencies','B2Frequencies','B3Frequencies','B4Frequencies',...
              'B5Frequencies','B6Frequencies','B7Frequencies','B8Frequencies',...
              'B9Frequencies','B10Frequencies'}
            if strcmp(props{idx}(3),'0')
              numStr = '10';
            else
              numStr = props{idx}(2);
            end
            if isShortName
              props{idx} = ['BandFrequencies' numStr];
            else
              props{idx} = ['Band frequencies' numStr];
            end
          case {'B1Amplitudes','B2Amplitudes','B3Amplitudes','B4Amplitudes',...
              'B5Amplitudes','B6Amplitudes','B7Amplitudes','B8Amplitudes',...
              'B9Amplitudes','B10Amplitudes'}
            if strcmp(props{idx}(3),'0')
              numStr = '10';
            else
              numStr = num2str(props{idx}(2));
            end
            if isShortName
              props{idx} = ['BandAmplitudes' numStr];
            else
              props{idx} = ['Band amplitudes' numStr];
            end
            % used by filterbuilder ------------------
          case 'B'
            if isShortName
              props{idx} = 'NumBands';
            else
              props{idx} = 'Number of bands';
            end
          case 'N'
            if isShortName
              props{idx} = 'FilterOrder';
            else
              props{idx} = 'Filter order';
            end
          case 'TW'
            if isShortName
              props{idx} = 'TransitionWidth';
            else
              props{idx} = 'Transition width';
            end
          case 'F'
            if isShortName
              props{idx} = 'Frequencies';
            else
              props{idx} = 'Frequencies';
            end
          case 'A'
            if isShortName
              props{idx} = 'Amplitudes';
            else
              props{idx} = 'Amplitudes';
            end
          case 'Nb'
            if isShortName
              props{idx} = 'NumeratorOrder';
            else
              props{idx} = 'Numerator order';
            end
          case 'Na'
            if isShortName
              props{idx} = 'DenominatorOrder';
            else
              props{idx} = 'Denominator order';
            end
          case 'F6dB'
            if isShortName
              props{idx} = 'CutoffFrequency';
            else
              props{idx} = 'Cutoff frequency';
            end            
          case 'F6dB1'
            if isShortName
              props{idx} = 'CutoffFrequency1';
            else
              props{idx} = 'First cutoff frequency';
            end
          case 'F6dB2'
            if isShortName
              props{idx} = 'CutoffFrequency2';
            else
              props{idx} = 'Second cutoff frequency';
            end

            
          case {'A1','A2','A3','A4','A5','A6','A7','A8','A9','A10'}
            numStr = props{idx}(2:end);
            if isShortName
              props{idx} = ['BandAmplitudes' numStr];
            else
              props{idx} = ['Band amplitudes' numStr];
            end
          case {'F1','F2','F3','F4','F5','F6','F7','F8','F9','F10'}
            numStr = props{idx}(2:end);
            if isShortName
              props{idx} = ['BandFrequencies' numStr];
            else
              props{idx} = ['Band frequencies' numStr];
            end
          otherwise
            props{idx} = props{idx};
        end
      end
    elseif strcmp(type,'designfiltToFdesign')
      % From filter designer to fspec names
      for idx = 1:length(props)
        switch props{idx}
          case 'SampleRate'
            props{idx} = 'Fs';
          case 'PassbandFrequency'
            props{idx} = 'Fpass';
          case 'PassbandFrequency1'
            props{idx} = 'Fpass1';
          case 'PassbandFrequency2'
            props{idx} = 'Fpass2';
          case 'StopbandFrequency'
            props{idx} = 'Fstop';
          case 'StopbandFrequency1'
            props{idx} = 'Fstop1';
          case 'StopbandFrequency2'
            props{idx} = 'Fstop2';
          case 'CutoffFrequency'
            props{idx} = 'Fcutoff';
          case 'CutoffFrequency1'
            props{idx} = 'Fcutoff1';
          case 'CutoffFrequency2'
            props{idx} = 'Fcutoff2';
          case 'HalfPowerFrequency'
            props{idx} = 'F3dB';
          case 'HalfPowerFrequency1'
            props{idx} = 'F3dB1';
          case 'HalfPowerFrequency2'
            props{idx} = 'F3dB2';
          case 'PassbandRipple'
            props{idx} = 'Apass';
          case 'PassbandRipple1'
            props{idx} = 'Apass1';
          case 'PassbandRipple2'
            props{idx} = 'Apass2';
          case 'StopbandAttenuation'
            props{idx} = 'Astop';
          case 'StopbandAttenuation1'
            props{idx} = 'Astop1';
          case 'StopbandAttenuation2'
            props{idx} = 'Astop2';
          case 'NumeratorOrder'
            props{idx} = 'NumOrder';
          case 'DenominatorOrder'
            props{idx} = 'DenOrder';
          case 'NumBands'
            props{idx} = 'NBands';
          case {'BandFrequencies1','BandFrequencies2','BandFrequencies3',...
              'BandFrequencies4','BandFrequencies5','BandFrequencies6',...
              'BandFrequencies7','BandFrequencies8','BandFrequencies9',...
              'BandFrequencies10'}
            numStr = props{idx}(end);
            if strcmp(numStr,'0')
              props{idx} = 'B10Frequencies';
            else
              props{idx} = ['B' numStr 'Frequencies'];
            end
          case {'BandAmplitudes1','BandAmplitudes2','BandAmplitudes3',...
              'BandAmplitudes4','BandAmplitudes5','BandAmplitudes6',...
              'BandAmplitudes7','BandAmplitudes8','BandAmplitudes9',...
              'BandAmplitudes10'}
            numStr = props{idx}(end);
            if strcmp(numStr,'0')
              props{idx} = 'B10Amplitudes';
            else
              props{idx} = ['B' numStr 'Amplitudes'];
            end
          otherwise
            props{idx} = props{idx};
        end
      end
    else
      error(message('signal:designfilt:InvalidInputArgument','TYPE'))
    end
    
    if notCell
      props = props{:};
    end
    
  end
  %------------------------------------------------------------------------
  function list = getWindowsList()
  %getWindowsList List of valid windows
    list = {'bartlett','blackman','blackmanharris','bohmanwin',...
      'chebwin','flattopwin','gausswin','hamming','hann','kaiser',...
      'nuttallwin','parzenwin','rectwin','taylorwin','triang','tukeywin'};
  end
  %------------------------------------------------------------------------
  function list = getDesignMethodsList()
  %getDesignMethodsList List of valid design methods
    
    list = {'equiripple','kaiserwin','window','ls',...
      'cls','freqsamp','butter','cheby1','cheby2','ellip','maxflat'};
  end
  %------------------------------------------------------------------------
  function [designfiltProps,designfiltPropValues] = getSpecPropertiesFromFdesign(fdesObj,fdesignFlag)
    %getSpecProperties Get spec property list from fdesign object
    %
    % This function looks at the fdesign object and strips down the
    % specification properties relevant to the designfilt function.
    % Inputs:
    %   fdesObj - FDESIGN object
    %   fdesignFlag - if set to true then we return property names in
    %   fdesign fortmat. Its default is false.     
    % Outputs:
    %   designfiltProps - cell array with property names in designfilt
    %   format or in fdesign format if fdesignFlag = true.
    %   designfiltPropValues - cell array with values for each property

    if nargin == 1
      fdesignFlag = false;
    end
    
    fdesStruct = get(fdesObj);
    
    fdesPropNames = fieldnames(fdesStruct);
    fdesPropValues = struct2cell(fdesStruct);
    [~, I] = setdiff(fdesPropNames,{'Response','Description','Specification','NormalizedFrequency','Fs'});
    fdesProps = fdesPropNames(sort(I));
    designfiltPropValues = fdesPropValues(sort(I));
    
    % Convert names from fdesign to designfilt
    if fdesignFlag
      designfiltProps = fdesProps;
    else
      designfiltProps = signal.internal.DesignfiltProcessCheck.convertpropnames(fdesProps,'fdesignToDesignfilt');                           
    end
  end
%--------------------------------------------------------------------------
function [testString,ambInputs,err,errType,errid] = validateStringList(listString,testString)
% Validates if the listString set contains the testString.

% errType = 0 means no errors
% errType = 1 means invalid input
% errType = 2 means ambiguous input

err = '';
errid = '';
ambInputs = {};
errType = 0;

if ~iscell(testString)
  testString = {testString};
end

for i = 1:length(testString)
  
  % Check for full match
  idx = strcmpi(listString,testString{i});
  if ~any(idx)
    % Check for partial match
    idx = strncmpi(listString,testString{i},length(testString{i}));
  end
  
  numMatches = sum(idx);
  if numMatches == 0
    errid = 'signal:designfilt:InvalidInputEntered';
    err = getString(message(errid));
    errType = 1;
  elseif numMatches > 1
    ambInputs = [ambInputs, testString{i}];
  else
    testString{i} = listString{idx};
  end  
end

if ~isempty(ambInputs)
  errType = 2;
  list =  cell2str(ambInputs,', ',true);
  if numel(ambInputs) > 1
    errid = 'signal:designfilt:AmbiguousParameterNames';
    err = getString(message(errid,list));    
  else
    errid = 'signal:designfilt:AmbiguousParameterName';
    err = getString(message(errid,list));    
  end
end

if nargout < 3 && ~isempty(err)
  % throw error here 
  error(errid,err)
end

end  
  
end % end of static methods
end % end of class definition
%--------------------------------------------------------------------------
% Helper functions
%--------------------------------------------------------------------------
function m = getDefaultDesignMethod(inputMethods)
%getDefaultDesignMethod Choose default method from a list

if ~iscell(inputMethods)
  inputMethods = {inputMethods};
end

if any(strcmpi(inputMethods,'freqsamp'))
  m = 'freqsamp';
elseif any(strcmpi(inputMethods,'equiripple'))
  m = 'equiripple';
elseif any(strcmpi(inputMethods,'butter'))
  m = 'butter';
else
  m = inputMethods{1};
end
end

%--------------------------------------------------------------------------
function c = csl2cell(str)
% Convert a comma-delimited string into a cell of strings.
c = textscan(str,'%s','delimiter',',');
c = c{1}; % puts obj one level down
end
%--------------------------------------------------------------------------
function [cachedInputs,cachedVals,cachedValsNames,remInputs,remVals,remValsNames] = cacheInputList(listOfInterest,allInputs,allVals,allValsNames)
% Compares allInputs with the listOfInterest and get the matching inputs
% and their values.

cachedInputs    = {};
cachedVals      = {};
cachedValsNames = {};
remInputs       = allInputs;
remVals         = allVals;
remValsNames    = allValsNames;

% Compare allInputs with the listOfInterest and get matching indices
strIdx = false(size(allInputs));
for idx = 1:length(listOfInterest)
  strIdx = strIdx | strcmpi(allInputs,listOfInterest{idx});
end

% Get the matching inputs and their values
if any(strIdx)
  cachedInputs    = allInputs(strIdx);
  cachedVals      = allVals(strIdx);
  cachedValsNames = allValsNames(strIdx);
  remInputs       = allInputs(~strIdx);
  remVals         = allVals(~strIdx);
  remValsNames    = allValsNames(~strIdx);
end
end
%--------------------------------------------------------------------------
function [s,N] = cell2str(c,d,addAndFlag)
% Return delimiter-separated string.
%   cell2str(C,D) returns a single string created by concatenating each
%   string in cell-string C, separating each string by the delimiter string
%   D.  If omitted, D is a comma.

if ~iscellstr(c)
  error(message('signal:designfilt:InputMustBeCellString'));    
end
if nargin < 2
  d = ',';
end
if nargin < 3
  addAndFlag = false;
end

N = numel(c);
s = '';
for i=1:N
  s = [s c{i}];
  if i<N
    if i == (N-1) && N>1 && addAndFlag
      s = [s d 'and '];
    else
      s = [s d];
    end
  end
end
end
