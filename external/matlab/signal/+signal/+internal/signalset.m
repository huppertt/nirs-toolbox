function varargout = signalset(obj,varargin)
%SET  Set object property values
%   SET(obj,'PropertyName',PropertyValue) sets the value of the
%   specified property for the object, obj.
%
%   SET(obj,'PropertyName1',Value1,'PropertyName2',Value2,...) sets
%   multiple property values with a single statement.
%
%   Given a structure S, whose field names are object property names,
%   SET(obj,S) sets the properties identified by each field name of S
%   with the values contained in the structure.
%
%   A = SET(obj, 'PropertyName') returns the possible values for the
%   specified property of the System object, obj. The returned array
%   is a cell array of possible value strings or an empty cell array
%   if the property does not have a finite set of possible string
%   values.
%
%   A = SET(obj) returns all property names and their possible values
%   for the object, obj. The return value is a structure whose
%   field names are the property names of obj, and whose values are
%   cell arrays of possible property value strings or empty cell
%   arrays.
%

% only support scalar (for now)
if numel(obj) ~= 1
  matlab.system.internal.error('signal:sigtools:nonScalarSetGet','set');
end

switch(nargin)
  case 1
    % S = set(obj)
    fns = fieldnames(obj);
    st = [];
    for ii = 1:length(fns)
      fn = fns{ii};
      fnprop = findprop(obj,fn);
      if ~isempty(fnprop) && strcmp(fnprop.SetAccess,'public')
        val = {set(obj,fn)};
        if isempty(val{1})
          st.(fn) = {};
        else
          st.(fn) = val;
        end
      end
    end
    varargout = {st};
  case 2
    if isstruct(varargin{1})
      % set(obj, struct)
      st = varargin{1};
      stfn = fieldnames(st);
      for ii = 1:length(stfn)
        prop = stfn{ii};
        obj.(prop) = st.(prop);
      end
      varargout = {{}};
    else
      prop = varargin{1};
      mp = findprop(obj,prop);
      if isempty(mp)
        matlab.system.internal.error(...
          'signal:sigtools:invalidProperty', prop, class(obj));
      elseif ~strcmp(mp.SetAccess, 'public')
        matlab.system.internal.error(...
          'signal:sigtools:propertyInvalidSetAccess', prop, class(obj));
      end
      varargout = {getAllowedStringValues(obj,prop)};
      if isempty(varargout)
        varargout = {{}};
      end
    end
  otherwise
    % set(obj, <PV Pairs>)
    if mod(length(varargin),2)      
      error(message('signal:sigtools:invalidPvp'));
    end
    for ii = 1:2:length(varargin)
      % we call off to a public function, otherwise this would allow
      % users to set protected/private properties
      obj.(varargin{ii}) =  varargin{ii+1};
    end
    varargout = {{}};
end
end