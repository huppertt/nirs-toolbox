function [y,msg,xuser_params] = parse_pv_pairs(candidate_pv, varargin)
% PARSE_PV_PAIRS Parse parameter/value pairs
%   PARSE_PV_PAIRS(param_list, args)
%   PARSE_PV_PAIRS(param_list, args, led_object)
%
% If the last argument of the input list is an LED object,
% it is copied to the return variable for use as default values.
%
% If the candidate_pv list includes a second column, then
% those values are taken as default values for all P/V pairs.
%
% If both an LED object and default values are passed,
% the default values are applied to the LED object.
%
% If the param list does not include default values,
% and no LED object is passed, then only the parameters
% that the user passed are returned.

% Author: D. Orofino
% Copyright 1988-2002 The MathWorks, Inc.

if nargin>1 & isa(varargin{end},'led7'),
   % LED object passed as last option
   % Copy it to output, and remove from varargin list
   y = varargin{end};
   varargin(end) = [];
else
   y = [];
end

msg = '';
xuser_params = {};

setmode = size(candidate_pv,2)>1;

if setmode,
   % Defaults included:
   
   % This is probably the first time we have seen this
   % param list, so check it carefully:
   
   % Check candidate pv list:
   if ~iscell(candidate_pv),
      msg = getString(message('signal:private:CandidateParametervalueListMustBeACellArray'));
      return;
   end
   
   % Get vector cell array of parameter names:
        candidate_params = lower(candidate_pv(:,1));

   % Get params from list of candidate (p,v) pairs
   % Check that candidate params are all strings:
   for i=1:length(candidate_params),
      if ~ischar(candidate_params{i}),
         msg = getString(message('signal:private:CandidateParametersMustBeStrings'));
         return;
      end
   end

   % If default PV pairs are to be retained,
   % uncomment the following line:
   %
   y.candidate_params = candidate_params;

   % Get vector cell array of default parameters:
        candidate_values = candidate_pv(:,2);

   % Record default candidate parameter values:
   for i=1:length(candidate_params), 
      y = setfield(y, candidate_params{i}, candidate_values{i});
   end

else
   % Get vector cell array of parameter names:
   % They should be lowercase, and in a vector.
        candidate_params = candidate_pv;
   
end

% If user passed no options, return:
if nargin<2, return; end

% If first user argument is numeric, we could assume
% that it is a handle to the parent.  Uncomment the
% following lines if this behavior is desired:
%
%if isnumeric(varargin{1}),
%   y.parent = varargin{1};
%   varargin(1)=[];
%end
%
% No need to continue if only a parent handle passed:
% if length(varargin)==0, return; end

% Get params from list of user's (p,v) pairs
passed_pv = (length(varargin)>1);
if passed_pv,
        user_params = lower(varargin(1:2:end));
   user_values = varargin(2:2:end);
else
   user_params = lower(varargin);
   user_values = [];
end

for i=1:length(user_params),
   idx = strmatch(user_params{i}, candidate_params);
   
   % Check for ambiguous property names:
   if length(idx)>1,
      msg = [getString(message('signal:private:AmbiguousProperty')) ': ''' user_params{i} ''''];
      return
      
   elseif isempty(idx),
      % Check for invalid property names:
      
      if ~ischar(user_params{i}),
         % Non-string passed as parameter name:
         msg = getString(message('signal:private:InvalidParametervaluePairArguments'));
      else
         msg = [getString(message('signal:private:InvalidProperty')) ': ''' user_params{i} ''''];
      end
      return
   end
   
   % Record "expanded" (matched) user_param's
   % Don't record duplicates
   if isempty(strmatch(candidate_params{idx}, xuser_params)),
           xuser_params{end+1} = candidate_params{idx};
   end
   
   % Record p,v pair:
   if passed_pv,
      y = setfield(y, candidate_params{idx},user_values{i});
   end
end

