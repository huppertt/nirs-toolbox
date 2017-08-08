function [options,msg,msgobj] = psdoptions(isreal_x,options,varargin)
%PSDOPTIONS   Parse the optional inputs to most psd related functions.
%
%
%  Inputs:
%   isreal_x             - flag indicating if the signal is real or complex
%   options              - the same structure that will be returned with the
%                          fields set to the default values
%   varargin             - optional input arguments to the calling function
%
%  Outputs:
%  PSDOPTIONS returns a structure, OPTIONS, with following fields:
%   options.nfft         - number of freq. points at which the psd is estimated
%   options.Fs           - sampling freq. if any
%   options.range        - 'onesided' or 'twosided'
%   options.centerdc     - true if 'centered' specified
%   options.conflevel    - confidence level between 0 and 1 ('omitted' when unspecified)
%   options.ConfInt      - this field only exists when called by PMTM 
%   options.MTMethod     - this field only exists when called by PMTM
%   options.NW           - this field only exists when called by PMUSIC (see PMUSIC for explanation)
%   options.Noverlap     - this field only exists when called by PMUSIC (see PMUSIC for explanation)
%   options.CorrFlag     - this field only exists when called by PMUSIC (see PMUSIC for explanation)
%   options.EVFlag       - this field only exists when called by PMUSIC (see PMUSIC for explanation)



%  Authors: D. Orofino and R. Losada
%  Copyright 1988-2012 The MathWorks, Inc.
   
[s,msg,msgobj] = psd_parse(varargin{:});
if ~isempty(msg),
   return
end

omit = 'omitted';

% Replace defaults when appropriate

% If empty or omitted, nfft = default nfft (default is calling function dependent)
if ~any([isempty(s.NFFT),strcmp(s.NFFT,omit)]),    
   options.nfft = s.NFFT;   
end

if ~strcmp(s.Fs,omit), % If omitted, Fs = [], work in rad/sample
   options.Fs = s.Fs;
   if isempty(options.Fs),
      options.Fs = 1;  % If Fs specified as [], use 1 Hz as default
   end
end

% Only PMTM has the field options.ConfInt; 
if isfield(options,'ConfInt'),
   if ~strcmp(s.ConfInt,omit) && ~isempty(s.ConfInt),
      % Override the default option ONLY if user specified a new value
      options.ConfInt = s.ConfInt;
   end
else
   % If a third scalar was specified, error out unless PMUSIC is the calling function
   % (hence options.NW exists)
   if ~strcmp(s.ConfInt,omit) && ~isfield(options,'nw'),  
      msgobj = message('signal:psdoptions:TooManyNumericOptions');
      msg = getString(msgobj);
      return
   end
end

% Only PMUSIC has the field options.nw; 
if isfield(options,'nw'),
   if ~strcmp(s.nw,omit) && ~isempty(s.nw),
      % Override the default option ONLY if user specified a new value
      options.nw = s.nw;
      if ~any(size(options.nw)==1),
         msgobj = message('signal:psdoptions:MustBeScalarOrVector','NW');
         msg = getString(msgobj);
         return
      elseif length(options.nw) > 1,
         options.window = options.nw;
         options.nw = length(options.nw);
      end
   end
else
   % If a third scalar was specified, error out unless PMTM is the calling function
   % (hence options.ConfInt exists)
   if ~strcmp(s.nw,omit) && ~isfield(options,'ConfInt'),  
      msgobj = message('signal:psdoptions:TooManyNumericOptions');
      msg = getString(msgobj);
      return
   end
end

% Only PMUSIC has the field options.noverlap; 
if isfield(options,'noverlap'),
   if ~strcmp(s.noverlap,omit) && ~isempty(s.noverlap),
      % Override the default option ONLY if user specified a new value
      options.noverlap = s.noverlap;
   else
      % Use default
      options.noverlap = options.nw -1;
   end
end
 
options.centerdc = ~strcmp(s.DCFlag,'omitted');
options.conflevel = s.ConfLevel;

if ~strcmp(s.Range,omit),
   options.range = s.Range;
end

if ~isreal_x & strcmpi(options.range,'onesided'), %#ok
   msgobj = message('signal:psdoptions:ComplexInputDoesNotHaveOnesidedPSD');
   msg = getString(msgobj);
   return
end

% Only PMTM has the field options.MTMethod
if ~isfield(options,'MTMethod'),
   if ~strcmp(s.MTMethod,omit),
      % A string particular to pmtm is not supported here
      msgobj = message('signal:psdoptions:UnrecognizedString');
      msg = getString(msgobj);
      return
   end
elseif ~strcmp(s.MTMethod,omit),
   % Override the default option ONLY if user specified a new value
   % psd_parse has already determined the validity of this string
   options.MTMethod = s.MTMethod; 
end

% Only PMUSIC has the field options.CorrFlag
if ~isfield(options,'CorrFlag'),
   if ~strcmp(s.CorrFlag,omit),
      % A string particular to pmusic is not supported here
      msgobj = message('signal:psdoptions:UnrecognizedString');
      msg = getString(msgobj);
      return
   end
elseif ~strcmp(s.CorrFlag,omit),
   % Override the default option ONLY if user specified a new value
   % psd_parse has already determined the validity of this string, we set the flag to one
   options.CorrFlag = 1; 
end

% Only PMUSIC has the field options.EVFlag
if ~isfield(options,'EVFlag'),
   if ~strcmp(s.EVFlag,omit),
      % A string particular to pmusic is not supported here
      msgobj = message('signal:psdoptions:UnrecognizedString');
      msg = getString(msgobj);
      return
   end
elseif ~strcmp(s.EVFlag,omit),
   % Override the default option ONLY if user specified a new value
   % psd_parse has already determined the validity of this string, we set the flag to one
   options.EVFlag = 1; 
end


%------------------------------------------------------------------------------
function [s,msg,msgobj] = psd_parse(varargin)
%PSD_PARSE Parse trailing inputs for PSD functions.
%  [S,MSG] = PSD_PARSE(varargin) parses the input argument list and
%  returns a structure, S, with fields corresponding to each of the
%  possible PSD options.  If an option is omitted, a default value
%  is returned such that the structure will always return a value
%  for every possible option field.
%
%  MSG is a error message returned from the parse operation.
%  It may be empty, in which case no parsing errors have occurred.
%
%  The structure fields are as follows:
%   S.NFFT    - FFT length (scalar).  If it appears in the argument
%               list, it must be the 1st scalar argument, but not
%               necessarily the 1st argument in the list.  A string option
%               could occur before, instead of, or after this option.
%               If omitted, 'omitted' is returned.
%
%   S.Fs      - Sample rate in Hz (scalar).  If it appears in the argument
%               list, it must be the 2nd scalar argument, but not
%               necessarily the 2nd argument in the list.  A string option
%               could occur before, instead of, or after this option.
%               If omitted, 'omitted' is returned.
%
%   S.ConfInt - Confidence interval (scalar).  If it appears in the argument
%   S.nw        list, it must be the 3rd scalar argument, but not
%   (synonyms)  necessarily the 3rd argument in the list.  A string option
%               could occur before, instead of, or after this option.
%               If omitted, 'omitted' is returned.
%               NOTE: Only pmusic expects S.NW at present, and S.ConfInt and
%               S.NW are simply synonyms of each other (both are set to exactly
%               the same values).
%
%   S.ConfLevel Confidence level property value (scalar) between 0 and 1.
%               If omitted, 'omitted' is returned.  This property does not
%               interfere with S.ConfInt which is used only with PMTM.
%
%   S.DCFlag  - CenterDC.  Returns true iff 'centered' is specified in the argument
%               list.  Recognized argument strings: 'centered'
%
%   S.noverlap  Sample overlap (scalar).  If it appears in the argument list,
%               it must be the 4th scalar argument, but not necessarily the
%               4th argument in the list.  A string option could occur before,
%               instead of, or after this option.  If omitted, 'omitted' is
%               returned.  NOTE: Only pmusic expects S.NOverlap at present.
%
%   S.Range   - X-axis scaling method (string).  This option is specified in
%               the argument list as one of the mutually exclusive strings
%               specified in the list below.  Note that any case (including
%               mixed) will be accepted, including partial unique string matching.
%               This option may appear anywhere in the argument list.  The
%               structure field is set to the full name of the option string
%               as provided in the list below.  If omitted, 'omitted' is returned.
%               Recognized argument strings:
%                       Input             Field Contents
%                       'half'            'one-sided'
%                       'whole'           'two-sided'
%                       'one-sided'       'one-sided'
%                       'two-sided'       'two-sided'
%               NOTE: The hyphens may be omitted or replaced with a single space.
%
%   S.MTMethod -MTM method name (string).  This option is specified in
%               the argument list as one of the mutually exclusive strings
%               specified in the list below.  Note that any case (including
%               mixed) will be accepted, including partial unique string matching.
%               This option may appear anywhere in the argument list.  The
%               structure field is set to the full name of the option string
%               as provided in the list below.  If omitted, 'omitted' is returned.
%               Recognized argument strings: 'adapt', 'eigen', 'unity'.
%
%   S.CorrFlag -Method for pmusic (string).  This option is specified in
%               the argument list as one of the mutually exclusive strings
%               specified in the list below.  Note that any case (including
%               mixed) will be accepted, including partial unique string matching.
%               This option may appear anywhere in the argument list.  The
%               structure field is set to the full name of the option string
%               as provided in the list below.  If omitted, 'omitted' is returned.
%               Recognized argument strings: 'corr'
%
%  s.EVFlag   - Eigenvector method for pmusic (string).  This option is specified in
%               the argument list as one of the mutually exclusive strings
%               specified in the list below.  Note that any case (including
%               mixed) will be accepted, including partial unique string matching.
%               This option may appear anywhere in the argument list.  The
%               structure field is set to the full name of the option string
%               as provided in the list below.  If omitted, 'omitted' is returned.
%               Recognized argument strings: 'EV'


% Initialize return arguments
msg        = '';
msgobj     = [];
omit       = 'omitted';
s.NFFT     = omit;
s.Fs       = omit;
s.ConfInt  = omit;
s.Range    = omit;
s.MTMethod = omit;
s.nw       = omit;
s.noverlap  = omit;
s.CorrFlag = omit;
s.EVFlag   = omit;
s.DCFlag   = omit;
s.ConfLevel = omit;

% look for conflevel pv pair and set flag if user-specified
matchIdx = strcmpi('ConfidenceLevel',varargin);
numMatches = sum(matchIdx(:));
if numMatches>1
   msgobj = message('signal:psdoptions:MultipleValues');
   msg = getString(msgobj);
   return
elseif numMatches == 1
   pIdx = find(matchIdx);
   if pIdx == numel(varargin)
      msgobj = message('signal:psdoptions:MissingConfLevelValue');
      msg = getString(msgobj);
      return
   end
   % obtain the property value.
   confLevel = varargin{pIdx+1};
   if isscalar(confLevel) && isreal(confLevel) && 0<confLevel && confLevel<1
      s.ConfLevel = confLevel;
   else
      msgobj = message('signal:psdoptions:InvalidConfLevelValue');
      msg = getString(msgobj);
      return
   end
   % remove the pv-pair from the argument list
   varargin([pIdx pIdx+1]) = [];    
end

% look for flags

% Define all possible string options
% Lower case, no punctuation:
strOpts = {'half','onesided','whole','twosided', ...
           'adapt','unity','eigen', ...
           'corr', 'ev', ...
           'centered'};

% Check for mutually exclusive options
exclusiveOpts = strOpts([1 3; 2 4; 5 6; 5 7; 6 7; 1 10; 2 10]);
for i=1:size(exclusiveOpts,1)
   if any(strcmpi(exclusiveOpts(i,1),varargin)) && ...
         any(strcmpi(exclusiveOpts(i,2),varargin))
     msgobj = message('signal:psdoptions:ConflictingOptions', ...
                      exclusiveOpts{i,1}, exclusiveOpts{i,2});
     msg = getString(msgobj);
     return
   end
end

% No options passed - return all defaults:
if numel(varargin)<1, return; end

numReals = 0;  % Initialize count of numeric options

for argNum=1:numel(varargin),
   arg = varargin{argNum};
   
   isStr     = ischar(arg);
   isRealVector = isnumeric(arg) & ~issparse(arg) & isreal(arg) & any(size(arg)<=1);
   isValid   = isStr | isRealVector;
   
   if ~isValid,
      msgobj = message('signal:psdoptions:NeedValidOptionsType');
      msg = getString(msgobj);
      return;
   end
   
   if isStr,
      % String option
      %
      % Convert to lowercase and remove special chars:
      arg = RemoveSpecialChar(arg);
      
      % Match option set:
      i = find(strncmp(arg, strOpts, length(arg)));
      if isempty(i) || length(i)>1,
         msgobj = message('signal:psdoptions:UnknStringOption');
         msg = getString(msgobj);
         return
      end
      
      % Determine option set to which this string applies:
      switch i
      case {1,2,3,4},
         field = 'Range';
      case {5,6,7},
         field = 'MTMethod';
      case 8,
         field = 'CorrFlag';
      case 9,
         field = 'EVFlag';
      case 10,
         field = 'DCFlag';
      otherwise,
         error(message('signal:psdoptions:IdxOutOfBound'));
      end
      
      % Verify that no other options from (exclusive) set have
      % already been parsed:
      existingOpt = s.(field);
      if ~strcmp(existingOpt,'omitted'),
         msgobj = message('signal:psdoptions:MultipleValues');
         msg = getString(msgobj);
         return
      end
      
      % Map half to one-sided (index 1 -> index 2)
      % Map whole to two-sided (index 3 -> index 4)
      if i==1 || i==3,
         i=i+1;
      end
      
      % Set full option string into appropriate field:
      s.(field) = strOpts{i};
      
   else
      % Non-string options
      %
      numReals = numReals + 1;
      switch numReals
      case 1, s.NFFT    = arg;
      case 2,
         if length(arg)<=1,
            s.Fs = arg;
         else
            msgobj = message('signal:psdoptions:FsMustBeScalar');
            msg = getString(msgobj);
            return
         end
      case 3,
            s.ConfInt = arg;
            s.nw = arg;  % synonym
            
            % NOTE: Cannot error out for non-scalars, as NW may be a vector!
            %       We cannot disambiguate ConfInt from NW strictly from context.
            %
            %if length(arg)<=1,
            %   s.ConfInt = arg;
            %else
            %   msg = 'The confidence interval must be a scalar.';
            %   return
            %end
      case 4,
         if length(arg)<=1,
            s.noverlap = arg;
         else
            msgobj = message('signal:psdoptions:OverlapMustBeScalar');
            msg = getString(msgobj);
            return
         end
      otherwise
         msgobj = message('signal:psdoptions:TooManyNumericOptions');
         msg = getString(msgobj);
         return
      end
   end
   
end


%--------------------------------------------------------------------------------
function y=RemoveSpecialChar(x)
% RemoveSpecialChar
%   Remove one space of hyphen from 4th position, but
%   only if first 3 chars are 'one' or 'two'

y = lower(x);

% If string is less than 4 chars, do nothing:
if length(y)<=3, return; end

% If first 3 chars are not 'one' or 'two', do nothing:
if ~strncmp(y,'one',3) && ~strncmp(y,'two',3), return; end

% If 4th char is space or hyphen, remove it
if y(4)==' ' || y(4) == '-', y(4)=''; end

% [EOF] psdoptions.m
