function options = statset(varargin)
%STATSET Create/alter STATS options structure.
%   OPTIONS = STATSET('PARAM1',VALUE1,'PARAM2',VALUE2,...) creates a
%   statistics options structure OPTIONS in which the named parameters have
%   the specified values.  Any unspecified parameters are set to [].  When
%   you pass OPTIONS to a statistics function, a parameter set to []
%   indicates that the function uses its default value for that parameter.
%   Case is ignored for parameter names, and unique partial matches are
%   allowed.  NOTE: For parameters that are string-valued, the complete
%   string is required for the value; if an invalid string is provided, the
%   default is used.
%
%   OPTIONS = STATSET(OLDOPTS,'PARAM1',VALUE1,...) creates a copy of
%   OLDOPTS with the named parameters altered with the specified values.
%
%   OPTIONS = STATSET(OLDOPTS,NEWOPTS) combines an existing options
%   structure OLDOPTS with a new options structure NEWOPTS.  Any parameters
%   in NEWOPTS with non-empty values overwrite the corresponding old
%   parameters in OLDOPTS.
%
%   STATSET with no input arguments and no output arguments displays all
%   parameter names and their possible values, with defaults shown in {}
%   when the default is the same for all functions that use that option.
%   Use STATSET(STATSFUNCTION) (see below) to see function-specific
%   defaults for a specific function.
%
%   OPTIONS = STATSET (with no input arguments) creates an options
%   structure OPTIONS where all the fields are set to [].
%
%   OPTIONS = STATSET(STATSFUNCTION) creates an options structure with all
%   the parameter names and default values relevant to the statistics
%   function named in STATSFUNCTION.  STATSET sets parameters in OPTIONS to
%   [] for parameters that are not valid for STATSFUNCTION.  For example,
%   statset('factoran') or statset(@factoran) returns an options structure
%   containing all the parameter names and default values relevant to the
%   function 'factoran'.
%
%   STATSET parameters:
%      Display     - Level of display.  'off', 'iter', or 'final'.
%      MaxFunEvals - Maximum number of objective function evaluations
%                    allowed.  A positive integer.
%      MaxIter     - Maximum number of iterations allowed.  A positive integer.
%      TolBnd      - Parameter bound tolerance.  A positive scalar.
%      TolFun      - Termination tolerance for the objective function
%                    value.  A positive scalar.
%      TolTypeFun  - Flag to indicate how to use 'TolFun'. 'abs' indicates
%                    using 'TolFun' as an absolute tolerance; 'rel'
%                    indicates using it as a relative tolerance.
%      TolX        - Termination tolerance for the parameters.  A positive scalar.
%      TolTypeX    - Flag to indicate how to use 'TolX'. 'abs' indicates
%                    using 'TolX' as an absolute tolerance; 'rel'
%                    indicates using it as a relative tolerance.
%      GradObj     - Flag to indicate whether the objective function can return
%                    a gradient vector as a second output.  'off' or 'on'.
%      Jacobian    - Flag to indicate whether the model function can return a
%                    Jacobian as a second output.  'off' or 'on'.
%      DerivStep   - Relative difference used in finite difference derivative
%                    calculations.  A positive scalar, or a vector of
%                    positive scalars the same size as the vector of model
%                    parameters being estimated.
%      FunValCheck - Check for invalid values, such as NaN or Inf, from
%                    the objective function.  'off' or 'on'.
%      Robust      - Robust will be removed in a future release. Use
%                    RobustWgtFun to invoke the robust fitting option.
%      RobustWgtFun- A weight function for robust fitting. '', 'bisquare', 
%                    'andrews', 'cauchy', 'fair', 'huber', 'logistic',
%                    'talwar', or 'welsch'. Can also be a function handle that
%                    accepts a normalized residual as input and returns the
%                    robust weights as output. The default is [] to indicate
%                    no robust fitting. 
%      WgtFun      - WgtFun will be removed in the future release. Use
%                    RobustWgtFun instead.
%      Tune        - The tuning constant used in robust fitting to normalize the
%                    residuals before applying the weight function.  A positive
%                    scalar.  The default value depends upon the weight function.
%                    This parameter is required if the weight function is
%                    specified as a function handle.
%      UseParallel - Flag to indicate whether eligible functions should use
%                    capabilities of the Parallel Computing Toolbox (PCT),
%                    if the capabilities are available. Valid values are
%                    FALSE (the default), to indicate serial computation,
%                    and TRUE, to request parallel computation.
%                    See PARALLELSTATS for more complete information about
%                    this parameter.
%      UseSubstreams-Flag to indicate whether the random number generator
%                    in eligible functions should use the substream feature
%                    of the random number stream. FALSE (default) or
%                    TRUE. If TRUE, the function will use the
%                    substream feature during iterative loops in such a way
%                    as to generate reproducible random number streams in
%                    parallel and/or serial mode computation.
%                    See PARALLELSTATS for more complete information about
%                    this parameter.
%      Streams     - A single random number stream, or a cell array of random
%                    number streams.  The Streams option is accepted by
%                    some Statistics and Machine Learning Toolbox functions
%                    to govern what stream(s) to use when generating random
%                    numbers within the function. Acceptable values for the
%                    Streams parameter vary depending on the Statistics and
%                    Machine Learning Toolbox function and on other
%                    factors. See PARALLELSTATS for more complete
%                    information about this parameter.
%      OutputFcn   - Function handle specified using @, a cell array with
%                    function handles or an empty array (default).  The
%                    solver calls all output functions after each
%                    iteration.
%
%   See also STATGET.

%   Copyright 1993-2014 The MathWorks, Inc.

% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    fprintf('                Display: [ {off} | final | iter ]\n');
    fprintf('            MaxFunEvals: [ positive integer ]\n');
    fprintf('                MaxIter: [ positive integer ]\n');
    fprintf('                 TolBnd: [ positive scalar ]\n');
    fprintf('                 TolFun: [ positive scalar ]\n');
    fprintf('             TolTypeFun: [''abs'' |''rel'']\n')
    fprintf('                   TolX: [ positive scalar ]\n')
    fprintf('               TolTypeX: [''abs'' |''rel'']\n')
    fprintf('                GradObj: [ {off} | on ]\n')
    fprintf('               Jacobian: [ {off} | on ]\n')
    fprintf('              DerivStep: [ positive scalar or vector ]\n')
    fprintf('            FunValCheck: [ off | {on} ]\n')
    fprintf('           RobustWgtFun: [ {[]} | bisquare | andrews | cauchy | fair | huber | logistic | talwar | welsch | function handle ]\n')
    fprintf('                   Tune: [ positive scalar ]\n')
    fprintf('            UseParallel: [ {false} | true ]\n')
    fprintf('          UseSubstreams: [ {false} | true ]\n')
    fprintf('                Streams: [ {} | RandStream or cell array ]\n')
    fprintf('              OutputFcn: [ {[]} | function handle or cell array ]\n')
    fprintf('\n');
    return;
end

options = struct('Display', [], 'MaxFunEvals', [], 'MaxIter', [], ...
    'TolBnd', [], 'TolFun', [], 'TolTypeFun',[],'TolX', [], 'TolTypeX',[], ...
    'GradObj', [], 'Jacobian', [], 'DerivStep', [], 'FunValCheck', [], ...
    'Robust',[], 'RobustWgtFun', [], 'WgtFun',[], 'Tune',[], ...
    'UseParallel',[], 'UseSubstreams',[], 'Streams', {{}}, 'OutputFcn', []);

% If a function name/handle was passed in, then return the defaults.
if nargin == 1
    arg = varargin{1};
    if (ischar(arg) || isa(arg,'function_handle'))
        if isa(arg,'function_handle')
            arg = func2str(arg);
        end
        % Display is off by default.  The individual fitters have their own
        % warning/error messages that can be controlled via IDs.  The
        % optimizers print out text when display is on, but do not generate
        % warnings or errors per se.
        options.Display = 'off';
        switch lower(arg)
            case 'factoran' % this uses statsfminbx
                options.MaxFunEvals = 400;
                options.MaxIter = 100;
                options.TolFun = 1e-8;
                options.TolX = 1e-8;
            case {'normfit' 'lognfit' 'gamfit' 'bisafit' 'invgfit' 'logifit'...
                    'loglfit' 'nakafit' 'coxphfit'} % these use statsfminbx
                options.MaxFunEvals = 200;
                options.MaxIter = 100;
                options.TolBnd = 1e-6;
                options.TolFun = 1e-8;
                options.TolX = 1e-8;
            case {'evfit' 'wblfit'} % these use fzero (gamfit sometimes does too)
                options.TolX = 1e-6;
            case 'copulafit' % this uses fminbnd
                options.MaxFunEvals = 200;
                options.MaxIter = 100;
                options.TolX = 1e-6;
                options.TolBnd = 1e-6;
            case {'gpfit' 'gevfit' 'nbinfit' 'ricefit' 'tlsfit'} % these use fminsearch
                options.MaxFunEvals = 400;
                options.MaxIter = 200;
                options.TolBnd = 1e-6;
                options.TolFun = 1e-6;
                options.TolX = 1e-6;
            case 'kmeans'
                options.MaxIter = 100;
                options.UseParallel = false;
                options.UseSubstreams = false;
                options.Streams = {};
            case 'kmedoids'
                options.MaxIter = 100;
                options.UseParallel = false;
                options.UseSubstreams = false;
                options.Streams = {};
            case 'svmtrain'
                % do nothing because QP and SMO have different MaxIter values.
            case {'nlinfit','fitnlm'}
                options.MaxIter = 200;
                options.TolFun = 1e-8;
                options.TolX = 1e-8;
                options.DerivStep = eps^(1/3);
                options.FunValCheck = 'on';
                options.Robust = 'off';
                options.WgtFun = 'bisquare';
                options.Tune = []; % default varies by WgtFun, must be supplied for user-defined WgtFun
            case 'nlmefit'
                options.MaxIter = 200;
                options.TolFun = 1e-4;
                options.TolX = 1e-4;
                options.DerivStep = eps^(1/3);
                options.Jacobian = 'off';
                options.FunValCheck = 'on';
                options.OutputFcn = '';
            case 'nlmefitsa'
                options.DerivStep = eps^(1/3);
                options.FunValCheck = 'on';
                options.Streams = {};
                options.OutputFcn = {@nlmefitoutputfcn};
            case 'mlecustom' % this uses fminsearch, or maybe fmincon
                options.MaxFunEvals = 400;
                options.MaxIter = 200;
                options.TolBnd = 1e-6;
                options.TolFun = 1e-6;
                options.TolX = 1e-6;
                options.GradObj = 'off';
                options.DerivStep = eps^(1/3);
                options.FunValCheck = 'on';
            case 'mlecov'
                options.GradObj = 'off';
                options.DerivStep = eps^(1/4);
            case 'mdscale'
                options.MaxIter = 200;
                options.TolFun = 1e-6;
                options.TolX = 1e-6;
            case {'mvncdf' 'mvtcdf'}
                options.MaxFunEvals = 1e7;
                options.TolFun = []; % 1e-8 for dim < 4, 1e-4 otherwise
            case {'gmdistribution'}
                options.MaxIter = 100;
                options.TolFun = 1e-6;
            case {'nnmf'}
                options.MaxIter = 100;
                options.TolFun = 1e-4;
                options.TolX = 1e-4;
                options.UseParallel = false;
                options.UseSubstreams = false;
                options.Streams = {};
            case {'sequentialfs'}
                options.MaxIter = Inf;
                options.TolTypeFun = 'rel';
                options.UseParallel = false;
                options.UseSubstreams = false;
                options.Streams = {};
            case {'bootci' 'bootstrp' 'candexch' 'cordexch' 'crossval' 'daugment' ...
                    'dcovary' 'plsregress' 'rowexch' 'treebagger' 'lasso' 'lassoglm' ...
                    'parallel'}
                options.UseParallel = false;
                options.UseSubstreams = false;
                options.Streams = {};
            case {'jackknife'}
                options.UseParallel = false;
            case {'alsmf' 'pca'}
                options.MaxIter = 1e3;
                options.TolFun = 1e-6;
                options.TolX = 1e-6;
            case {'ppca'}
                options.MaxIter = 1e3;
                options.TolFun = 1e-6;
                options.TolX = 1e-6;
            case {'linearmixedmodel','lme','fitlme','fitlmematrix'}
                options.MaxIter = 10000;
                options.TolFun = 1e-6;
                options.TolX = 1e-12;
                options.Display = 'off';
            case {'generalizedlinearmixedmodel','fitglme'}
                options.MaxIter = 10000;
                options.TolFun = 1e-6;
                options.TolX = 1e-12;
                options.Display = 'off';                
            otherwise
                error(message('stats:statset:BadFunctionName', arg));
        end
        return
    end
end

names = fieldnames(options);
lowNames = lower(names);
numNames = numel(names);

% Process OLDOPTS and NEWOPTS, if it's there.
i = 1;
while i <= nargin
    arg = varargin{i};
    % Check if we're into the param name/value pairs yet.
    if ischar(arg), break; end
    
    if ~isempty(arg) % [] is a valid options argument
        if ~isa(arg,'struct')
            error(message('stats:statset:BadOptions', i));
        end
        argNames = fieldnames(arg);
        for j = 1:numNames
            name = names{j};
            if any(strcmp(name,argNames))
                val = arg.(name);
                if ~isempty(val)
                    if ischar(val)
                        val = lower(deblank(val));
                    end
                    val = checkparam(name,val);
                    options.(name) = val;
                end
            end
        end
    end
    i = i + 1;
end

% Done with OLDOPTS and NEWOPTS, now parse parameter name/value pairs.
if rem(nargin-i+1,2) ~= 0
    error(message('stats:statset:BadInput'));
end
expectval = false; % start expecting a name, not a value
while i <= nargin
    arg = varargin{i};
    
    if ~expectval
        % Process a parameter name.
        [~,j] = getParamVal(arg,lowNames,'parameter name');
        name = names{j};
        expectval = true; % expect a value next
    else
        % Process a parameter value.
        if ischar(arg)
            arg = lower(deblank(arg));
        end
        arg = checkparam(name,arg);
        options.(name) = arg;
        expectval = false; % expect a name next
    end
    i = i + 1;
end

% The default wgt function for robust fit is bisquare.
if (strcmp(options.Robust, 'on') && isempty(options.RobustWgtFun))
    if isempty(options.WgtFun)
        options.RobustWgtFun = 'bisquare';
    else
        options.RobustWgtFun = options.WgtFun;
    end
end

if ~isempty(options.RobustWgtFun)
    if ischar(options.RobustWgtFun)
        [~,options.Tune] = statrobustwfun(options.RobustWgtFun,[]);
    end
    if isempty(options.Tune) || ~isnumeric(options.Tune)
        error(message('stats:statset:BadWgtFun'));
    end
end

%-------------------------------------------------
function value = checkparam(name,value)
%CHECKPARAM Validate a STATSET parameter value.
%   VALUE=CHECKPARAM('name',VALUE) checks that the specified value VALUE is valid
%   for the parameter 'name'.
%   The return value allows coercion of multiple forms of specification
%   of a parameter to a single canonical form.

% Empty is always a valid parameter value.
if isempty(value)
    return
end

switch name
    case {'TolBnd','TolX'} % positive real scalar
        if ~isfloat(value) || ~isreal(value) || ~isscalar(value) || value <= 0
            m = message('stats:statset:BadTolerance',name);
            throwAsCaller(MException(m.Identifier,'%s',getString(m)));
        end
    case  {'TolFun'}
        if ~isfloat(value) || ~isreal(value) || ~isscalar(value) || value < 0
            m = message('stats:statset:BadTolerance',name);
            throwAsCaller(MException(m.Identifier,'%s',getString(m)));
        end
    case {'TolTypeFun', 'TolTypeX'}
        values = {'abs' 'rel'};
        if ~ischar(value) || ~any(strcmpi(value,values))
            m = message('stats:statset:BadTolType',name);
            throwAsCaller(MException(m.Identifier,'%s',getString(m)));
        end
    case {'Display'} % off,final,iter; we accept notify, but don't advertise it
        values = {'off' 'notify' 'final' 'iter'};
        if ~ischar(value) || ~any(strcmpi(value,values))
            m = message('stats:statset:BadDisplay',name);
            throwAsCaller(MException(m.Identifier,'%s',getString(m)));
        end
    case {'MaxIter' 'MaxFunEvals'} % non-negative integer, possibly inf
        if ~isfloat(value) || ~isreal(value) || ~isscalar(value) || any(value < 0)
            m = message('stats:statset:BadMaxValue',name);
            throwAsCaller(MException(m.Identifier,'%s',getString(m)));
        end
    case {'GradObj' 'Jacobian' 'FunValCheck' 'Robust'}
        values = {'off' 'on'};
        if ~ischar(value) || ~any(strcmpi(value,values))
            m = message('stats:statset:BadFlagValue',name);
            throwAsCaller(MException(m.Identifier,'%s',getString(m)));
        end
    case {'WgtFun'}
        values = {'andrews' 'bisquare' 'cauchy' 'fair' 'huber' 'logistic' 'talwar' 'welsch'};
        % allow custom wgt function
        if ~isa(value, 'function_handle')&&(~ischar(value) || ~any(strcmpi(value,values)))
            m = message('stats:statset:BadWeightFunction',name);
            throwAsCaller(MException(m.Identifier,'%s',getString(m)));
        end
    case {'RobustWgtFun'}
        values = {[] 'andrews' 'bisquare' 'cauchy' 'fair' 'huber' 'logistic' 'talwar' 'welsch'};
        % allow custom wgt function
        if ~isa(value, 'function_handle')&&(~ischar(value) || ~any(strcmpi(value,values)))
            m = message('stats:statset:BadWeightFunction',name);
            throwAsCaller(MException(m.Identifier,'%s',getString(m)));
        end
    case 'DerivStep'
        if ~isfloat(value) || ~isreal(value) || any(value <= 0) || ~isvector(value)
            m = message('stats:statset:BadDifference',name);
            throwAsCaller(MException(m.Identifier,'%s',getString(m)));
        end
    case 'Tune'
        if ~isfloat(value) || ~isreal(value) || any(value <= 0)
            m = message('stats:statset:BadTune',name);
            throwAsCaller(MException(m.Identifier,'%s',getString(m)));
        end
    case {'UseParallel' 'UseSubstreams'}
        % Prefer logical true/false (or scalar 1/0) 
        % but allow deprecated 'always'/'never'.
        deprecated_values = {'always' 'never'};
        if  ~( (isscalar(value) && ( islogical(value) || (isnumeric(value) && (value ==0 || value==1)) )) || ...
                (ischar(value) && any(strcmpi(value,deprecated_values))) );
            m = message('stats:statset:BadParallelParameter',name);
            throwAsCaller(MException(m.Identifier,'%s',getString(m)));
        end
        if isnumeric(value)
            value = logical(value);
        end
        if ischar(value)
            if strcmpi(value,'always')
                value = true;
            else
                value = false;
            end
        end
    case {'Streams'}
        if ~( isempty(value) || isa(value,'RandStream') || (iscell(value) && all(cellfun(@(x)isa(x,'RandStream'),value))) )
            m = message('stats:statset:BadStreams',name);
            throwAsCaller(MException(m.Identifier,'%s',getString(m)));
        end
    case {'OutputFcn'}
        if ~((iscell(value) && all(cellfun(@(x) isa(x,'function_handle'),value))) || isa(value,'function_handle'))
            m = message('stats:statset:BadOutputFcn',name);
            throwAsCaller(MException(m.Identifier,'%s',getString(m)));
        end
        
    otherwise
        m = message('stats:statset:BadParameter',name);
        throwAsCaller(MException(m.Identifier,'%s',getString(m)));
end
