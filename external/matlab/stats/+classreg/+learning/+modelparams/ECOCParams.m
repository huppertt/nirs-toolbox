classdef ECOCParams < classreg.learning.modelparams.ModelParams
%ECOCParams Parameters for Error-Correcting Output Code model.
%
%   ECOCParams properties:
%       BinaryLearners         - Cell array of binary learners.
%       Coding                 - String for one of the built-in codings or
%                                custom matrix.
%       FitPosterior           - Logical flag for fitting posterior
%                                probabilities.
%       Options                - Options for parallel computation.
%       VerbosityLevel         - Verbosity level.
    
%   Copyright 2013-2014 The MathWorks, Inc.    
    
    properties
        BinaryLearners = '';
        Coding = [];
        FitPosterior = [];
        Options = [];
        VerbosityLevel = [];
    end

    methods(Access=protected)
        function this = ECOCParams(learners,coding,doposterior,paropts,verbose)
            this = this@classreg.learning.modelparams.ModelParams('ECOC','classification');
            this.BinaryLearners              = learners;
            this.Coding                      = coding;
            this.FitPosterior                = doposterior;
            this.Options                     = paropts;
            this.VerbosityLevel              = verbose;
        end
    end

    methods(Static,Hidden)
        function [holder,extraArgs] = make(type,varargin) %#ok<INUSL>
            args = {'learners' 'coding' 'fitposterior' 'options' 'verbose'};
            defs = {        ''       ''             []        []        []};
            [learners,M,doposterior,paropts,verbose,~,extraArgs] = ...
                internal.stats.parseArgs(args,defs,varargin{:});
            
            if ~isempty(learners)
                if ~ischar(learners) ...
                        && ~isa(learners,'classreg.learning.FitTemplate') ...
                        && ~iscell(learners)
                    error(message('stats:classreg:learning:modelparams:ECOCParams:make:BadLearners'));
                end
                if ischar(learners)
                    learners = classreg.learning.FitTemplate.make(learners,'type','classification');
                end
                if iscell(learners)
                    f = @(x) isa(x,'classreg.learning.FitTemplate');
                    isgood = cellfun(f,learners);
                    if ~all(isgood)
                        error(message('stats:classreg:learning:modelparams:ECOCParams:make:BadCellArrayLearners'));
                    end
                end
            end
            
            if ~isempty(M)
                if ischar(M)                    
                    allowedVals = {'onevsone' 'allpairs' 'onevsall' 'binarycomplete' 'ternarycomplete' ...
                        'ordinal' 'sparserandom' 'denserandom'};
                    tf = strncmpi(M,allowedVals,length(M));
                    Nfound = sum(tf);
                    if Nfound~=1
                        error(message('stats:classreg:learning:modelparams:ECOCParams:make:BadCodingName'));
                    end
                    M = allowedVals{tf};
                    if strcmp(M,'allpairs')
                        M = 'onevsone';
                    end
                else
                    if ~isnumeric(M) || ~ismatrix(M)
                        error(message('stats:classreg:learning:modelparams:ECOCParams:make:BadCodingType'));
                    end
                    
                    vals = unique(M(:));
                    if numel(vals)>3 || any(vals~=-1 & vals~=0 & vals~=1)
                        error(message('stats:classreg:learning:modelparams:ECOCParams:make:BadCodingElements'));
                    end
                    
                    L = size(M,2);                    
                    for l=1:L
                        if ~any(M(:,l)==-1)
                            error(message('stats:classreg:learning:modelparams:ECOCParams:make:CodingColumnWithoutNegOne',l));
                        end
                        if ~any(M(:,l)==1)
                            error(message('stats:classreg:learning:modelparams:ECOCParams:make:CodingColumnWithoutPosOne',l));
                        end
                    end
                    
                    for i=1:L-1
                        for j=i+1:L
                            if all(M(:,i)==M(:,j)) || all(M(:,i)==-M(:,j))
                                error(message('stats:classreg:learning:modelparams:ECOCParams:make:CodingHasIdenticalColumns',i,j));
                            end
                        end
                    end
                    
                    [tf,i,j] = isconnected(M);
                    if ~tf
                        fprintf('%s\n',getString(message('stats:classreg:learning:modelparams:ECOCParams:make:CodingHasInseparableRows')));
                        for n=1:numel(i)
                            fprintf('     %5i     %5i\n',i(n),j(n));
                        end
                        error(message('stats:classreg:learning:modelparams:ECOCParams:make:DisconnectedCoding'));
                    end
                end
            end
            
            if ~isempty(doposterior)
                doposterior = internal.stats.parseOnOff(doposterior,'FitPosterior');
            end
            
            if ~isempty(paropts) && ~isstruct(paropts)
                error(message('stats:classreg:learning:modelparams:ECOCParams:make:BadOptions'));
            end
            
            if ~isempty(verbose) && ...
                    (~isscalar(verbose) || ~isfloat(verbose) ...
                    || verbose<0 || isnan(verbose))
                error(message('stats:classreg:learning:modelparams:ECOCParams:make:BadVerbose'));
            end
            
            holder = classreg.learning.modelparams.ECOCParams(...
                learners,M,doposterior,paropts,verbose);
        end
    end

    methods(Hidden)
        function this = fillDefaultParams(this,X,Y,W,dataSummary,classSummary) %#ok<INUSL>

            %
            % Fill binary learners
            %
            if isempty(this.BinaryLearners)
                templates = templateSVM;
            else
                templates = this.BinaryLearners;
            end
            
            if isscalar(templates) && ~iscell(templates)
                templates = {templates};
            end
            
            for l=1:numel(templates)
                learner = templates{l};
                
                learner = fillIfNeeded(learner,'classification');
                
                learner = setBaseArg(learner,'predictornames',dataSummary.PredictorNames);
                learner = setBaseArg(learner,'categoricalpredictors',dataSummary.CategoricalPredictors);
                learner = setBaseArg(learner,'responsename',dataSummary.ResponseName);

                switch learner.Method
                    case 'Discriminant'
                        if isempty(learner.ModelParams.FillCoeffs)
                            learner.ModelParams.FillCoeffs = false;
                        end
                        if isempty(learner.ModelParams.DiscrimType)
                            learner.ModelParams.DiscrimType = 'pseudoLinear';
                        end
                        
                    case 'SVM'
                        if isempty(learner.ModelParams.SaveSupportVectors) && ...
                                (isempty(learner.ModelParams.KernelFunction) || ...
                                strcmp(learner.ModelParams.KernelFunction,'linear'))
                            learner.ModelParams.SaveSupportVectors = false;
                        end
                end
                
                templates{l} = learner;
            end
            
            if isscalar(templates)
                templates = templates{1};
            end
            
            this.BinaryLearners = templates;
            
            %
            % Fill coding
            %
            if     isempty(this.Coding)
                this.Coding = 'onevsone';
            elseif isnumeric(this.Coding)
                M = this.Coding;
                
                K = numel(classSummary.ClassNames);
                if size(M,1)~=K
                    error(message('stats:classreg:learning:modelparams:ECOCParams:fillDefaultParams:BadNumRowsInM',K));
                end
                
                zeroclass = sum(M~=0,2)==0;
                if any(zeroclass)
                    k = find(zeroclass,1,'first');
                    s = cellstr(classSummary.ClassNames(k));
                    error(message('stats:classreg:learning:modelparams:ECOCParams:fillDefaultParams:ZeroRowInM',s{:},k));
                end
                
                for i=1:K-1
                    for j=i+1:K
                        if all(M(i,:)==M(j,:))
                            s1 = cellstr(classSummary.ClassNames(i));
                            s2 = cellstr(classSummary.ClassNames(j));
                            error(message('stats:classreg:learning:modelparams:ECOCParams:fillDefaultParams:IdenticalRowsInM',i,j,s1{:},s2{:}));
                        end
                    end
                end
            end
            
            if isempty(this.FitPosterior)
                this.FitPosterior = false;
            end
            
            if isempty(this.Options)
                this.Options = statset('parallel');
            end
            
            if isempty(this.VerbosityLevel)
                this.VerbosityLevel = 0;
            end
        end
    end
    
end


% Test coding matrix M for connectivity. Loop over classes. For each class,
% find all classes connected to this one either directly or indirectly. Two
% classes (rows of M) are connected directly if they have opposite-sign
% values in the same column of M. Two classes are connected indirectly if
% they are both connected to a 3rd class.
function [tf,i,j] = isconnected(M)

K = size(M,1); % number of classes

F = eye(K)==1; % Connectivity matrix. F(I,J)=1 if class I is connected to class J.

for k=1:K
    lookAt = k; % list of classes to inspect
    lookedAt = []; % list of inspected classes
    
    Mreduced = M; % coding matrix with inspected columns dropped
    
    % This loop finds classes directly connected to class lookAt(1) in each
    % iteration. The loop continues until no direct connection can be
    % found. The final list includes all classes directly and indirectly
    % connected to class K.
    while ~isempty(lookAt)
        % Get a list of classes directly connected to class lookAt(1) and
        % the appropriately reduced coding matrix.
        [newk,Mreduced] = findConnectedClasses(Mreduced,lookAt(1));
        
        F(k,newk) = true; % update connectivity matrix
        
        % Update lists of inspected classes and classes to inspect
        lookedAt = unique([lookedAt lookAt(1)]);
        lookAt(1) = [];
        
        tf = ismember(newk,lookedAt);
        newk(tf) = [];
        
        lookAt = unique([lookAt newk]);
    end
end

tf = all(all(F|F')); % Are all classes connected?

if nargout>1
    if tf
        i = [];
        j = [];
    else % Find indices for disconnected classes.
        idx = find(~(F|F'));
        [i,j] = ind2sub([K K],idx);
        keep = i<j;
        i = i(keep);
        j = j(keep);
    end
end

end


% Find list of classes KLIST directly connected to class K in the coding
% matrix M. Class J is directly connected to class K if there exists a
% column of M such that Kth and Jth elements have opposite sign. Classes
% that are directly connected to K are then dropped from the coding matrix
% so we would not have to re-analyze them.
function [klist,M] = findConnectedClasses(M,k)
klist = [];

idxnonzero = find(M(k,:));

for n=1:numel(idxnonzero)
    idxcol = idxnonzero(n);
    elem = M(k,idxcol);
    newk = find( M(:,idxcol) == -elem );
    klist = unique([klist newk(:)']);
end

M(:,idxnonzero) = [];

end
