classdef SubjLevelConnectivityStats < nirs.modules.AbstractModule
    % SubjLevelConnectivityStats
    %
    % Combines multiple nirs.core.sFCStats objects for each participant
    % into one subject-level connectivity-statistics object.
    %
    % Example:
    %   job = nirs.modules.SubjLevelConnectivityStats();
    %   job.sortfield = 'Subject';
    %   job.formula = 'R ~ -1 + cond';
    %   job.robust = true;
    %   SubjStats = job.run(Stats);

    properties
        sortfield = 'Subject';
        formula = 'R ~ -1 + cond';
        robust = false;
        include_diagnostics = false;
    end

    methods

        function obj = SubjLevelConnectivityStats(prevJob)

            obj.name = 'Subject-Level Connectivity Effects Model';

            if nargin > 0
                obj.prevJob = prevJob;
            end
        end

        function G = runThis(obj, S)

            %% Validate input

            if isempty(S)
                G = S;
                return
            end

            if ~isa(S(1), 'nirs.core.sFCStats')
                error('SubjLevelConnectivityStats:InvalidInput', ...
                    ['Input must contain nirs.core.sFCStats objects. ' ...
                     'The supplied input contains %s objects.'], ...
                    class(S(1)));
            end

            %% Obtain demographics

            demo = nirs.createDemographicsTable(S);

            if isempty(demo)
                warning('SubjLevelConnectivityStats:NoDemographics', ...
                    ['No demographics were found. The input will be ' ...
                     'returned without subject-level averaging.']);
                G = S;
                return
            end

            variableNames = demo.Properties.VariableNames;

            % Find the grouping variable without requiring matching case.
            fieldIndex = find( ...
                strcmpi(variableNames, obj.sortfield), 1);

            if isempty(fieldIndex)
                warning('SubjLevelConnectivityStats:MissingSortField', ...
                    ['The demographic variable "%s" was not found. ' ...
                     'The input will be returned unchanged.'], ...
                    obj.sortfield);

                G = S;
                return
            end

            sortfieldActual = variableNames{fieldIndex};

            % Convert participant identifiers to strings so numeric,
            % categorical, character, and cell identifiers are supported.
            subjectValues = string(demo.(sortfieldActual));
            subjects = unique(subjectValues, 'stable');

            numberOfSubjects = numel(subjects);

            %% Generate one connectivity object per participant

            G(numberOfSubjects,1) = nirs.core.sFCStats();

            for subjectIndex = 1:numberOfSubjects

                currentSubject = subjects(subjectIndex);

                subjectRows = find( ...
                    subjectValues == currentSubject);

                subjectStats = S(subjectRows);

                fprintf('Running subject %s (%d files)\n', ...
                    currentSubject, numel(subjectRows));

                % If only one connectivity result exists for the
                % participant, no within-subject model is necessary.
                if numel(subjectStats) == 1

                    subjectResult = subjectStats;

                else
                    job = nirs.modules.MixedEffectsConnectivity();

                    job.formula = obj.formula;
                    job.robust = obj.robust;
                    job.include_diagnostics = ...
                        obj.include_diagnostics;

                    subjectResult = job.run(subjectStats);
                end

                % Preserve the participant's demographics.
                subjectResult.demographics = ...
                    subjectStats(1).demographics;

                subjectResult.description = sprintf( ...
                    'Subject-Level Connectivity: %s', ...
                    currentSubject);

                G(subjectIndex,1) = subjectResult;
            end

            fprintf('Completed %d subjects.\n', numberOfSubjects);
        end
    end
end