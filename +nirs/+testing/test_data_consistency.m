function stats = test_data_consistency(data,pipeline,niter)
%% test_data_consistency - Tests the reliability/consistency of a data
% analysis pipeline by repeatedly splitting the data in half (randomly
% subdividing the events of each stimulus condition into two equal-count
% groups), rerunning the pipeline on each half, and comparing the results
% to the pipeline run on the full data.
%
% Args:
%     data     - the nirs.core.Data (or similar) array to test
%     pipeline - an analysis pipeline (chain of nirs.modules objects) that
%                must include a nirs.modules.MixedEffects (or subclass)
%                module.  The include_diagnostics flag of that module will
%                be forced to true (a warning is issued if it was false).
%                An error is thrown if no mixed-effects module is found in
%                the pipeline.
%     niter    - (optional) number of random split iterations to run
%                (default = 10)
%
% Returns:
%     stats - a table summarizing the consistency (reliability) of the two
%             data halves for each channel/condition, including a
%             partialed F-statistic and p-value comparing the fit of the
%             full model to the fit obtained when the model is refit
%             separately to each half, along with the intraclass
%             correlation (ICC) between the two halves.  The table is also
%             displayed indicating whether the two halves of the data
%             appear to be statistically consistent.

if(nargin<3 || isempty(niter))
    niter=10;
end


%% Run the pipeline on the full (unsplit) data
sFull = pipeline.run(data);

tblFull = sFull.table;

%% Iterate over random splits
ResultsFull=[];
ResultsSplitA=[];
ResultsSplitB=[];

ResultsTFull=[];
ResultsTSplitA=[];
ResultsTSplitB=[];

for iter=1:niter
    
    dataS = nirs.testing.split_data_conditions(data);
    
    sSplit = pipeline.run(dataS);
    
    tblSplit = sSplit.table;
    
    % match each condition in the full model to its two halves in the
    % split model
    condFull = tblFull.cond;
        
    condsAll = unique(condFull);


    for cidx=1:length(condsAll)
        
        
        condA = [condsAll{cidx} '_halfA'];
        condB = [condsAll{cidx} '_halfB'];
        
        lstF = find(strcmp(condFull,condsAll{cidx}));
        lstA = find(strcmp(tblSplit.cond,condA));
        lstB = find(strcmp(tblSplit.cond,condB));
        
        if(isempty(lstA) || isempty(lstB))
            continue
        end
       
       
        
        %% Intraclass correlation between the beta estimates of the two
        % halves (a simple two-observation ICC)
        betaF = tblFull(lstF,:).beta;
        betaA = tblSplit(lstA,:).beta;
        betaB = tblSplit(lstB,:).beta;

        tstatF = tblFull(lstF,:).tstat;
        tstatA = tblSplit(lstA,:).tstat;
        tstatB = tblSplit(lstB,:).tstat;

       
        ResultsFull=[ResultsFull; betaF];
        ResultsSplitA=[ResultsSplitA; betaA];
        ResultsSplitB=[ResultsSplitB; betaB];

        ResultsTFull=[ResultsTFull; tstatF];
        ResultsTSplitA=[ResultsTSplitA; tstatA];
        ResultsTSplitB=[ResultsTSplitB; tstatB];

    end
   
end

[r,p]=corrcoef((ResultsTSplitA+ResultsTSplitA)/2,ResultsTFull);
disp(['ICC (Full vs Split): ' num2str(r(1,2)) ' (p<' num2str(p(1,2)) ')']);

[r,p]=corrcoef(ResultsTSplitA,ResultsTFull);
disp(['ICC (Full vs SplitA): ' num2str(r(1,2)) ' (p<' num2str(p(1,2)) ')']);

[r,p]=corrcoef(ResultsTSplitB,ResultsTFull);
disp(['ICC (Full vs SplitB): ' num2str(r(1,2)) ' (p<' num2str(p(1,2)) ')']);

[r,p]=corrcoef(ResultsTSplitB,ResultsTSplitA);
disp(['ICC (SplitA vs SplitB): ' num2str(r(1,2)) ' (p<' num2str(p(1,2)) ')']);