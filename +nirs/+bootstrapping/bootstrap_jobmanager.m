classdef bootstrap_jobmanager
    % Job manager to run permutation and bootstrapping on pipelines

    properties
        pipeline = nirs.modules.default_modules.single_subject;  % The pipeline to run.  Can be single or multiple (cell) pipelines
        data_split_function = @(data)nirs.bootstrapping.data_splitting.randomize_event_timing(data);
        data_pooling_function = @(data)nirs.bootstrapping.data_pooling.pool_channelStats(data,'beta');
    
        max_iterations = 100;
        use_parallel_computing=true;
    end

    methods
        function varargout = run(obj,data)
            result=cell(obj.max_iterations,1);
            if(iscell(obj.pipeline))

                truth=cell(length(obj.pipeline),1);
                for j=1:length(obj.pipeline)
                    for j=1:length(obj.pipeline)
                        tmpTruth = obj.pipeline{j}.run(data);
                        truth{j}=feval(obj.data_pooling_function,tmpTruth);
                    end
                end
            else
                tmpTruth = obj.pipeline.run(data);
                truth{1}=feval(obj.data_pooling_function,tmpTruth);

            end

            if(obj.use_parallel_computing)
                parfor iter=1:obj.max_iterations
                    result{iter}=run_iteration(obj,data);
                end
            else
                for iter=1:obj.max_iterations
                    result{iter}=run_iteration(obj,data);
                end
            end

            ResultBS=nirs.bootstrapping.bootstrap_result;
            ResultBS.truth=horzcat(truth{:});


            % Collapse the results
            for j=1:length(result{1})
                R=zeros(size(result{1}{j},1),length(result));
                for iter=1:length(result)
                    R(:,iter)=result{iter}{j};
                end
                
                for i=1:size(R,1)
                    [ResultBS.pcdf_estimate(:,i,j),...
                        ResultBS.value_bins(:,i,j),...
                        ResultBS.pcdf_lower_bounds(:,i,j),...
                        ResultBS.pcdf_upper_bounds(:,i,j)]=ecdf(R(i,:));
                end
                
            end
            
            if(nargout>0)
                tmpTruth.pvalue_fixed=ResultBS;
                varargout{1}=tmpTruth;
            end
            if(nargout>1)
                varargout{2}=ResultBS;
            end

        
        end
    
    end
end

function result=run_iteration(obj,data)
    datatmp = feval(obj.data_split_function,data);
    if(iscell(obj.pipeline))
        for j=1:length(obj.pipeline)
            tmp = obj.pipeline{j}.run(datatmp);
            result{j}=feval(obj.data_pooling_function,tmp);
        end
    else
        tmp = obj.pipeline.run(datatmp);
        result{1}=feval(obj.data_pooling_function,tmp);
    end

end