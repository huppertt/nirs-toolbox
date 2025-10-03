classdef bootstrap_jobmanager
    % Job manager to run permutation and bootstrapping on pipelines

    properties
        pipeline = nirs.modules.default_modules.single_subject;  % The pipeline to run.  Can be single or multiple (cell) pipelines
        data_split_function = @(data)nirs.bootstrapping.data_splitting.randomize_event_timing(data);
        data_pooling_function = @(data)nirs.bootstrapping.data_pooling.pool_channelStats(data,'beta');

        max_iterations = 10000;
        use_parallel_computing=true;
        save_temp_file=true;
        save_temp_file_frequency=5000;
        temp_file_name='temp_bootstrap_jobmanager.mat';

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
                queue = createParallelProgressBar(obj.max_iterations);
            end
            if(obj.save_temp_file)
                numsaves=ceil(obj.max_iterations/obj.save_temp_file_frequency);
                iterrunning=1;
                for saveIter=1:numsaves
                    if(obj.use_parallel_computing)
                        iit=min(iterrunning+obj.save_temp_file_frequency-1,obj.max_iterations);
                        parfor iter=iterrunning:iit
                            send(queue,iter);
                            result{iter}=obj.run_iteration(data);
                        end
                    else
                        for iter=1:iit
                            disp(iter);
                            result{iter}=obj.run_iteration(data);
                        end
                    end
                    disp(['Iteration ' num2str(iit) ' of ' num2str(obj.max_iterations) ': Saving temp file as ' obj.temp_file_name])
                    save(obj.temp_file_name,"result","iterrunning",'-mat');
                    iterrunning=iit+1;
                end


            else

                if(obj.use_parallel_computing)
                    parfor iter=1:obj.max_iterations
                        send(queue,iter);
                        result{iter}=obj.run_iteration(data);
                    end
                else
                    for iter=1:obj.max_iterations
                        disp(iter);
                        result{iter}=obj.run_iteration(data);
                    end
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
                    
                    if(all(isnan(R(i,:))))
                        ResultBS.value_bins(:,i,j)=nan(size(R,2)+1,1);
                        ResultBS.ecdf(:,i,j)=nan(size(R,2)+1,1);
                    else
                        [ResultBS.ecdf(:,i,j),...
                            ResultBS.value_bins(:,i,j)]=ecdf(R(i,:));
                        lst=find(ResultBS.value_bins(:,i,j)>mean(R(i,:)));
                        ResultBS.ecdf(lst,i,j)=1-ResultBS.ecdf(lst,i,j);
                       

                    end
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

    end
end

