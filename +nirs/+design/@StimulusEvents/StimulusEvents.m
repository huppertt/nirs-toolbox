classdef StimulusEvents
    properties
        name
        onset
        dur
        amp
        regressor_no_interest=false;
        metadata=table;
        
    end

    properties( Dependent = true )
        count % number of stimuli
    end
    
    methods
        function obj = StimulusEvents( name, onset, dur, amp )
           if nargin > 0, obj.name  = name; end
           if nargin > 1, obj.onset = onset; end
           if nargin > 2, obj.dur   = dur; end
           if nargin > 3, obj.amp   = amp; end
           
        end

        function count = get.count( obj )
            %% count - returns the number of items in dictionary
            count = length(obj.onset);
        end
        
%         function disp(obj)
%             disp(['nirs.design.StimulusEvents [' num2str(length(obj.onset)) ' events']);
%         end
%         
        
        function vec = getStimVector( obj, time )
            assert( isvector( time ) )
            
            vec = zeros(size(time(:)));
            
            for i = 1:length( obj.onset )
                % list of points for this onset
                lst = time >= obj.onset(i) & time < obj.onset(i)+obj.dur(i);
                
                if(isempty(find(lst)))
                    % for really short durations this can happen because of
                    % the >= vs < issue
                    lst=min(find(time >= obj.onset(i)));
                end
                
                % set them to correct amplitude
                vec(lst) = vec(lst) + obj.amp(i);
                    
            end
            
        end
        
        function draw(obj)
            time=[0:.1:max(obj.onset+obj.dur)+30];
            vec=obj.getStimVector(time);
            figure;
            plot(time,vec);
            xlabel('time (sec)')
            ylabel(obj.name);
        end

    end

    methods (Hidden = true)


        % assignment, i.e. dict('hello') = 1234
        function obj = subsasgn(obj,s,b)
            
            numSubRef=length(s);
            if strcmp(s(1).type,'()')
                % Assignment of dictionary item
                
                if(numSubRef>1&&~strcmp(s(2).type,'()')) % assigning stim(1).amp vs stim.amp(1) should behave similarly
                    obj= builtin('subsasgn',obj,s([2,1,3:end]),b);
                elseif(isempty(b)&&numSubRef==1)
                    key = s(1).subs;
                    idx=key{1};
                    log_idx=true(1,obj.count);
                    log_idx(idx)=false;
                    s(1).subs={log_idx};
                    obj=subsref(obj,s);
                else
                   obj = builtin('subsasgn',obj,s,b); 
                end
            else
               	obj = builtin('subsasgn',obj,s,b);
            end
        end
        
        % retrieval;
        function [varargout] = subsref(obj,s)
            varargout=cell(1,1);
            out=cell(1,1);

            numSubRef=length(s);
            if(numSubRef>=1)
                switch s(1).type
                    case '()'
                        key = s(1).subs;
                        
                        idx=key{1};
                        
                        temp=obj;

                        if(length(obj)>1)
                            temp=temp(idx);
                        else
                            try
                                temp.onset=obj.onset(idx); 
                            catch ex
                                warning('Invalid index for stimuli');
                                throw(ex);
                            end
                                
                            if(length(obj.onset)==length(obj.dur))
                                temp.dur=obj.dur(idx);
                            end
                            
                            if(length(obj.onset)==length(obj.amp))
                                temp.amp=obj.amp(idx);
                            end
                        end

                        if(numSubRef>1)
                            temp=builtin('subsref',temp,s(2:end));
                        end
                        
                        out={temp};
                    case '.'
                        key = s.subs;
                        if isprop(obj,key)
                            out{:} = builtin('subsref',obj,s);
                        else
                            error('%s is not a property of StimulusEvents object',key);
                        end
                    otherwise
                        out{:} = builtin('subsref',obj,s);
                end
            else
                out{:} = builtin('subsref',obj,s);
            end


            if(nargout==0)
                varargout=out;
            else
                for k=1:max(nargout,length(out)) % At least assign the first output
                    varargout{k}=out{k};
                end
            end

            
        end

        function obj = plus(obj,obj2)
            if(strcmp(obj.name,obj2.name))
               % retain same name 
            else
                newName=sprintf('%s_%s',obj.name,obj2.name);
                fprintf('Merging "%s" and "%s" as "%s"\n',obj.name,obj2.name,newName);
                obj.name=newName;
            end

            obj.onset=[obj.onset,obj2.onset];
            obj.dur=[obj.dur,obj2.dur];
            obj.amp=[obj.amp,obj2.amp];

            % Resort in time
            [obj.onset,b_idx]=sort(obj.onset);
            obj.dur=obj.dur(b_idx);
            obj.amp=obj.amp(b_idx);

            if(obj.regressor_no_interest~=obj2.regressor_no_interest)
                warning('Mismatch in regressor_no_interest, using first value');
            end

            try
                obj.metadata=[obj.metadata;obj2.metadata];
            catch
                error('Unable to merge metadata');
            end
        
        end
        
    end
end