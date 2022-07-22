classdef Snirfwriter < nirs.realtime.modules.AbstractModule
    % real-time implementation of a snirf data file write.  This uses the
    % infinite data chunk notation to allow dynamic writing of the data
    % field as the data is being collected
    
    properties
       filename;
       overwrite=false;
    end
    
    properties(Hidden=true);
        fieldname='/nirs/data1';
        cnt=1;
    end
    
    methods
        function obj=Snirfwriter(prevJob)
            obj.name='RT-SNIRF data format writer';
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function obj=resetThis(obj)
        end
        
        function set.filename(obj,filename)
            obj.filename=filename;

            
        end



        function [d,t,probe,stimulus] = updateThis(obj,d,t,probe,stimulus)
            
            if(obj.cnt==1)
                if(~exist(obj.filename,'file') || obj.overwrite)
                    fid   = H5F.create(obj.filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
                    H5F.close(fid);
                    obj.fieldname='/nirs/data1';
                else
                    %TODO
                    obj.fieldname='/nirs/data1';
                end
                h5create(obj.filename,obj.fieldname,[size(d,2),Inf]);
            end
            h5write(obj.filename,obj.fieldname,d,[cnt 1],[size(d,1) size(d,2)]);
            cnt=cnt+size(d,1);

        end
            
    end
    
end