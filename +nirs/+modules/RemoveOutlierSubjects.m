classdef RemoveOutlierSubjects < nirs.modules.AbstractModule
%% RemoveOutlierSubjects - Flags and removes outlier subjects based on leverage at the group level
% 
% Options:
%     formula - formula intended to use for group level model
%     allow_partial_removal- if false, then all files belonging to a subject will be removed if ANY file is bad.
%     cutoff - p-val cuttoff
    
    properties
        formula = 'beta ~ -1 + cond';
        allow_partial_removal=true;  % if false, then all files belonging to a subject will be removed if ANY file is bad.
        cutoff=0.05;
    end
    
    methods
        function obj = RemoveOutlierSubjects( prevJob )
           obj.name = 'Remove outlier subjects';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
           tbl=nirs.util.grouplevelleveragestats(data,obj.formula);
           
           lst=find(tbl.pval_subject<obj.cutoff);
           if(length(lst)==0)
               disp('no subjects removed');
               return
           end
           
           lst=unique(tbl.file_idx(lst));
           if(obj.allow_partial_removal)
               disp(['Removing ' num2str(length(lst)) ' entries']);
               disp([table(lst,'VariableNames',{'FileIndex'}) nirs.createDemographicsTable(data(lst))]);
               data(lst)=[];
           else
               demo=nirs.createDemographicsTable(data); 
               lst=find(ismember(demo,unique(demo(lst,:))));
               disp(['Removing ' num2str(length(lst)) ' entries']);
               disp([table(lst,'VariableNames',{'FileIndex'}) nirs.createDemographicsTable(data(lst))]);
               data(lst)=[];
               
           end
           
           
            
        end
    end
    
end

