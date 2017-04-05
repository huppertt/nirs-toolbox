classdef Run_HOMER2 < nirs.modules.AbstractModule
    %% This is a generic wrapper to run a HOMER-2 job
    %
    % Options:
    %     fcn - the homer function to run (needs to be file name to run)
    %           - should autopopulate the fields upin invoke of the fcn set
    %           command
    
    % Known issues:  
    %    aux variables are not held in my toolbox at the moment.  
    %    none of the image recon HOMER2 code is handled
    %    My toolbox does not use time masks (e.g. tInc)
    %    This function will not change data classes
    %    This function does not update the probe.link fields
    
    properties
        fcn=[];
        inputs={};
        outputs={};
        vars=[];
    end
    
    methods
        function obj = Run_HOMER2( prevJob )
            obj.name = 'Run generic HOMER2 function';
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function obj = set.fcn(obj,fcn)
            [I,O,vars] = getHOMER2inputs(fcn);
            obj.inputs=I;
            obj.outputs=O;
            obj.vars=vars;
            obj.fcn=fcn;
        end
        
        
        function data = runThis( obj, data )
            
            for i = 1:length(data)
                varargin={};
                for j=1:size(obj.inputs,1)
                    varargin{j}=obj.inputs{j,2}(obj,data(i));
                end
                
                vargout=cell(size(obj.outputs,1),1);
                [vargout{:}]=feval(obj.fcn,varargin{:});
                
                for j=1:size(obj.outputs,1)
                    data(i)=obj.outputs{j,2}(data(i),vargout{j});
                end
                
            end
        end
    end
    
end


function [I,O,vars] = getHOMER2inputs(fcn);

fcn=which(fcn);
if(isempty(fcn))
    error('file not found: check that the file is on the MatlabPath');
end

fid=fopen(fcn,'r');
line=fgetl(fid);
fclose(fid);

line(1:min(find(double(line)~=double('%')))-1)=[];

if(~isempty(strfind(line,'function')))
    line(1:strfind(line,'function')+length('function'))=[];
end
    
outputs=line(1:strfind(line,'=')-1);
inputs=line(strfind(line,'(')+1:strfind(line,')')-1);

outputs(strfind(outputs,' '))=[];
inputs(strfind(inputs,' '))=[];
outputs(strfind(outputs,'['))=[];
outputs(strfind(outputs,']'))=[];

outputs(end+1)=',';
inputs(end+1)=',';

I={};
while(1)
    lst=min(strfind(inputs,','))-1;
    if(isempty(lst))
        break
    end
    I{end+1,1}=inputs(1:lst);
    inputs(1:lst+1)=[];
end


O={};
while(1)
    lst=min(strfind(outputs,','))-1;
    if(isempty(lst))
        break
    end
    O{end+1,1}=outputs(1:lst);
    outputs(1:lst+1)=[];
end

% Some of the key words can translate to my object classes.  However, HOMER
% seems to use a variety of input/output names with little consistency, so
% I have only the ones that exist as of 3/5/2017

% any term that is not identified needs to be specified in the "job.vars"
% struct

keywordsIn={'d', @(obj,data)data.data;...
    'd2', @(obj,data)data.data;...
    'y', @(obj,data)data.data;...
    'yo', @(obj,data)data.data;...
    'y2', @(obj,data)data.data;...
    'dod', @(obj,data)data.data;...
    'dod2', @(obj,data)data.data;...
    'dc',@(obj,data)data.data;...
    'dc2',@(obj,data)data.data;...
    'fs',@(obj,data)data.Fs;...
    't',@(obj,data)data.time;...
    'time',@(obj,data)data.time;...
    'SD',@(obj,data)nirs.util.probe2sd(data.probe);...
    'ml',@(obj,data)getfield(nirs.util.probe2sd(data.probe),'MeasList');...
    's',@(obj,data)localgetstim(obj,data)};

keywordsOut={'d', @(data,d)setfield(data,'data',d);...
    'd2',@(data,d)setfield(data,'data',d);...
    'y', @(data,d)setfield(data,'data',d);...
    'y2', @(data,d)setfield(data,'data',d);...
    'dc',@(data,d)setfield(data,'data',d);...
    'dc2',@(data,d)setfield(data,'data',d);...
    'dod',@(data,d)setfield(data,'data',d);...
    'dcCbsi',@(data,d)setfield(data,'data',d);...
    'dodSpline',@(data,d)setfield(data,'data',d);...
    'dodWavelet',@(data,d)setfield(data,'data',d);...
    'yc',@(data,d)setfield(data,'data',d);...
    'yhrf',@(data,d)setfield(data,'data',d);...
    'dN',@(data,d)setfield(data,'data',d);...
    'yavg',@(data,d)setfield(data,'data',d);...
    'tHRF',@(data,t)setfield(data,'time',t);...
    'dodWavelet',@(data,d)setfield(data,'data',d);...
    'dod2',@(data,d)setfield(data,'data',d);...
    'SD', @(data,SD)setfield(data,'probe',nirs.util.sd2probe(SD))};


vars=struct;

for i=1:length(I)
    [found,a]=ismember(I{i,1},{keywordsIn{:,1}});
    if(found)
        I{i,2}=keywordsIn{a,2};
    else
        vars=setfield(vars,I{i,1},[]);
        I{i,2}=@(obj,data)getfield(obj.vars,I{i,1});
    end
end


for i=1:length(O)
    [found,a]=ismember(O{i,1},{keywordsOut{:,1}});
    if(found)
        O{i,2}=keywordsOut{a,2};
    else
        O{i,2}=@(data,tmp)data;
    end
end

end


function X = localgetstim(obj,data)

b=nirs.design.basis.FIR;
b.binwidth=1;
b.nbins=1;
b.isIRF=false;

basis=Dictionary;
basis('default')=b;
[X, names] = nirs.design.createDesignMatrix( data.stimulus,data.time, basis);
end
