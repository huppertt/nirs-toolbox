classdef synthetic_measurement < nirs.modules.AbstractModule
%% Combines ChannelStats data into a common space
% 
    properties
        commonprobe = 'combine'; 
    end
    methods
        function obj = synthetic_measurement( prevJob )
           obj.name = 'synthetic measurement';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            if(~isa(obj.commonprobe,'nirs.core.Probe') && ...
                    ~isa(obj.commonprobe,'nirs.core.Probe1020'))
                probe=combineprobes(data);
            else
                probe=obj.commonprobe;
            end
            for i = 1:length(data)
                data(i)=synthetic_meas(data(i),probe);
            end
        end
    end
    
end


function probe = combineprobes(data)





ndet=0;
nsrc=0;
link=[];
optodes=[];
for i=1:length(data)
    thislink=data(i).probe.link;
    thislink.source=thislink.source+nsrc;
    thislink.detector=thislink.detector+ndet;
    
    thisoptodes=data(i).probe.optodes;
    for j=1:size(data(i).probe.detPos,1)
        dI=['0000' num2str(j)];
        dI=dI(end-3:end);
        dI2=['0000' num2str(j+ndet)];
        dI2=dI2(end-3:end);
        
        lst=find(ismember(thisoptodes.Name,['Detector-' dI]));
        thisoptodes.Name{lst}=['Detector-' dI2];
    end
    
    for j=1:size(data(i).probe.srcPos,1)
        dI=['0000' num2str(j)];
        dI=dI(end-3:end);
        dI2=['0000' num2str(j+nsrc)];
        dI2=dI2(end-3:end);
        
        lst=find(ismember(thisoptodes.Name,['Source-' dI]));
        thisoptodes.Name{lst}=['Source-' dI2];
    end
    link=[link; thislink];
    optodes=[optodes; thisoptodes];
    ndet=ndet+size(data(i).probe.detPos,1);
    nsrc=nsrc+size(data(i).probe.srcPos,1);
    
end
  
    
probe=data(1).probe;
probe.optodes=optodes;
probe.link=link;

if(isa(data(1).probe,'nirs.core.Probe1020'))
    warning('this will not work if the underlying probe mesh differ');
    
    optodes=[];
    for i=1:length(data)
        thislink=data(i).probe.link;
        thislink.source=thislink.source+nsrc;
        thislink.detector=thislink.detector+ndet;
        
        thisoptodes=data(i).probe.swap_reg.optodes;
        for j=1:size(data(i).probe.swap_reg.detPos,1)
            dI=['0000' num2str(j)];
            dI=dI(end-3:end);
            dI2=['0000' num2str(j+ndet)];
            dI2=dI2(end-3:end);
            
            lst=find(ismember(thisoptodes.Name,['Detector-' dI]));
            thisoptodes.Name{lst}=['Detector-' dI2];
        end
        
        for j=1:size(data(i).probe.swap_reg.srcPos,1)
            dI=['0000' num2str(j)];
            dI=dI(end-3:end);
            dI2=['0000' num2str(j+nsrc)];
            dI2=dI2(end-3:end);
            
            lst=find(ismember(thisoptodes.Name,['Source-' dI]));
            thisoptodes.Name{lst}=['Source-' dI2];
        end
        link=[link; thislink];
        optodes=[optodes; thisoptodes];
        ndet=ndet+size(data(i).probe.detPos,1);
        nsrc=nsrc+size(data(i).probe.srcPos,1);
        
    end
    probe.optodes_registered=optodes;    
end


end


function ChanStatsNew = synthetic_meas(ChanStats,NewProbe,FwdModel)
% This function converts the ChannelStats info from the imput into the
% space of the NewProbe.  If the forward model is not provided, the
% semi-infinte slab model with default parameters is used.
%
%
%% Example:
% data=nirs.testing.simData;
% job=nirs.modules.Resample();
% job=nirs.modules.OpticalDensity(job);
% job=nirs.modules.AR_IRLS(job);
% SubjStats = job.run(data);
% newprobe=data.probe;
% newprobe.optodes.X=newprobe.optodes.X+15*randn(height(newprobe.optodes),1);
% newprobe.optodes.Y=newprobe.optodes.Y+15*randn(height(newprobe.optodes),1);
% NewSubjStats = synthetic_meas(SubjStats,newprobe);

%% Construct the forward model for the old and new probe

if(nargin<3)
    
    if(~isa(NewProbe,'nirs.core.Probe1020'))
    % No forward model provided, then build the slab version
    minX = min(min(ChanStats.probe.optodes.X),min(NewProbe.optodes.X));
    maxX = max(max(ChanStats.probe.optodes.X),max(NewProbe.optodes.X));
    dX = (maxX-minX)/3;
    minY = min(min(ChanStats.probe.optodes.Y),min(NewProbe.optodes.Y));
    maxY = max(max(ChanStats.probe.optodes.Y),max(NewProbe.optodes.Y));
    dY = (maxY-minY)/3;
    
    [X,Y,Z]=meshgrid([minX-dX:dX/30:maxX+dX],[minY-dY:dY/30:maxY+dY],[-5 -10 -15]);
    
    mesh=nirs.core.Mesh;
    mesh.nodes=[X(:) Y(:) Z(:)];
    else
        mesh=NewProbe.getmesh;
        %mesh=mesh(end);
    end
    lambda=unique(NewProbe.link.type);
    
    FwdModel=nirs.forward.ApproxSlab;
    FwdModel.mesh=mesh;
    FwdModel.prop=nirs.media.tissues.brain(.7,50,lambda);
    FwdModel.Fm=0;
end

% Do some FwdModel checking here


ChanStats=sorted(ChanStats,{'type','source','detector','cond'});
NewProbe.link=sortrows(NewProbe.link,{'type','source','detector'});

 if(~isa(ChanStats.probe,'nirs.core.Probe1020'))
     FwdModel.probe=ChanStats.probe;
 else
     FwdModel.probe=ChanStats.probe.swap_reg;
 end
 
[Lold]=FwdModel.jacobian;
Lold=Lold.mua;

 if(~isa(NewProbe,'nirs.core.Probe1020'))
     FwdModel.probe=NewProbe;
 else
     FwdModel.probe=NewProbe.swap_reg;
 end
[Lnew]=FwdModel.jacobian;
Lnew=Lnew.mua;

%[Unew,Snew,Vnew]=nirs.math.mysvd(Lnew);
% [Uold,Sold,Vold]=nirs.math.mysvd(Lold);
% s=diag(Sold);
% lst=find(s>s(1)/1E4);
% S=diag(s(lst));
% H = Lnew*Vold(:,lst)*S*Uold(:,lst)';

[c,~,i]=unique(ChanStats.probe.link.type);
H=[];
for j=1:length(unique(i))
    lst=find(i==j);
    lst2=find(ismember(NewProbe.link.type,c(j)));
    H=blkdiag(H,Lnew(lst2,:)*pinv(Lold(lst,:)));
end

if(isa(ChanStats,'nirs.core.Data'))
    ChanStatsNew=ChanStats;
    ChanStatsNew.probe=NewProbe;
    lam=.1;
    [U,S,V]=nirs.math.mysvd(Lold);
    H = Lnew*V*inv(S+eye(size(S))*lam)*U';
    ChanStatsNew.data=(H*ChanStats.data')';
elseif(isa(ChanStats,'nirs.core.ChannelStats'))
    ChanStatsNew=ChanStats;
    ChanStatsNew.probe=NewProbe;
    
    cond=ChanStats.conditions;
    H2=[];
    c={};
    for id=1:length(cond)
        H3=zeros(size(H,1)*length(cond),size(H,2)*length(cond));
        lst=find(ismember(ChanStats.variables.cond,cond{id}));
        H3((id-1)*size(H,1)+1:id*size(H,1),lst)=H;
        H2=blkdiag(H2,H3);
        c=[c; repmat({cond{id}},size(H,1),1)];
    end
    
    if(length(cond)>1)
        error('not fully tested')
    end
    cond=c;
    
    ChanStatsNew.beta=H2*ChanStats.beta;
    covb=H2*ChanStats.covb*H2';
    covb=covb+1E-4*max(diag(covb))*eye(size(covb));
    
    c=(triu(covb,1)+tril(covb,-1)')/2;
    c=c+c'+diag(diag(covb));
    ChanStatsNew.covb=covb;
    
    ChanStatsNew.variables=[ChanStatsNew.probe.link table(cond)];
    
end
end
