function cfg = getConfig( obj, idx, idxType )
%GETCONFIG Summary of this function goes here
%   Detailed explanation goes here

%             node: [201113x3 double]
%             elem: [1225347x4 double]
%          nphoton: 100000
%         elemprop: [1225347x1 double]
%           srcpos: [30 30 0]
%           srcdir: [0 0 1]
%         srcspace: [5124x1 double]
%            media: {1x3 cell}
%             prop: [4x4 double]
%           tstart: 0
%             tend: 5.0000e-09
%            tstep: 5.0000e-10
%       debuglevel: 'TP'
%     isreoriented: 1
%           facenb: [1225347x4 double]
%             evol: [1225347x1 double]
%               e0: 159801

iLambda=idx(2);
idx=idx(1);

cfg.nphoton = obj.nPhotons;

cfg.node=obj.mesh.nodes;
cfg.elem=obj.mesh.elems;


prop(:,1)=obj.mesh.regions(obj.mesh.elems(:,1));
prop(:,2)=obj.mesh.regions(obj.mesh.elems(:,2));
prop(:,3)=obj.mesh.regions(obj.mesh.elems(:,3));
prop(:,4)=obj.mesh.regions(obj.mesh.elems(:,4));


cfg.elemprop=mode(prop,2);




    %Assume optodes point in the direction of the center of mass
  
    CenterMass=mean(cfg.node,1);             
                  
                  
    if strcmp( idxType,'source' )
        cfg.srcpos = obj.probe.srcPos(idx,:);
           cfg.srcdir = (cfg.srcpos-CenterMass.*ones(1,size(cfg.srcpos,1)))./...
            sqrt(sum((cfg.srcpos-CenterMass.*ones(1,size(cfg.srcpos,1))).^2,2)*ones(1,3));
        cfg.detpos = [obj.probe.detPos ones(size(obj.probe.detPos,1),1)];
        lst = obj.probe.link.source == idx;
    elseif strcmp( idxType, 'detector' )
        cfg.srcpos = obj.probe.detPos(idx,:);
         cfg.srcdir = (cfg.srcpos-CenterMass.*ones(1,size(cfg.srcpos,1)))./...
            sqrt(sum((cfg.srcpos-CenterMass.*ones(1,size(cfg.srcpos,1))).^2,2)*ones(1,3));
        cfg.detpos = [obj.probe.srcPos ones(size(obj.probe.srcPos,1),1)];
        lst = obj.probe.detector == idx;
    else
        error( 'Must specify "source" or "detector".' )
    end
    
    cfg.autopilot = 1;
    cfg.isnormalized = 0;
    
    cfg.prop = [0 0 1 1];
    
    for iLayer = 1:length( obj.prop )
        [~,iW]=ismember(iLambda,obj.prop{iLayer}.lambda);
       
        % SUPER DUMB BUG; this will not work if cfg.prop is not a double
            % CANNOT be a single or it will return 0 flux everywhere
            if iscell(obj.prop)
                cfg.prop = double( [cfg.prop;...
                obj.prop{iLayer}.mua(iW)...
                obj.prop{iLayer}.mus(iW)...
                obj.prop{iLayer}.g(1)...
                obj.prop{iLayer}.ri] );
            else
                cfg.prop = double( [cfg.prop;...
                obj.prop(iLayer).mua(iW)...
                obj.prop(iLayer).mus(iW)...
                obj.prop(iLayer).g(iW)...
                obj.prop(iLayer).ri] );
            end
    end
    
    cfg.tstart = 0;
    cfg.tstep = obj.timeStep;
    cfg.tend = obj.nTimeGates * obj.timeStep;
    
    
    cfg.isreflect = 1; % external boundary reflection 
    cfg.isrefint = 1; % internal boundary reflection
    
  
    
    % seed rng
    % cfg.seed = uint32( randi(2^32-1) );
    
    % not functioning properly as of MCX 0.89
	% cfg.respin = obj.nRepetitions;
end

