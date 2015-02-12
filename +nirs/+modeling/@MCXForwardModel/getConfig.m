function cfg = getConfig( obj, idx, idxType )
%GETCONFIG Summary of this function goes here
%   Detailed explanation goes here
    assert( all(obj.image.dim == obj.image.dim(1)) )
       
	cfg.unitinmm = obj.image.dim(1);
    
    cfg.nphoton = obj.nPhotons;
    cfg.vol = obj.image.vol;
    
    if strcmp( idxType,'source' )
        cfg.srcpos = obj.probe.srcPos(idx,:);
        cfg.srcdir = obj.probe.srcDir(idx,:);
        cfg.detpos = [obj.probe.detPos ones(size(obj.probe.detPos,1),1)];
        lst = obj.probe.link(:,1) == idx;
    elseif strcmp( idxType, 'detector' )
        cfg.srcpos = obj.probe.detPos(idx,:);
        cfg.srcdir = obj.probe.detDir(idx,:);
        cfg.detpos = [obj.probe.srcPos ones(size(obj.probe.srcPos,1),1)];
        lst = obj.probe.link(:,2) == idx;
    else
        error( 'Must specifiy "source" or "detector".' )
    end
    
    cfg.srcpos = cfg.srcpos / obj.image.dim(1) + obj.image.origin; % convert to grid units
    cfg.detpos = cfg.detpos / obj.image.dim(1) ...
        + repmat([obj.image.origin 0],[size(cfg.detpos,1) 1]); % convert to grid units
    
    cfg.gpuid = obj.gpuId;
    
    cfg.autopilot = 1;
    cfg.isnormalized = 0;
    
    cfg.prop = [0 0 1 1];
    
    g = 0;
    
    for iLayer = 1:length( obj.prop )
        iLambda = unique(obj.probe.link(lst,3));
        
            % SUPER DUMB BUG; this will not work if cfg.prop is not a double
            % CANNOT be a single or it will return 0 flux everywhere
            if iscell(obj.prop)
                cfg.prop = double( [cfg.prop;...
                obj.prop{iLayer}.mua(iLambda)...
                obj.prop{iLayer}.mus(iLambda)/(1-g)...
                g obj.prop{iLayer}.ri] );
            else
                cfg.prop = double( [cfg.prop;...
                    obj.prop(iLayer).mua(iLambda)...
                    obj.prop(iLayer).mus(iLambda)/(1-g)...
                    g obj.prop(iLayer).ri] );
            end
    end
    
    cfg.tstart = 0;
    cfg.tstep = obj.timeStep;
    cfg.tend = obj.nTimeGates * obj.timeStep;
    
    cfg.isgpuinfo = 0;
    
    cfg.isreflect = 1; % external boundary reflection 
    cfg.isrefint = 1; % internal boundary reflection
    
    % seed rng
    % cfg.seed = uint32( randi(2^32-1) );
    
    % not functioning properly as of MCX 0.89
	% cfg.respin = obj.nRepetitions;
end

