function cfg = getConfig( obj, idx, idxType )
%GETCONFIG Summary of this function goes here
%   Detailed explanation goes here
    if any( obj.image.dim ~= obj.image.dim(1) )
        error( 'Image must be isotropic dimensions.' );
    else
        cfg.unitinmm = obj.image.dim(1);
    end
    
    cfg.nphoton = obj.nPhotons;
    cfg.vol = obj.image.volume;
    
    if strcmp( idxType,'source' )
        cfg.srcpos = obj.probe.srcPos(idx,:);
        cfg.srcdir = obj.probe.srcDir(idx,:);
        cfg.detpos = [obj.probe.detPos 2*ones(size(obj.probe.detPos,1),1)];
        lst = obj.probe.link(:,1) == idx;
    elseif strcmp( idxType, 'detector' )
        cfg.srcpos = obj.probe.detPos(idx,:);
        cfg.srcdir = obj.probe.detDir(idx,:);
        cfg.detpos = [obj.probe.srcPos 2*ones(size(obj.probe.srcPos,1),1)];
        lst = obj.probe.link(:,2) == idx;
    else
        error( 'Must specifiy "source" or "detector".' )
    end
    
    cfg.srcpos = cfg.srcpos / obj.image.dim(1); % convert to grid units
%     cfg.detpos = cfg.detpos / obj.image.dim(1); 
    
    cfg.gpuid = obj.gpuId;
    cfg.autopilot = 1;
    cfg.isnormalized = 0;
    
    cfg.prop = [0 0 1 1];
    
    g = 0;
    
    for iLayer = 1:length( obj.optProp )
        iLambda = unique(obj.probe.link(lst,3));
%         iLambda = find( obj.optProp(iLayer).lambda == lambda );

        % SUPER DUMB BUG; this will not work if cfg.prop is not a double
        % CANNOT be a single or it will return 0 fluence everywhere
%         cfg.prop = double( [cfg.prop;...
%             obj.optProp(iLayer).mua(iLambda)...
%             obj.optProp(iLayer).mus(iLambda)...
%             0 obj.optProp(iLayer).ri] );
        
        cfg.prop = double( [cfg.prop;...
            obj.optProp(iLayer).mua(iLambda)...
            obj.optProp(iLayer).mus(iLambda)/(1-g)...
            g obj.optProp(iLayer).ri] );
    end
    
    cfg.tstart = 0;
    cfg.tstep = obj.timeStep;
    cfg.tend = obj.nTimeGates * obj.timeStep;
    
    cfg.isgpuinfo = 0;
    cfg.isreflect = 1;
    cfg.isrefint = 1;
    cfg.seed = uint32( randi(2^32-1) );
    
    cfg.respin = obj.nRepetitions;
end

