function cfg = getConfig( obj, idx, idxType )
%GETCONFIG Summary of this function goes here
%   Detailed explanation goes here

iLambda=idx(2);
idx=idx(1);

    assert( all(obj.image.dim == obj.image.dim(1)) )
       
	cfg.unitinmm = obj.image.dim(1);
    
    cfg.nphoton = obj.nPhotons;
    cfg.vol = uint8(obj.image.vol);
    
    %Assume optodes point in the direction of the center of mass
    [x,y,z]=meshgrid([0:size(obj.image.vol,1)-1]*obj.image.dim(1)+obj.image.origin(1),...
                      [0:size(obj.image.vol,2)-1]*obj.image.dim(2)+obj.image.origin(2),...
                      [0:size(obj.image.vol,3)-1]*obj.image.dim(3)+obj.image.origin(3));
     CenterMass=[sum(obj.image.vol(:).*x(:)),sum(obj.image.vol(:).*y(:)),sum(obj.image.vol(:).*z(:))]./sum(obj.image.vol(:));             
                  
                  
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
    
    cfg.srcpos = cfg.srcpos / obj.image.dim(1) + obj.image.origin; % convert to grid units
    cfg.detpos = cfg.detpos / obj.image.dim(1) ...
        + repmat([obj.image.origin 0],[size(cfg.detpos,1) 1]); % convert to grid units
    
    cfg.gpuid = obj.gpuId;
    
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
                obj.prop{iLayer}.g...
                obj.prop{iLayer}.ri] );
            else
                cfg.prop = double( [cfg.prop;...
                obj.prop(iLayer).mua(iW)...
                obj.prop(iLayer).mus(iW)...
                obj.prop(iLayer).g...
                obj.prop(iLayer).ri] );
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

