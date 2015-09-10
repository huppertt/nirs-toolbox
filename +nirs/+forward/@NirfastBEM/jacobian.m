function [J, meas] = jacobian( obj, type )

    if nargin < 2
        isSpectral = false;
    elseif strcmpi( type,'standard' )
        isSpectral = false;
    elseif strcmpi( type,'spectral' );
        isSpectral = true;
    else
        error('Jacobian can either be ''standard'' or ''spectral''.')
    end
    
    mesh = obj.getNirfastMeshes();
    d_src = zeros(1,size(obj.probe.link,1));
	d_det = zeros(1,size(obj.probe.link,1));
    
%     phi_src = zeros( size(obj.probe.link,1), size(obj.mesh{1}.nodes,1) );
%     phi_det = zeros( size(obj.probe.link,1), size(obj.mesh{1}.nodes,1) );
    
    types = unique( obj.probe.link.type );
    assert( isnumeric( types ) );
    
    % sources
    for i = 1:length( types ) 

        lst = obj.probe.link.type == types(i);
        data = bemdata_stnd( mesh{i},obj.Fm );
        
        if obj.Fm == 0
            thisD = data.paa(:,1);
            thisD( thisD < 0 ) = 0;
        else
            thisD = data.paa(:,1) .* exp(1i * pi * data.paa(:,2)/180);
        end
        
        
        d_src(lst) = thisD;
        phi_src{i} = data.phi;
        
    end
    
    thisObj = obj;
    thisObj.probe = thisObj.probe.swapSD();
    mesh = thisObj.getNirfastMeshes();
    
    % detectors
    paa=zeros(height(obj.probe.link),1);
    for i = 1:length( types )

        lst = obj.probe.link.type == types(i);
        data = bemdata_stnd( mesh{i},obj.Fm );
        
        if obj.Fm == 0
            thisD = data.paa(:,1);
            thisD( thisD < 0 ) = 0;
        else
            thisD = data.paa(:,1) .* exp(1i * pi * data.paa(:,2)/180);
        end
        
        d_det(lst) = thisD;
        phi_det{i} = data.phi;
        paa(lst)=data.paa(:,1);
        
    end
    
    meas = nirs.core.Data( d_det/2 + d_src/2,0,obj.probe,obj.Fm );
    
    % assemble jacobian
    Jmua = zeros( size(obj.probe.link,1), size(mesh{1}.nodes,1) );

    [~,~,iType] = unique( obj.probe.link.type );
    for i = 1:length(obj.probe.link.type)
        iSrc = obj.probe.link.source(i);
        iDet = obj.probe.link.detector(i);
        
        Jmua(i,:) = ( phi_src{iType(i)}(:,iSrc) .* phi_det{iType(i)}(iDet) ).';
    end
    Jmua = abs(Jmua);
    Jmua = Jmua./(sum(Jmua,2)*ones(1,size(Jmua,2)));
    Jmua = Jmua.*(paa*ones(1,size(Jmua,2)));

    
    if ~isSpectral
        J.mua = Jmua;
    else
        % convert jacobian to conc
        ext = nirs.media.getspectra( types );
                
        ehbo = ext(iType,1);
        ehbr = ext(iType,2);
        
        J.hbo = bsxfun(@times,ehbo,Jmua);
        J.hbr = bsxfun(@times,ehbr,Jmua);
        
    end
    
    
    
end

