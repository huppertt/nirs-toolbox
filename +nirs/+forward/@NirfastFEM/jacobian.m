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
    
    d = zeros(1,size(obj.probe.link,1));
    Jmua = zeros(size(obj.probe.link,1),size(obj.mesh.nodes,1));
    Jkappa = zeros(size(obj.probe.link,1),size(obj.mesh.nodes,1));
    
    types = unique( obj.probe.link.type );
    
    
    
    
    for i = 1:length( types )
        lst = obj.probe.link.type == types(i);
        [jac,data] = jacobian(mesh{i},obj.Fm);
        
        d(lst) = data.complex';
        
        if obj.Fm > 0
            Jkappa(lst,:) = conj( jac.complex(:,1:end/2) );
            Jmua(lst,:) = conj( jac.complex(:,end/2+1:end) );
        else
            Jmua(lst,:) = jac.complex;
        end
    end
    
    Jmua = abs(Jmua);
    

    
    if ~isSpectral
        J.mua = Jmua;
        
        if obj.Fm > 0
            J.kappa = Jkappa;
        end
    else
        ext = nirs.media.getspectra( types );

        [~,~,iType] = unique( obj.probe.link.type );
        
        ehbo = ext(iType,1);
        ehbr = ext(iType,2);
        
        J.hbo = bsxfun(@times,ehbo,Jmua);
        J.hbr = bsxfun(@times,ehbr,Jmua);

        if obj.Fm > 0
            J.kappa = Jkappa; 
        end            
%             % preallocation
%             J.mus = J.kappa;
%             J.a = J.kappa;
%             J.b = J.kappa;
%             for i = 1:length(obj.prop)
%                 lambda = obj.prop{i}.lambda(iLambda);
%                 lambda = lambda(:);
% 
%                 mus = obj.prop{i}.mus(iLambda);
%                 mus = mus(:);
% 
%                 J.mus(:,i) = -1/3./mus.^2 .* J.kappa(:,i);
% 
%                 if isprop( obj.prop{i},'b' )
%                     b = obj.prop{i}.b;
%                     J.a(:,i) = (lambda/500).^-b .* J.mus(:,i);
%                     J.b(:,i) = -mus .* log( lambda/500 ) .* J.mus(:,i);
%                 else
%                     J.a(:,i) = NaN;
%                     J.b(:,i) = NaN;
%                 end
%             end
    end
    
    meas = nirs.core.Data( d,0,obj.probe,obj.Fm );
    
end

%     for i = 1:max( obj.mesh.region(:) )
%         lst = obj.mesh.region == i;
%         J.mua(:,i) = sum( fullJ.mua(:,lst),2 ) ./ d.';
%         J.kappa(:,i) = sum( fullJ.kappa(:,lst),2 ) ./ d.';
%     end