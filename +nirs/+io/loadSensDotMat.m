function [J, mesh, probe] = loadSensDotMat( filename )

    tmp = load(filename, 'sensitivity','probe');
    
    %% mesh
    mesh=nirs.core.Mesh;
    mesh.nodes = tmp.sensitivity.mesh.vertices;
    mesh.faces = tmp.sensitivity.mesh.faces;
    mesh.nodes(:,[2 3])=mesh.nodes(:,[3 2]);
    mesh.nodes(:,[3])=-mesh.nodes(:,[3]);
    
    %% probe
    % can't find wavelength info in file so making big assumptions here
    if length( unique( tmp.probe.ml(:,4) ) ) == 1
        lambda = 808; 
    else
        lambda = [690 830]';
    end
    
	lambda = lambda( tmp.probe.ml(:,4) );
    
    link = table(tmp.probe.ml(:,1),tmp.probe.ml(:,2),lambda, ...
        'VariableNames',{'source','detector','type'});
    
    srcPos = tmp.probe.srcpos;
    detPos = tmp.probe.detpos;
    
    [link, idx] = sortrows(link,{'type','source','detector'});
    
    probe = nirs.core.Probe( srcPos, detPos, link );
    
    
    %% jacobian
    A = double( tmp.sensitivity.Adot );
    A = A(idx,:)*1e-7; %/ max(abs(A(:))) * 1; % scaling?
    
    ext = nirs.media.getspectra( link.type );
    
    J.hbo = bsxfun(@times, ext(:,1), A);
    J.hbr = bsxfun(@times, ext(:,2), A);

    J.hbo( abs(J.hbo(:)) < 1e-16*max(abs(J.hbo(:))) ) = 0;
    J.hbr( abs(J.hbr(:)) < 1e-16*max(abs(J.hbr(:))) ) = 0;

    J.hbo = sparse(J.hbo);
    J.hbr = sparse(J.hbr);

end

