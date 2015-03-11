function [J, mesh, probe] = loadSensDotMat( filename )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    tmp = load(filename, 'sensitivity','probe');
    
    %% mesh
    mesh.nodes = tmp.sensitivity.mesh.vertices;
    mesh.faces = tmp.sensitivity.mesh.faces;
    
    
    %% probe
    lambda = [690 830]';
    lambda = lambda( tmp.probe.ml(:,4) );
    
    link = table(tmp.probe.ml(:,1),tmp.probe.ml(:,2),lambda, ...
        'VariableNames',{'source','detector','type'});
    
    srcPos = tmp.probe.srcpos;
    detPos = tmp.probe.detpos;
    
    [link, idx] = sortrows(link,{'type','source','detector'});
    
    probe = nirs.Probe( srcPos, detPos, link );
    
    
    %% jacobian
    A = double( tmp.sensitivity.Adot );
    A = -A(idx,:) / max(abs(A(:))) * 1e-3; % scaling?
    
    ext = nirs.getSpectra( link.type );
    
    J.hbo = bsxfun(@times, ext(:,1), A);
    J.hbr = bsxfun(@times, ext(:,2), A);

    J.hbo( abs(J.hbo(:)) < 1e-12*max(abs(J.hbo(:))) ) = 0;
    J.hbr( abs(J.hbr(:)) < 1e-12*max(abs(J.hbr(:))) ) = 0;

    J.hbo = sparse(J.hbo);
    J.hbr = sparse(J.hbr);

%     %% output
%     fwdModel.J = J;
%     fwdModel.mesh = mesh;
%     fwdModel.probe = probe;

end

