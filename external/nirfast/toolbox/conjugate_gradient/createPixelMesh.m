function [ mesh, recon_mesh ] = createPixelMesh( mesh,pixelBasis )
%CREATEPIXELBASIS Creates a pixel basis and adds interpolation data to
%original mesh.
%   
[mesh.fine2coarse,recon_mesh] = pixel_basis(pixelBasis,mesh);
    recon_mesh.type = mesh.type;
    recon_mesh.link = mesh.link;
    recon_mesh.source = mesh.source;
    recon_mesh.meas = mesh.meas;
    recon_mesh.dimension = mesh.dimension;
    if recon_mesh.dimension == 2
        recon_mesh.element_area = ele_area_c(recon_mesh.nodes(:,1:2),...
                                             recon_mesh.elements);
    else
        recon_mesh.element_area = ele_area_c(recon_mesh.nodes,...
                                             recon_mesh.elements);
    end

end

