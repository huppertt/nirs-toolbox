function [Jacobian,datafl] = jacobian_multispectral(mesh,frequency)

% [J,datafl] = jacobian_multispectral(mesh,frequency,datax)
%
% Calculates jacobian multispectral.
%
% mesh is the input mesh (variable or filename)
% frequency is the modulation frequency (MHz)



if frequency < 0
    errordlg('Frequency must be nonnegative','NIRFAST Error');
    error('Frequency must be nonnegative');
end

% If not a workspace variable, load mesh
if ischar(mesh)== 1
    mesh = load_mesh(mesh);
end

[Data, mesh] = femdata_fl(mesh, frequency);

% create adjoint mesh
AdjMesh = mesh;
AdjMesh.source = mesh.meas;
AdjMesh.source.fwhm = ones(size(mesh.source.coord,1),1).*mesh.source.fwhm(1);
AdjMesh.meas = mesh.source;
if mesh.dimension == 2
  [ind,int_func] = mytsearchn(mesh,AdjMesh.meas.coord(:,1:2));
elseif mesh.dimension == 3
  [ind,int_func] = mytsearchn(mesh,AdjMesh.meas.coord);
end
AdjMesh.meas.int_func = [ind int_func];

[AdjData, AdjMesh] = femdata_fl(AdjMesh, frequency);
Phix = real(Data.phix');
AdjPhix = real(AdjData.phix');
Link = mesh.link;
NumberNodes = size(mesh.muax, 1);
NumberMeasurements = size(find(Link ~= 0), 1);

JacLineCount = 1;
Jacobian = zeros(NumberMeasurements, NumberNodes);
for i = 1:size(mesh.source.coord)
    Field = Phix(i, :);
    for j = 1:length(mesh.link(i,:))
        RankDetector = Link(i, j);
        if (RankDetector ~= 0)
            AdjField = AdjPhix(RankDetector, :); 
            Jacobian(JacLineCount, :) = Field .* AdjField;
            JacLineCount = JacLineCount + 1;
        end
    end
end

datafl = Data;