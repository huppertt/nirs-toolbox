function [ R ] = generateRMatrix( mesh, frequency, recParam )
%GENERATERMATRIX Summary of this function goes here
%   Detailed explanation goes here
if nargin == 2
    recParam.useWavelengthDependentR = 1;
end
if strcmp(mesh.type,'stnd')
    if isfield(mesh,'R') == 0 %This may have already been calculated.
        if length(mesh.nodes) >= 3800
            R = cholinc(generateMassMatrix(mesh,frequency),1e-3);
        elseif length(mesh.nodes) < 3800
            R = [];
        end
    else
        R = mesh.R;
    end
elseif strcmp(mesh.type,'spec') || strcmp(mesh.type,'specPenn')
    if recParam.useWavelengthDependentR
        if isfield(mesh,'R') == 0 %This may have already been calculated.
            if length(mesh.nodes) >= 3800
                R = cell(length(mesh.wv),1);
                for ii = 1:length(mesh.wv)
                    fprintf('Generating preconditioner for wavelength %i\n',ii);
                    R{ii} = cholinc(generateMassMatrix(mesh,frequency,mesh.wv(ii)),1e-3);
                end
            elseif length(mesh.nodes) < 3800
                R = [];
            end
        else
            R = mesh.R;
        end
    else
        if isfield(mesh,'R') == 0 %This may have already been calculated.
            if length(mesh.nodes) >= 3800
                R = cholinc(generateMassMatrix(mesh,frequency,mesh.wv(1)),1e-3);
            elseif length(mesh.nodes) < 3800
                R = [];
            end
        else
            R = mesh.R;
        end
    end
end

end

