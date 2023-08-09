function [ z,data,z_amp,z_phase ] = adjointGradientStnd( mesh, measurements,frequency,regParam)
%ADJOINTGRADIENT Calculates gradient for standard mesh using adjoint method.
%   Less memory and computationally expensive than calculating gradient
%   from Jacobian. Returns z of form [z_kappa; z_mua].

%% If no regularisation parameters have been provided, assume none.
if nargin < 4
    regParam.lambda = 0;
end
if ~isfield(regParam, 'P')
    regParam.P = 1;
end
if ~isfield(mesh,'cgFormatLink')
    mesh = convertLinkFormat(mesh);
end

%% Calculate MASS matrix (code from jacobian_stnd.m)
% MASS is N_n x N_n but sparse
MASS = generateMassMatrix(mesh,frequency);

mesh.R = generateRMatrix(mesh,MASS);

%% Calculate source vectors (code from jacobian_stnd.m)
% qvec is N_n x N_s but sparse
[nnodes]=size(mesh.nodes,1);
[nsource]=size(mesh.source.coord,1);

%Presort elements to save time.
sortedElements = sort(mesh.elements,2);

%Create matrix transposed as faster for sparse matrices.
qvec = spalloc(nsource,nnodes,nsource*100);
if mesh.dimension == 2
    for i = 1 : nsource
        if mesh.source.fwhm(i) == 0
            qvec(i,:) = gen_source_point(mesh,mesh.source.coord(i,1:2));
        else
            qvec(i,:) = gen_source(mesh.nodes(:,1:2),...
                sortedElements,...
                mesh.dimension,...
                mesh.source.coord(i,1:2),...
                mesh.source.fwhm(i));
        end
    end
elseif mesh.dimension == 3
    for i = 1 : nsource
        if mesh.source.fwhm(i) == 0
            qvec(i,:) = gen_source_point(mesh,mesh.source.coord(i,1:3));
        else
            qvec(i,:) = gen_source(mesh.nodes,...
                sortedElements,...
                mesh.dimension,...
                mesh.source.coord(i,:),...
                mesh.source.fwhm(i));
        end
    end
end
qvec = conj(qvec');
clear junk i w;

% Catch zero frequency (CW) here
if frequency == 0
    MASS = real(MASS);
    qvec = real(qvec);
end

% catch error in source vector
junk = sum(qvec);
junk = find(junk==0);
if ~isempty(junk)
    display(['WARNING...Check the FWHM of Sources ' num2str(junk)]);
end
clear junk
%% Calculate forward fields due to sources.
% Field is N_n x N_s
[data.phi] = get_field(MASS,mesh,qvec);

%% Calculate boundary data. (Copied from femdata_stnd.m)
% Data is k x N_m

% Calculate boundary data
[data.complex]=get_boundary_data(mesh,data.phi);

% Map complex data to amplitude and phase
data.amplitude = abs(data.complex);

data.phase = atan2(imag(data.complex),...
    real(data.complex));
data.phase = mod(data.phase,2*pi);
data.phase = data.phase*180/pi;
if frequency == 0
    data.phase(:) = 0;
end

data.paa = [data.amplitude data.phase];
%% Create residual matrices of log measurements and phases.
lnb = log(measurements.paa((mesh.cgFormatLink(:,3) == 1),1)./data.paa((mesh.cgFormatLink(:,3) == 1),1));
phaseb = pi*(measurements.paa((mesh.cgFormatLink(:,3) == 1),2) - data.paa((mesh.cgFormatLink(:,3) == 1),2))/180;
% Make phase difference periodic, so it belongs to the interval [-pi, pi].
phaseb(phaseb > pi) = phaseb(phaseb > pi) - 2*pi;
phaseb(phaseb < -pi) = phaseb(phaseb < -pi)+ 2*pi;

b = [lnb;phaseb];
%Create measurement scaling matrix and apply it to residual.
if isscalar(regParam.P)
    pMat = regParam.P;
else
    pMat = sparse(regParam.P(:,1),regParam.P(:,2),regParam.P(:,3));
end
b = (pMat.^0.5)*b;
lnb = b(1:length(lnb));
phaseb = b((length(lnb)+1):end);
tmplnb = zeros(size(mesh.cgFormatLink,1),1);
tmplnb((mesh.cgFormatLink(:,3) == 1)) = lnb;
tmpphaseb = zeros(size(mesh.cgFormatLink,1),1);
tmpphaseb((mesh.cgFormatLink(:,3) == 1)) = phaseb;
lnb = tmplnb;
phaseb = tmpphaseb;
clear b
%% Calculate adjoint fields.
eta0 = zeros(nnodes,nsource);
if frequency > 0
    eta1 = zeros(nnodes,nsource);
end

numberOfSourcesUsed = length(unique(mesh.cgFormatLink(:,1)));
numberOfDetectorsUsed = length(unique(mesh.cgFormatLink(:,2)));

if ((frequency == 0) && (numberOfDetectorsUsed < numberOfSourcesUsed)) ||...
        ((frequency > 0) && (numberOfDetectorsUsed < (2*numberOfSourcesUsed)))
    
    fprintf('Using adjoint field calculation method optimised for detectors.\n');
    
    % To avoid unnecessary for loops pre-calculate indexing. We are calculating
    % adjoint fields, so want to calculate an adjoint field from a detector and
    % deal with the all the measurements arising from sources paired with the
    % specific detector.
    
    %%
    
    detectors = unique(mesh.cgFormatLink(:,2));
    for ii = 1:length(detectors)
        currentDetector = detectors(ii);
        tmpI = mesh.elements(mesh.meas.int_func(currentDetector,1),:);
        tmpS = mesh.meas.int_func(currentDetector,2:end)' .* ...
            complex(cos(0.15),sin(0.15));
        qvec = sparse(tmpI,ones(size(tmpI)),tmpS,nnodes,1);
        if frequency == 0
            qvec = real(qvec);
        end
        aphi = get_field(conj(MASS),mesh,conj(qvec));
        clear tmpI tempS qvec
        measurementIndices = find(mesh.cgFormatLink(:,2) == currentDetector);
        for jj = 1:length(measurementIndices)
            currentMeasurement = measurementIndices(jj);
            if (mesh.cgFormatLink(currentMeasurement,3) == 0)
                continue;
            end
            eta0Delta = lnb(currentMeasurement)*aphi/conj(data.complex(currentMeasurement));
            eta0(:,mesh.cgFormatLink(currentMeasurement,1)) = eta0(:,mesh.cgFormatLink(currentMeasurement,1)) + eta0Delta;
            if frequency > 0
                eta1Delta = phaseb(currentMeasurement)*aphi/conj(data.complex(currentMeasurement));
                eta1(:,mesh.cgFormatLink(currentMeasurement,1)) = eta1(:,mesh.cgFormatLink(currentMeasurement,1)) + eta1Delta;
            end
        end
    end
    
else
    fprintf('Using adjoint field calculation method optimised for sources.\n');
    
    % Limit number of matrix inversions by summing up adjoint sources
    % first, and then solving to obtain total adjoint field. Requires
    % noSources matrix inversions in CW mode, and 2 x noSources matrix
    % inversions in frequency mode.
    %%
    sources = unique(mesh.cgFormatLink(:,1));
    for ii = 1:length(sources)
        currentSource = sources(ii);
        qvec = zeros(size(mesh.nodes,1),1);
        if frequency > 0
            qvec1 = zeros(size(qvec));
        end
        measurementIndices = find(mesh.cgFormatLink(:,1) == currentSource);
        for jj = 1:length(measurementIndices)
            currentMeasurement = measurementIndices(jj);
            if (mesh.cgFormatLink(currentMeasurement,3) == 0)
                continue;
            end
            currentDetector = mesh.cgFormatLink(currentMeasurement,2);
            tmpI = mesh.elements(mesh.meas.int_func(currentDetector,1),:);
            tmpS = mesh.meas.int_func(currentDetector,2:end)' .* ...
                complex(cos(0.15),sin(0.15));
            dqvec = zeros(size(qvec));
            dqvec(tmpI) = tmpS;
            if frequency == 0
                dqvec = real(dqvec);
            end
            qvec = qvec + lnb(currentMeasurement)*conj(dqvec)/conj(data.complex(currentMeasurement));
            if frequency > 0
                qvec1 = qvec1 + phaseb(currentMeasurement)*conj(dqvec)/conj(data.complex(currentMeasurement));
            end
        end
        eta0(:,ii) = get_field(conj(MASS),mesh,qvec);
        if frequency > 0
            eta1(:,ii) = get_field(conj(MASS),mesh,qvec1);
        end
    end    
    
end
%% Preallocate gradient vectors.
z_kappa_amp = zeros(length(mesh.nodes),1);
z_mua_amp = zeros(length(mesh.nodes),1);
if frequency > 0
    z_kappa_phase = zeros(length(mesh.nodes),1);
    z_mua_phase = zeros(length(mesh.nodes),1);
end

%% Sum gradients resulting from adjoint fields.
for jj = 1:nsource
    dphi = full(data.phi(:,jj));
    aphi0 = full(eta0(:,jj));
    if frequency > 0
        aphi1 = full(eta1(:,jj));
    end
    
    %% Calculate mua gradient.
    if mesh.dimension == 2
        tmp = IntFG(mesh.nodes(:,1:2),sortedElements,...
            mesh.element_area,dphi,conj(aphi0));
        if frequency > 0
            tmp2 = IntFG(mesh.nodes(:,1:2),sortedElements,...
                mesh.element_area,dphi,conj(aphi1));
        end
    elseif mesh.dimension == 3
        tmp = IntFG_tet4(mesh.nodes,sortedElements,...
            mesh.element_area,dphi,conj(aphi0));
        if frequency > 0
            tmp2 = IntFG_tet4(mesh.nodes,sortedElements,...
                mesh.element_area,dphi,conj(aphi1));
        end
    end
    z_mua_amp(:,1) = z_mua_amp(:,1) + real(tmp);
    if frequency > 0
        z_mua_phase(:,1) = z_mua_phase(:,1) + imag(tmp2);
    end
    %% Calculate kappa gradient (sorting must be removed for this
    %% gradient).
    if mesh.dimension == 2
        tmp = IntgradFgradG(mesh.nodes(:,1:2),mesh.elements,...
            mesh.element_area,dphi,conj(aphi0));
        if frequency > 0
            tmp2 = IntgradFgradG(mesh.nodes(:,1:2),mesh.elements,...
                mesh.element_area,dphi,conj(aphi1));
        end
    elseif mesh.dimension == 3
        tmp = IntgradFgradG_tet4(mesh.nodes,mesh.elements,...
            mesh.element_area,dphi,conj(aphi0));
        if frequency > 0
            tmp2 = IntgradFgradG_tet4(mesh.nodes,mesh.elements,...
                mesh.element_area,dphi,conj(aphi1));
        end
    end
    z_kappa_amp(:,1) = z_kappa_amp(:,1) + real(tmp);
    if frequency > 0
        z_kappa_phase(:,1) = z_kappa_phase(:,1) + imag(tmp2);
    end
end

%% Construct gradient from amplitude and phase components.
if frequency > 0
    fprintf('kappa gradient amp/phase ratio is %f\n',norm(z_kappa_amp)/norm(z_kappa_phase));
    fprintf('mua gradient amp/phase ratio is %f\n',norm(z_mua_amp)/norm(z_mua_phase));
    z = [z_kappa_amp+z_kappa_phase;...
        z_mua_amp+z_mua_phase];
    z_amp = [z_kappa_amp;z_mua_amp];
    z_phase = [z_kappa_phase;z_mua_phase];
else
    z = [z_kappa_amp;z_mua_amp];
    z_amp = z;
    z_phase = zeros(size(z));
end
end


