function [data,mesh]=femdata_fl(mesh,frequency)

% [data,mesh]=femdata_fl(fn,frequency)
%
% Calculates fluorescence data (phase and amplitude) for a given
% problem (fn) at a given frequency (MHz).
% outputs phase and amplitude in structure data
% and mesh information in mesh

% 2011 May 2: Modified program to allow user to choose what forward models to
% calculate.  The programs currently calculates the excitation field, an
% intrinsic field at the emission wavelength (usually not used) and a
% fluroescence emission field.  In the reconstruct program, we often use
% femdata_fl to get the excitation field and don't need the other fields.
% Also, the intrinsic emission field is almost never used.  So, a flagging
% system based on additional fields in the mesh has been introduced.  The
% user can add mesh.mm (for intrinsic emission) or mesh.fl (for fluorescence
% emission) to signal what fields to calculate.
% - If neither field (mesh.mm or mesh. fl) exists, the fields
%           WILL be calculated by default.
% - if a field (mesh.mm or mesh. fl) exists and ~= 1,
%           that field WILL NOT be calculated.
% - if a field (mesh.mm or mesh. fl) exists and == 1,
%           that field WILL be calculated
% - if mesh.phix exists, skip calulating the excitation field and use the
% input field mesh.phix.

% If not a workspace variable, load mesh
if ischar(mesh)== 1
    mesh = load_mesh(mesh);
end

% Use fields to flag what forward models to calculate
if isfield(mesh,'phix') == 1
    data.phix = mesh.phix;
    xflag = 0;  % Excitation field phix is a field in the mesh
    % structured variable, so skip this calculation.
else
    xflag = 1;
end

% Use fields to flag what forward models to calculate
if isfield(mesh,'mm') == 0
    mmflag = 1; % Run forward model for intrinsic emission field
elseif isfield(mesh,'mm') == 1
    if mesh.mm == 1
        mmflag = 1; % Run forward model for intrinsic emission field
    elseif mesh.mm ~= 1
        mmflag = 0; % Skip forward model for intrinsic emission field
    end
end

if isfield(mesh,'fl') == 0
    flflag = 1; % Run forward model for fluorescence emission field
elseif isfield(mesh,'fl') == 1
    if mesh.fl == 1
        flflag = 1; % Run forward model for fluorescence emission field
    elseif mesh.fl ~= 1
        flflag = 0; % Skip forward model for fluorescence emission field
        % (used to calculate excitation field only in recons)
    end
end

% error checking
if frequency < 0
    errordlg('Frequency must be nonnegative','NIRFAST Error');
    error('Frequency must be nonnegative');
end

% modulation frequency
omega = 2*pi*frequency*1e6;

% set fluorescence variables
mesh.gamma = (mesh.eta.*mesh.muaf)./(1+(omega.*mesh.tau).^2);

% Create FEM matricex
if mesh.dimension == 2
    if xflag == 1
        % Excitation FEM matrix
        [i,j,s] = gen_matrices_2d(mesh.nodes(:,1:2),...
            sort(mesh.elements')', ...
            mesh.bndvtx,...
            mesh.muax,...
            mesh.kappax,...
            mesh.ksi,...
            mesh.c,...
            omega);
        junk = length(find(i==0));
        MASS_x = sparse(i(1:end-junk),j(1:end-junk),s(1:end-junk));
        clear junk i j s;
    end
    
    if mmflag == 1 || flflag == 1
        % Emission FEM matrix
        [i,j,s] = gen_matrices_2d(mesh.nodes(:,1:2),...
            sort(mesh.elements')', ...
            mesh.bndvtx,...
            mesh.muam,...
            mesh.kappam,...
            mesh.ksi,...
            mesh.c,...
            omega);
        junk = length(find(i==0));
        MASS_m = sparse(i(1:end-junk),j(1:end-junk),s(1:end-junk));
        clear junk i j s;
    end
elseif mesh.dimension ==3
    
    if xflag == 1
        % Excitation FEM matrix
        [i,j,s] = gen_matrices_3d(mesh.nodes,...
            sort(mesh.elements')', ...
            mesh.bndvtx,...
            mesh.muax,...
            mesh.kappax,...
            mesh.ksi,...
            mesh.c,...
            omega);
        junk = length(find(i==0));
        MASS_x = sparse(i(1:end-junk),j(1:end-junk),s(1:end-junk));
        clear junk i j s;
    end
    
    if mmflag == 1 || flflag == 1
        % Emission FEM matrix
        [i,j,s] = gen_matrices_3d(mesh.nodes,...
            sort(mesh.elements')', ...
            mesh.bndvtx,...
            mesh.muam,...
            mesh.kappam,...
            mesh.ksi,...
            mesh.c,...
            omega);
        junk = length(find(i==0));
        MASS_m = sparse(i(1:end-junk),j(1:end-junk),s(1:end-junk));
        clear junk i j s;
    end
end

% If the fn.ident exists, then we must modify the FEM matrices to
% account for refractive index mismatch within internal boundaries
% if isfield(mesh,'ident') == 1
%     disp('Modifying for refractive index')
%     M = bound_int(MASS,mesh);
%     MASS = M;
%     clear M
% end

% Calculate the RHS (the source vectors. For simplicity, we are
% just going to use a Gaussian Source, The width of the Gaussian is
% changeable (last argument). The source is assumed to have a
% complex amplitude of complex(cos(0.15),sin(0.15));

% Now calculate excitation source vector
% NOTE last term in mex file 'qvec' is the source FWHM
% First, reconcile link and source variables
source = unique(mesh.link(:,1));
[nnodes,junk]=size(mesh.nodes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory
ind = mesh.link(:,3)==0;
foo = mesh.link;
foo(ind,:)=[]; clear ind
source = unique(foo(:,1));
nsource = length(source);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if xflag == 1 || mmflag == 1
    qvec = spalloc(nnodes,nsource,nsource*100);
    if mesh.dimension == 2
        for i = 1 : nsource
            s_ind = mesh.source.num == source(i);
            if mesh.source.fwhm(s_ind) == 0
                qvec(:,i) = gen_source_point(mesh,mesh.source.coord(s_ind,1:2));
            else
                qvec(:,i) = gen_source(mesh.nodes(:,1:2),...
                    sort(mesh.elements')',...
                    mesh.dimension,...
                    mesh.source.coord(s_ind,1:2),...
                    mesh.source.fwhm(s_ind));
            end
        end
    elseif mesh.dimension == 3
        for i = 1 : nsource
            s_ind = mesh.source.num == source(i);
            if mesh.source.fwhm(s_ind) == 0
                qvec(:,i) = gen_source_point(mesh,mesh.source.coord(s_ind,1:3));
            else
                qvec(:,i) = gen_source(mesh.nodes,...
                    sort(mesh.elements')',...
                    mesh.dimension,...
                    mesh.source.coord(s_ind,:),...
                    mesh.source.fwhm(s_ind));
            end
        end
    end
    
    clear junk i;
    
    % catch error in source vector
    junk = sum(qvec);
    junk = find(junk==0);
    if ~isempty(junk)
        display(['WARNING...Check the FWHM of Sources ' num2str(junk)]);
    end
    clear junk
    
    % Catch zero frequency (CW) here
    if frequency == 0
        if xflag == 1
            MASS_x = real(MASS_x);
        end
        if mmflag == 1 || flflag == 1
            MASS_m = real(MASS_m);
        end
        qvec = real(qvec);
        %   [i,j,s]=find(qvec);
        %   [m,n] = size(s);
        %   s = complex(s,1e-20);
        %   qvec = sparse(i,j,s,nnodes,nsource);
        %   clear i j s m n
    end
end

% Check for distributed source
if mesh.source.distributed == 1
    qvec = sum(qvec,2);
end

%*******************************************************
% Calculate INTRINSIC FIELDS:

% Calculate INTRINSIC EXCITATION field for all sources
if xflag == 1
    [data.phix]=get_field(MASS_x,mesh,qvec);
end

% Calculate INTRINSIC field at EMISSION WAVELENGTH laser source
if mmflag == 1
    [data.phimm]=get_field(MASS_m,mesh,qvec);
end

clear qvec;
%********************************************************
% FLUORESCENCE EMISSION FIELDS
if flflag == 1
    % Calculate the RHS (the source vectors) for the FLUORESCENCE EMISSION.
    qvec = zeros(nnodes,nsource);
    % Simplify the RHS of emission equation
    beta = mesh.gamma.*(1-(sqrt(-1).*omega.*mesh.tau));
    % get rid of any zeros!
    if frequency == 0
        beta(beta==0) = 1e-20;
    else
        beta(beta==0) = complex(1e-20,1e-20);
    end
    
    if mesh.dimension == 2
        for i = 1 : nsource
            val = beta.*data.phix(:,i);
            qvec(:,i) = gen_source_fl(mesh.nodes(:,1:2),...
                sort(mesh.elements')',...
                mesh.dimension,...
                val);
        end
    elseif mesh.dimension == 3
        for i = 1 : nsource
            val = beta.*data.phix(:,i);
            qvec(:,i) = gen_source_fl(mesh.nodes,...
                sort(mesh.elements')',...
                mesh.dimension,...
                val);
        end
    end
    clear junk i nnodes nsource val beta;
    
    % Catch zero frequency (CW) here
    if frequency == 0
        qvec = real(qvec);
    end
    
    %*********************************************************
    % Calculate FLUORESCENCE EMISSION field for all sources
    [data.phifl]=get_field(MASS_m,mesh,qvec);
    clear qvec;
end

%*********************************************************
% EXTRACT DATA
data.link = mesh.link;

% Extract Excitation data
[data.complexx]=get_boundary_data(mesh,data.phix);
% Map complex data to amplitude
data.amplitudex = abs(data.complexx);
% Map complex data to phase
data.phasex = atan2(imag(data.complexx),...
    real(data.complexx));
% Calculate phase in degrees
data.phasex(data.phasex<0) = data.phasex(data.phasex<0) + (2*pi);
data.phasex = data.phasex*180/pi;
% Build data format
data.paax = [data.amplitudex data.phasex];


% Extract Fluorescence Emission data
if flflag == 1
    [data.complexfl]=get_boundary_data(mesh,data.phifl);
    % Map complex data to amplitude
    data.amplitudefl = abs(data.complexfl);
    % Map complex data to phase
    data.phasefl = atan2(imag(data.complexfl),...
        real(data.complexfl));
    % Calculate phase in degrees
    data.phasefl(data.phasefl<0) = data.phasefl(data.phasefl<0) + (2*pi);
    data.phasefl = data.phasefl*180/pi;
    % Build data format
    data.paafl = [data.amplitudefl data.phasefl];
    data.paaxfl = [data.amplitudex data.phasex data.amplitudefl data.phasefl];
end

% Exrtact intrinsic emssion field data
if mmflag == 1
    [data.complexmm]=get_boundary_data(mesh,data.phimm);
    % Map complex data to amplitude
    data.amplitudemm = abs(data.complexmm);
    % Map complex data to phase
    data.phasemm = atan2(imag(data.complexmm),...
        real(data.complexmm));
    % Calculate phase in degrees
    data.phasemm(data.phasemm<0) = data.phasemm(data.phasemm<0) + (2*pi);
    data.phasemm = data.phasemm*180/pi;
    % Build data format
    data.paamm = [data.amplitudemm data.phasemm];
    if flflag == 1
        data.paaxflmm = [data.amplitudex data.phasex data.amplitudefl data.phasefl data.amplitudemm data.phasemm];
    end
end