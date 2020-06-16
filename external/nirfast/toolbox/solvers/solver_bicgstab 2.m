function [phi,R]=solver_bicgstab(Mass,mesh,qvec)

% [phi,R]=solver_bicgstab(Mass,mesh,qvec)
%
% Used by femdata and jacobian
% Calculates the field, given the Mass matrix and RHS source
% vector.
% Uses the biconjugate gradient stabilized method
%
% Mass is the mass matrix
% mesh is the input mesh
% qvec is the RHS source vector
% R is the preconditioner
% phi is the field

[nnodes,nsource]=size(qvec);
msg=[];
flag = 0;
phi=zeros(nnodes,nsource);

if isfield(mesh,'R') == 0
    if exist('ichol')
        % keeping this out until         
        R = ichol(Mass);
        %R = cholinc(Mass,1e-3);
    else
        R = cholinc(Mass,1e-3);
    end
    mesh.R = R;
else
    R = mesh.R;
end

% if ~isfield(mesh,'L')
%     [L,U] = ilu(Mass,struct('type','ilutp','droptol',1e-2));
%     mesh.L = L;
%     mesh.U = U;
% else
%     L = mesh.L;
%     U = mesh.U;
% end

for i = 1 : nsource
%     [x,flag] = bicgstabl(Mass,qvec(:,i),1e-12,100,L,U);
    [x,flag] = bicgstabl(Mass,qvec(:,i),1e-12,100,R',R);
    msg = [msg flag];
    phi(:,i) = x;
end

if any(msg==1)
    disp('some solutions did not converge')
    errordlg('Some solutions did not converge; this could be caused by noisy/bad data or mesh issues','NIRFAST Error');
    error('Some solutions did not converge; this could be caused by noisy/bad data or mesh issues');
elseif any(msg==2)
    disp('some solutions are unusable')
    errordlg('Some solutions are unusable; this could be caused by noisy/bad data or mesh issues','NIRFAST Error');
    error('Some solutions are unusable; this could be caused by noisy/bad data or mesh issues');
elseif any(msg==3)
    disp('some solutions from stagnated iterations')
    errordlg('Some solutions are unusable; this could be caused by noisy/bad data or mesh issues','NIRFAST Error');
    error('Some solutions are unusable; this could be caused by noisy/bad data or mesh issues');
end
