function [phi,R]=solver_pardiso(MASS,mesh,qvec)

% [phi,R]=solver_pardiso(Mass,mesh,qvec)
%
% Used by femdata and jacobian
% Calculates the field, given the Mass matrix and RHS source
% vector.
% Uses the Pardiso solver (you must have it configured on
% your computer/cluster)
% solvertype should 
% be 0 for direct solve and 1 for iterative solve.
% Visit www.pardiso-project.org for more details on the solver.
%
% Mass is the mass matrix
% mesh is the input mesh
% qvec is the RHS source vector
% R is the preconditioner
% phi is the field

[nnodes,nsource]=size(qvec);
phi=zeros(nnodes,nsource);
solvertype = 0;

fullqvec = full(qvec); %NIRFAST stores the RHS matrix as sparse but PARDISO wants it as full
info = pardisoinit(13,solvertype);

% set num of comp threads explicitly to workaround some situations where
% OMP_NUM_THREADS is not honored
info.iparm(3) = 8;
 
info = pardisoreorder(MASS,info,false);
[fullphi, info] = pardisosolve(MASS,fullqvec,info,false);
phi = sparse(fullphi); %PARDISO will return a full solution matrix but NIRFAST wants it sparse

pardisofree(info);
clear info;
R = [];