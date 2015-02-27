function [phi,R]=solver_matlab(Mass,mesh,qvec)

% [phi,R]=solver_matlab(Mass,mesh,qvec)
%
% Used by femdata and jacobian
% Calculates the field, given the Mass matrix and RHS source
% vector.
% Uses the Matlab solver, which uses a method dependent on the
% matrix structure
%
% Mass is the mass matrix
% mesh is the input mesh
% qvec is the RHS source vector
% R is the preconditioner
% phi is the field

[nnodes,nsource]=size(qvec);
phi=zeros(nnodes,nsource);
R = [];
phi = Mass\qvec;