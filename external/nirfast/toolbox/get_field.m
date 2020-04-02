function [phi,R]=get_field(Mass,mesh,qvec)

% [phi,R]=get_field(Mass,mesh,qvec,R)
%
% Used by femdata and jacobian
% Calculates the field, given the Mass matrix and RHS source
% vector.
%
% Mass is the mass matrix
% mesh is the input mesh
% qvec is the RHS source vector
% R is the preconditioner
% phi is the field

if ispref('nirfast','solver')
    
    solver = getpref('nirfast','solver');
    
    % if pardiso is chosen, make sure it's available
    if strcmp(solver,'pardiso')
        env = getenv('LD_LIBRARY_PATH');
        if isempty(findstr(env,'pardiso'))
            pardiso = 0;
        else
            pardiso = 1;
        end
        if exist('pardisolist')
            eval('pardisolist')
            hostname = getComputerName();
            hostname = cellstr(repmat(hostname,length(pardisohost),1));
        else
            hostname = 'a';
            pardisohost = 'b';
        end

        if sum(strcmp(pardisohost,hostname)) == 0 || pardiso == 0
            disp('Pardiso is not available, using Matlab solver');
            solver = 'matlab';
        end
    end
    
else
    
    % determine which solver is best
    env = getenv('LD_LIBRARY_PATH');
    if isempty(findstr(env,'pardiso'))
        pardiso = 0;
    else
        pardiso = 1;
    end
    if exist('pardisolist')
        eval('pardisolist')
        hostname = getComputerName();
        hostname = cellstr(repmat(hostname,length(pardisohost),1));
    else
        hostname = 'a';
        pardisohost = 'b';
    end
    
    if sum(strcmp(pardisohost,hostname)) ~= 0 && pardiso == 1
        solver = 'pardiso';
    elseif length(mesh.nodes) >= 3800
        solver = 'bicgstab';
    else
        solver = 'matlab';
    end
    
end

eval(['[phi,R]=solver_' solver '(Mass,mesh,qvec);']);
