function [X,names]=getdesign_matrix(data,basis)

if(nargin<2)
    basis=Dictionary;
    basis('default')=nirs.design.basis.Canonical;
end

t       = data.time;
stims   = data.stimulus;

[X, names] = nirs.design. ...
    createDesignMatrix( stims, t, basis );

return