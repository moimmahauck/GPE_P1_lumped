function P = computeDG0toDG1(mesh)
    % get prolongation from dg0 to dg1 in nodal basis
    %
    % Input:
    %     mesh:  current mesh 
    %   
    % Output: 
    %        P:  prolongation matrix
    %
    % M. Hauck, Y. Liang, D. Peterseim

    d = size(mesh.p,2);
    nt = size(mesh.t,1);
    P = kron(speye(nt),ones(d+1,1));   
end % function