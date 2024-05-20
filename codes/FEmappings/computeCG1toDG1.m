function P = computeCG1toDG1(mesh)
    % get prolongation from cg1 to dg1 in nodal basis
    %
    % Input:
    %     mesh:  current mesh 
    %   
    % Output: 
    %        P:  prolongation matrix
    %
    % M. Hauck, Y. Liang, D. Peterseim

    [np,d] = size(mesh.p);
    nt = size(mesh.t,1);
    P = sparse(1:nt*(d+1),reshape(mesh.t.',1,[]),1,nt*(d+1),np);
end % function