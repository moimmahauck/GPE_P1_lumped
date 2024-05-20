function P = computeCGtoDG(mesh,ord)
    % get prolongation/projection from cg1 to dg0/dg1
    %
    % Input:
    %     mesh:  current mesh 
    %      ord:  order of dg space of projection/prolongation
    %   
    % Output: 
    %        P:  prolongation/projection matrix
    %
    % M. Hauck, Y. Liang, D. Peterseim

    switch ord
        case 0
            d = size(mesh.p,2);
            P = computeDG0toDG1(mesh).'*computeCG1toDG1(mesh)./(d+1);
        case 1
            P = computeCG1toDG1(mesh);
    end % switch
end % function