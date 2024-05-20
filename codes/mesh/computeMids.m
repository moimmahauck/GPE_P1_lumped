function mids = computeMids(mesh)
    % compute midpoints of simplices of mesh
    %
    % Input:
    %     mesh:  simplicial mesh 
    %   
    % Output: 
    %     mids: midpoints of simplices in t
    %
    % M. Hauck, Y. Liang, D. Peterseim

    d = size(mesh.t,2);

    mids = mesh.p(mesh.t(:,1),:);
    for k=2:d
        mids = mids + mesh.p(mesh.t(:,k),:);
    end % for
    
    mids = mids./d;
end % function