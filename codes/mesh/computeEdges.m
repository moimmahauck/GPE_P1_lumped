function e = computeEdges(mesh)
    % Computes edges of mesh
    %
    % Input:
    %     mesh:  simplicial mesh with points p and elements t
    %   
    % Output: 
    %        e:  list of edges of mesh
    %     nmbe:  index vector how often edges appears. if index is one then 
    %            edge is bondary edge
    %
    % M. Hauck, Y. Liang, D. Peterseim

    [np,d] = size(mesh.p); 
    
    switch d
        case 1
            e = mesh.t;
        case 2
            e = [mesh.t(:,[1,2]); mesh.t(:,[1,3]); mesh.t(:,[2,3])];
        case 3
            e = [mesh.t(:,[1,2]); mesh.t(:,[1,3]); mesh.t(:,[2,3]);...
                 mesh.t(:,[1,4]); mesh.t(:,[2,4]); mesh.t(:,[3,4])];
        otherwise
            error('dimension error')
    end % switch
    
    e = sort(e,2);
    d2p = sparse(e(:,1),e(:,2),1,np,np);
    [e1,e2] = find(d2p);
    e = [e1,e2];
end % function