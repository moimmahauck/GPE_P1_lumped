function [mesh,P,P0,P1dg] = refine(mesh)
    % refines the current triangulation by dividing each element
    % into 2^dim congruent elements (red refinement)
    %
    % Input:
    %     mesh:  current mesh 
    % 
    % Output: 
    %     mesh:  new mesh after refinement
    %        P:  transformation matrix
    %            extrapolates nodal grid functions (continuous) from the old to the refined mesh
    %       P0:  transformation matrix
    %            extrapolates element grid functions from the old to the refined mesh
    %     P1dg:  transformation matrix
    %            extrapolates element P1-functions (discontinuous) from the old to the refined mesh
    %
    % M. Hauck, Y. Liang, D. Peterseim

    % Construct data structure
    [np,d] = size(mesh.p);
    nt = size(mesh.t,1);
    e = computeEdges(mesh);
    ne = size(e,1);
    d2p = sparse(e,e(:,[2,1]),(1:ne)'*[1 1],np,np);

    % New nodes from the mid points of each edge
    newnode = 0.5.*(mesh.p(e(:,1),:)+mesh.p(e(:,2),:)); 
    P = sparse((1:ne)'*[1 1],e,.5,ne,np);
    mesh.p = [mesh.p; newnode]; 
    P = [speye(np);P];
    emarker = (np+1:np+ne)';
    p = zeros(nt,d+1+sum(1:d));

    switch d
        case 1
            p(:,1:2) = mesh.t; 
            p(:,3) = emarker(d2p(mesh.t(:,1)+np*(mesh.t(:,2)-1)));
            mesh.t = [p(:,[1,3])
                      p(:,[3,2])];
            tnew2t = repmat(1:nt,1,2)';
            if nargout > 3
                P1dg = [sparse(1:1:2*nt,1:1:2*nt,repmat([1;.5],nt,1),2*nt,2*nt)+sparse(2:2:2*nt,1:2:2*nt,.5,2*nt,2*nt)
                        sparse(1:1:2*nt,1:1:2*nt,repmat([.5;1],nt,1),2*nt,2*nt)+sparse(1:2:2*nt,2:2:2*nt,.5,2*nt,2*nt)];
            end % if
        case 2
            p(:,1:3) = mesh.t; 
            p(:,4) = emarker(d2p(mesh.t(:,1)+np*(mesh.t(:,2)-1)));
            p(:,5) = emarker(d2p(mesh.t(:,2)+np*(mesh.t(:,3)-1)));
            p(:,6) = emarker(d2p(mesh.t(:,3)+np*(mesh.t(:,1)-1)));
            mesh.t = [p(:,[1,4,6])
                      p(:,[4,2,5])
                      p(:,[6,5,3])
                      p(:,[4,5,6])];
            tnew2t = repmat(1:nt,1,4)';
            if nargout > 3
                E = speye(nt);     
                P1dg = [kron(E,[1 0 0;.5 .5 0;.5 0 .5])
                        kron(E,[.5 .5 0;0 1 0;0 .5 .5])
                        kron(E,[.5 0 .5;0 .5 .5;0 0 1])
                        kron(E,[.5 .5 0;0 .5 .5;.5 0 .5])];
            end % if
        case 3
            p(:,1:4) = mesh.t(:,[1 2 3 4]);
            p(:,5) = emarker(d2p(mesh.t(:,1)+np*(mesh.t(:,2)-1)));
            p(:,6) = emarker(d2p(mesh.t(:,1)+np*(mesh.t(:,3)-1)));
            p(:,7) = emarker(d2p(mesh.t(:,2)+np*(mesh.t(:,3)-1)));
            p(:,8) = emarker(d2p(mesh.t(:,1)+np*(mesh.t(:,4)-1)));
            p(:,9) = emarker(d2p(mesh.t(:,2)+np*(mesh.t(:,4)-1)));
            p(:,10) = emarker(d2p(mesh.t(:,3)+np*(mesh.t(:,4)-1)));
            mesh.t = [p(:,[1 5 6 8])
                      p(:,[5 2 7 9]) 
                      p(:,[6 7 3 10])
                      p(:,[8 9 10 4])
                      p(:,[5 8 9 7])
                      p(:,[5 7 6 8])
                      p(:,[6 8 7 10])
                      p(:,[8 10 9 7])];
            tnew2t = repmat(1:nt,1,8)';
            if nargout > 3
                E = speye(nt);     
                P1dg = [kron(E,[1 0 0 0;.5 .5 0 0;.5 0 .5 0;.5 0 0 .5])
                        kron(E,[.5 .5 0 0;0 1 0 0;0 .5 .5 0;0 .5 0 .5])
                        kron(E,[.5 0 .5 0;0 .5 .5 0;0 0 0 1;0 0 .5 .5])
                        kron(E,[.5 0 0 .5;0 .5 0 .5;0 0 .5 .5;0 0 0 1])
                        kron(E,[.5 .5 0 0;.5 0 0 .5;0 .5 0 .5;0 .5 .5 0])
                        kron(E,[.5 .5 0 0;0 .5 .5 0;.5 0 .5 0;.5 0 0 .5])
                        kron(E,[.5 0 .5 0;.5 0 0 .5;0 .5 .5 0;0 0 .5 .5])
                        kron(E,[.5 0 0 .5;0 0 .5 .5;0 .5 0 .5;0 .5 .5 0])];
            end % if
        otherwise
            error('dimension error')
    end % switch

    P0 = sparse((1:length(tnew2t))',tnew2t,1,length(tnew2t),nt);
end % function
