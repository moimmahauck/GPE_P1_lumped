function S = assembleStiffness(mesh,coeff,type)
    % calculates possibly weighted stiffness matrix without quadrature error
    %
    % Input:
    %     mesh:  simplicial mesh with points p and elements t
    %    coeff:  piecewise constant scalar coefficient vector
    %     type:  type of stiffness matrix: cont. Galerkin 
    %            (type = 'cg', default), discont. Galerkin (type = 'dg')
    %   
    % Output: 
    %        S:  stiffness matrix
    %
    % M. Hauck, Y. Liang, D. Peterseim

    % set default coefficient to 1
    if nargin < 2 || isempty(coeff)
        coeff = ones(size(mesh.t,1),1);
    end % if

    % set default type to cg
    if nargin < 3
        type = 'cg';
    end % if

    % number of vertices, physical dimension
    [np,d] = size(mesh.p);
    nt = size(mesh.t,1);

    % local stiffness matrices
    g = computeGrads(mesh);
    c = coeff.*simpvol(mesh);
    locstiffness = zeros(nt,(d+1)^2);
    for j = 1:d+1
        for k = 1:d+1
            locstiffness(:,k+(d+1)*(j-1)) = c.*dot(g(:,:,j),g(:,:,k),2);
        end % for
    end % for

    % create global matrix
    I = (1:(d+1))'*ones(1,d+1);
    switch type
        case 'dg'
            ndof = (d+1)*nt;
            idx = reshape((1:ndof),d+1,[])';
            S = sparse(idx(:,reshape(I',1,(d+1)^2)),...
                        idx(:,reshape(I,1,(d+1)^2)),locstiffness,ndof,ndof);
        case 'cg'
            S = sparse(mesh.t(:,reshape(I',1,(d+1)^2)),...
                        mesh.t(:,reshape(I,1,(d+1)^2)),locstiffness,np,np);
    end % switch
end % function