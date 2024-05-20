function M = assembleMass(mesh,coeff,type)
    % calculates possibly weighted mass matrix without quadrature error
    %
    % Input:
    %     mesh:  simplicial mesh with points p and elements t
    %    coeff:  piecewise constant scalar coefficient vector
    %     type:  type of stiffness matrix: cont. Galerkin 
    %            (type = 'cg', default), discont. Galerkin (type = 'dg')
    %   
    % Output: 
    %        M:  mass matrix
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
    np = size(mesh.p,1);
    nt = size(mesh.t,1);
    d = size(mesh.t,2)-1;


    % local mass matrix
    locmass = (ones(d+1)+eye(d+1))./((d+1)*(d+2));

    % create global matrix
    switch type
        case 'dg'
            M = kron(spdiags(simpvol(mesh),0,nt,nt),locmass);
        case 'cg'
            c = coeff.*simpvol(mesh);
            I = (1:(d+1))'*ones(1,d+1);
            M = sparse(mesh.t(:,reshape(I',1,(d+1)^2)),...
                  mesh.t(:,reshape(I,1,(d+1)^2)),c.*reshape(locmass,1,[]),np,np);
    end % switch
end % function