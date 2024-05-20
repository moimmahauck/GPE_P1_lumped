function plotDiscFun(mesh,x,varargin)
    % plots dg0/dg1 functions as surface plot, only implemented for 2d!
    %
    % Input:
    %     mesh:  simplicial mesh with points p and elements t
    %        x:  vector of nodal values to plot (nt x 1 or nt x d+1)
    % varargin:  optional input arguments for plotting
    %   
    % no Output
    %
    % M. Hauck, Y. Liang, D. Peterseim

    % throw error if dimension is not two
    d = size(mesh.p,2);
    assert(d == 2,'dimension must be two.');

    % P0 input
    if size(x,2) == 1
        x = repmat(x,1,3);
    end % if

    % coordinates for all nodes of sides
    X1 = mesh.p(mesh.t(:,1), 1);
    X2 = mesh.p(mesh.t(:,2), 1);
    X3 = mesh.p(mesh.t(:,3), 1);
    Y1 = mesh.p(mesh.t(:,1), 2);
    Y2 = mesh.p(mesh.t(:,2), 2);
    Y3 = mesh.p(mesh.t(:,3), 2);
    
    % quads for the function in patch-style
    valX = [X1';X2';X3'];
    valY = [Y1';Y2';Y3'];
    
    % nodal values at each element
    valZ = [x(:,1)';x(:,2)';x(:,3)'];

    patch(valX,valY,valZ,valZ,varargin{:});
end