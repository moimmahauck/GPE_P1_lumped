function plotMeshFun(mesh,x,varargin)
    % plots cg1 functions as surface plot, only implemented for 2d!
    %
    % Input:
    %     mesh:  simplicial mesh with points p and elements t
    %        x:  vector of nodal values to plot (np x 1)
    % varargin:  optional input arguments for plotting
    %   
    % no Output
    %
    % M. Hauck, Y. Liang, D. Peterseim

    % throw error if dimension is not two
    d = size(mesh.p,2);
    assert(d == 2,'dimension must be two.');

    trisurf(mesh.t,mesh.p(:,1),mesh.p(:,2),x,x,varargin{:});
end % function