function  mesh = scaleMesh(mesh,L)
    % center mesh at origin and scale it by 2*L
    %
    % Input:
    %     mesh:  mesh of unit cube [0,1]^d 
    %   
    % Output: 
    %     mesh:  scaled mesh of the cube [-L,L]^d
    %
    % M. Hauck, Y. Liang, D. Peterseim

    centerandscale = @(x) (x-.5).*2*L; 
    mesh.p = centerandscale(mesh.p);
end % function
