function mesh = getMesh(geom)
    % get initial mesh
    %
    % Input:
    %     geom:  geometry of mesh to load 
    %   
    % Output: 
    %     mesh:  initial mesh
    %
    % M. Hauck, Y. Liang, D. Peterseim

    switch geom
        case 'Interval'
            p = [0;1];
            t = [1 2];
        case 'Square'
            p = [0 0; 1 0; 1 1; 0 1];
            t = [1 2 3; 1 3 4];
        case 'Cube'
            p = [0 0 0; 0 1 0; 1 0 0; 1 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 1];
            t = [6 2 3 1; 6 3 5 1; 6 2 4 3; 6 4 8 3; 6 7 5 3; 6 8 7 3];
        otherwise
            error('geometry not found.');
    end % switch

    % create mesh
    mesh = struct('p',p,'t',t,'np',size(p,1),'nt',size(t,1));
end % function