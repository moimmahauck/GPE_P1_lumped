function v=simpvol(mesh)
    % computes volume of the mesh elements of mesh
    %
    % Input:
    %     mesh:  simplicial mesh 
    %   
    % Output: 
    %      vlo:  volume of mesh elements
    %
    % M. Hauck, Y. Liang, D. Peterseim

    p = mesh.p;
    t = mesh.t;
    switch (size(t,2)-1)
        case 1 % interval in R1
            d12=p(t(:,2),:)-p(t(:,1),:);
            if size(p,2)==1
                v=d12;
            else  % interval in R2 or R3
                v=sqrt(sum(d12.^2,2));
            end % if
        case 2
            d12=p(t(:,2),:)-p(t(:,1),:);
            d13=p(t(:,3),:)-p(t(:,1),:);
            if size(p,2)==2 % triangle in R2
                v=(d12(:,1).*d13(:,2)-d12(:,2).*d13(:,1))/2;
            else % triangle in R3
                n = cross(d12,d13,2);
                v = sqrt(sum(n.^2,2))./2;
            end % if
        case 3
            d12=p(t(:,2),:)-p(t(:,1),:);
            d13=p(t(:,3),:)-p(t(:,1),:);
            d14=p(t(:,4),:)-p(t(:,1),:);
            v=dot(cross(d12,d13,2),d14,2)/6;
        otherwise
            error('dimension error')
    end % switch
end % function