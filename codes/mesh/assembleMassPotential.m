function M = assembleMassPotential(mesh,V,type)
    % calculates weighted mass matrix with a high order quadrature rule
    %
    % Input:
    %     mesh:  simplicial mesh with points p and elements t
    %        V:  potential function as function handle
    %   
    % Output: 
    %        M:  weighted mass matrix
    %
    % M. Hauck, Y. Liang, D. Peterseim
    
    % number of vertices, physical dimension
    [np,d] = size(mesh.p);
    nt = size(mesh.t,1);
    
    % set default type to 'exact'
    if nargin < 3 || isempty(type)
        type = 'exact';
    end % if

    % set default weighting functions to 1
    if nargin < 2 || isempty(V)
        V = @(x) ones(size(x,1),1);
    end % if

    % load quadrature rule
    switch d
        case 1
            lam1 = @(X) 1-X;
            lam2 = @(X) X;

            evalbasisP1 = @(X) [lam1(X) lam2(X)];
            
            if strcmp(type,'lumped')
                X = [0; 1];
                W = [.5; .5];
            else
                X = [.5-.5*sqrt(3/5); .5; .5+.5*sqrt(3/5)];
                W = [5/18; 8/18; 5/18];
            end % if
        case 2
            lam1 = @(X) 1-X(:,1)-X(:,2);
            lam2 = @(X) X(:,1);
            lam3 = @(X) X(:,2);

            evalbasisP1 = @(X) [lam1(X) lam2(X) lam3(X)];
            
            c1 = [0 0]; c2 = [1 0]; c3 = [0 1];

            if strcmp(type,'lumped')
                X = [c1; c2; c3];
                W = [1/3; 1/3; 1/3];
            else
                a1 = 0.445948490915965;
                a2 = 0.091576213509771;

                X = [a1*c1(1)+a1*c2(1)+(1-2*a1)*c3(1) a1*c1(2)+a1*c2(2)+(1-2*a1)*c3(2); ...
                     a1*c1(1)+(1-2*a1)*c2(1)+a1*c3(1) a1*c1(2)+(1-2*a1)*c2(2)+a1*c3(2); ...
                     (1-2*a1)*c1(1)+a1*c2(1)+a1*c3(1) (1-2*a1)*c1(2)+a1*c2(2)+a1*c3(2); ...
                     a2*c1(1)+a2*c2(1)+(1-2*a2)*c3(1) a2*c1(2)+a2*c2(2)+(1-2*a2)*c3(2); ...
                     a2*c1(1)+(1-2*a2)*c2(1)+a2*c3(1) a2*c1(2)+(1-2*a2)*c2(2)+a2*c3(2); ...
                     (1-2*a2)*c1(1)+a2*c2(1)+a2*c3(1) (1-2*a2)*c1(2)+a2*c2(2)+a2*c3(2);];

                w1 = 0.223381589678010;
                w2 = 0.109951743655322;

                W = [w1; w1; w1; w2; w2; w2];
            end % if
        case 3
            error('not implemented')
    end % switch
    
    % local mass matrix
    locmass = zeros(nt,(d+1).^2);

    w = zeros(mesh.nt,length(W));
    
    % figure; plotMeshFun(mesh,zeros(mesh.np,1),'facecolor','none');
    
    for q = 1:length(W)
        qx = (mesh.p(mesh.t(:,2),:)-mesh.p(mesh.t(:,1),:)).*X(q,1) + ...
             (mesh.p(mesh.t(:,3),:)-mesh.p(mesh.t(:,1),:)).*X(q,2) + ...
              mesh.p(mesh.t(:,1),:);
        % hold on; plot(qx(:,1),qx(:,2),'.r');
        w(:,q) = V(qx);
    end % for

    % evaluate basis and weight functions in gauss points
    F = evalbasisP1(X).';
    
    % loop over local pairs of dofs (parallel in t)
    vol = simpvol(mesh);
    for k1 = 1:d+1
        for k2 = 1:d+1
            locmass(:,k1+(d+1)*(k2-1)) = vol.*(w*(W.*(F(k1,:).*F(k2,:)).'));
        end % for
    end % for

    % create global matrix
    I = (1:(d+1))'*ones(1,d+1);       
    M = sparse(mesh.t(:,reshape(I',1,(d+1)^2)),...
               mesh.t(:,reshape(I,1,(d+1)^2)),locmass,np,np);
end % function