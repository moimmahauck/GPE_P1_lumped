function G = computeGrads(mesh,u)
    % compute piecewise gradient of the grid function u 
    %
    % Input:
    %     mesh:  simplicial mesh
    %        u:  grid function (size(u) = [size(p,1),1])
    %   
    % Output: 
    %        G: gradients (row-wise) according to the simplices in t
    %
    % M. Hauck, Y. Liang, D. Peterseim

    t = mesh.t;
    p = mesh.p;
    nt=size(t,1);
    [~,dim]=size(p);

    switch dim
        case 1
            % Compute vedge: edge as a vector and area of each element
            ve = p(t(:,2),:)-p(t(:,1),:);
            area = abs(ve);

            G(:,:,1)=-sign(ve)./area;
            G(:,:,2)=-G(:,:,1);

        case 2
            % Compute vedge: edge as a vector and area of each element
            ve1 = p(t(:,3),:)-p(t(:,2),:);
            ve2 = p(t(:,1),:)-p(t(:,3),:);
            ve3 = p(t(:,2),:)-p(t(:,1),:);
            area2 = abs(-ve3(:,1).*ve2(:,2)+ve3(:,2).*ve2(:,1));

            G(:,:,1)=[-ve1(:,2) ve1(:,1)]./(area2*[1 1]);
            G(:,:,2)=[-ve2(:,2) ve2(:,1)]./(area2*[1 1]);
            G(:,:,3)=[-ve3(:,2) ve3(:,1)]./(area2*[1 1]);

        case 3
            % Compute ve: edge as a vector and area of each element
            ve1(:,:) = p(t(:,3),:)-p(t(:,2),:);
            ve2(:,:) = p(t(:,4),:)-p(t(:,3),:);
            ve3(:,:) = p(t(:,2),:)-p(t(:,4),:);
            ve4(:,:) = p(t(:,2),:)-p(t(:,1),:);
            ve5(:,:) = p(t(:,3),:)-p(t(:,1),:);
            ve6(:,:) = p(t(:,4),:)-p(t(:,1),:);

            G=zeros(nt,3,4);
            G(:,:,1) = cross(ve2,ve3,2);
            G(:,:,2) = cross(ve5,ve6,2);
            G(:,:,3) = cross(ve4,ve6,2);
            G(:,:,4) = cross(ve4,ve5,2);

            vol = simpvol(mesh);

            G(:,:,1)=G(:,:,1)./(6.*sign(dot(-ve4,G(:,:,1),2)).*vol*ones(1,3));
            G(:,:,2)=G(:,:,2)./(6.*sign(dot(-ve1,G(:,:,2),2)).*vol*ones(1,3));
            G(:,:,3)=G(:,:,3)./(6.*sign(dot(ve5,G(:,:,3),2)).*vol*ones(1,3));
            G(:,:,4)=G(:,:,4)./(6.*sign(dot(ve6,G(:,:,4),2)).*vol*ones(1,3));
        otherwise
            error('dimension error')
    end % switch

    if nargin==3
        for k=1:dim+1
            G(:,:,k)=(u(t(:,k))*ones(1,dim)).*G(:,:,k);
        end % for
        G=sum(G,3);
    end % if
end % function