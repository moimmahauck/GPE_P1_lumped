function [f,f2t,bf]=getBoundaryFaces(mesh)
    % Compute boundary edges/faces of mesh t
    %
    % Input:
    %     mesh:  simplicial mesh 
    %   
    % Output: 
    %        f:  list of faces of t
    %      f2t:  ???
    %       bf:  index vector, contains indizes of faces in f that lie on the 
    %            boundary
    %
    % M. Hauck, Y. Liang, D. Peterseim

    t = mesh.t;
    [nt,dim]=size(t);
    dim=dim-1;

    switch dim
        case 1
            f = [t(:,1); t(:,2)];
            f2t = [(1:nt)' ones(nt,1)*1;
                   (1:nt)' ones(nt,1)*2];
        case 2
            f = [t(:,[1,2]); t(:,[1,3]); t(:,[2,3])];
            f2t = [(1:nt)' ones(nt,1)*[1,2];
                   (1:nt)' ones(nt,1)*[1,3];
                   (1:nt)' ones(nt,1)*[2,3]];
        case 3
            f = [t(:,[1,2,3]); t(:,[1,2,4]); t(:,[1,3,4]); t(:,[2,3,4])];
            f2t = [(1:nt)' ones(nt,1)*[1,2,3];
                   (1:nt)' ones(nt,1)*[1,2,4];
                   (1:nt)' ones(nt,1)*[1,3,4];
                   (1:nt)' ones(nt,1)*[2,3,4]];    
       otherwise
            error('dimension error')
    end % switch
    
    f1=sort(f,2);
    [~,ix1,~]=unique(f1,'rows','first');
    f=f(ix1,:);
    [~,ix2,jx]=unique(f1,'rows','last');
    f2t=[f2t(ix1,:) f2t(ix2,:)];
    vec=histc(jx,1:max(jx));
    [bf,~]=find(vec==1);   
end % function