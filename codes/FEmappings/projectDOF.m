function E = projectDOF(N,freenodes)
    % get matrix which projects a vector of values at dof to a larger
    % vector with zeros where Dirichlet nodes are
    %
    % Input:
    %         N:  number of all nodes
    % freenodes:  vector of dof 
    %   
    % Output: 
    %         E:  projection matrix
    %
    % M. Hauck, Y. Liang, D. Peterseim

    % handle logical input vector
    if islogical(freenodes)
        freenodes = find(freenodes);
    end % if
    
    Ndof = length(freenodes);
    E = sparse(freenodes,1:Ndof,1,N,Ndof);
end % function