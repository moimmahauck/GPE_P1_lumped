function [Th,P,P0,P1dg] = refineMesh(TH,nref)
    % refines a coarse triangulation TH nref times by applying the function
    % refine nref times. Suitable restriction and prolongation operator
    % between corresponding mesh functions are computed
    %
    % Input:
    %       TH:  coarse mesh (structure array consisting of array of vertices p
    %            and array of elements represented by vertex indices
    %     nref:  number of refinments [default = 1]
    % 
    % Output: 
    %       Th:  fine mesh after refinement
    %        P:  Prolongation matrix P1
    %            extrapolates P1 grid functions from the coarse to the fine mesh
    %       P0:  Prolongation matrix P0
    %            extrapolates P0 mesh functions from the coarse to the fine mesh
    %     P1dg:  Prolongation matrix P1dg
    %            extrapolates P1dg mesh functions from the coarse to the fine mesh
    %
    % M. Hauck, Y. Liang, D. Peterseim

    % Construct data structure
    if nargin < 2
        nref = 1;
    end % if
    
    % compute fine reference mesh
    P = speye(TH.np);
    P0 = speye(TH.nt);
    if nargout > 3
        P1dg = speye(TH.nt*size(TH.t,2));
    end % if
    
    Th = TH;

    for k = 1:nref
        if nargout > 3
            [Th,p,p0,p1dg] = refine(Th);
            P1dg = p1dg*P1dg;
        else
            [Th,p,p0] = refine(Th);
        end % if
        P = p*P;
        P0 = p0*P0;
    end % for
    
    Th.nt = size(Th.t,1);
    Th.np = size(Th.p,1);
end % function