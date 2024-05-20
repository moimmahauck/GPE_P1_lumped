function [u,flag,lam,erg,res] = Amethod4GPE(T,A0,M0,E,kappa,u0,t,tol,maxit,type)
    % performs Amethod which is a Sobolev gradient flow method with a
    % weighted metric (see [HP20]) to practically solve the discrete 
    % minimization for the discrete ground state
    %
    % Input:
    %        T:  simplicial mesh with points p and elements t
    %       A0:  stiffness matrix restricted to interior dof
    %       M0:  mass matrix restricted to interior dof
    %        E:  dof restriction matrix
    %    kappa:  particle interaction parameter of GPE
    %       u0:  initial iterate (needs to be chosen non-negative if 
    %            one seeks the ground state)
    %        t:  range of step sizes (needs to be chosen [0 1] to guarantee
    %            the non-negativity of the iterates)
    %      tol:  relative tolerance for termination criterion 
    %    maxit:  maximum number of iterations
    %     type:  'lumped' for the lumped P1FEM and 'exact' for the P1FEM
    %   
    % Output: 
    %        u:  discrete ground state
    %     flag:  convergence flag
    %      lam:  vector of eigenvalue approximations
    %      erg:  vector of energy approximations
    %      res:  vector of relative residuals of the approximations
    %
    % M. Hauck, Y. Liang, D. Peterseim

    d = size(T.p,2);
    
    % initial guess
    if isempty(u0)
        u0 = E*rand(size(A0,1),1);
    end
    u0 = u0/sqrt((E'*u0)'*(M0*(E'*u0)));

    % initialize output
    lam = zeros(maxit,1);
    normrelres = zeros(maxit,1);
    erg = zeros(maxit,1);
    flag = 1;

    % assemble initial system matrix
    u = E'*u0;
    PCGtoDG = computeCGtoDG(T,1);
    Muu0 = E'*(assembleMassNonlinear(T,reshape(PCGtoDG*u0,d+1,[]).',reshape(PCGtoDG*conj(u0),d+1,[]).',type)*E);

    erg(1) = real(0.5*u'*(A0+0.5*kappa*Muu0)*u);
    lam(1) = real(u'*(A0+kappa*Muu0)*u);
    res = (A0+kappa*Muu0)*u-lam(1)*M0*u;
    normrelres(1) = sqrt(real(res'*M0*res))/sqrt(real(u'*M0*u)); 
    GRAD = u*0;          % old search direction

    % parameter
    autot = false;
    if isempty(t)
        tmin = 0.1;
        tmax = 2;
        t = 1;
    else
        tmin = t(1);
        tmax = t(2);
        t = 0.5*(tmin+tmax);
        if tmax-tmin>tol
            autot = true;
        end
        alpha    = 0.95;
        sigma     = 1e-4;
        delta    = 0.5;
        gamma0   = 0.01;
        gammamin = 1e-4;
        gammamax = 1.8;
        q = 1;
        cErg = erg(1);
        maxnmlsiter = 20; 
    end

    % display
    k = 1;
    fprintf('\n----------------------------------------------------------------------------\n A-Method, delta = %g, tau = [%g,%g] \n----------------------------------------------------------------------------\n   k  t  \t  energy\t lambda\t     relresidual\n----------------------------------------------------------------------------\n',kappa,tmin,tmax);
    fprintf('%4i  %7.5e  %7.5e  %7.5e  %7.5e\n',k-1,0,erg(k),lam(k),normrelres(k))

    % iteration
    for k = 2:maxit+1
        uold = u;
        GRADold = GRAD;

        % compute u = A\uold
        A = A0+kappa*Muu0;
        rhs = (M0*uold);   
        GRAD = A\rhs;
        GRAD = (uold - GRAD./(GRAD'*rhs)); 
        dGRAD = GRAD - GRADold;

        % adaptive time step (quadratic line search)
        if autot
            %% compute step size
            if k == 2
                gamma = gamma0;
            else
                if ~mod(k,2) 
                    gamma = (dPHI'*dPHI)./abs(dPHI'*dGRAD);
                else
                    gamma = abs((dPHI'*dGRAD))./(dGRAD'*dGRAD);
               end
            end
            gamma = max([gammamin,min([gamma,gammamax])]);
            for kk = 0:maxnmlsiter
                t = gamma.*delta^kk;
                dPHI = -GRAD*t; % difference of the iterates
                %% retraction
                u1 = uold + dPHI;
                u1 = u1./sqrt(u1'*(M0*u1));

                c = E*u1;
                Muu0 = E'*(assembleMassNonlinear(T,reshape(PCGtoDG*c,d+1,[]).',reshape(PCGtoDG*conj(c),d+1,[]).',type)*E);
                % Muu0 = E'*(assembleMass(Th,(PCGtoDG*c).^2)*E);

                Erg = real(0.5*u1'*(A0+0.5*kappa*Muu0)*u1);

                if Erg <= cErg - sigma.*t.*GRAD'*(A0*GRAD)
                    break
                end
            end
        else
            dPHI = -GRAD*t; % difference of the iterates
            %% retraction
            u1 = uold + dPHI;
            u1 = u1./sqrt(u1'*(M0*u1));

            c = E*u1;
            Muu0 = E'*(assembleMassNonlinear(T,reshape(PCGtoDG*c,d+1,[]).',reshape(PCGtoDG*conj(c),d+1,[]).',type)*E);
            % Muu0 = E'*(assembleMass(Th,(PCGtoDG*c).^2)*E);

            Erg = real(0.5*u1'*(A0+0.5*kappa*Muu0)*u1);
        end

        % update u and compute lam and erg
        u = u1;
        lam(k) = real(u'*(A0+kappa*Muu0)*u);
        erg(k) = Erg;

        % compute L^2-norm of residual
        res = (A0+kappa*Muu0)*u-lam(k)*M0*u;
        normrelres(k) = sqrt(real(res'*M0*res))/sqrt(real(u'*M0*u)); 

        if autot
            q = alpha*q+1;
            cErg = (1-1/q)*cErg + 1/q*erg(k);
        end

        % display
        fprintf('%4i  %7.5e %7.5e  %7.5e  %7.5e\n',k-1,t,erg(k),lam(k),normrelres(k))
        
        if k >= 2 && (erg(k-1)-erg(k))/erg(k) < tol && normrelres(k) < tol
            flag = 0;
        end

        if ~flag
            break;
        end
    end
    if sum(u)<0
        u = -u;
    end
    u = E*u;
    lam = lam(1:k);
    erg = erg(1:k);
    res = normrelres(1:k);
end % function