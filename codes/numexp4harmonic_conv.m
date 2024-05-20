%% This code reproduces the first numerical experiment in the paper
%
% Note: the below parameter configuration is suitable for a laptop.
%
% M. Hauck, Y. Liang, D. Peterseim

init
close all
clear all

%% paramters
% domain
L = 8; % domain is [-L,L]^2

% physical parameters
kappa = 1000; % particle itneraction
potential = @(x) .5*(x(:,1).^2 + x(:,2).^2); % potential function

% discretization and solver parameters
lvlhs = 1:8; % refinement level for discretization
lvlhh = lvlhs(end)+2; % refinement level for reference solution
tol = 1e-12; % tolerance for nonlinear sovler
maxit = 1000; % maximal number of iterations for nonlinear solver

%% compute reference solution
fprintf('computing reference solution.\n');

% create meshes
T0 = scaleMesh(getMesh('Square'),L);
Thh = refineMesh(T0,lvlhs(end)+2);

% assemble matrices
Shh = assembleStiffness(Thh);
Mhh = assembleMassPotential(Thh,[],'exact');
MVhh = assembleMassPotential(Thh,potential,'exact');
Ahh = Shh + MVhh;

% get Dirichlet nodes and dof 
Dnodes = max(abs(Thh.p),[],2) > L - 1e-10;
dof = ~Dnodes;

% restrict matrices to dof
Shh0 = Ahh(dof,dof); 
Mhh0 = Mhh(dof,dof);

% assemble dof projection matrix
Ehh = projectDOF(Thh.np,dof); 

% apply energy-adapted Riemannian gradient descent
[uhh,flag,lamhh,erghh,reshh] = Amethod4GPE(Thh,Shh0,Mhh0,Ehh,kappa,ones(size(Mhh,1),1),[0 1],tol,maxit,'exact');

%% compute errors for exact approximation
% initialize error vectors
err_uh_l2_exact = zeros(1,length(lvlhs));
err_uh_H1_exact = zeros(1,length(lvlhs));
err_lam_exact = zeros(1,length(lvlhs));
err_erg_exact = zeros(1,length(lvlhs));

for indlvlh = 1:length(lvlhs)
    lvlh = lvlhs(indlvlh);
    Th = refineMesh(T0,lvlh);

    fprintf(['computing solution for lvlh = ' num2str(lvlh) '.\n']);

    % assemble matrices
    Sh = assembleStiffness(Th);
    Mh = assembleMassPotential(Th,[],'exact');
    MVh = assembleMassPotential(Th,potential,'exact');
    Ah = Sh + MVh;

    % get Dirichlet nodes and dof 
    Dnodes = max(abs(Th.p),[],2) > L - 1e-10;
    dof = ~Dnodes;

    % restrict matrices to dof
    Sh0 = Ah(dof,dof); 
    Mh0 = Mh(dof,dof);

    % assemble dof projection matrix
    Eh = projectDOF(Th.np,dof); 

    % apply energy-adapted Riemannian gradient descent
    [uh,flag,lamh,ergh,resh] = Amethod4GPE(Th,Sh0,Mh0,Eh,kappa,ones(size(Mh,1),1),[0 1],tol,maxit,'exact');
    
    % prolongation to Thh
    [~,P] = refineMesh(Th,lvlhh-lvlh);
    e = P*uh - uhh;
    
    % compute errors
    err_uh_l2_exact(indlvlh) = sqrt(e'*(Mhh*e))/sqrt(uhh'*(Mhh*uhh));
    err_uh_H1_exact(indlvlh) = sqrt(e'*(Shh*e))/sqrt(uhh'*(Shh*uhh));
    err_lam_exact(indlvlh) = abs(lamh(end) - lamhh(end))/lamhh(end);
    err_erg_exact(indlvlh) = abs(ergh(end) - erghh(end))/erghh(end);
end % for

%% compute errors for lumped approximation
% initialize error vectors
err_uh_l2_lumped = zeros(1,length(lvlhs));
err_uh_H1_lumped = zeros(1,length(lvlhs));
err_lam_lumped = zeros(1,length(lvlhs));
err_erg_lumped = zeros(1,length(lvlhs));

for indlvlh = 1:length(lvlhs)
    lvlh = lvlhs(indlvlh);
    Th = refineMesh(T0,lvlh);

    fprintf(['computing solution for lvlh = ' num2str(lvlh) '.\n']);

    % assemble matrices
    Sh = assembleStiffness(Th);
    Mh = assembleMassPotential(Th,[],'lumped');
    MVh = assembleMassPotential(Th,potential,'lumped');
    Ah = Sh + MVh;

    % get Dirichlet nodes and dof 
    Dnodes = max(abs(Th.p),[],2) > L - 1e-10;
    dof = ~Dnodes;

    % restrict matrices to dof
    Sh0 = Ah(dof,dof); 
    Mh0 = Mh(dof,dof);

    % assemble dof projection matrix
    Eh = projectDOF(Th.np,dof); 

    % apply energy-adapted Riemannian gradient descent
    [uh,flag,lamh,ergh,resh] = Amethod4GPE(Th,Sh0,Mh0,Eh,kappa,ones(size(Mh,1),1),[0 1],tol,maxit,'lumped');
    
    % prolongation to Thh
    [~,P] = refineMesh(Th,lvlhh-lvlh);
    e = P*uh - uhh;
    
    % compute errors
    err_uh_l2_lumped(indlvlh) = sqrt(e'*(Mhh*e))/sqrt(uhh'*(Mhh*uhh));
    err_uh_H1_lumped(indlvlh) = sqrt(e'*(Shh*e))/sqrt(uhh'*(Shh*uhh));
    err_lam_lumped(indlvlh) = abs(lamh(end) - lamhh(end))/lamhh(end);
    err_erg_lumped(indlvlh) = abs(ergh(end) - erghh(end))/erghh(end);
end % for

%% visualize errors
f1 = figure('position',[100,100,450,450]);
loglog(2.^(-1:-1:-length(lvlhs)),real(err_uh_H1_exact),'.-','LineWidth',1,'MarkerSize',10,'DisplayName', '$h\mapsto \|\nabla(u-u_h^{\mathcal{P}^1})\|_{L^2}/\|\nabla u\|_{L^2}$')
hold on
loglog(2.^(-1:-1:-length(lvlhs)),real(err_uh_H1_lumped),'.-','LineWidth',1,'MarkerSize',10,'DisplayName', '$h\mapsto \|\nabla(u-u_h^{\mathcal{P}^1\mathrm{-lumped}})\|_{L^2}/\|\nabla u\|_{L^2}$')
loglog(2.^(-1:-1:-length(lvlhs)),real(err_uh_l2_exact),'.-','LineWidth',1,'MarkerSize',10,'DisplayName', '$h\mapsto \|u-u_h^{\mathcal{P}^1}\|_{L^2}/\|u\|_{L^2}$')
loglog(2.^(-1:-1:-length(lvlhs)),real(err_uh_l2_lumped),'.-','LineWidth',1,'MarkerSize',10,'DisplayName', '$h\mapsto \|u-u_h^{\mathcal{P}^1\mathrm{-lumped}}\|_{L^2}/\|u\|_{L^2}$')
loglog(2.^(-1:-1:-length(lvlhs)),2.^(-1:-1:-length(lvlhs)),'--k','DisplayName','$h\mapsto h$')
loglog(2.^(-1:-1:-length(lvlhs)),(2.^(-1:-1:-length(lvlhs))).^2,'-.k','DisplayName','$h\mapsto h^2$')
xticks(flip(2.^(-1:-1:-length(lvlhs))));
xlim([2^-8,2^-1]);
xticklabels({'2^{-8}','2^{-7}','2^{-6}','2^{-5}','2^{-4}','2^{-3}','2^{-2}','2^{-1}'});
axis square
hl = legend('show');
set(hl, 'Interpreter','latex','Location','SouthEast')
% filename = 'figures/harmonic_err_u.png';
% exportgraphics(f1,filename);
% filename = 'figures/harmonic_err_u.tex';
% matlab2tikz(filename,'standalone',true);

f2 = figure('position',[100,100,450,450]);
loglog(2.^(-1:-1:-length(lvlhs)),real(err_erg_exact),'.-','LineWidth',1,'MarkerSize',10,'DisplayName', '$h\mapsto |E-E_h^{\mathcal{P}^1}|/|E|$')
hold on 
loglog(2.^(-1:-1:-length(lvlhs)),real(err_erg_lumped),'.-','LineWidth',1,'MarkerSize',10,'DisplayName', '$h\mapsto |E-E_h^{\mathcal{P}^1\mathrm{-lumped}}|/|E|$')
loglog(2.^(-1:-1:-length(lvlhs)),real(err_lam_exact),'.-','LineWidth',1,'MarkerSize',10,'DisplayName', '$h\mapsto |\lambda-\lambda_h^{\mathcal{P}^1})|/|\lambda|$')
loglog(2.^(-1:-1:-length(lvlhs)),real(err_lam_lumped),'.-','LineWidth',1,'MarkerSize',10,'DisplayName', '$h\mapsto |\lambda-\lambda_h^{\mathcal{P}^1\mathrm{-lumped}}|/|\lambda|$')
loglog(2.^(-1:-1:-length(lvlhs)),(2.^(-1:-1:-length(lvlhs))).^2,'-.k','DisplayName','$h\mapsto h^2$')
xticks(flip(2.^(-1:-1:-length(lvlhs))));
xlim([2^-8,2^-1]);
xticklabels({'2^{-8}','2^{-7}','2^{-6}','2^{-5}','2^{-4}','2^{-3}','2^{-2}','2^{-1}'});
axis square
hl = legend('show');
set(hl, 'Interpreter','latex','Location','SouthEast')
% filename = 'figures/harmonic_err_lamerg.png';
% exportgraphics(f1,filename);
% filename = 'figures/harmonic_err_lamerg.tex';
% matlab2tikz(filename,'standalone',true);

% plot reference solution
% figure; plotMeshFun(Thh,uhh,'edgecolor','none')