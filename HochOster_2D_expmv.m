%Hochburck-Osterman Equation in 2D+time with expmv%
mkdir(pwd,'Results_HochOster_Expmv')
mypath = fullfile(pwd,'Results_HochOster_Expmv');

rmax = 8; %Maximum elements in each space direction: 2^r
mmax = 9; %Maximum time step size: 2^m

for r = 8:rmax 
    
    %Initialize vectors for the results
    Errors1 = zeros(mmax,1);
    Errors2 = zeros(mmax,1);
    Errors3 = zeros(mmax,1);
    times_Euler = zeros(mmax,1);
    times_RK2 = zeros(mmax,1);
    times_RK3 = zeros(mmax,1);
    taum = zeros(mmax,1);
    nsteps = zeros(mmax,1);
    
    P = 3; %Hochburck-Osterman equation
    %Import the data of the problem
    [eps,beta,u0,x1,x2,y1,y2,T] = data(P);
    uexact = @(x,y,t) exp(t).*x.*(1-x).*y.*(1-y);
    source = @(x,y,t) uexact(x,y,t)+2*exp(t).*(x.*(1-x)+y.*(1-y))-1./(1+exp(2*t).*(x.*(1-x).*y.*(1-y)).^2);
    
    %Meshes and Parameters
    [xsol,ysol,nelx,nely] = mesh2D(r,x1,x2,y1,y2,P);
    [nx,ny,nel,nnode,coord,nodes] = parameters(xsol,ysol);
    nquad = 2; %quadrature points (Lobatto quadrature)
    nquadf = 2; %number of quadrature points to integrating the source term
    
    %Full sparse matrix
    dir = BC(nx,ny,P);
    [A,M] = FEM_matrices_sparse(eps,beta,dir,nel,coord,nodes,nquad);
    
    %Initial condition vector
    [xmat,ymat] = meshgrid(xsol,ysol);
    U0 = u0(xmat,ymat);
    U0 = reshape(U0',[],1);
    U0(dir) = [];
    dimx = size(U0,1);
    
    UexactT = uexact(xmat,ymat,T);
    UexactT = reshape(UexactT',[],1);
    UexactT(dir) = [];
    
    for m = 1:mmax %tau=T/2^m time step size  
    
        %Time grid and step size
        steps = 2^m;
        nsteps(m) = steps;
        tau = T/steps;
        taum(m) = tau;
        t = 0:tau:T;
    
        %Assembly the linear source term
        fprintf('Assembly of the RHS \n\n')
        flin = zeros(dimx,steps+1);
        flin12 = zeros(dimx,steps);
        flin13 = zeros(dimx,steps);
        flin23 = zeros(dimx,steps);
        for i = 1:steps
            flin(:,i) = source_term(source,t(i),M,nel,nnode,coord,nodes,nquadf,dir);
            flin12(:,i) = source_term(source,t(i)+tau/2,M,nel,nnode,coord,nodes,nquadf,dir);
            flin13(:,i) = source_term(source,t(i)+tau/3,M,nel,nnode,coord,nodes,nquadf,dir);
            flin23(:,i) = source_term(source,t(i)+2*tau/3,M,nel,nnode,coord,nodes,nquadf,dir);
        end
        flin(:,end) = source_term(source,t(end),M,nel,nnode,coord,nodes,nquadf,dir);
            
        %%%Exponential Euler method%%%
        fprintf('Running Exponential Euler for r=%d and %d time steps\n\n',r,steps)
        U = zeros(dimx,steps+1);
        U(:,1) = U0;
        tic
        for i = 1:steps
            f1 = 1./(1+U(:,i).^2)+flin(:,i);
            U(:,i+1) = U(:,i)+tau*phiB(-tau*A,f1-A*U(:,i));
        end
        times_Euler(m) = toc;
        Errors1(m) = norm(UexactT-U(:,end),'inf');
    
        %RK2 method%%%
        fprintf('Running RK2 for r=%d and %d time steps\n\n',r,steps)
        U = zeros(dimx,steps+1);
        U(:,1) = U0;
        tic
        for i = 1:steps
            %stage 1
            f1 = 1./(1+U(:,i).^2)+flin(:,i);
            %stage 2
            Us1 = U(:,i)+(tau/2)*phiB(-(tau/2)*A,f1-A*U(:,i));
            f2 = 1./(1+Us1.^2)+flin12(:,i);
            %compute solution
            U(:,i+1) = U(:,i)+tau*phiB(-tau*A,f2-A*U(:,i));
        end
        times_RK2(m) = toc;
        Errors2(m) = norm(UexactT-U(:,end),'inf');
     
        %third-order method%%%
        fprintf('Running RK3 for r=%d and %d time steps\n\n',r,steps)
        U = zeros(dimx,steps+1);
        U(:,1) = U0;
        tic
        for i = 1:steps
            %stage 1
            f1 = 1./(1+U(:,i).^2)+flin(:,i);
            %stage 2
            Us1 = U(:,i)+(tau/3)*phiB(-tau*A/3,f1-A*U(:,i));
            f2 = 1./(1+Us1.^2)+flin13(:,i);
            %stage 3
            Us2 = U(:,i)+tau*phiB(-2*tau*A/3,[(4/3)*(-f1+f2) (2/3)*(f1-A*U(:,i))]);
            f3 = 1./(1+Us2.^2)+flin23(:,i);
            %compute solution
            U(:,i+1) = U(:,i)+tau*phiB(-tau*A,[(3/2)*(-f1+f3) f1-A*U(:,i)]);
        end
        times_RK3(m) = toc;
        Errors3(m) = norm(UexactT-U(:,end),'inf');
    end
    
    %Save errors and computational times
    Tab_times = table(nsteps,times_Euler,times_RK2,times_RK3,'VariableNames',{'Nsteps','Times_EulerExpmv','Times_RK2Expmv','Times_RK3Expmv'});
    writetable(Tab_times,fullfile(mypath,['TimesHochOster_Expmv_r' num2str(r) '_m' num2str(m) '.txt']),'Delimiter',' ');
    
    Tab_errors = table(taum,Errors1,Errors2,Errors3,'VariableNames',{'tau','Err_EulerExpmv','Err_RK2Expmv','Err_RK3Expmv'});
    writetable(Tab_errors,fullfile(mypath,['ErrHochOster_Expmv_r' num2str(r) '_m' num2str(m) '.txt']),'Delimiter',' ');

end


%Function to compute the linear combination of actions with expmv
function out = phiB(A,B)
    q=size(B,2);
    if q == 1
        At = [A, B; zeros(q, length(A)) 0];
        if issparse(A)
            At = sparse(At);
        end  
    else
        if issparse(A)
            At = sparse([A, B; sparse(zeros(q,length(A))), sparse(diag(ones(q-1,1),1))]);
        else
            At = [A, B; (zeros(q,length(A))), (diag(ones(q-1,1),1))];
        end
    end
    e = zeros(length(At),1); e(length(A)+q) = 1;
    out = expmv(1,At,e);
    out = out(1:length(A));
end
