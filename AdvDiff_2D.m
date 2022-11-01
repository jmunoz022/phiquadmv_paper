%Eriksson-Johnson problem 2D+time%
mkdir(pwd,'Results_AdvDiff')
mypath=fullfile(pwd,'Results_AdvDiff');

tol = 1.0e-14;
As = cell(2,1);

P = 2; %select problem P=1 (Heat) P=2 (Eriksson-Johnson)

for m = 3:3 %tau=T/2^m time step size
    for r = 5:6 %2^r elements in each space direction
    
        %Compute 1D and 2D matrices
        [As{1},As{2},A,b] = compute_matrices(P,r,2^m);

        sizeA = size(A,1);
        normA = norm(A,'inf');
        
        q = 20;
       
        times_expmv = zeros(1,q);
        error_gauss = zeros (1,q);
        error_cheb = zeros(1,q);

        tic
        [phis_gauss,npts_gauss,err_gauss,cost_gauss,l_gauss] = phiquadmv(1:q,As,b,normA,-1,'gauss',tol); 
        times_gauss = toc;
        fprintf('Gauss, npts: %d, nevals: %d, scaling: %d\n\n', npts_gauss, cost_gauss, l_gauss)
        
        tic
        [phis_cheb,npts_cheb,err_cheb,cost_cheb,l_cheb] = phiquadmv(1:q,As,b,normA,-1,'chebadaptive',tol); 
        times_cheb = toc;
        fprintf('Chebyshev, npts: %d, nevals: %d, scaling: %d\n\n', npts_cheb, cost_cheb, l_cheb)
        
        for i = 1:q 
            tic
            ex = phi(i,A,b);
            times_expmv(i)=toc;
        
            %compute relative errors
            error_gauss(i) = norm(ex-phis_gauss(:,i),'inf')/norm(ex,'inf');
            error_cheb(i) = norm(ex-phis_cheb(:,i),'inf')/norm(ex,'inf');
            fprintf('phi_%d - error gauss: %e, error cheb: %e\n\n', i, error_gauss(i), error_cheb(i))
            %fprintf('# of time steps=%d, # of elements=%d, varphi order=%d\n\n', 2^m,2^r,i)
        end
        
        %Save errors and computational times
        Tab_timesExpmv=table((1:q)',times_expmv','VariableNames',{'q','Times_expmv'});
        writetable(Tab_timesExpmv,fullfile(mypath,['TimesExpmvAdvDiff_r' num2str(r) '_m' num2str(m) '.txt']),'Delimiter',' ');
        
        Tab_times = table(sizeA,times_gauss,times_cheb,sum(times_expmv),'VariableNames',{'sizeA','Time_gauss','Time_cheb','Time_expmv'});
        writetable(Tab_times,fullfile(mypath,['TimesAdvDiff_r' num2str(r) '_m' num2str(m) '.txt']),'Delimiter',' ');
        
        Tab_errors = table((1:q)',error_gauss',error_cheb','VariableNames',{'q','Error_gauss','Error_cheb'});
        writetable(Tab_errors,fullfile(mypath,['ErrAdvDiff_r' num2str(r) '_m' num2str(m) '.txt']),'Delimiter',' ');

        Tab_npts = table(sizeA,normA,npts_gauss,cost_gauss,l_gauss,npts_cheb,cost_cheb,l_cheb,'VariableNames',{'sizeA','normA','N_gauss','Cost_gauss','Sc_gauss','N_cheb','Cost_cheb','Sc_cheb'});
        writetable(Tab_npts,fullfile(mypath,['NptsAdvDiff_r' num2str(r) '_m' num2str(m) '.txt']),'Delimiter',' ');
    end
end

%Function to compute the actions with expmv
function out = phi(q,A,b)
    if q == 0
        out = expmv(1,A,b);
        return
    end
    if q == 1
        At = [A, b; zeros(q, length(A)) 0];
        if issparse(A)
            At = sparse(At);
        end
            
    else
        if issparse(A)
            At = sparse([A, b, sparse(zeros(length(A),q-1)); sparse(zeros(q,length(A))), sparse(diag(ones(q-1,1),1))]);
        else
            At = [A, b, (zeros(length(A),q-1)); (zeros(q,length(A))), (diag(ones(q-1,1),1))];
        end
    end
    e = zeros(length(At),1); e(length(A)+q) = 1;
    out = expmv(1,At,e);
    out = out(1:length(A));
end

