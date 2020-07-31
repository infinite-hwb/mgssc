function [Af,C] = centroid_MGSSC(n, K, mu,lambda)

num_iter = 5;
max_mu = 1e6;
rho = 1.5;
err_thr = 10^-4;
num_views = 3;

C1 = repmat({zeros(n,n)}, 1, num_views);
C2 = repmat({zeros(n,n)}, 1, num_views);
A = repmat({zeros(n,n)}, 1, num_views);
A_prev = repmat({zeros(n,n)}, 1, num_views);
Lambda1 = repmat({zeros(n,n)}, 1, num_views);
Lambda2 = repmat({zeros(n,n)}, 1, num_views);
C_centroid = zeros(n,n);

mu1 = mu;
mu2 = mu;
mu = [mu1 mu2];

iter = 0;
converged = false;

%for iter = 1:20
while iter < num_iter && ~converged

    
    for v = 1:num_views
        A_prev{v} = A{v}; % save previous value
        [C1{v}, C2{v}, Lambda1{v}, Lambda2{v}, A{v}] = centroid_MGSSC_1view...
            (n, K{v}, C1{v}, C2{v}, C_centroid, Lambda1{v}, Lambda2{v}, ...
            lambda, mu);
    end
    
    % update centroid
    for v = 1:num_views
        C_centroid = C_centroid+lambda(2)*C1{v};
    end
    C_centroid = C_centroid/(num_views*lambda(2));
    
    % check convergence
    converged = true;
    for v=1 : num_views
        err1 = max(max(abs(A_prev{v}-A{v})));
        err2 = max(max(abs(A{v}-C1{v}+diag(diag(C1{v})))));
        err3 = max(max(abs(A{v}-C2{v})));
        xx = err1 + err2 + err3; 
        if err1>err_thr || err2>err_thr || err3>err_thr
            converged = false;
            break
        end    
    end
    iter = iter + 1;
    mu = min(rho*mu,max_mu);
end

C = C_centroid;
Af = abs(C)+abs(C');
end
