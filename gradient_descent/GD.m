function [z_c,ctr] = GD(z_start, stepsize, maxiter, tol, n,m,I1,I2,J1,J2,xi,mean_a0,diff_mean_a0)
    D = inf;
    discretization = 1e-8;
    z_c = z_start;
    ctr=0;
    while D>tol && ctr<maxiter
        [normK, grad] = numeric_gradient(z_c, discretization, n,m,I1,I2,J1,J2,xi,mean_a0,diff_mean_a0);
        
        z_c = z_c - grad(1)*stepsize - 1i*grad(2)*stepsize;
        if normK>D
            stepsize = 0.5*stepsize;
        end
        ctr = ctr+1;
        D = normK;
    end
    if ctr==maxiter
        warning('Maximum number of GD steps reached. Results may be inaccurate.')
    end
end