function [normK, grad] = numeric_gradient(z_0, discretization, n,m,I1,I2,J1,J2,xi,mean_a0,diff_mean_a0)
% Computes the gradient of ||(I-K)|| numerically
    stencil = [z_0, z_0+discretization, z_0+1i*discretization];
    norms = zeros(1,length(stencil));

    for j=1:length(stencil)
        z = stencil(j);
        K = build_K(z,m,I1,I2,J1,J2,xi,mean_a0,diff_mean_a0);
        norms(j) = svds(eye(n+1)-K,1,'smallest');
    end
    grad = [(norms(2)-norms(1)), (norms(3)-norms(1))];
    grad = grad./discretization;
    normK = norms(1);
end