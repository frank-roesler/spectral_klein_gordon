function [integral1, integral2, corrector1, corrector2] = potential_integrals(R,n,N,xi,V)
%     Computes the integrals of e_k*e_l*V^2 and e_k*e_l*(2V) and first order 
%     correction terms for potential V
    x = linspace(-R,R,N+1);
    volQ = (max(xi)-min(xi))/(length(xi)-1);
    r_n = max(xi);
    e = sqrt(volQ/(2*pi)) * exp(1i*xi'*x) .* sinc(x/pi*r_n/n);                % Integral of e^(i*xi*x)
    cor = -1i * sqrt(volQ/(2*pi)) * volQ/2 * exp(1i*xi'*x) .* dsinc(x*r_n/n); % Integral of (xi-xi_m)*e^(i*xi*x)

    W1 = V.^2;
    W2 = 2*V;
    integral1 = zeros(n+1);
    integral2 = zeros(n+1);
    corrector1 = zeros(n+1);
    corrector2 = zeros(n+1);
    for k=1:n+1
        for l=1:n+1
            integral1(k,l)  = simpson_integral( conj(e(k,:)).*e(l,:).*W1, x);
            integral2(k,l)  = simpson_integral( conj(e(k,:)).*e(l,:).*W2, x);
            corrector1(k,l) = simpson_integral( conj(cor(k,:)).*cor(l,:).*W1, x);
            corrector2(k,l) = simpson_integral( conj(cor(k,:)).*cor(l,:).*W2, x);
        end
    end
end






