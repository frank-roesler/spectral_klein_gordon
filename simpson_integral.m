function integral = simpson_integral(f,x)
% f is function, x is linspace vector describing domain of integration.
    x_coarse = x(1:2:end);
    f_coarse = f(1:2:end);
    a = x_coarse(1:end-1);
    b = x_coarse(2:end);
    fa = f_coarse(1:end-1);
    fb = f_coarse(2:end);
    fm = f(2:2:end-1);
    integral = sum((b-a)/6 .* (fa + 4*fm + fb));
end