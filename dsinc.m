function y = dsinc(x)
%   derivative of sin(x)/x;
    tol = 1e-6;
    y = zeros(size(x));
    where_x_big = abs(x)>tol;
    where_x_small = abs(x)<=tol;
    
    y(where_x_big) = (cos(x(where_x_big)) - sinc(x(where_x_big)/pi))./x(where_x_big);
    y(where_x_small) = -x(where_x_small)/3 + x(where_x_small).^3/30 - x(where_x_small).^5/840;
end