function K = build_K(z,m,I1,I2,J1,J2,xi,mean_a0,diff_mean_a0)
    % Taylor coefficiants for a_0 and a_lambda:
    mean_az = sqrt(xi.^2+m^2)./(xi.^2+m^2-z^2);
    diff_mean_az = -xi./sqrt(xi.^2+m^2) .* (xi.^2+m^2+z^2)/(xi.^2+m^2-z^2).^2;
    
    K = mean_a0.*mean_az.'.*(I1-z*I2) + diff_mean_a0.*diff_mean_az.'.*(J1-z*J2);
end


