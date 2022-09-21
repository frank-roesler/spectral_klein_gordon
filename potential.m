function V = potential(x)
    % Sauter-like potential:
    v0 = 3.7;
    D = 3.2;
    W = 0.3;
    V = -v0*(tanh((x+D/2)/W) - tanh((x-D/2)/W))/2;
end