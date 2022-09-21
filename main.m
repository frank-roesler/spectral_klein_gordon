%% Setup:
clear
addpath('gradient_descent')

% Equation parameters:
m  = 1;

% Lattice in Fourier space:
r_n = 10;
n = 2*r_n^2;
xi = linspace(-r_n,r_n,n+1);

% W-integral:
% Choose R large enough to ensure |V|<1e-9 outside[-R,R]:
R = 0;
crit = 1;
while crit>1e-8
    R = R+1;
    crit = max(abs(potential(R:0.01:R+1)));
end
s = 2*round(n*max(1,R));  
x = linspace(-R,R,s+1);
V = potential(x);

% Precompute all integrals independent of z:
[I1,I2,J1,J2] = potential_integrals(R,n,s,xi,V);
% Precompute a0 mean values:
mean_a0 = 1./sqrt(xi.^2+m^2);
diff_mean_a0 = -xi./(xi.^2+m^2).^(3/2);

% Lattice in complex plane (defines the region where we look for eigenvalues):
h_res = 200;
bound = max(m, sqrt(abs(max(abs(V))^2-m^2)));
z1 = -2.5 - 0.2i;
z2 =  1 + 0.2i;
[Ln,dist_L,adj] = build_lattice(z1, z2, h_res);

%% Main loop:
norms = zeros(size(Ln));
parfor i=1:length(Ln(:))
    z = Ln(i);
    K = build_K(z,m,I1,I2,J1,J2,xi,mean_a0,diff_mean_a0);
    try
        norms(i) = svds(eye(n+1)-K,1,'smallest');
    catch
        disp(['svds failed for z = ',num2str(z)])
        norms(i) = 1./norm(inv(eye(n+1) - K));
    end
end

% Find local minima on lattice:
minima = islocalmin(abs(norms),1) & islocalmin(abs(norms),2);
spectrum_coarse = Ln(minima);
spectrum_coarse = spectrum_coarse(abs(real(spectrum_coarse))<m|...
    abs(imag(spectrum_coarse))>dist_L);  % only take points away from essential spectrum
spectrum_coarse = cast(spectrum_coarse,'like',1+1i);

%% GD minimization:
stepsize = 5e-4;
maxiter = 1000;
tol = 1e-8;
spectrum=zeros(size(spectrum_coarse));
parfor i=1:length(spectrum_coarse)
    z = spectrum_coarse(i);
    [z_end,ctr_gd] = GD(z,stepsize, maxiter, tol,n,m,I1,I2,J1,J2,xi,mean_a0,diff_mean_a0);
    spectrum(i) = z_end;
end
spectrum

%% Plot results:
figure('Position',[100,300,900,700])
subplot(3,1,1)
plot(x,V,'LineWidth',1)
title('Potential','fontsize',16,'fontweight','normal')
subplot(3,1,2)
contour(real(Ln),imag(Ln),log(1./abs(norms)),20,'LineWidth',0.6)
colormap winter
colorbar;
title('Logarithmic contour plot of ||(I-K(z))^{-1}||','fontsize',16,'fontweight','normal')
subplot(3,1,3)
plot(spectrum,'.','MarkerSize',15,'MarkerEdgeColor',[0.8,0,0])
title('Eigenvalues','fontsize',16,'fontweight','normal')
xlim([min(real(Ln(:))),max(real(Ln(:)))])
ylim([min(imag(Ln(:))),max(imag(Ln(:)))])

