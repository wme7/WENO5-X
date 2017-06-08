%% Mapping function for WENO5-M
clear, clc, close all;

% Coded by Manuel A. Diaz, NHRI, 2017.05.20

%% Reproduction of Figure 5 in:
% Henrick, Andrew K., Tariq D. Aslam, and Joseph M. Powers. "Mapped
% weighted essentially non-oscillatory schemes: achieving optimal order
% near critical points." JCP 207.2 (2005): 542-567.

% Initialize discrete space
w = 0:0.01:1;

% The linear weigths
d0 = 1/10;
d1 = 6/10;
d2 = 3/10;

% The mapped functions for WENO5-M
g0 = w.*(d0+d0^2-3*d0*w+w.^2) ./ (d0^2 + w*(1-2*d0));
g1 = w.*(d1+d1^2-3*d1*w+w.^2) ./ (d1^2 + w*(1-2*d1));
g2 = w.*(d2+d2^2-3*d2*w+w.^2) ./ (d2^2 + w*(1-2*d2));

% Plot the mapped functions
figure(1); plot(w,w,'--',w,g0,w,g1,w,g2); axis([0,1,0,1]);
xlabel('$\omega$','interpreter','latex'); 
ylabel('$\alpha_{k}^{*}$','interpreter','latex');
legend('Identity mapping','k=0','k=1','k=2','Location','Northwest'); 
legend boxoff;