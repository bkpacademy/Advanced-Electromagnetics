% Matlab Programme for E-H filed plot for circular waveguide TMmn mode
% OCT 30, 2018
% AUTHOR: AJEET KUMAR

clc;
close all;

% % Zeros Xmn_d of derivatives of Bessel fumction Jm_d(Xmn_d)
% % m starts from 0, n starts form 1
% Xmn_d = [3.8318 1.8412 3.0542 4.2012 5.3175 6.4155 7.5013 8.5777 9.6474 10.7114 11.7708 12.8264
% 7.0156 5.3315 6.7062 8.0153 9.2824 10.5199 11.7349 12.9324 14.1155 15.2867 16.4479 17.6003
% 10.1735 8.5363 9.9695 11.3459 12.6819 13.9872 15.2682 16.5294 17.7740 19.0046 20.2230 21.4309
% 13.3237 11.7060 13.1704 14.5859 15.9641 17.3129 18.6375 19.9419 21.2291 22.5014 23.7607 25.0085
% 16.4706 14.8636 16.3475 17.7888 19.1960 20.5755 21.9317 23.2681 24.5872 25.8913 27.1820 28.4609];


% Zeros Xmn of Bessel fumction Jm(Xmn)
% m starts from 0, n starts form 1
Xmn = [2.4049 3.8318 5.1357 6.3802 7.5884 8.7715 9.9361 11.0864 12.2251 13.3543 14.4755 15.5898;
5.5201 7.0156 8.4173 9.7610 11.0647 12.3386 13.5893 14.8213 16.0378 17.2412 18.4335 19.6160;
8.6537 10.1735 11.6199 13.0152 14.3726 15.7002 17.0038 18.2876 19.5545 20.8071 22.0470 23.2759;
11.7915 13.3237 14.7960 16.2235 17.6160 18.9801 20.3208 21.6415 22.9452 24.2339 25.5095 26.7733;
14.9309 16.4706 17.9598 19.4094 20.8269 22.2178 23.5861 24.9349 26.2668 27.5838 28.8874 30.1791];



% Waveguide dimensions
a = 1;  % radius in cm in x-direction
f = 15e+9;   % Frequency of operation 45GHz
c = 3e+8;    % Velocity of light

% m = 0;    % Mode number in rho-Direction
% n = 1;    % Mode number in phi-Direction
m = input('Enter TM mode value m:');
n = input('Enter TM mode value n:');

% Xmn_dTE = Xmn_d(n,m+1);             % TEmn zero for derivative of Bessel function
Xmn_TM = Xmn(n,m+1);             % TMmn zeros for Bessel function

% Wave propagation in Z-Direction
%********************************%


epsilon = 8.8540e-12;           % Permittivity constant
epsilon_r = 1;                  % Relative Permittivity constant

mu1 = 4*pi*10e-7;               % Permeability constant
mu1_r = 1;                      % Relative Permeability constant
fc = c*100*Xmn_TM/(2*pi*a*1e9*sqrt(epsilon_r*mu1_r));      % Cutoff frequency calculation in GHz
lambda = c*100/(fc*1e9);        % Wavelength in cm
omega = 2*pi*f;                 % Frequency of operation in rad/s
M = 100;                         % Number of points to be poltted

beta = omega*(sqrt(mu1*mu1_r*epsilon*epsilon_r));   % Propagation constant
beta_rho = Xmn_TM/a;                   % beta in rho direction
beta_z = sqrt(beta.^2-beta_rho.^2);         % beta in z - direction

% Front View field plot
z = 0;
rho = linspace(0,a,M);
phi = linspace(0,2*pi,M);
% z = linspace(0,2*lambda,M);
[rho,phi] = meshgrid(rho,phi);
k = (beta_rho.*rho);
% Jm = besselj(m,k);
% Jm_d = diff(Jm);
% Jm1 = repmat(Jm,M,1);
% [p, q] = size(Jm1);
% Jm_d = repmat(Jm_d,M,1);
% Jm_d1 = padarray(Jm_d,[0,1],Jm_d(p-1,q-1),'post');      % Pad one column with last values
% [rho,phi] = meshgrid(rho,phi);
Jm1 = besselj(m,k);            % Bessel Function of 1st kind of order m
% Jm_d1 = -besselj(m+1,k)+(m./(k)).*besselj(m,k);       % Diff. of Bessel function
Jm_d1 = (besselj(m-1,k)-besselj(m+1,k))/2;       % Diff. of Bessel function

% Constants Assumptions in the filed equation
Bmn = 1;
C2 = 1;
D2 = 1;

%%
%Field Equations for TMmn in circular waveguide
E_rho = -Bmn.*(beta_rho*beta_z/(omega*mu1*epsilon)).*Jm_d1.*(C2.*cos(m.*phi)+D2.*sin(m.*phi)).*exp(-1i*beta_z*z);
E_phi = Bmn.*(m*beta_z./(omega*mu1*epsilon.*rho)).*Jm1.*(-C2.*sin(m.*phi)+D2.*cos(m.*phi)).*exp(-1i*beta_z*z);
Ez = -1j*Bmn.*(beta_z.^2./(omega*mu1*epsilon)).*Jm1.*(C2.*cos(m.*phi)+D2.*sin(m.*phi)).*exp(-1i*beta_z*z);

H_rho = -Bmn.*(m./(mu1.*rho)).*Jm1.*(-C2.*sin(m.*phi)+D2.*cos(m.*phi)).*exp(-1i*beta_z*z);
H_phi = -Bmn.*(beta_rho/(mu1)).*Jm_d1.*(C2.*cos(m.*phi)+D2.*sin(m.*phi)).*exp(-1i*beta_z*z);
Hz = 0;



%%
% Cylindrical to cartesian conversion
[x,y]=pol2cart(phi,rho);
% Carteisian Field Matrix Initialisation
Ex = zeros(M,M);
Ey = zeros(M,M);
E_cart_net = zeros(M,M);
Hx = zeros(M,M);
Hy = zeros(M,M);
H_cart_net = zeros(M,M);

%E-H-Field conversion from cylindrical to cartesian
% [Ex, Ey] = pol2cart(E_phi, E_rho);
for ii = 1:M
    for jj = 1:M
Ex(ii,jj) = E_rho(ii,jj).*cos(phi(ii,jj)) - E_phi(ii,jj).*sin(phi(ii,jj));
Ey(ii,jj) = E_rho(ii,jj).*sin(phi(ii,jj)) + E_phi(ii,jj).*cos(phi(ii,jj));
E_cart_net(ii,jj) = sqrt(Ex(ii,jj).^2+Ey(ii,jj).^2);

Hx(ii,jj) = H_rho(ii,jj).*cos(phi(ii,jj)) - H_phi(ii,jj).*sin(phi(ii,jj));
Hy(ii,jj) = H_rho(ii,jj).*sin(phi(ii,jj)) + H_phi(ii,jj).*cos(phi(ii,jj));
H_cart_net(ii,jj) = sqrt(Hx(ii,jj).^2+Hy(ii,jj).^2);
    end
end

vect_length = 1.5;
figure();
quiver(x,y,Ex,Ey,vect_length);
title(['Plot of front view for TM_',num2str(m),'_',num2str(n),' E-Field']);
legend('E-Field');
xlabel('x-axis');
ylabel('y-asix');
% axis([-0.1 0.1 -0.1 0.1]);

figure();
% set(H,'PropertyName',PropertyValue)  % For objects H

quiver(x,y,Hx,Hy,vect_length);
% quiver(x,y,real(Hx),real(Hy),1.5,'b','LinEwidth',1);
% contour(x,y,real(H_cart_net));
title(['Plot of front view for TM_',num2str(m),'_',num2str(n),' H-Field']);

legend('H-Field');
xlabel('x-axis');
ylabel('y-asix');
% axis([-0.1 0.1 -0.1 0.1]);

figure()

quiver(x,y,Ex,Ey);
% set(gca, 'color', 'blue');          % Set background color
hold on
quiver(x,y,Hx,Hy);
title(['E-H Filed of TM_',num2str(m),'_',num2str(n)]);
legend('E-Field','H-Field');

% %
% % E-H-Field Plot together in one window
% figure()
% subplot(2,3,1);
% quiver(x,y,real(Ex),real(Ey));
% % contour(x,y,real(E_cart_net));
% title(['TM_',num2str(m),'_',num2str(n),' E-Field']);
% legend('E-Field');
% subplot(2,3,2);
% quiver(x,y,real(Hx),real(Hy));
% % contour(x,y,real(H_cart_net));
% title(['TM_',num2str(m),'_',num2str(n),' H-Field']);
% legend('H-Field');
% subplot(2,3,3);
% quiver(x,y,real(Ex),real(Ey));
% % contour(x,y,real(E_cart_net));
% hold on
% quiver(x,y,real(Hx),real(Hy));
% % contour(x,y,real(H_cart_net));
% title(['TM_',num2str(m),'_',num2str(n)]);
% legend('E-Field','H-Field');
%%