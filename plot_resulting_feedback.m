%xlim([7000 10000])
%ylim([-6e-3 6e-3])

% xlim([9000 10000])
% ylim([0.999993 1])
meana = a0;
% waveform:
a1=al(1); a2=al(2); a3=al(3); a4=al(4); a5=al(5);
a_1=conj(a1); a_2=conj(a2); a_3=conj(a3); a_4=conj(a4); a_5=conj(5);

% phase response curve:
z_1=conj(zl(1)); z_2=conj(zl(2)); z_3=conj(zl(3)); z_4=conj(zl(4));

omega = 0.242622176836821;

x_wave = @(theta) meana;
prc = @(theta) meanz;

for n = 1:5
      x_wave = @(theta) x_wave(theta) + 2*real(al(n))*cos(n*omega*theta) + 2*imag(al(n))*sin(n*omega*theta);
      prc = @(theta) prc(theta) + 2*real(zl(n))*cos(n*omega*theta) + 2*imag(zl(n))*sin(n*omega*theta);
end

%%

figure
theta = linspace(0,2*pi/omega,100);
plot(theta,x_wave(theta),'o')
hold on;
plot(theta, prc(theta), 'o')

%%
tau100 = 14.3647942499046;
k100 = 0.000648419285917711;
tau300 = 1.67599439739307;
k300 = 0.00189831681865462;

h_sync_unstable = @(psi1,psi2) k100*x_wave(psi1 - omega*tau100)...
  + k300*x_wave(psi1 - omega*tau300).^3;

tau100 = 12.1774667204947;
k100 = 0.00120246591503002;
tau300 = 2.16777143200637;
k300 = 0.00332274026924014;
tau120 = 19.5579367631832;
k120 = 0.00428601684294361;
tau210 = 6.24120179118049;
k210 = 0.0031692533558537;
tau110 = 22.9888818752266;
k110 = 0.00130253560476841;

h_sync_stable = @(psi1,psi2) k100*x_wave(psi1 - omega*tau100)...
  + k300*x_wave(psi1 - omega*tau300).^3 ...
  + k120*x_wave(psi1 - omega*tau120) .* x_wave(psi2 - omega*tau120).^2 ...
  + k210*x_wave(psi1 - omega*tau210).^2 .* x_wave(psi2 - omega*tau210) ...
  + k110*x_wave(psi1 - omega*tau210) .* x_wave(psi2 - omega*tau210);

%%
[psi1,psi2]  = meshgrid( linspace(0,2*pi/omega,100));

figure
surf(psi1,psi2,h_sync_stable(psi1,psi2))
hold on;
surf(psi1,psi2,h_sync_unstable(psi1,psi2))

xlabel('\psi_1')
ylabel('\psi_2')
zlabel('h')

%%
xi1s = 0.5;  chi1s = -pi/2;
xi2s = 0.5;  chi2s = -pi/4;
xi3s = 0.5;  chi3s = -pi/2;

H100 = 0.5 * xi1s * exp(-1i * chi1s);
H110 = 0.5 * xi3s * exp(-1i * chi3s);
H200 = 0.5 * xi2s * exp(-2 * 1i * chi2s);

H_sync_stable = @(psi1,psi2) H100 * exp( -1i*(1*psi1) ) ...
          + H110 * exp( -1i*(1*psi1 + 1*psi2)) ...
          + H200 * exp( -1i*(2*psi1)  );

xi3us = 0;
H100 = 0.5 * xi1s * exp(-1i * chi1s);
H110 = 0.5 * xi3us * exp(-1i * chi3s);
H200 = 0.5 * xi2s * exp(-2 * 1i * chi2s);

H_sync_unstable = @(psi1,psi2) H100 * exp( -1i*(1*psi1) ) ...
          + H110 * exp( -1i*(1*psi1 + 1*psi2)) ...
          + H200 * exp( -1i*(2*psi1)  );
%%
figure
[psi1,psi2]  = meshgrid( linspace(0,2*pi,100));
surf(psi1,psi2,real(H_sync_stable(psi1,psi2)))
hold on;
surf(psi1,psi2,real(H_sync_unstable(psi1,psi2)))

xlabel('\psi_1')
ylabel('\psi_2')
zlabel('H')


