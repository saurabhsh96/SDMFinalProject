%Final Project SDM
%Aim:
%1) Reflection coeff model and Tx Line Model 
%2) Quarter wave transformer with normal dielectric
%3) Multilayer-normal dielectric (Binomial/Tchebyshev)
%4) ADL-Slab (one slab-multilayer)
%5) Multiple ADL-slabs (many multilayer slabs)
%6) Shifted ADL (one slab-multilayer)
%7) Shifted ADL (many multilayer slabs)
%8) Analysis of Bandwidth in each case!

%1) This code deals with reflection coeff model and Tx Line Model
clear;
close all;

%% Defining inputs
er = 20;

%Some frequency range
freq = 5e9;
c = 3e8;
lam = c./freq;
k0 = 2*pi./lam;
ks = k0.*sqrt(er);

eps_0 = 8.854187817e-12;
mu_0 = 1.2566370614e-6;

zeta0 = 120*pi;
zetad = zeta0./sqrt(er);

%Refractive index
n0 = 1;
nd = sqrt(er);

drad = pi/180;
thC = asin((nd/n0))/drad;
if(isreal(thC))
    th = (0:.5:thC)*drad;
else
    th = (0:.5:90)*drad;
end 
ph = 0;

%% Reflection and Transmission using simple boundary condition problem
%E = E0 e^(-jkz) x
%H = E0/zeta0 e^(-jkz) y
%Transmit angle
thT = asin((n0./nd).*sin(th));

%For TE Incidence
rTE = (n0.*cos(th) - nd.*cos(thT))./(n0.*cos(th) + nd.*cos(thT));
RTE = abs(rTE).^2;
TTE = 1 - RTE;
%tCoeffTE = (2.*n0.*cos(th))./(n0.*cos(th) + nd.*cos(thT));

%For TM Incidence
rTM = (nd.*cos(th) - n0.*cos(thT))./(nd.*cos(th) + n0.*cos(thT));
RTM = abs(rTM).^2;
TTM = 1 - RTM;
%tCoeffTM = (2.*n0.*cos(th))./(n0.*cos(thT) + nd.*cos(th));

%Plotting
%TE
figure(1);
%plot(th./drad, abs(rTE), 'LineWidth', 1.5, 'DisplayName', '|\Gamma|'); hold on
%plot(th./drad, abs(tCoeffTE), 'LineWidth', 1.5, 'DisplayName', '|T|');
plot(th./drad, abs(RTE), 'LineWidth', 1.5, 'DisplayName', '|\Gamma|^2 BC'); hold on
plot(th./drad, abs(TTE), 'LineWidth', 1.5, 'DisplayName', '|T|^2 BC');
%plot(th./drad, abs(tCoeffTE) + abs(RTE), 'LineWidth', 1.5, 'DisplayName', '|Total|');
xlabel('\theta (degree)');
ylabel('|\Gamma|^2 and |T|^2');
title('Reflection Coeff of TE (Boundary Condition)');
grid on;
legend show;
hold off;

%TM
figure(2);
% plot(th./drad, abs(rTM), 'LineWidth', 1.5, 'DisplayName', '|\Gamma|'); hold on
% plot(th./drad, abs(tCoeffTM), 'LineWidth', 1.5, 'DisplayName', '|T|');
plot(th./drad, abs(RTM), 'LineWidth', 1.5, 'DisplayName', '|\Gamma|^2'); hold on
plot(th./drad, abs(TTM), 'LineWidth', 1.5, 'DisplayName', '|T|^2');
%plot(th./drad, abs(tCoeffTM) + abs(rTM), 'LineWidth', 1.5, 'DisplayName', '|Total|');
xlabel('\theta (degree)');
ylabel('|\Gamma|^2 and |T|^2');
title('Reflection Coeff of TM (Boundary Condition)');
grid on;
legend show;
hold off;

%% The Tx line model
% TE and TM

%Prop const
[~, ~, ~, kz0] = propConst(k0, th, ph);
[~, ~, ~, kzs] = propConst(ks, thT, ph);

%Air TE TM impedance
[Z0TE, Z0TM] = imped(zeta0, k0, kz0);

%Dielectric TE TM impedance
[ZdTE, ZdTM] = imped(zetad, ks, kzs);

%Zup and Zdown
ZupTE = Z0TE;
ZupTM = Z0TM;

ZdownTE = ZdTE;
ZdownTM = ZdTM;

[gammaTE_tx, tCoeffTE_tx] = refCoeff(ZupTE, ZdownTE);
[gammaTM_tx, tCoeffTM_tx] = refCoeff(ZupTM, ZdownTM);

%Plotting
%TE
figure(3);
plot(th./drad, abs(gammaTE_tx).^2, 'LineWidth', 1.5, 'DisplayName', '|\Gamma|^2 TL'); hold on
plot(th./drad, (abs(tCoeffTE_tx)), 'LineWidth', 1.5, 'DisplayName', '|T|^2 TL');

plot(th./drad, abs(RTE), 'LineWidth', 1.5, 'DisplayName', '|\Gamma|^2 BC');
plot(th./drad, abs(TTE), 'LineWidth', 1.5, 'DisplayName', '|T|^2 BC');

%plot(th./drad, abs(tCoeffTE_tx) + abs(gammaTE_tx), 'LineWidth', 1.5, 'DisplayName', '|Total|');
xlabel('\theta (degree)');
ylabel('|\Gamma|^2 and |T|^2');
title('Reflection and Transmission Coeff of TE');
grid on;
legend show;
hold off;

%TM
figure(4);
plot(th./drad, abs(gammaTM_tx).^2, 'LineWidth', 1.5, 'DisplayName', '|\Gamma|^2 TL'); hold on
plot(th./drad, (abs(tCoeffTM_tx)), 'LineWidth', 1.5, 'DisplayName', '|T|^2 TL');

plot(th./drad, abs(RTM), 'LineWidth', 1.5, 'DisplayName', '|\Gamma|^2 BC'); 
plot(th./drad, abs(TTM), 'LineWidth', 1.5, 'DisplayName', '|T|^2 BC');

%plot(th./drad, abs(tCoeffTM_tx) + abs(gammaTM_tx), 'LineWidth', 1.5, 'DisplayName', '|Total|');
xlabel('\theta (degree)');
ylabel('|\Gamma|^2 and |T|^2');
title('Reflection and Transmission Coeff of TM');
grid on;
legend show;
hold off;

%Note: It all makes sense, see, if you solve the equations of fields in a
%media in spectral domain, they result in VI relation with k and kz and at
%the same time the HE relation is with cos theta_i and with cos theta_t.
%Point to be noted here, is that both of these things are one and the same
%and the transmission line model does consider all the angles of
%incidences! Brilliant! So, first quesiton solved!! 

%% Tx Line Model freq variation

%Frequency variation
freq = 2e9:0.1e9:10e9;
c = 3e8;
lam = c./freq;
k0 = 2*pi./lam;
ks = k0.*sqrt(er);

th = 0;
ph = 0;

%Prop const
[~, ~, ~, kz0] = propConst(k0, th, ph);
[~, ~, ~, kzs] = propConst(ks, th, ph);

%Air TE TM impedance
[Z0TE, Z0TM] = imped(zeta0, k0, kz0);

%Dielectric TE TM impedance
[ZdTE, ZdTM] = imped(zetad, ks, kzs);

%Zup and Zdown
ZupTE = Z0TE;
ZupTM = Z0TM;

ZdownTE = ZdTE;
ZdownTM = ZdTM;

[gammaTE_tx, tCoeffTE_tx] = refCoeff(ZupTE, ZdownTE);
[gammaTM_tx, tCoeffTM_tx] = refCoeff(ZupTM, ZdownTM);

%Plotting
%TE
figure(5);
plot(freq./10^9, abs(gammaTE_tx).^2, 'LineWidth', 1.5, 'DisplayName', '|\Gamma|^2'); hold on
plot(freq./10^9, (abs(tCoeffTE_tx)), 'LineWidth', 1.5, 'DisplayName', '|T|^2');
%plot(th./drad, abs(tCoeffTE_tx) + abs(gammaTE_tx), 'LineWidth', 1.5, 'DisplayName', '|Total|');
xlabel('Frequency (GHz)');
ylabel('|\Gamma|^2 and |T|^2');
title('Reflection and Transmission Coeff Vs. Frequency (Tx-Line)');
grid on;
legend show;
hold off;

%TM
figure(6);
plot(freq./10^9, pow2db(abs(gammaTM_tx).^2), 'LineWidth', 1.5, 'DisplayName', 'TM'); hold on;
plot(freq./10^9, pow2db(abs(gammaTE_tx).^2), 'LineWidth', 1.5, 'DisplayName', 'TE');
%plot(freq./10^9, (abs(tCoeffTM_tx)), 'LineWidth', 1.5, 'DisplayName', '|T|^2');
%plot(th./drad, abs(tCoeffTM_tx) + abs(gammaTM_tx), 'LineWidth', 1.5, 'DisplayName', '|Total|');
xlabel('Frequency (GHz)');
ylabel('|\Gamma| (in dB)');
title('Reflection and Transmission Coeff in dB (Tx-Line)');
grid on;
legend show;
hold off;
ylim([-40, 0]);