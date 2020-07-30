%2) Quarter wave transformer (Normal dielectric)
clear;
close all;

%% Defining inputs
er = 20;
%er1 = sqrt(er);
%Some frequency range
freq = .01e9:0.1e9:10e9;
c = 3e8;
lam = c./freq;
k = 2*pi./lam;
ks = k.*sqrt(er);

%Fr0 To be used later, Center freq for BW calculations
fr0 = 5.01e9;
lam0 = c./fr0;
k0 = 2*pi./lam0;

drad = pi/180;
th = 0;
ph = 0;

eps_0 = 8.854187817e-12;
mu_0 = 1.2566370614e-6;

zeta0 = 120*pi;
zetad = zeta0./sqrt(er);

%% Matching with QWT
%Prop const
[~, ~, ~, kz] = propConst(k, th, ph);
[~, ~, ~, kzs] = propConst(ks, th, ph);

%Air TE TM impedance
[Z0TE, Z0TM] = imped(zeta0, k, kz);

%Dielectric TE TM impedance
[ZdTE, ZdTM] = imped(zetad, ks, kzs);

%Zup and Zdown
ZupTE = Z0TE;
ZupTM = Z0TM;

%QWT
%Impedance of QWT Tx Line
Z1TE = sqrt(ZupTE.*ZdTE);
Z1TM = sqrt(ZupTM.*ZdTM);

%zeta1
%Er value of the medium that is going to be coated on top of the
%semi-inifinite medium is also important. As that will decide the length
%and impedance of the line ultimately. Below calculation is done for the
%broadside case. Similar calculation was done for Lens Antennas as well.
zeta1 = Z1TE;
er1 = sqrt(er);
ks1 = k.*sqrt(er1);
h = lam0./(4.*sqrt(er1));
[~, ~, ~, kzs1] = propConst(ks1, th, ph);

%Impedance at point A, Z0, Z1 boundary
ZATE = findZ(Z1TE, ZdTE, kzs1, h);
ZATM = findZ(Z1TM, ZdTM, kzs1, h);

%RefCoeff Calculation
[gammaTE_tx, tCoeffTE_tx] = refCoeff(ZupTE, ZATE);
[gammaTM_tx, tCoeffTM_tx] = refCoeff(ZupTM, ZATM);

%Plotting
%TE
% figure(1);
% plot(freq./10^9, abs(gammaTE_tx).^2, 'LineWidth', 1.5, 'DisplayName', '|\Gamma|^2'); hold on
% plot(freq./10^9, (abs(tCoeffTE_tx)), 'LineWidth', 1.5, 'DisplayName', '|T|^2');
% %plot(th./drad, abs(tCoeffTE_tx) + abs(gammaTE_tx), 'LineWidth', 1.5, 'DisplayName', '|Total|');
% xlabel('Frequency (GHz)');
% ylabel('|\Gamma|^2 and |T|^2');
% title('Reflection Coeff of TE (with QWT)');
% grid on;
% legend show;
% hold off;

%TM
figure(2);
plot(freq./10^9, abs(gammaTM_tx).^2, 'LineWidth', 1.5, 'DisplayName', '|\Gamma|^2'); hold on
plot(freq./10^9, (abs(tCoeffTM_tx)), 'LineWidth', 1.5, 'DisplayName', '|T|^2');
%plot(th./drad, abs(tCoeffTM_tx) + abs(gammaTM_tx), 'LineWidth', 1.5, 'DisplayName', '|Total|');
xlabel('Frequency (GHz)');
ylabel('|\Gamma|^2 and |T|^2');
title('Reflection and Transmission Coeff (with QWT)');
grid on;
legend show;
hold off;

%TM and TE Ref Coeff DB
figure(3);
%plot(freq./10^9, pow2db(abs(gammaTM_tx).^2), 'LineWidth', 1.5, 'DisplayName', 'TM'); hold on
plot(freq./10^9, pow2db(abs(gammaTE_tx).^2), 'LineWidth', 1.5, 'DisplayName', 'TE');
xlabel('Frequency (GHz)');
ylabel('Reflection Coefficient (in dB)');
title('Reflection Coeff of QWT (in dB)');
grid on;
%legend show;
hold off;
ylim([-40, 0]);

%Just reflection coefficient
figure(4);
plot(freq./10^9, abs(gammaTE_tx), 'LineWidth', 1.5, 'DisplayName', '|\Gamma TE|'); hold on
plot(freq./10^9, abs(gammaTM_tx), 'LineWidth', 1.5, 'DisplayName', '|\Gamma TM|');
%plot(th./drad, abs(tCoeffTE_tx) + abs(gammaTE_tx), 'LineWidth', 1.5, 'DisplayName', '|Total|');
xlabel('Frequency (GHz)');
ylabel('|\Gamma|');
title('Reflection Coeff of structure with QWT');
grid on;
legend show;
hold off;

%% BW Calculation

BWTM = BWCalc(freq, pow2db(abs(gammaTM_tx).^2));
BWTE = BWCalc(freq, pow2db(abs(gammaTE_tx).^2));

%% Saving important data
save('singleQWT.mat', 'er1', 'h', 'Z1TE', 'Z1TM', 'freq', 'fr0', 'er1');