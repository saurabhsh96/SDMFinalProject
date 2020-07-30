%3) Multilayer-normal dielectric (Binomial/Tchebyshev)
clear;
close all;

%% Defining Inputs
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

%% Multiple QWT Binomial distribution
%Number of stratifications 2, 3, and 4 
N = 4;

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
% [ZA, Z, erZ, h] = binomialImped(ZL, Z0, N, er, th, ph, freq, fr0)
[ZATE, ZTE, erZTE, hTE, ZpTE] = binomialImped(ZdTE, ZupTE, N, er, th, ph, freq, fr0);
[ZATM, ZTM, erZTM, hTM, ZpTM] = binomialImped(ZdTM, ZupTM, N, er, th, ph, freq, fr0);

%RefCoeff Calculation
[gammaTE_tx, tCoeffTE_tx] = refCoeff(ZupTE, ZATE);
[gammaTM_tx, tCoeffTM_tx] = refCoeff(ZupTM, ZATM);

%Plotting
%TE
figure(1);
plot(freq./10^9, abs(gammaTE_tx).^2, 'LineWidth', 1.5, 'DisplayName', '|\Gamma|^2'); hold on
plot(freq./10^9, (abs(tCoeffTE_tx)), 'LineWidth', 1.5, 'DisplayName', '|T|^2');
%plot(th./drad, abs(tCoeffTE_tx) + abs(gammaTE_tx), 'LineWidth', 1.5, 'DisplayName', '|Total|');
xlabel('Frequency (GHz)');
ylabel('|\Gamma|^2 and |T|^2');
title(['Reflection and Transmission Coeff ', num2str(N), '-section', ' QWT']);
grid on;
legend show;
hold off;

% %TM
% figure(2);
% plot(freq./10^9, abs(gammaTM_tx).^2, 'LineWidth', 1.5, 'DisplayName', '|\Gamma|^2'); hold on
% plot(freq./10^9, (abs(tCoeffTM_tx)), 'LineWidth', 1.5, 'DisplayName', '|T|^2');
% %plot(th./drad, abs(tCoeffTM_tx) + abs(gammaTM_tx), 'LineWidth', 1.5, 'DisplayName', '|Total|');
% xlabel('Frequency (GHz)');
% ylabel('|\Gamma|^2 and |T|^2');
% title('Reflection Coeff of TM (with QWT)');
% grid on;
% legend show;
% hold off;

%TM and TE Ref Coeff DB
figure(3);
plot(freq./10^9, pow2db(abs(gammaTM_tx).^2), 'LineWidth', 1.5, 'DisplayName', 'TM'); hold on
plot(freq./10^9, pow2db(abs(gammaTE_tx).^2), 'LineWidth', 1.5, 'DisplayName', 'TE');
xlabel('Frequency (GHz)');
ylabel('Reflection Coefficient (in dB)');
title(['Reflection and Transmission Coeff ', num2str(N), '-section', ' QWT']);grid on;
legend show;
hold off;
ylim([-40, 0]);

%Just reflection coefficient
% figure(4);
% plot(freq./10^9, abs(gammaTE_tx), 'LineWidth', 1.5, 'DisplayName', '|\Gamma TE|'); hold on
% plot(freq./10^9, abs(gammaTM_tx), 'LineWidth', 1.5, 'DisplayName', '|\Gamma TM|');
% %plot(th./drad, abs(tCoeffTE_tx) + abs(gammaTE_tx), 'LineWidth', 1.5, 'DisplayName', '|Total|');
% xlabel('Frequency (GHz)');
% ylabel('|\Gamma|');
% title('Reflection Coeff of TE (with QWT)');
% grid on;
% legend show;
% hold off;

%BW Calculation

BWTM = BWCalc(freq, pow2db(abs(gammaTM_tx).^2));
BWTE = BWCalc(freq, pow2db(abs(gammaTE_tx).^2));

%% Saving important data
save('multipleQWT.mat', 'erZTM', 'erZTE', 'hTM', 'hTE', 'ZTE', 'ZTM', 'ZpTE',...
    'ZpTM', 'N', 'freq', 'fr0');

%% Multiple QWT Chebyshev distribution
%Maybe later if time permits :P! LoL!! :P
%Advantage of Chebyshev is that it can have better BW than that of binomial
%for the same number of sections however, it will have ripples instead of
%flat passband, and that is a disadvantage, if a smooth curve is expected!!