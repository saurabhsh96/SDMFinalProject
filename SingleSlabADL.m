%4) ADL-Slab (one slab-multilayer)
clear;
close all;

%% Defining Inputs
er = 20;
dielectricData = load('singleQWT.mat');

%Some frequency range
freq = dielectricData.freq;
c = 3e8;
lam = c./freq;
k = 2*pi./lam;
ks = k.*sqrt(er);

%Fr0 To be used later, Center freq for BW calculations
fr0 = dielectricData.fr0;
lam0 = c./fr0;
k0 = 2*pi./lam0;

drad = pi/180;
th = 0;
ph = 0;

eps_0 = 8.854187817e-12;
mu_0 = 1.2566370614e-6;

zeta0 = 120*pi;
zetad = zeta0./sqrt(er);

%% ADL Slab Inputs
%Number of layers
N = 6;

%Distance between the patches
%W = 0.01.*lam0;
dy = 0.2.*lam(end);

%m - slightly different because of the function
m = [-30:-1, 1:30];

%Shift
%shift = 0;
shift = 0*dy;

%% ADL Slab Design
%Determining required Zin for the structure
%ZinTE = dielectricData.Z1TE(find(ismember(freq,fr0)));
%ZinTM = dielectricData.Z1TM(find(ismember(freq,fr0)));

%Calculation must be done only for the center frequency as for other set of
%frequencies of the equations change drastically. 
erHost = 1; %Air 
h = dielectricData.h; %Required height of the slab
erR = dielectricData.er1;
[W, dz, S, dzud, Zprev] = ADLSlabDesign3(zeta0, zetad, dy, h, N, erHost, m, fr0, th, ph, 0.022, "TE", shift);
%[epsilon, mu, eta, n, zeta] = StoEpsilon(S, th, h, k0, "TE");

for indF = 1:length(freq)
    %Calculating required impedances for the particular frequency
    omega = 2.*pi.*freq(indF);
    lam = c./(freq(indF).*sqrt(erHost));
    k = 2.*pi./lam;
        
    %Propagation constant
    kx = k.*sin(th).*cos(ph);
    ky = k.*sin(th).*sin(ph);
    kRho = sqrt(kx.^2 + ky.^2);
    kz0 = -1j*sqrt(-((k.^2)-(kRho.^2)));

    %Tx Line Impedance
%   ZTx = (zeta0.*k)./kz0;
    ZTx = (zeta0.*kz0)./k;
     
    %Susceptance
    if(shift > 0)
        B_SI = suscpetance_SIShifted(omega, dy, dz, m, W, shift);
        B_I = suscpetance_IShifted(omega, dy, dz, m, W, shift);
    else
        B_SI = suscpetance_SI(omega, dy, dz, m, W);
        B_I = suscpetance_I(omega, dy, dz, m, W);
    end
    
    Zedge = -1j./B_SI.*(1./(1-((sin(th).^2)./2)));
    %Middle impedance
    Zinf = -1j./B_I.*(1./(1-((sin(th).^2)./2)));

    Zprev1 = effectiveImp(N, zetad, Zedge, Zinf, ZTx, dzud, dz, k);
    %Zreq(indF) = Zprev1;
    [S11(indF), S12(indF)] = refCoeff(zeta0, Zprev1);    
end
%Plotting
figure(1);
plot(freq./10^9, abs(S11).^2, 'LineWidth', 1.5, 'DisplayName', '|\Gamma|^2'); hold on;
%plot(freq./10^9, abs(S11te).^2, 'LineWidth', 1.5, 'DisplayName', 'S11(TE)');
%plot(freq./10^9, abs(S12tm).^2, 'LineWidth', 1.5, 'DisplayName', 'S12(TE)');
plot(freq./10^9, abs(S12), 'LineWidth', 1.5, 'DisplayName', '|T|^2');
title(['|S11|^2 and |S12|^2; ADL with N =', num2str(N)]);
xlabel('Frequency (in GHz)');
ylabel('|\Gamma|^2 and |T|^2');
legend show;
grid on;
hold off;

figure(2);
plot(freq./10^9, pow2db(abs(S11).^2), 'LineWidth', 1.5, 'DisplayName', 'TM'); hold on
%plot(freq./10^9, pow2db(abs(gammaTE_tx).^2), 'LineWidth', 1.5, 'DisplayName', 'TE');
xlabel('Frequency (GHz)');
ylabel('Reflection Coefficient (in dB)');
title('Reflection Coeff of ADL slab');
grid on;
% legend show;
hold off;
%ylim([-40, 0]);

figure(3);
plot(freq./10^9, abs(S11), 'LineWidth', 1.5, 'DisplayName', '|\Gamma TE|'); hold on
%plot(freq./10^9, abs(gammaTM_tx), 'LineWidth', 1.5, 'DisplayName', '|\Gamma TM|');
%plot(th./drad, abs(tCoeffTE_tx) + abs(gammaTE_tx), 'LineWidth', 1.5, 'DisplayName', '|Total|');
xlabel('Frequency (GHz)');
ylabel('|\Gamma|');
title('Reflection Coeff of ADL slab');
grid on;
legend show;
hold off;

% figure(4);
% plot(freq./10^9, abs(Zreq), 'LineWidth', 1.5, 'DisplayName', '|\Gamma TE|'); hold on
% %plot(freq./10^9, abs(gammaTM_tx), 'LineWidth', 1.5, 'DisplayName', '|\Gamma TM|');
% %plot(th./drad, abs(tCoeffTE_tx) + abs(gammaTE_tx), 'LineWidth', 1.5, 'DisplayName', '|Total|');
% xlabel('Frequency (GHz)');
% ylabel('|\Gamma|');
% title('Reflection Coeff of TE (with QWT)');
% grid on;
% legend show;
% hold off;

%% BW Calculation
BW = BWCalc(freq, pow2db(abs(S11).^2));