%4) ADL-Slab (multiple slab-multilayer)
clear;
close all;

%% Defining Inputs
er = 20;
dielectricData = load('multipleQWT.mat');

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
%Number of slabs
M = dielectricData.N;

%Number of layers
%N = [6 5 3];
% N = [4 4];
N = [4 3 3 2];

%Distance between the patches
%W = 0.01.*lam0;
%dy = [0.2.*lam(end) 0.2.*lam(end) 0.4.*lam(end)];
%dy = [0.2.*lam(end) 0.2.*lam(end)];
dy = [0.2.*lam(end) 0.2.*lam(end) 0.4.*lam(end) 0.4.*lam(end)];

%m - slightly different because of the function
m = [-30:-1, 1:30];

%Shift
%shift = [0 0 0].*dy;
shift = [0.5 0.5 0.5 0.5].*dy;

%Tolerances
%tol = [0.01 0.01 0.03];
tol = [0.02 0.02 0.03 0.03];

%% ADL Slab Design
%Determining required Zin for the structure
%ZinTE = dielectricData.Z1TE(find(ismember(freq,fr0)));
%ZinTM = dielectricData.Z1TM(find(ismember(freq,fr0)));

%Calculation must be done only for the center frequency as for other set of
%frequencies of the equations change drastically. 
%Now as TE = TM in this case, we will straightforward implement fot TE mode
%and later will make changes in the code to accomodate the TM. Now, there
%is no time!
erHost = 1; %Air 
h = dielectricData.hTE; %Required height of the slab
erR = dielectricData.erZTE;    
for ind = 1:M
    Zin = dielectricData.ZpTE(ind);
    if(ind == M)
        ZL = zetad;
    else
        ZL = dielectricData.ZpTE(ind+1);
%       Zin = dielectricData.ZpTE(ind);
    end
    [W(ind), dz(ind), S(ind, :, :), dzud(ind), Zprev(ind)] = ...
        ADLSlabDesignMS(Zin, ZL, dy(ind), h(ind), N(ind), erHost, m, fr0,...
        th, ph, tol(ind), "TE", M, ind, shift(ind));
end
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
    Zprev1 = zetad;
    for iter = M:-1:1
    %Susceptance
        if(shift(iter) > 0)
            B_SI = suscpetance_SIShifted(omega, dy(iter), dz(iter), m, W(iter), shift(iter));
            B_I = suscpetance_IShifted(omega, dy(iter), dz(iter), m, W(iter), shift(iter));
        else
            B_SI = suscpetance_SI(omega, dy(iter), dz(iter), m, W(iter));
            B_I = suscpetance_I(omega, dy(iter), dz(iter), m, W(iter));
        end

        Zedge = -1j./B_SI.*(1./(1-((sin(th).^2)./2)));
        %Middle impedance
        Zinf = -1j./B_I.*(1./(1-((sin(th).^2)./2)));

        Zprev1 = effectiveImpMS(N(iter), Zprev1, Zedge, Zinf, ZTx,...
            dzud(iter), dz(iter), k, M, iter);
        %Zreq(indF) = Zprev1;
    end 
    [S11(indF), S12(indF)] = refCoeff(zeta0, Zprev1);    
end
%Plotting
figure(1);
plot(freq./10^9, abs(S11).^2, 'LineWidth', 1.5, 'DisplayName', '|\Gamma|^2'); hold on;
%plot(freq./10^9, abs(S11te).^2, 'LineWidth', 1.5, 'DisplayName', 'S11(TE)');
%plot(freq./10^9, abs(S12tm).^2, 'LineWidth', 1.5, 'DisplayName', 'S12(TE)');
plot(freq./10^9, abs(S12), 'LineWidth', 1.5, 'DisplayName', '|T|^2');
title(['|S11|^2 and |S12|^2; ADL with slabs =', num2str(M)]);
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
title('Reflection Coeff of ADL slabs');
grid on;
% legend show;
hold off;
ylim([-40, 0]);

% figure(3);
% plot(freq./10^9, abs(S11), 'LineWidth', 1.5, 'DisplayName', '|\Gamma TE|'); hold on
% %plot(freq./10^9, abs(gammaTM_tx), 'LineWidth', 1.5, 'DisplayName', '|\Gamma TM|');
% %plot(th./drad, abs(tCoeffTE_tx) + abs(gammaTE_tx), 'LineWidth', 1.5, 'DisplayName', '|Total|');
% xlabel('Frequency (GHz)');
% ylabel('|\Gamma|');
% title('Reflection Coeff of ADL slab');
% grid on;
% legend show;
% hold off;

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