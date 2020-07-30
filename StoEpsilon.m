%Function implementing Cohen formulae
%Formulae are for all angles and in this particular scenario, we are
%implementing the system only for th = 0;
function [epsilon, mu, eta, n, zeta] = StoEpsilon(S, th, h, k0, type)
    eta1 = sqrt(((1+S(1,1))^2 - (S(2,1))^2)./((1-S(1,1))^2 - (S(2,1))^2)).*sec(th);
    eta2 = -sqrt(((1+S(1,1))^2 - (S(2,1))^2)./((1-S(1,1))^2 - (S(2,1))^2)).*sec(th);
    zeta1 = S(2,1)./(1 - S(1,1).*(eta1.*cos(th) - 1)./(eta1.*cos(th) + 1));
    zeta2 = S(2,1)./(1 - S(1,1).*(eta2.*cos(th) - 1)./(eta2.*cos(th) + 1));
    if abs(real(eta1))>=0.005 & real(eta1)>=0
        zeta=zeta1;
        Z2=eta1;
    end
    if abs(real(eta1))>=0.005 & real(eta1)<0
        zeta=zeta2;
        Z2=eta2;
    end
    if abs(real(eta1))<0.005 & abs(zeta1)<=1
        zeta=zeta1;
        Z2=eta1;
    end
    if abs(real(eta1))<0.005 & abs(zeta1)>1
        zeta=zeta2;
        Z2=eta2;
    end
    ni=-1/(k0*h)*1j*real(log(zeta)); % imag of n
    nr=1/(k0*h)*(imag(log(zeta))+2*(0)*pi); %real of n         n2(it,itt)=nr(it,itt)+ni;
    %n2 = -sqrt(((log(abs(zeta)) + 1j*angle(zeta))./(-1j.*k0.*h))^2 + (sin(th))^2);
    n2 = nr + ni;
    Er2=n2/Z2;
    Mr2=n2*Z2;
    
    n = n2;
    eta = Z2;
    epsilon = Er2;
    mu = Mr2;
%     if(type=="TE")
%         zeta = S(2,1)./(1 - S(1,1).*(eta.*cos(th) - 1)./(eta.*cos(th) + 1));
%     else
%         zeta = S(2,1)./(1 - S(1,1).*(eta./cos(th) - 1)./(eta./cos(th) + 1));
%     end 
    %In ntm, zetatm is used and in nte zetate is used, but here as both the
    %formulae are same for th = 0, I am defining only one!!
%     n = sqrt(((log(abs(zeta)) + 1j*angle(zeta))./(-1j.*k0.*h))^2 + (sin(th))^2);
    %Here, as TE = TM as th = 0, eps(1) = eps(2) otherwise, eps(1) is TM
    %and eps(2) is TE components of given formulae.
%     n0 = n; %In other cases, n for theta and n for th = 0 will differ.
%     epsilon(1) = n./eta;
%     epsilon(2) = epsilon(1);
%     epsilon(3) = epsilon(1).*(sin(th))^2./((sin(th))^2 - n^2 + n0^2); %TM values are used
    
%     mu(1) = n.*eta; %otherwise TE;
%     mu(2) = mu(1);
%     mu(3) = mu(1).*(sin(th))^2./((sin(th))^2 - n^2 + n0^2); %TE values are used otherwise
end