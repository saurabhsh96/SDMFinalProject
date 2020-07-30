%For multiple slabs
%Function to calculation given parameters when required Zin is given
%Varying W
%This function is essentially same as ADLSlabDesign3 a minor difference in
%calculating effective impedance
function [W, dz, S, dzud, Zprev] = ADLSlabDesignMS(Zin, ZL, dy, h, N, erHost,...
    m, fr0, th, ph, tol, type, TotSlabs, SlabNum, shift)
    %Initial setups
    omega = 2.*pi.*fr0;
    c = 3e8;
    lam = c./(fr0.*sqrt(erHost));
    k = 2.*pi./lam;
        
    %Propagation constant
    kx = k.*sin(th).*cos(ph);
    ky = k.*sin(th).*sin(ph);
    kRho = sqrt(kx.^2 + ky.^2);
    kz0 = -1j*sqrt(-((k.^2)-(kRho.^2)));

    %Tx Line Impedance
    zeta0 = 120.*pi./sqrt(erHost);
    if(type == "TE")
        ZTx = (zeta0.*k)./kz0;
    else
        ZTx = (zeta0.*kz0)./k;
    end
    
    %ADL Parameters
    dz = h./(N); %Infinitesimly thin line
    dzud = dz./2;
    
    %dz = W;
    %Variation of dy (It is strictly sub-lam0)
    %A question why is is sub-lambda?
    W_iter = 0.0001*lam:0.0001*lam:0.08*lam;
    
    for ind = W_iter
        %Susceptance
        if(shift > 0)
            B_SI = suscpetance_SIShifted(omega, dy, dz, m, ind, shift);
            B_I = suscpetance_IShifted(omega, dy, dz, m, ind, shift);
        else
            B_SI = suscpetance_SI(omega, dy, dz, m, ind);
            B_I = suscpetance_I(omega, dy, dz, m, ind);
        end
        
        if(type == "TE")
            %Edge impedance
            Zedge = -1j./B_SI.*(1./(1-((sin(th).^2)./2)));
            %Middle impedance
            Zinf = -1j./B_I.*(1./(1-((sin(th).^2)./2)));
        else
            Zedge = -1j./B_SI;
            Zinf = -1j./B_I;
        end
        
        %Matrices
        %Semi infinite
        matEdge = ABCD_Z(Zedge); 
        
        %Infinite
        matBet = ABCD_Z(Zinf);

        %Air in between
        matTx = ABCD_TxLine(ZTx, k, dz);

        %Edge air matrix
        matTxud = ABCD_TxLine(ZTx, k, dzud);

        %Final Matrix
        %N+1 represents that the structure is entirely surrounded by the Tx
        %lines. Here, a doubt is whether to take the air side also
        %srrounded if the host is going to be the air itself, let's solve
        %that later! First let's focus on the structure surrounded
        %entirely!
        Fmat = (matEdge^2)*(matBet^(N-2))*(matTx^(N-1))*(matTxud^2); %*matEnd;
        %Converting ABCD parameters to S parameters
        S = ABCDtoS(Fmat, ZTx);
        
        % This part we calculate the impedance overall at point A, let's
        % see if it works?!
        % Method of finding effective impedance
        [Zprev] = effectiveImpMS(N, ZL, Zedge, Zinf, ZTx, dzud, dz, k, TotSlabs, SlabNum);
        %[Zprev] = effectiveImp(N, ZL, Zedge, Zinf, ZTx, dzud, dz, k);
        %[Zprev1] = effectiveImp1(N, ZL, Zedge, Zinf, ZTx, dzud, dz, k);
        %[epsilon, mu, eta, n, zeta] = StoEpsilon(S, th, h, k, "TE");
        %disp(epsilon(1));
        if(abs(imag(Zprev))<tol*100*5 & abs(abs(Zin) - real(Zprev))<tol*abs(Zin)) 
            W = ind;
            break;
        else
            continue;
        end
    end
end