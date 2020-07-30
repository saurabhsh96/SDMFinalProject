%Function to find the effective impedance
function [Zprev] = effectiveImp(N, ZL, Zedge, Zinf, ZTx, dzud, dz, k)
    Imp_Iter = 1:2*N+2;
    Z_iter = zeros(size(Imp_Iter));
    %ZL = Zin;
    Zprev = 0;
    for iter = Imp_Iter
        if(mod(iter, 2)==1)
            if(iter == 1)
                Z_iter(iter) = ZL;
                Zprev = Z_iter(iter);
            elseif(iter == 3 | iter == length(Imp_Iter)-1)
                Z_iter(iter) = Zedge;
                Ztemp = Zprev*Z_iter(iter)./(Zprev+Z_iter(iter));
                Zprev = Ztemp;
            else
                Z_iter(iter) = Zinf;
                Ztemp = Zprev*Z_iter(iter)./(Zprev+Z_iter(iter));
                Zprev = Ztemp;
            end     
        else
            Z_iter(iter) = ZTx;
            if(iter == 2 | iter == length(Imp_Iter))
                Ztemp = findZ(ZTx, Zprev, k, dzud);
                Zprev = Ztemp;
            else
                Ztemp = findZ(ZTx, Zprev, k, dz);
                Zprev = Ztemp;
            end
        end
    end
end