%Impedences
function [Zte, Ztm] = imped(zeta, k0, kz)
    Zte = (zeta.*k0)./(kz);
    Ztm = (zeta.*kz)./k0;
end