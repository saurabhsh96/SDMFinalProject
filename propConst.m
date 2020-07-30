%Fining propagation constants
function [kx, ky, krho, kz] = propConst(k0, th, ph)
    kx = k0.*sin(th).*cos(ph);
    ky = k0.*sin(th).*sin(ph);
    
    krho = sqrt(kx.^2 + ky.^2); 
    
    kz = -1j*sqrt(-((k0.^2)-(krho.^2)));
end