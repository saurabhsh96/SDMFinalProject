%RefCoeff
function [R, T] = refCoeff(Za, Zb)
    R = (Zb - Za)./(Zb + Za);
    T = 1 - abs(R).^2;
end