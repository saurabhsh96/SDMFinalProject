%Function to find the impedance
function z = findZ(z0, zd, kz, h)
    z = z0.*(zd + 1j.*z0.*tan(kz.*h))./(z0 + 1j.*zd.*tan(kz.*h));
end