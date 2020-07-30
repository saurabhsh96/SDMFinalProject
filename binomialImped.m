%Binomial QWT design
%Works only for single value of th and ph (I guess, not defining Z values
%considering th, ph)
function [ZA, Z, erZ, h, Zp] = binomialImped(ZL, Z0, N, er, th, ph, freq, fr0)
    c = 3e8;
    lam = c./freq;
    k = 2.*pi./lam;
    lam0 = c./fr0;
%     k0 = 2.*pi.*fr0./c;
    
    Z = zeros(N, length(freq));
    Ztemp = zeros(N, length(freq));
    h = zeros(N, 1);
    erZ = zeros(N, 1);
    const = (2.^-N).*(log(ZL./Z0));
    %Finding binomial impedances
    for ind = 1:N
        C = nchoosek(N,ind-1);
        if(ind == 1)
            Z(ind, :) = exp(log(Z0) + C.*const);
            %erZ(ind) = sqrt(er);
        else
            Z(ind, :) = exp(log(Z(ind-1)) + C.*const);
%             erZ(ind) = sqrt(erZ(ind-1));
        end
        %erZ step is pretty much not robust but it is the need of the time
        %can't help! Let's see, later we can write a more robust code!!
        erZ(ind) = (120*pi/Z(ind, find(freq == fr0)))^2;
        h(ind) = lam0./(4.*sqrt(erZ(ind)));
    end
    %Figuring out impedances at each point of connection
    for ind = N:-1:1
        ks1 = k.*sqrt(erZ(ind));
        [~, ~, ~, kzs1] = propConst(ks1, th, ph);
        if(ind == N)
            Ztemp(ind, :) = findZ(Z(ind, :), ZL, kzs1, h(ind));
        else
            Ztemp(ind, :) = findZ(Z(ind, :), Ztemp(ind+1, :), kzs1, h(ind));
        end
        %Point impedance
        Zp(ind) = Ztemp(ind, find(freq == fr0));
    end
    ZA = Ztemp(1, :);
end