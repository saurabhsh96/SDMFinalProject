%BW Calculation
%Take refCoeff in dB
function BW = BWCalc(freq, refCoeff)
    distance = abs(refCoeff(2) - refCoeff(1));
    %Max = max(refCoeff);
    %refCoeff10dB = -10;
    index = find(refCoeff < -10);

    %LowerFreq
    Fl = freq(index(1));

    %HigherFreq
    Fh = freq(index(end));

    %BW (-10 dB BW is being taken) 
    BW = (Fh - Fl);
end