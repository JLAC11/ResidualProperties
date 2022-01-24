function H = enthalpyIdeal(DeltaH, T, T0, Cp, y)
    % Returns enthalpy of an ideal gas, given formation enthalpy, initial and
    % final temperatures, mole fraction and heat capacity given as a matrix,
    % where each row corresponds to each component and each column to the power
    % of the polynomial. It only receives heat capacities given as polynomials.
    % The polynomials are given in the following fashion:
    % Cp(i,T) = Cp(i,1)*T^n + Cp(i,2)*T^(n-1) + ... + Cp(i,n)*T + Cp(i,n+1)
    H = 0;
    dCP = deltaCP(Cp,y)
    f = @(T) polyval(dCP,T);
    H = integral(f,T0,T);

end

function dCP = deltaCP(Cp,y)
    dCP = y*Cp; % Sumados con operacion matricial
end
