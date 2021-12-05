function deltaS = entropyIdeal_PT(S, T, T0, P, P0, Cp, y,R)
    % Returns change in entropy in J/mol K; requires heat capacity to be in J/mol K.
    % Assumes heat capacity is in polynomial form as a matrix, where each row
    % corresponds to each term, and each column to the exponent
    % in the following fashion:
    % Cp(i,T) = Cp(i,1)*T^n + Cp(i,2)*T^(n-1) + ... + Cp(i,n)*T + Cp(i,n+1)
    % T: end temperature, in absolute units
    % T0: initial temperature, in absolute units
    % P: end pressure
    % P0: initial pressure
    % Cp: heat capacity in J/mol K.
    % S: entropy of component i, listed as a row vector
    % y: mole fraction of component i, listed as a row vector
    if nargin < 8
       R = 8.31447; 
    end
    for i = 1:length(y)
        dS(i) = entropy_noP(Cp(i, :), T, T0, S(i), y(i),R);
    end

    deltaS = sum(dS) -R * log(P / P0);

end

function dS = entropy_noP(Cp, T, T0, S, y,R)
    CpentreT = @(T) polyval(Cp, T) ./ T;
    dS = S - R * y .* log(y) + integral(CpentreT, T0, T);
end
