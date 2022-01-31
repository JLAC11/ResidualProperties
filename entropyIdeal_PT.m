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
    % Cp: heat capacity in units of R.
    % S: entropy of component i, listed as a row vector
    % y: mole fraction of component i, listed as a row vector
    if nargin < 8
       R = 8.31447; % J/mol K
    end
    dCP = deltaCP(Cp,y)
    dS = entropy_noP(dCP,T,T0,S,y,R)
    pmix = - R * log(P / P0)
    deltaS = dS + pmix
    

end

function dS = entropy_noP(Cp, T, T0, S, y,R)
    CpentreT = @(T) polyval(Cp, T) ./ T;
    sMix = - R * y * log(y)'
    intgr = integral(CpentreT, T0, T)
    dS = sum(S) + sMix + intgr;
end

function dCP = deltaCP(Cp,y)
    dCP = y*Cp; % Sumados con operacion matricial
end
