function [uu, hh, ss, theta, phi, A, B] = residualRK(Pr, Tr)
    % residualRK: calculates residual properties for the Redlich-Kwong equation of
    % state, given the reduced temperature and pressure
    %
    % Syntax: [uu, hh, ss, theta, phi] = residualRK(Pr,Tr)
    %
    % @param {float} Pr - reduced pressure
    % @param {float} Tr - reduced temperature
    %
    % Returns
    %   uu {float} - residual internal energy
    %   hh {float} - residual enthalpy
    %   ss {float} - residual entropy
    %   theta {float} - !!!add description here!!
    %   phi {float} - fugacity coefficient at reduced temperature and pressure
    %
    A = 0.42748 * (Pr ./ Tr.^(2.5))
    B = 0.08664 * (Pr ./ Tr)
    p = -1
    q = A - B - B.^2
    r = -A .* B

    for i = 1:length(p)

        zeta = roots([1 p q r])
        z(i) = zeta(imag(zeta) == 0)
    end

    uu = -3 * A ./ (2 * B) .* log(1 + B ./ z);
    hh = z - 1 + uu;
    ss = log(z - B) - A ./ (2 * B) * log(1 + B ./ z);
    theta = (A ./ B) .* log(1 + B ./ z);
    phi = exp(z - 1 - log(z - B) - theta)
end
