function [uu, hh, ss, theta, phi, A, B] = residualSRK(Pr, Tr, w)
    % residualSRK - Calculates residual properties for the Soave-Redlich-Kwong
    % equation of state, given the reduced temperature and pressure.
    %
    % Syntax: [uu,hh,ss,theta,phi] = residualSRK(Pr,Tr)
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
    alpha = (1 + (0.48 + 1.574 * w - 0.176 * w.^2) * (1 - sqrt(Tr))).^2;
    A = 0.42748 * (Pr ./ Tr.^(2.5)) * alpha;
    B = 0.08664 * (Pr / Tr);
    p = -1;
    q = A - B - B.^2;
    r = -A .* B;

    for i = 1:length(p)
        zeta = roots([1 p q r]);
        z = zeta(imag(zeta) == 0);
    end

    uu = -3 * A ./ (2 * B) .* log(1 + B ./ z);
    hh = z - 1 + uu;
    ss = log(z - B) - A ./ (2 * B) .* log(1 + B ./ z);
    theta = (A ./ B) * log(1 + B ./ z);
    phi = exp(z - 1 - log(z - B) - theta);
end
