function [uu, hh, ss, theta, phi, A, B] = residualPR(Pr, Tr, w)
    % residualPR - Calculates residual properties for the Peng Robinson equation
    % of state, given the reduced temperature and pressure, and acentric factor.
    %
    % Syntax: [uu,hh,ss,theta,phi] = residualPR(Pr,Tr,w)
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
    coef_acentric = (0.37464 + 1.54226 * w - 0.26992 * w.^2)
    alpha = (1 + coef_acentric * (1 - sqrt(Tr))).^2
    A = 0.45724 * (Pr ./ Tr.^2) * alpha
    B = 0.07780 * (Pr ./ Tr)
    p = -1 + B
    q = A - 2 * B - 3 * B.^2
    r = -A * B + B.^2 + B.^3

    for i = 1:length(p)
        zeta = roots([1 p(i) q(i) r(i)]);
        z(i) = zeta(imag(zeta) == 0);
    end

    gamma = coef_acentric .* (sqrt(Tr ./ alpha));

    uu = -A * (1 + gamma) ./ (B * sqrt(8)) .* log((z + B * (1 + sqrt(2))) ./ (z + B * (1 - sqrt(2))));
    hh = z - 1 + uu;
    ss = log(z - B) - A * gamma ./ (B * sqrt(8)) .* log((z + B * (1 + sqrt(2))) ./ (z + B * (1 - sqrt(2))));
    theta = (A ./ (B * sqrt(8)) .* log((z + B * (1 + sqrt(2))) ./ (z + B * (1 - sqrt(2)))));
    phi = exp(z - 1 - log(z - B) - theta);
end
