function [uu, hh, ss, theta, phi] = PR(Pr, Tr, x, w, kij)
    % residualPR - Calculates residual properties for the Peng Robinson equation
    % of state, given the reduced temperature and pressure, and acentric factor.
    %
    % Syntax: [uu,hh,ss,theta,phi] = PR(Pr,Tr,x,w,kij)
    %
    % @param {float} Pr - reduced pressure of each component
    % @param {float} Tr - reduced temperature of each component
    % @param {float} x - molar fractions of each component
    % @param {float} w - acentric factor
    % @param {float} kij - interaction parameter matrix, size N*N (N components)
    %                      if not given, automatically zeros matrix
    %
    % Returns
    %   uu {float} - residual internal energy
    %   hh {float} - residual enthalpy
    %   ss {float} - residual entropy
    %   theta {float} - !!!add description here!!
    %   phi {float} - fugacity coefficient at reduced temperature and pressure
    %
    if nargin = 4
        kij = zeros(length(Pr))
    end

    % Pre mixing
    coef_acentric = (0.37464 + 1.54226 * w - 0.26992 * w.^2)
    alpha = (1 + coef_acentric * (1 - sqrt(Tr))).^2
    A = 0.45724 * (Pr ./ Tr.^2) * alpha
    Aii = sqrt(A' * A) .* (1 - kij)
    B = 0.07780 * (Pr ./ Tr)
    % Applying kay's mixing rule
    A = x * Aii * x'
    B = x * B'
    % Solving equations
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
