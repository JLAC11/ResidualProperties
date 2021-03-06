function [uu, hh, ss, theta, phi, A, B] = SRK(Pr, Tr, x, w, kij)
    % residualSRK - Calculates residual properties for the Soave-Redlich-Kwong
    % equation of state, given the reduced temperature and pressure.
    %
    % Syntax: [uu,hh,ss,theta,phi] = SRK(Pr,Tr,x,w,kij)
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

    if nargin == 4
        kij = zeros(length(Pr));
    end

    % Individual coefficients
    alpha = (1 + (0.48 + 1.574 * w - 0.176 * w.^2) .* (1 - sqrt(Tr))).^2;
    A = 0.42748 * (Pr ./ Tr.^(2.5)) .* alpha;
    B = 0.08664 * (Pr ./ Tr);

    % Applying kay's mixing rule
    Aii = sqrt(A' * A) .* (1 - kij);
    A = x * Aii * x';
    B = x * B';

    p = -1;
    q = A - B - B.^2;
    r = -A .* B;

    for i = 1:length(p)
        zeta = roots([1 p q r]);
        z = zeta(imag(zeta) == 0);
        z = z(real(z)>0);
    end

    uu = -3 * A ./ (2 * B) .* log(1 + B ./ z);
    hh = z - 1 + uu;
    ss = log(z - B) - A ./ (2 * B) .* log(1 + B ./ z);
    theta = (A ./ B) * log(1 + B ./ z);
    phi = exp(z - 1 - log(z - B) - theta);
end
