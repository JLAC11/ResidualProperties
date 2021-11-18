function [uu, hh, ss, theta, phi, A, B] = residualVdW(Pr, Tr)
    % residualVdW - Calculates residual properites for Van der Waals equation of
    % state, given the reduced temperature and pressure
    %
    % Syntax: [uu,hh,ss,theta,phi] = residualVdW(Pr,Tr)
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
    A = (27/64) * Pr / (Tr^2)
    B = Pr / (8 * Tr)
    p = -1 - B
    q = A
    r = -A * B
    z = roots([1 p q r])
    z = z(imag(z) == 0)

    uu = -A / z;
    hh = z - 1 - A / z;
    ss = ln(z - B);
    theta = A / z;
    phi = exp(z - 1 - log(z - B) - theta)
end
