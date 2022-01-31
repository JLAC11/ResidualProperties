function ex = exergySRK(T, T0, P, P0, params, R)
    % Calculates exergy for SRK model
    % T:  final temperature, in absolute units
    % T0: initial temperature, in absolute units
    % P:  final pressure
    % P0: initial pressure
    % params: struct with the following fields:
    %     1. Pc: critical pressures of each component, as a row vector
    %     2. Tc: critical temperatures of each component, as a row vector.
    %            Must be in absolute units.
    %     3. y:  molar fractions of components in the mixture
    %     4. w:  acentric factor of each component
    %     5. Cp: heat capacity in J/mol K, which is in polynomial form as a
    %           matrix, where each row corresponds to each term, and each
    %           column to the exponent in the following fashion:
    %       Cp(i,T) = Cp(i,1)*T^n + Cp(i,2)*T^(n-1) + ... + Cp(i,n)*T + Cp(i,n+1)
    %     6. DH: Enthalpy of formation of each component at T0
    %     7. DS: Enthalpy of formation of each component at T0
    %
    % Returns exergy for the system at given conditions.
    if nargin < 6
        R = 8.31447;
    end

    % At initial conditions:
    [DH0, DS0, DG0] = SRKpoints(P0, params.Pc, T0, params.Tc, params.y, params.w);
    DH0 = DH0 * R * T0;
    DS0 = DS0 * R;
    DG0 = DG0 * R * T0;
    % At final conditions:
    [DHF, DSF, DGF] = SRKpoints(P, params.Pc, T, params.Tc, params.y, params.w);
    DHF = DHF * R * T;
    DSF = DSF * R;
    DGF = DGF * R * T;
    % Ideal trajectory
    H = enthalpyIdeal(params.DH, T, T0, params.Cp, params.y);
    S = entropyIdeal_PT(params.DS, T, T0, P, P0, params.Cp, params.y, R);
    % Delta M = -M0' + M* + MF'; from real to ideal, trajectory, from ideal to real
    DeltaH_t =- DH0 + H + DHF;
    DeltaS_t = (S + DSF - DS0);
    ex = DeltaH_t - T * DeltaS_t;
    chem_exergy = DGF - DG0 + params.y * params.ExCh';
    ex = chem_exergy + ex;
end

function [DeltaH, DeltaS, DeltaG] = SRKpoints(P, Pc, T, Tc, y, w)
    Pr = P ./ Pc;
    Tr = T ./ Tc;
    [~, DeltaH, DeltaS, DeltaG] = SRK(Pr, Tr, y, w);
end
