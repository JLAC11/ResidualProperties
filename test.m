clc, clear, close all, format compact
%% Parameters
P0 = 1; % bar
T0 = 298; % K
%  CO    CO2     CH4     H2 from Smith Van Ness
y = [60.175 34.575 3.375 1.875]; % molar fraction
Pc = [34.99 73.83 45.99 13.13]; % bar
Tc = [132.9 304.2 190.6 33.19]; % K
w = [0.048 0.224 0.012 -0.216]; % Acentric factor
% Delta Gf = Delta Hf - Tf Delta Sf => Delta S = ( Delta H - Delta G ) / T
%   CO     CO2   CH4   H2  from Smith Van Ness
DH = [-121.0 -413.8 -89.0 0.0];
DG = [-119.9 -386.0 -34.3 0.0];
DS = (DH - DG) / T0;
% From Cengel, thermodynamics
Cp = [%AT^3   +  BT^2   +   CT    + D; In J/mol K
    -2.222e-10 5.372e-6 1.675e-3 28.16 % CO
    7.469e-9 -3.501e-5 5.981e-2 22.26 % CO2
    -11.01e-9 1.269e-5 5.024e-2 19.89 % CH4
    -8.704e-10 4.003e-6 -1.916e-3 29.11 % H2
    ];

%% Initial state values:
Pr = P0 ./ Pc;
Tr = T0 ./ Tc;
[~, DH0, DS0, ~, phi0] = SRK(Pr, Tr, y, w);

%% First case:
P = 1; % bar
T = 425; % K
Pr = P ./ Pc;
Tr = T ./ Tc;
[~, DHF, DSF, ~, phiF] = SRK(Pr, Tr, y, w);
H = enthalpyIdeal(DH, T, T0, Cp, y);
S = entropyIdeal_PT(DS, T, T0, P, P0, Cp, y);
% Delta M = -M0' + M* + MF'; from real to ideal, trajectory, from ideal to real
Exergy1 = H + DHF - DH0 - T0 * (S + DSF - DS0)

%% Second case:
P = 5; % bar
T = 698.15; % K
Pr = P ./ Pc;
Tr = T ./ Tc;
[~, DHF, DSF, ~, phi] = SRK(Pr, Tr, y, w);
H = enthalpyIdeal(DH, T, T0, Cp, y);
S = entropyIdeal_PT(DS, T, T0, P, P0, Cp, y);
% Delta M = -M0' + M* + MF'; from real to ideal, trajectory, from ideal to real
Exergy2 = H + DHF - DH0 - T0 * (S + DSF - DS0)
