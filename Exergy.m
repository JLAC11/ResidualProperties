clc, clear, close all, format compact
T0 = 298.15; % K
P0 = 1; % bar
R = 8.31447;
data = readtable('comps.csv')
parameters = struct_wrapper(data);
% %% First case:
% P = 1; % bar
% T = 698.15; % K
% Exergy1 = exergySRK(T, T0, P, P0, parameters) %NOTA: Sale negativo porque el sistema recibe energía

%% Second case:
P = 5; % bar
T = 698.15; % K
Exergy2 = exergySRK(T, T0, P, P0, parameters);

function parameters = struct_wrapper(Data)
    parameters = struct;
    parameters.Pc = Data.Pc';
    parameters.Tc = Data.Tc';
    parameters.w = Data.w';
    parameters.DH = Data.DH';
    parameters.DS = Data.DS';
    parameters.y = Data.y';
    parameters.Cp = [Data.CpA Data.CpB Data.CpC Data.CpD];
    parameters.ExCh = Data.ExCh';
end

function ValidateData(data,tol)
    if nargin == 1
        tol = 1e-5;
    end
    if abs(sum(data.y)-1)>tol
        error('ERROR: Datos de fracción molar no consistentes.')
    end
    Names = {{'Pc'};{'Tc'};{'w'};{'DH'};{'DS'};{'y'}};
    Varnames = data.Properties.VariableNames;
%     if sum(strcmp(Names,Varnames))<6 || sum(contains(Varnames,'Cp'))==0
%         error(sprintf(['Error: no se encuentran las columnas requeridas.\n'...
%             'Las columnas requeridas son las siguientes:\n'...
%             '\t1. Pc - Presion critica\n'...
%             '\t2. Tc - Temperatura critica\n'...
%             '\t3. w  - Factor acentrico\n'...
%             '\t4. DH - Entalpia de formacion\n'...
%             '\t5. DS - Entropia de formacion\n'...
%             '\t6. y  - Fraccion molar\n'...
%             '\t7. Cp - Columnas que describan al polinomio caracteristico.\n'...
%             '\t        Diseñado sin requerir modificaciones para CpA, CpB, CpC y CpD.\n'
%             ])
%     end
end
