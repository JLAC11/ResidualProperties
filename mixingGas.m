function [Amix, Bmix] = mixingGas(y, Pc, Tc, P, T, w)
    Pr = P ./ Pc
    Tr = T ./ Tc

    [~, ~, ~, ~, ~, A, B] = residualSRK(Pr, Tr, w)

end
