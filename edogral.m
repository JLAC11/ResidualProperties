function [A, B, z] = edogral(T, Tc, P, Pc, w, ee)

    Tr = T / Tc;
    Pr = P / Pc;

    switch ee
        case 'VdW'
            [uu, hh, ss, theta, phi] = residualVdW(Pr, Tr)

            disp('H-H* / RT = ')
            disp(hh)
            disp('U-U* / RT = ')
            disp(uu)
            disp('S-S* / R = ')
            disp(ss)
            disp('Θ= ')
            disp(theta)
            disp('φ= ')
            disp(phi)

        case 'RK'
            [uu, hh, ss, theta, phi] = residualRK(Pr, Tr)
            disp('H-H* / RT = ')
            disp(hh)
            disp('U-U* / RT = ')
            disp(uu)
            disp('S-S* / R = ')
            disp(ss)
            disp('Θ= ')
            disp(theta)
            disp('φ= ')
            disp(phi)
        case 'SRK'
            [uu, hh, ss, theta, phi] = residualSRK(Pr, Tr)

            disp('H-H* / RT = ')
            disp(hh)
            disp('U-U* / RT = ')
            disp(uu)
            disp('S-S* / R = ')
            disp(ss)
            disp('Θ= ')
            disp(theta)
            disp('Θ= ')
            disp(theta)
            disp('φ= ')
            disp(phi)

        case 'PR'
            [uu, hh, ss, theta, phi] = residualPR(Pr, Tr, w)

            disp('H-H* / RT = ')
            disp(hh)
            disp('U-U* / RT = ')
            disp(uu)
            disp('S-S* / R = ')
            disp(ss)
            disp('Θ= ')
            disp(theta)
            disp('φ= ')
            disp(phi)

        case 'Virial'
            phi = residualFOVirial(Pr, Tr, w)
        otherwise
            disp("Option does not exist")
    end

end
