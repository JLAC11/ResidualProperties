function phi = residualFOVirial(Pr, Tr, w)

    % ! TODO, still not finished
    B0 = 0.083 - 0.0422 / (Tr).^(1.6);
    B1 = 0.139 - 0.172 / (Tr).^(4.2);
    phi = exp(Pr ./ Tr .* (B0 + w * B1));
    disp('φ= ')
    disp(phi)
end
