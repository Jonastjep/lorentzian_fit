function res = S11_residual(params, f, S11_measured)
    S11_fit = S11_complex_MPL(params, f);
    res = [real(S11_fit - S11_measured); imag(S11_fit - S11_measured)];
end