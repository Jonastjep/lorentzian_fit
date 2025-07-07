function S11 = S11_complex_SPL(params, f)
%S11_complex_SPL stands for Single Lorentzian Peak (single resonance) complex S11 curve. It has as parameters the frequency range f, and the parameters p, ordered as p = [f0, k_int, k_ext], respectively the resonant frequency f0, the internal dissipation k_int and the external coupling rate (load loss)
    f0 = params(1);
    k_int = params(2);
    k_ext = params(3);
    S11 = -1 - (2 * k_ext) ./ (1i * (f - f0) - (k_int + k_ext));
end