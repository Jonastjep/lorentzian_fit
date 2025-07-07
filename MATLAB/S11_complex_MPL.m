function S11 = S11_complex_MPL(params, f)
% params: [f0_1, k_int_1, k_ext_1, ..., f0_N, k_int_N, k_ext_N]
    S11 = -1 * ones(size(f)); % baseline reflection
    N = length(params)/3;
    for i = 1:N
        f0     = params(3*(i-1)+1);
        k_int  = params(3*(i-1)+2);
        k_ext  = params(3*(i-1)+3);

        p = [f0,k_int,k_ext];

        S11 = S11 - S11_complex_SPL(p,f);
    end
end