function  pn=ProbN_muEquals0(h, n, s, x)
    %% coefficients
    [alpha_ns, beta_ns] = coefficients_muEuqals0(h, n, s);
    
    %% monomi
    p = 0:(n-1);
    monomi = transpose((x).^p);

    %% probability
    pn = 1/(2^n*s^(n-1))*( mtimes(alpha_ns, monomi)*exp(x/s) + mtimes(beta_ns, monomi)*exp(-x/s) );
end