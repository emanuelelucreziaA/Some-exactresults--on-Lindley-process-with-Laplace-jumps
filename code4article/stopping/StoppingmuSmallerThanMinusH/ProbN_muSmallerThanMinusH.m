function  pn=ProbN_muSmallerThanMinusH(h, n, mu, s, x)
if n==1
    %% formula (xx)
    pn = exp((mu-h+x)/s)/2;
else
    %% formula (xx)
    [alpha_n, eta_n] = coefficients_muSmallerThanMinusH(h, n, mu, s);
    pn = alpha_n*exp(x/s)+eta_n;
end