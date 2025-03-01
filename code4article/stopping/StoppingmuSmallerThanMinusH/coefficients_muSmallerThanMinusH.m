function [alpha_n, eta_n] = coefficients_muSmallerThanMinusH(h,n,mu,s)
    if n==2
        %% formula (xx)
        alpha_n = h*exp((2*mu-h)/s)/(4*s)-exp((2*mu-h)/s)/4;
        eta_n = exp((mu-h)/s)/2;
    else
        [alpha_nMeno1, eta_nMeno1] = coefficients_muSmallerThanMinusH(h,n-1,mu,s);
        %% formula (xx)
        alpha_n = (h*exp(mu/s)/(2*s)-exp(mu/s)/2)*alpha_nMeno1 + (exp(mu/s)/2*(1-exp(-h/s))-exp(mu/s)/2)*eta_nMeno1;
        eta_n = alpha_nMeno1 + eta_nMeno1;
    end
end