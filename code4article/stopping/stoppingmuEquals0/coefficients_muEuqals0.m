function [alpha_ns, beta_ns] = coefficients_muEuqals0(h,n,s)
    alpha_ns = zeros(1,n);
    beta_ns = zeros(1,n);
    if n==1
        alpha_ns(1,1) = exp(-h/s);
    else
        [alpha_nMeno1s, beta_nMeno1s] = coefficients_muEuqals0(h,n-1,s);
        for j = 1:n-1
           a=0;
           b=0;
           for k = (j-1):(n-2)
               %% formula (5.95)
               a = a + alpha_nMeno1s(k+1)*(-s/2)^(k-j+1)*factorial(k)/factorial(j);
               %% formula (5.96)
               b = b + beta_nMeno1s(k+1)*(s/2)^(k-j+1)*factorial(k)/factorial(j);
           end
           alpha_ns(1,j+1) = -a;
           beta_ns(1,j+1) = b;
        end
    a = 0;
    b = 0;
    for k = 0:(n-2)
        gamma_bk = 0;
        for i = 0:k
            gamma_bk = gamma_bk - h^i*factorial(k)/factorial(i)*exp(-2*h/s)*(s/2)^(k-i+1);
        end
        a = a + alpha_nMeno1s(k+1)*( (-1)^k*(s/2)^(k+1)*factorial(k) + h^(k+1)/(k+1) ) + beta_nMeno1s(k+1)*gamma_bk;
        b = b +-alpha_nMeno1s(k+1)*(-1)^k*(s/2)^(k+1)*factorial(k) + beta_nMeno1s(k+1)*(s/2)^(k+1)*factorial(k);
    end
    alpha_ns(1,1) = a;
    K0n = s*(alpha_nMeno1s(1,1)+beta_nMeno1s(1,1));
    %K0n = ProbN_muEquals0(h, n-1, s, 0);
    beta_ns(1,1) = b+K0n;
    end
end