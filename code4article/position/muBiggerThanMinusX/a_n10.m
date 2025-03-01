function a10 = a_n10(n,mu,x,s,A_nMeno1, B_nMeno1)
    %% fromula (3.61a) (3.61b)
    A2 = 0;
    for k = 0:(n-2)
        fun = @(x,k,n,mu,s) (x-(n-1)*mu).^k.*exp(-2.*x./s);
        gamma_k = integral(@(x) fun(x,k,n,mu,s), x+(n-1)*mu, inf);
        % equivalently, but with some numerical issue
        %gamma_k = gammainc(2*x/s,k+1,'upper')*exp(-2*(n-1)*mu/s)*(s/2)^(k+1);
        A2 = A2 + B_nMeno1(2,k+1)*gamma_k;
    end

    a10 = A2;
    for k = 0:(n-2)
        gamma_k = 0;
           for d=0:k
               gamma_k = gamma_k - (x)^(k-d)*exp(-2*(x+(n-1)*mu)/s)*factorial(k)/factorial(k-d)*(s/2)^(d+1);
           end   
        a10 = a10 + A_nMeno1(1,k+1)*( (x)^(k+1)/(k+1) + (-1)^k*factorial(k)*(s/2)^(k+1) )...
                  + B_nMeno1(1,k+1)*gamma_k;
    end
    a10 = exp(-mu/s)*a10;

end