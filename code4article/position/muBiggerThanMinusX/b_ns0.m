function [b10,b20] = b_ns0(n,mu,x,s,A_nMeno1, B_nMeno1,C_nmeno1)
    
    %% fromula (3.61c)-(3.61d)
    b10 = 0;
    for k = 0:(n-2)
        gamma_k=0;
        for d=0:k
            gamma_k = gamma_k + (-1)^d*(s/2)^(d+1)*factorial(k)/factorial(k-d)*(-(n-1)*mu)^(k-d);
        end
       %b10 = b10 + B_nMeno1(1,k+1)*((-(n-1)*mu)^(k+1)/(k+1)+(-1)^k*factorial(k)*(s/2)^(k+1) )...
        b10 = b10 + B_nMeno1(1,k+1)*(-(-(n-1)*mu)^(k+1)/(k+1)+      factorial(k)*(s/2)^(k+1) )...
                  - A_nMeno1(1,k+1)*gamma_k;
    end
    bin_index = (x >= -(n-1)*mu)+1; 
    b10 = exp(mu/s)*(b10 + (2*s)^(n-1)*C_nmeno1(bin_index));


    %% fromula (3.62b) (3.62c)
    b20 = 0;
    if x+(n-1)*mu>0
        A1 = 0;
        for k = 0:(n-2)
            fun = @(z,k,n,mu,s) (z-(n-1)*mu).^k.*exp(2.*z./s);
            gamma_k = integral(@(z) fun(z,k,n,mu,s), 0, x+(n-1)*mu);
            A1 = A1 + A_nMeno1(1,k+1)*gamma_k...
                    + B_nMeno1(1,k+1)*( x^(k+1)/(k+1)-(-(n-1)*mu)^(k+1)/(k+1));
        end
        b20 = A1;
        for k = 0:(n-2)
            %{
            for d=0:k
                gamma_k = gamma_k + (-1)^d*(s/2)^(d+1)*factorial(k)/factorial(k-d)*(x)^(k-d);
            end
            %}
            b20 = b20 + B_nMeno1(2,k+1)*(-(x)^(k+1)/(k+1)+factorial(k)*(s/2)^(k+1) ); ...
                     % - A_nMeno1(2,k)*gamma_k;
        end
        b20 = exp(mu/s)*(b20 + (2*s)^(n-1)*C_nmeno1(bin_index));
    else
        for k = 0:(n-2)
            %{
            for d=0:k
                gamma_k = gamma_k + (-1)^d*(s/2)^(d+1)*factorial(k)/factorial(k-d)*(x)^(k-d);
            end
            %}
            b20 = b20 + B_nMeno1(2,k+1)*(-(-(n-1)*mu)^(k+1)/(k+1)+factorial(k)*(s/2)^(k+1) ); ...
                     % - A_nMeno1(2,k)*gamma_k;
        end
        b20 = exp(mu/s)*(b20 + (2*s)^(n-1)*C_nmeno1(bin_index));
    end

end