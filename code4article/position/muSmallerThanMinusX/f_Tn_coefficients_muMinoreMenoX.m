function [B_n, c_n] = f_Tn_coefficients_muMinoreMenoX(n,mu,x,s)
    
    if n == 1
        B_n = exp((mu+x)/s);
        c_n = 1- exp((mu+x)/s)/2;
    else
        B_n = zeros([1,n]); % the max degree is ´n-1´ but we include also the constant term
        [B_nMeno1, c_nMeno1] = f_Tn_coefficients_muMinoreMenoX(n-1,mu,x,s);
    
        %% formula (3.126b)
        for j = 1:(n-1)
            b=0;
            for k = (j-1):(n-2)
                b = b + B_nMeno1(k+1)*(s/2)^(k-j+1)*factorial(k)/factorial(j);
            end
            B_n(1+j) = exp(mu/s)*b;
        end
        
        %% formula (3.126d)
        b=0;
        for k = 0:(n-2)
            b = b + B_nMeno1(k+1)*(-(-(n-2)*mu)^(k+1)/(k+1) + factorial(k)*(s/2)^(k+1));
        end
        B_n(1) = exp(mu/s)*(b + (2*s)^(n-1)*c_nMeno1);

        %% massa in 0
        c = 0;
        for k = 0:(n-2)
            % formula (3.133): this gives numerical error, please also consider using the function ´gammaMia´
            %gamma_k_upper = exp(-2*(n-2)*mu/s)*(s/2)^(k+1)*gammainc(k+1, -2*(n-1)*mu/s, 'upper');
            % formula (3.130)
            fun = @(x,p,n,m,s) (x-(n-2)*m).^p.*exp(-2.*x./s);
            gamma_k_upper = integral(@(x) fun(x,k,n,mu,s), -mu, inf);
            c = c + B_nMeno1(k+1)*exp(-mu/s)/(2^(n)*s^(n-1))*gamma_k_upper;
            % formula (3.134): this gives numerical error, please also consider using the function ´gammaMia´
            %gamma_k_lower = exp(-(n-2)*mu/s)*s^(k+1)*(gammainc(k+1, -(n-1)*mu/s, 'lower')-gammainc(k+1, -(n-2)*mu/s, 'lower'));
            % formula (3.31)
            fun = @(x,p,n,m,s) (x-(n-2)*m).^p.*exp(-x./s);
            gamma_k_lower = integral(@(x) fun(x,k,n,mu,s), 0, -mu);
            c = c + 1/(2^(n-1)*s^(n-1))*B_nMeno1(k+1)*gamma_k_lower;
            % formula (3.129)
            c = c - B_nMeno1(k+1)*exp(mu/s)/(2^(n)*s^(n-1))*( (-(n-1)*mu)^(k+1)/(k+1) - (-(n-2)*mu)^(k+1)/(k+1));
        end
        % formula (3.130)
        c = c + (1-exp(mu/s)/2)*c_nMeno1;
        c_n = c;
    end
end