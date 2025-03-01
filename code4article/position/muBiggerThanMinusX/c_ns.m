function [c1,c2] = c_ns(n,mu,x,s,A_nMeno1, B_nMeno1,c_nMeno1)
    bin_index = (x >= -(n-1)*mu)+1;

    c1=0;

    if bin_index == 2
        for k = 0:(n-2)
            %% (3.117e)
            fun = @(z,k,n,mu,s) (z-(n-1)*mu).^k.*exp(z./s);
            gamma_ak1 = integral(@(z) fun(z,k,n,mu,s), 0, x+(n-1)*mu);
            fun = @(z,k,n,mu,s) (z-(n-1)*mu).^k.*exp(-z./s);
            gamma_bk1 = integral(@(z) fun(z,k,n,mu,s), 0, x+(n-1)*mu);
            c1 = c1 + A_nMeno1(1,k+1)*gamma_ak1 +B_nMeno1(1,k+1)*gamma_bk1;
            %% (3.117f)
            fun = @(z,k,n,mu,s) (z-(n-1)*mu).^k.*exp(2.*z./s);
            gamma_ak2 = integral(@(z) fun(z,k,n,mu,s), 0, x+(n-1)*mu);
            power_bk2 = ((x)^(k+1) - (-(n-1)*mu)^(k+1) )/(k+1);
            c1 = c1 - exp(mu/s)/2*(A_nMeno1(1,k+1)*gamma_ak2+ B_nMeno1(1,k+1)*power_bk2);
            %% (3.117g)
            fun = @(z,k,n,mu,s) (z-(n-1)*mu).^k.*exp(z./s);
            gamma_ak3 = integral(@(z) fun(z,k,n,mu,s), x+(n-1)*mu,-mu);
            fun = @(z,k,n,mu,s) (z-(n-1)*mu).^k.*exp(-z./s);
            gamma_bk3 = integral(@(z) fun(z,k,n,mu,s), x+(n-1)*mu,-mu);
            c1 = c1 + A_nMeno1(2,k+1)*gamma_ak3 +B_nMeno1(2,k+1)*gamma_bk3;
            %% (3.117h)
            fun = @(z,k,n,mu,s) (z-(n-1)*mu).^k.*exp(2.*z./s);
            gamma_ak4 = integral(@(z) fun(z,k,n,mu,s), x+(n-1)*mu,-mu);
            power_bk4 = ((-n*mu)^(k+1) - (x)^(k+1) )/(k+1);
            c1 = c1 - exp(mu/s)/2*(A_nMeno1(2,k+1)*gamma_ak4+ B_nMeno1(2,k+1)*power_bk4);
            %% (3.117i)
            power_ak5 = 0;
            fun = @(z,k,n,mu,s) (z-(n-1)*mu).^k.*exp(-2.*z./s);
            gamma_bk5 = integral(@(z) fun(z,k,n,mu,s),-mu,inf);
            c1 = c1 + exp(-mu/s)/2*(A_nMeno1(2,k+1)*power_ak5+ B_nMeno1(2,k+1)*gamma_bk5);
        end
        %% (3.117j)
        c1 = c1/(2*s)^(n-1);
        c1 = c1 + (1-exp(mu/s)/2)*c_nMeno1(bin_index);
    
    else
    
        for k = 0:(n-2)
            %% (3.117g)
            fun = @(z,k,n,mu,s) (z-(n-1)*mu).^k.*exp(z./s);
            gamma_ak3 = integral(@(z) fun(z,k,n,mu,s), 0,-mu);
            fun = @(z,k,n,mu,s) (z-(n-1)*mu).^k.*exp(-z./s);
            gamma_bk3 = integral(@(z) fun(z,k,n,mu,s), 0,-mu);
            c1 = c1 + A_nMeno1(2,k+1)*gamma_ak3 +B_nMeno1(2,k+1)*gamma_bk3;
            %% (3.117h)
            fun = @(z,k,n,mu,s) (z-(n-1)*mu).^k.*exp(2.*z./s);
            gamma_ak4 = integral(@(z) fun(z,k,n,mu,s), 0,-mu);
            power_bk4 = ((-n*mu)^(k+1) - (-(n-1)*mu)^(k+1) )/(k+1);
            c1 = c1 - exp(mu/s)/2*(A_nMeno1(2,k+1)*gamma_ak4+ B_nMeno1(2,k+1)*power_bk4);
            %% (3.117i)
            power_ak5 = 0;
            fun = @(z,k,n,mu,s) (z-(n-1)*mu).^k.*exp(-2.*z./s);
            gamma_bk5 = integral(@(z) fun(z,k,n,mu,s),-mu,inf);
            c1 = c1 + exp(-mu/s)/2*(A_nMeno1(2,k+1)*power_ak5+ B_nMeno1(2,k+1)*gamma_bk5);
        end
        %% (3.117j)
        c1 = c1/(2*s)^(n-1);
        c1 = c1 + (1-exp(mu/s)/2)*c_nMeno1(bin_index);
    end   


    c2=0;
    for k = 0:(n-2)
         %% (3.119e)
        fun = @(z,k,n,mu,s) (z-(n-1)*mu).^k.*exp(z./s);
        gamma_ak1 = integral(@(z) fun(z,k,n,mu,s), 0, -mu);
        fun = @(z,k,n,mu,s) (z-(n-1)*mu).^k.*exp(-z./s);
        gamma_bk1 = integral(@(z) fun(z,k,n,mu,s), 0, -mu);
        c2 = c2 + A_nMeno1(1,k+1)*gamma_ak1 +B_nMeno1(1,k+1)*gamma_bk1;
        %% (3.119f)
        fun = @(z,k,n,mu,s) (z-(n-1)*mu).^k.*exp(2.*z./s);
        gamma_ak2 = integral(@(z) fun(z,k,n,mu,s), 0, -mu);
        power_bk2 = ((-n*mu)^(k+1) - (-(n-1)*mu)^(k+1) )/(k+1);
        c2 = c2 - exp(mu/s)/2*(A_nMeno1(1,k+1)*gamma_ak2+ B_nMeno1(1,k+1)*power_bk2);
        %% (3.119g)
        power_ak3 = ((x)^(k+1) - (-n*mu)^(k+1) )/(k+1);
        fun = @(z,k,n,mu,s) (z-(n-1)*mu).^k.*exp(-2*z./s);
        gamma_bk3 = integral(@(z) fun(z,k,n,mu,s),-mu, x+(n-1)*mu);
        c2 = c2 + exp(-mu/s)/2*(A_nMeno1(1,k+1)*power_ak3 +B_nMeno1(1,k+1)*gamma_bk3);
        %% (3.119h)
        power_ak4 = 0;
        fun = @(z,k,n,mu,s) (z-(n-1)*mu).^k.*exp(-2*z./s);
        gamma_bk4 = integral(@(z) fun(z,k,n,mu,s),x+(n-1)*mu,inf);
        c2 = c2 + exp(-mu/s)/2*(A_nMeno1(2,k+1)*power_ak4 + B_nMeno1(2,k+1)*gamma_bk4);
    end
    %% (3.117i)
    c2 = c2/(2*s)^(n-1);
    c2 = c2 + (1-exp(mu/s)/2)*c_nMeno1(bin_index);

    
end