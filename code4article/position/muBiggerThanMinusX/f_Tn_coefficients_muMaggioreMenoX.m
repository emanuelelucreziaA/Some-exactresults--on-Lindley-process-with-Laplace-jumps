function [A_n, B_n, c_n] = f_Tn_coefficients_muMaggioreMenoX(n,mu,x,s)
    

    % initialize matrix of coefficients
    A_n = zeros([2,n]); % the mex degree is ´n-1´ but we include also the constant term
    B_n = zeros([2,n]);
    c_n = zeros([1,2]);

    if n == 1
        %% formule (3.63)-(3.67)
        A_n(1,1) = exp(-(mu+x)/s);
        B_n(2,1) = exp((mu+x)/s);
        c_n(2) = exp(-(mu+x)/s)/2;
    else
        [A_nMeno1, B_nMeno1, c_nMeno1] = f_Tn_coefficients_muMaggioreMenoX(n-1,mu,x,s);
    
        for i = 1:2
            for j = 1:(n-1)
                a=0;
                b=0;
                for k = (j-1):(n-2)
                    %% formula (3.60b)
                    a = a + A_nMeno1(i,k+1)*(-1)^(k+j)*(s/2)^(k-j+1)*factorial(k)/factorial(j);
                    %% formula (3.60d)
                    b = b + B_nMeno1(i,k+1)*(s/2)^(k-j+1)*factorial(k)/factorial(j);
                end
                A_n(i,j+1) = exp(-mu/s)*a;
                B_n(i,1+j) = exp(mu/s)*b;
            end
        end
        
        % Constant terms
        A_n(1,1) = a_n10(n,mu,x,s,A_nMeno1, B_nMeno1);
        [b10,b20] = b_ns0(n,mu,x,s,A_nMeno1, B_nMeno1,c_nMeno1);
        B_n(:,1) = [b10,b20]; 

        %% massa in 0
        [c_n1,c_n2] = c_ns(n,mu,x,s,A_nMeno1, B_nMeno1,c_nMeno1);
        c_n = [c_n1,c_n2];
    end
end