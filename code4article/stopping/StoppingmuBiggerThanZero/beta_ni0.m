function b_ni0s = beta_ni0(h,n,L,mu,s, A_n_pre,B_n_pre,C_n_pre)
    %{ ==============================================================
    %{  Compute the coefficients a_ni0 that P_i[N=n] multiply e^{-u}
    %{ ==============================================================
    b_ni0s = zeros(L,1);
    %% calcolo dei 'b_{n,i,0}'
    % prob di fermarsi in ´n-1´ passi partendo da 0.
    %p0_nmeno1 = ProbN(h,n-1,mu,s,0);
    % equivalently
    m_n_1 =  min(n-1,L-1+1)-1;
    p = 0:m_n_1;
    monomi = transpose((0+(n-1)*mu).^p);
    p0_nmeno1 = ( mtimes(A_n_pre(1,:), monomi) + mtimes(B_n_pre(1,:), monomi) ) + C_n_pre(1);
    
    for i = 1:L
        i_star = i+1;
        B = 0;
        for j = 1:(i_star-1)
            % extreme of integration
            Ij_l = max(h-(L-j+1)*mu,0);   
            Ij_r =     h-(L-j)*mu;  
            % max degree in the polynomial
            m_nmeno1_j = min(n-1,L-j+1)-1;
            % coefficienti in P_j[N=n+1]
            b_nj = B_n_pre(j,:);
            a_nj = A_n_pre(j,:);
            c_nj = C_n_pre(j,:);
            for p = 0:m_nmeno1_j
                % gamma_jp = (gammainc(-2/s*(Ij_l+(n-1)*mu),p+1) - gammainc(-2/s*(Ij_r+(n-1)*mu),p+1))*gamma(p+1)*(s/2)^(p+1)*(-1)^p*exp(-2*(n-1)*mu/s)
                %y = sym('y');
                %gamma_jp = int((y+(n-1)*mu)^p*exp(2*y/s), y,Ij_l, Ij_r);
                % equivalently
                fun = @(y,p,n,k,s) (y+(n-1)*k).^p.*exp(2.*y./s);
                gamma_jp = integral(@(y) fun(y,p,n,mu,s), Ij_l, Ij_r); 
                power_jp = ((Ij_r+(n-1)*mu)^(p+1)-(Ij_l+(n-1)*mu)^(p+1))/(p+1);
                B = B + a_nj(p+1)*gamma_jp + b_nj(p+1)*power_jp;   
            end
            B = B + c_nj*s*(exp(Ij_r/s)-exp(Ij_l/s));
        end
        if L > 1 && i+1<=L
            % extreme of integration
            Iistar_l = max(h-(L-i_star+1)*mu,0);         
            Iistar_r =     h-(L-i_star)*mu;  
            m_istar_nmeno1 = min((n-1),L-i_star+1)-1;
            % coefficients in P_{i+1}[N=n-1]
            b_nistar = B_n_pre(i_star,:);
            a_nistar = A_n_pre(i_star,:);
            c_nistar = C_n_pre(i_star);
            for p = 0:m_istar_nmeno1                                                       
               power_istarp = ((Iistar_l+(n-1)*mu)^(p+1))/(p+1);  
               % gamma_istarp = -(-1)^p*(s/2)^(p+1)*exp(-2*(n-1)*mu/s)*(
               % gammainc(-2/s*(Iipiu1_l+(n-1)*mu),p+1)*gamma(p+1)) % numerical errors
               gamma_istarp = 0;
               for d=0:p               
                   gamma_istarp = gamma_istarp + (Iistar_l+(n-1)*mu)^(p-d)*exp(2*(Iistar_l)/s)*factorial(p)/factorial(p-d)*(s/2)^(d+1)*(-1)^d;
               end
               B = B - b_nistar(p+1)*(power_istarp - (s/2)^(p+1)*factorial(p)) - a_nistar(p+1)*gamma_istarp;                                                                                                
            end
            B = B - c_nistar*s*exp(Iistar_l/s);
        end
        B = B*exp(-mu/s)/(2*s) + p0_nmeno1*exp(-mu/s)/2;
        b_ni0s(i) = B;
    end
end