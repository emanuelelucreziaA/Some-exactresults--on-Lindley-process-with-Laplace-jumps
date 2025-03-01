function a_ni0s = alpha_ni0(h,n,L,mu,s, A_n_pre,B_n_pre,C_n_pre)
    %{ ==========================================================
    %{ Compute the coefficients a_ni0 that P_i[N=n] multiply e^u
    %{ ==========================================================
    a_ni0s = zeros(L,1);
    %% calcolo degli 'a_{n,i,0}'
    for i = 1:L
      A = 0;
      i_star = i+1;
      for j = (i_star+1):L
         % extreme of integration
         Ij_l = max(h-(L-j+1)*mu,0);    
         Ij_r =     h-(L-j)*mu;     
         % max degree at the previous step
         m_nmeno1_j = min(n-1,L-j+1)-1;
         % coefficients in P_j[N=n-1]
         b_nj = B_n_pre(j,:);
         a_nj = A_n_pre(j,:);
         c_nj = C_n_pre(j);
         for p = 0:m_nmeno1_j 
               %integrand = generalizedGamma(p,s/2,(n-1)*mu);
               %gamma_jp = integrand.evaluate_integral(Ij_l,Ij_r);
               % or equivalently
               fun = @(y,p,n,k,s) (y+(n-1)*k).^p.*exp(-2.*y./s);
               gamma_jp = integral(@(y) fun(y,p,n,mu,s), Ij_l, Ij_r);
               power_jp = ((Ij_r+(n-1)*mu)^(p+1)-(Ij_l+(n-1)*mu)^(p+1))/(p+1);
               A = A + a_nj(p+1)*power_jp + b_nj(p+1)*gamma_jp;   
         end
         A = A - c_nj*s*(exp(-Ij_r/s)-exp(-Ij_l/s));
      end
       
     if L > 1 && i+1 <= L
          %extreme of integration
          Iistar_l = max(h-(L-i_star+1)*mu,0);  
          Iistar_r =     h-(L-i_star)*mu;     
          % degree of polynomial
          m_ipiu1_nmeno1 = min((n-1),L-i_star+1)-1; 
          % coefficienti in P_j[N=n-1]
          b_nistar = B_n_pre(i_star,:);
          a_nistar = A_n_pre(i_star,:);
          c_nistar = C_n_pre(i_star);  
          for p = 0:m_ipiu1_nmeno1                                                                   
              power_istarp = ((Iistar_r+(n-1)*mu)^(p+1))/(p+1);
              gamma_istarp = 0;
              for d=0:p
                    gamma_istarp = gamma_istarp - (Iistar_r+(n-1)*mu)^(p-d)*exp(-2*(Iistar_r)/s)*factorial(p)/factorial(p-d)*(s/2)^(d+1);
              end                   
              A = A + a_nistar(p+1)*( power_istarp + (-1)^(p)*factorial(p)*(s/2)^(p+1) ) + b_nistar(p+1)*gamma_istarp;
          end
          A = A - c_nistar*s*exp(-Iistar_r/s);
     end
    A = A*exp(mu/s)/(2*s);
    a_ni0s(i) = A;
    end
end