function termini_noti = termini_noti_A_s(n,k,x,s, A_n_pre,B_n_pre,c_n_pre)
    %{ ====================================================================
    %{ Computation of the coefficient `a_ni0` that in f_n^i(u) multiply e^u
    %{ ====================================================================
   termini_noti = [];
   b_n = B_n_pre(:,n+1);
   %% calcolo dei 'a_{n,i,0}'
   % the cycle arrives untill `i=n+1` because `a_{n,n+2,0}=0`
   for i = 1:(n+1)
       A=0;
       % max degree of monomials
       m_ni = min(i-1,n-1)-1;
       if i==1
          a_ni = zeros(n+1,1);
          b_ni = zeros(n+1,1);
       else
          a_ni = A_n_pre(:,i-1);
          b_ni = B_n_pre(:,i-1);
       end
   
       %% terms from the integration with extreme 'u-k'
       if i==(n+1)
            Ii_r = (n-1)*k+x;
       else
            Ii_r = (i-1)*k;
       end
       
       for p=0:m_ni
           power_ip = (Ii_r-(n-1)*k)^(p+1)/(p+1);
           gamma_ip = 0;
           for d=0:p
               gamma_ip = gamma_ip - (Ii_r-(n-1)*k)^(p-d)*exp(-2*Ii_r/s)*factorial(p)/factorial(p-d)*(s/2)^(d+1);
           end
       % dall'integrale valutato nel nodo seguente 'u-k', estremo dx
       A = A + a_ni(p+1)*power_ip + b_ni(p+1)*gamma_ip;
       % dall'integrale valutato in 'u-k', come estremo dx
       A = A + a_ni(p+1)*(s/2)^(p+1)*factorial(p)*(-1)^p;
       end
       

       %% 1° summation - cicle on domain of integration
       for j = i:(n)
            a_nj = A_n_pre(:,j);
            b_nj = B_n_pre(:,j);
            % max degree of polynomial
            m_nj = min(j,n-1)-1;
            % extremes of integration
           if ( j==n ) 
                Ij_l = (n-1)*k;
                Ij_r = (n-1)*k+x;
           else
                Ij_l = (j-1)*k;
                Ij_r = (j*k);
           end
          
           for p = 0:m_nj
               power_jp = ((Ij_r-(n-1)*k)^(p+1)-(Ij_l-(n-1)*k)^(p+1))/(p+1);
               % change x ---> y/2+k, 2*(x-(n-1)*k)=y. Con nota: la
               % 'gammainc(c,p)' considera l'integrale tra 0 e c, anche se c<0 ---> quindi c'è un meno 
               fun = @(x,p,n,k,s) (x-(n-1)*k).^p.*exp(-2.*x./s);
               gamma_jp = integral(@(x) fun(x,p,n,k,s), Ij_l, Ij_r);
               % this below gives numerical problems, please consider also using the equivalent function ´gammaMia´
               %gamma_jp = (gammainc(2*(Ij_r-(n-1)*k)/s,p+1) - gammainc(2*(Ij_l-(n-1)*k)/s,p+1))*gamma(p+1)*(s/2)^(p+1)*exp((-2*(n-1)*k)/s);
               A = A + a_nj(p+1)*power_jp + b_nj(p+1)*gamma_jp;
           end
       end


       %% 2° sommation - cicle on exponents
       b_n1 = B_n_pre(:,n+1);
       for p = 0:(n-2)
           % Nota: '2x = 2*((n-1)*k+x-(n-1)*k)'
           Gamma_n1p = (gammainc(2*x/s,p+1,'upper')*gamma(p+1)*(s/2)^(p+1))*exp((-2*(n-1)*k)/s);
           A = A + b_n1(p+1)*Gamma_n1p;
       end
       a_ni0 = A*exp(-k/s)/(2*s);
       
       % when at the previous step `y=0` 
       if i == 1
            a_ni0 = a_ni0 + exp(-k/s)/(2*s)*c_n_pre;
       end
       termini_noti = [termini_noti, a_ni0];
   end
   
   termini_noti = [termini_noti, 0];
end
