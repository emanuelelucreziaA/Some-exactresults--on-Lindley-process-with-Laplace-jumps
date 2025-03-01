function termini_noti = termini_noti_B_s(n,k,x,s, A_n_pre,B_n_pre,c_n_pre)
   %{======================================================================
   %{ Computation of the coefficient `b_ni0` that in f_n^i(u) multiply e^{-u}
   %{======================================================================
    
   termini_noti = [];
   % Compuation of 'b_{n,i,0}'. We start from 2 since 'b_n10=0'
   for i = 2:(n+2)
       B = 0;
       m_ni = min(i-1,n-1)-1;
       %% 1° sommation - cicle on domain of integration
       for j=1:(i-2)
           m_nj = min(j,n-1)-1;
           a_nj_pre = A_n_pre(:,j);
           b_nj_pre = B_n_pre(:,j);
           % estremi di integrazione
           if ((i==(n+2)) && (j==(i-2)))
                Ij_l = (n-1)*k;
                Ij_r = (n-1)*k+x;
           else
                Ij_l = (j-1)*k;
                Ij_r = (j*k);
           end
           % 1° somma annidiata in B - ciclo sui gradi del monomio (y-kn)
           for p=0:m_nj
               potenza_jp = ( (Ij_r-(n-1)*k)^(p+1) - (Ij_l-(n-1)*k)^(p+1) )/(p+1);
               % cambio x ---> -y/2+k, -2*(x-(n-1)*k) = y, si invertono quindi gli estremi di intregrazione 
               fun = @(x,p,n,k,s) (x-(n-1)*k).^p.*exp(2.*x./s);
               gamma_jp = integral(@(x) fun(x,p,n,k,s), Ij_l, Ij_r);
               %gamma_jp = ( gammainc(-2/s*(Ij_l-(n-1)*k),p+1) - gammainc(-2/s*(Ij_r-(n-1)*k),p+1) )*gamma(p+1)*(-1)^p*exp(2*(n-1)*k/s)*(s/2)^(p+1);
               B = B + a_nj_pre(p+1)*gamma_jp + b_nj_pre(p+1)*potenza_jp;
           end
       end    
       %% %% termini dall'integrale con 'u-k'
       a_ni_pre  = A_n_pre(:,(i-1)); 
       b_ni_pre  = B_n_pre(:,(i-1));
       % estremo di integrazione sinistro precedente u-k
       if i==(n+2)
           Ii_l = (n-1)*k+x;
       else 
           Ii_l = (i-2)*k;
       end
       % ciclo sugli esponenti
       for p=0:m_ni
           potenza_i = (Ii_l-(n-1)*k)^(p+1)/(p+1);
           
           gamma_ip = 0;
           for d=0:p
               gamma_ip = gamma_ip + (Ii_l-(n-1)*k)^(p-d)*exp(2*Ii_l/s)*factorial(p)/factorial(p-d)*(s/2)^(d+1)*(-1)^d;
           end
           %{
           gamma_ip = 0;%(Ii_l-(n-1)*k)^(p)*exp(2*Ii_l/s)*(s/2); 
           for d=0:(p)
                % dall'integrazione per parti
                gamma_ip = gamma_ip + (Ii_l-(n-1)*k)^(d)*exp(2*Ii_l/s)*factorial(p)/factorial(d)*(s/2)^(p-d+1)*(-1)^(p-d);         
           end
           %}

           % dall'integrale valutato nel nodo precedente 'u-k', estremo sx
           B = B - a_ni_pre(p+1)*gamma_ip - b_ni_pre(p+1)*potenza_i;
           % dall'integrale valutato in 'u-k', come estremo sx 
           B = B + b_ni_pre(p+1)*(s/2)^(p+1)*factorial(p);
       end
       b_ni0 = B*exp(k/s)/(2*s);
       % quando al passo prima 'y=0' 
       b_ni0 = b_ni0 + exp(k/s)/(2*s)*c_n_pre;
       termini_noti = [termini_noti, b_ni0];
   end

   % 'b_n10=0'. Anche 'b_{n,n+1,0}=b_{n,n,0}'             
   termini_noti = [0, termini_noti];
end