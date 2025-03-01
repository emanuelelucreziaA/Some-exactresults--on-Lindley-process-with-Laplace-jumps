function [A_n,B_n,C_n] = coefficients_muSmallerThan0(h,n,mu,s)
   

   %% Matrix of coefficients:
   % -------------------------
   %       i° row : interval ´i´ in the partition
   %    p° colonna : coefficients of (y-nk)^p
   %    dimension : (L X \bar{m}_{n,1})
    
   % smaller integer such that ´L*k>h´. Number of bins in the partition
   H = ceil(h/-mu);
    
   % max degree of polynomial among the functions P_i(n,w) e P_i(n-1,w)
   m_n_1 = (n-1); 
   bar_m_n_1 = m_n_1+1; % we add ´1´ for the know term
   
   if n == 1
       %% Base case
       A_n = transpose(repelem(exp((-h+mu)/s)/2,H));
       B_n = transpose(repelem(0,               H));
       C_n = transpose(repelem(0,               H));
       
   else
       %% Recursion - formula (xx)
       % make coefficients at the previos step
       [A_n_pre,B_n_pre,C_n_pre] = coefficients_muSmallerThan0(h,n-1,mu,s);
       
       %% Probability to stop at `n-1 starting from `0`
       % grado massimo
       m_nMeno1_1 = (n-2); 
       %%ì monomi
       p = 0:m_nMeno1_1;
       monomi = transpose((0+(n-2)*mu).^p);
       K0 = ( mtimes(A_n_pre(1,:), monomi) + mtimes(B_n_pre(1,:), monomi) ) + C_n_pre(1);
      
       % initialize matrix of coefficients
       A_n = zeros([H,bar_m_n_1]);
       B_n = zeros([H,bar_m_n_1]);
       C_n = zeros([H,1]);
       for i = 2:H
           i_star = i-1;
           m_nMeno1_i_star = n-1;
           for j = 1:bar_m_n_1-1
               a=0;
               b=0;
               for k = (j-1):m_nMeno1_i_star-1
                   %% formula (4.5)
                   a = a + A_n_pre(i_star,k+1)*(-s/2)^(k-j+1)*factorial(k)/factorial(j);
                   %% formula (4.7)
                   b = b + B_n_pre(i_star,k+1)*(s/2)^(k-j+1)*factorial(k)/factorial(j);
               end
               A_n(i,j+1) = -exp(mu/s)*a/(2*s);
               B_n(i,j+1) = exp(-mu/s)*b/(2*s);
           end
           C_n(i) = C_n_pre(i_star);
          
       end
       %% The case ´i=1´ is treated separately 
       C_n(1) = K0;
        
       % compute the coefficients for j=0 (first column)
       a_ni0s = alpha_ni_muSmallerThan0(h,n,H,mu,s, A_n_pre,B_n_pre,C_n_pre);
       b_ni0s = beta_ni_muSmallerThan0(h,n,H,mu,s, A_n_pre,B_n_pre,C_n_pre);
       A_n(:,1) = a_ni0s;
       B_n(:,1) = b_ni0s;
   end
end




   