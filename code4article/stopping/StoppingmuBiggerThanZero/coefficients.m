function [A_n,B_n,C_n] = coefficients(h,n,mu,s)
   

   %% Matrix coefficients:
   % -------------------------
   %       i° row : ´i´ intervall in the partition
   %    p° row : coefficien of (y-nk)^p
   %    dimension : (L X \bar{m}_{n,1})
    
   % smaller integer such that ´L*k>h´. Number of bins in the partition
   L = ceil(h/mu);
    
   % max degree of polynomial among the functions P_i(n,w) 
   m_n_1 = min(n,L)-1;  
   bar_m_n_1 = m_n_1+1;


   %% Base case
   if n == 1
       A_n = transpose([repelem(exp((-h+mu)/s),L-1)/2,                    0]);
       B_n = transpose([repelem(0,                  L-1), -exp((+h-mu)/s)/2]);
       C_n = transpose([repelem(0,                  L-1),                 1]);
   %% Recursion
   else
       % make coefficients of the previous steo
       [A_n_pre,B_n_pre,C_n_pre] = coefficients(h,n-1,mu,s);
       % initialize matrix of coefficients
       A_n = zeros([L,bar_m_n_1]);
       B_n = zeros([L,bar_m_n_1]);
       C_n = zeros([L,1]);
       for i = 1:L
           i_star = i+1;
           m_nMeno1_i_star = min(n-1,L-i_star+1);
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
           if i < L
               C_n(i) = C_n_pre(i_star);
           end
       end

       %% Compute coefficients for j=0 (first column)
       a_ni0s = alpha_ni0(h,n,L,mu,s, A_n_pre,B_n_pre,C_n_pre);
       b_ni0s = beta_ni0(h,n,L,mu,s, A_n_pre,B_n_pre,C_n_pre);
       A_n(:,1) = a_ni0s;
       B_n(:,1) = b_ni0s;
   end
end




   