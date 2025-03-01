function y = f_Tn_muMaggioreMenoX(n,u,mu,x,s, A_n, B_n, c_n)
   
   %% formulae (3.60)-(3.62)
   if nargin == 5
        [A_n, B_n, c_n] = f_Tn_coefficients_muMaggioreMenoX(n,mu,x,s);
   end

   if u < 0
       y = 0;
   elseif u == 0
       if x < -n*mu
           %% formula (3.118)
           y = c_n(1);
       else
            %% formula (3.120)
           y = c_n(2);
       end
   else
       p = 0:(n-1);
       monomi = transpose((u-n*mu).^p);
       
       % find the bin of the partition where ´u´ belongs to
       if u < x+n*mu
           bin_index = 1;
       else
           bin_index = 2;
       end

       %% formula (3.58)
       y = (1/(2*s)^n)*( exp(-u/s)*mtimes(B_n(bin_index,:), monomi) + exp(u/s)*mtimes(A_n(bin_index,:), monomi));
   end
end