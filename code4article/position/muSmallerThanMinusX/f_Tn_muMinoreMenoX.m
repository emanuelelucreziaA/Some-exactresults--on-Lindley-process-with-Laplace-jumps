function y = f_Tn_muMinoreMenoX(n,u,mu,x,s,B_n,c_n)
   
   %% formulae (3.84a,3.85a,3.84c)
   if nargin == 5
        [B_n, c_n] = f_Tn_coefficients_muMinoreMenoX(n,mu,x,s);
   end

   if u < 0
       y = 0;
   elseif u == 0
       %% formula (3.83)
       y = c_n;
   else
       p = 0:(n-1);
       monomi = transpose((u-(n-1)*mu).^p);
       %% formula (3.82)
       y = (1/(2*s)^n)*exp(-u/s)*mtimes(B_n, monomi);
   end
end