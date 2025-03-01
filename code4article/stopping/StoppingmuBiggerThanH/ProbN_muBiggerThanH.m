function [p,b_n] = ProbN_muBiggerThanH(h,n,mu,s,w)
       if n==1
            %% formula (xx)
            p = 1-exp((h-mu-w)/s)/2;
            b_n=0;
       elseif n==2
            %% formula (xx)
            b_n = exp((h-mu)/s)/2 - ( h/(4*s) + 1/4 )*exp((h-2*mu)/s);
            p = exp(-w/s)*b_n;
       else 
            %% formula (xx)
            b_2 = exp((h-mu)/s)/2 - ( h/(4*s) + 1/4 )*exp((h-2*mu)/s);
            b_n = b_2*exp(-(n-2)*mu/s)*( h/(2*s)+1/2 )^(n-2);
            p = b_n*exp(-w/s);
       end
    
end