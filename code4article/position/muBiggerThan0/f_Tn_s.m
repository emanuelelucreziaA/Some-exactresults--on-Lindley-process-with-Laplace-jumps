function y = f_Tn_s(u,n,k,x,s,A_n,B_n,c_n) 

    %% create coefficients
    if nargin == 5
        [A_n,B_n,c_n] = f_Tn_coefficients_s(n,k,x,s);
    end

    %% identify the bin of partition at which `u` belongs to
    found=false;
    if u==0
        y = c_n;
    else
        for i=1:n
            if ((i-1)*k<u) && (u<=i*k)
                found=true;
                break
            end
        end
        if not(found)
            if ((n-1)*k<u) && (u<=n*k+x)
                i = n+1;
            else
                i = n+2;
            end
        end
        
        %% compute output
        % max degree of monomila
        m = min(n,i)-1;
        % create monomial functions
        monomi = [];
        for p=0:m
            monomi = [monomi, (u-n*k)^p];
        end
        if n>i
            monomi = [monomi, repelem(0,n-i)];
        end
        % output
        y = dot(monomi, A_n(:,i))*exp(u/s) + dot(monomi, B_n(:,i))*exp(-u/s);
    end
end

