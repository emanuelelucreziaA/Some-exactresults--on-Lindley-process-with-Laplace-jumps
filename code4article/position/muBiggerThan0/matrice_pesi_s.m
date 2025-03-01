function [Weight_a, Weight_b] = matrice_pesi_s(n,k,s)   
    %{ ==================================================
    %{ Compute the weights of the coefficients a_{n-1,i,j} 
    %{ ==================================================

    %% Create matrix of weights
    % With dimension (n-1)*(n-1): in entrance 'i,p' we found the weight of a_{n-1,i-1,p-1}
    % for a_{n,i,j}, i fixed.
    Weight_a = zeros(n-1,n-1);
    Weight_b = zeros(n-1,n-1);
    for d = 1:(n-1) 
        for p = (d):(n-1) 
            Weight_a(d,p) = [(-1)^(p-1-d)*(s/2)^(p-d)*factorial(p-1)/factorial(d)];
            Weight_b(d,p) = [(s/2)^(p-d)*factorial(p-1)/factorial(d)];
        end
    end
    Weight_a = Weight_a*exp(-k/s)/(2*s);
    Weight_b = Weight_b*exp(k/s)/(2*s);
end