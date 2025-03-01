function  [pn,A_n,B_n,C_n]=ProbN(h, n, k, s, w)
    
    % check initial position of the walk
    if w<0 || w>=h
        fprintf('Error, w is out of range')
    end

    %% Identify bin
    % smaller integer such that ´L*k>h´. Number of bins in the partition
    L = ceil(h/k);
    % index of the bin where the initial position ´w´ lies in
    for i=0:(L-1)
        if h-(i+1)*k <= w && w < h-i*k
            bin_index = L-i;
            break
        end
    end
    % ..equivalently
    %bin_index = L-floor((h-w)/k);    

    
    %% monomi
    % max degree
    m_n_1 =  min(n,L-1+1)-1; 
    p = 0:m_n_1;
    monomi = transpose((w+(n)*k).^p);


    %% Coefficients
    [A_n, B_n, C_n] = coefficients(h, n, k, s); 

    
    %% probabilità di arresto in n - formula (4.2)
    pn = ( mtimes(A_n(bin_index,:), monomi)*exp(w/s) + mtimes(B_n(bin_index,:), monomi)*exp(-w/s) ) + C_n(bin_index);
end