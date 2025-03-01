function  [pn,A_n,B_n,C_n]=ProbN_muSmallerThan0(h, n, k, s, w)
    
    % check initial position of the walk
    if w<0 || w>=h
        fprintf('Error, w is out of range')

    %% Find bin index of the interval containing ´w´
    end
    % smaller integer such that ´H*k>h´. Number of bins in the partition
    H = ceil(h/-k);
    % index of the bin where the initial position ´w´ lies in
    for i=1:(H)
        if -(i-1)*k <= w && w < -(i)*k
            bin_index = i;
            break
        end
    end
    % ..equivalently
    %bin_index = L-floor((h-w)/k);    

    
    %% monomi
    % max degree
    m_n_1 = (n-1);
    p = 0:m_n_1;
    monomi = transpose((w+(n-1)*k).^p);

    %% Coefficienti
    [A_n, B_n, C_n] = coefficients_muSmallerThan0(h, n, k, s); 

    
    %% probabilità di arresto in n - formula (4.2)
    pn = ( mtimes(A_n(bin_index,:), monomi)*exp(w/s) + mtimes(B_n(bin_index,:), monomi)*exp(-w/s) ) + C_n(bin_index);
end