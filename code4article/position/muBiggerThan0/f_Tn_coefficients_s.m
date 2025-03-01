function [A_n,B_n,c_n] = f_Tn_coefficients_s(n,k,x,s)
    % Matrix of coefficients:
    % ------------------------
    %       p° row : coefficient of (y-nk)^p
    %    i° column : i interval in the partition
    A_n = []; % in e^u
    B_n = []; % in e^-u
    
    if n == 1
        A_n = transpose([1/(2*s)*exp((-k-x)/s);1/(2*s)*exp((-k-x)/s);0]);
        B_n = transpose([0;0;1/(2*s)*exp((k+x)/s)]);
        c_n = 1/2*exp((-x-k)/s);
    else
        % create ocefficients at the previous step
        if ~exist('coefficients','var')
            [A_n_pre,B_n_pre,c_n_pre] = f_Tn_coefficients_s(n-1,k,x,s);
        end
       % matrix of weight of coefficient `a(b)_{n-1,i,j}` 
       [Wa,Wb] = matrice_pesi_s(n,k,s);
       % computation of coefficients `a_nij`; `i,j>0` (pag.7)
       A_n = mtimes(Wa, A_n_pre);
       A_n = [transpose(repelem(0, n-1)),A_n];
       % computation of coefficients `a_{n-1,i,0}`
       a_ni0s = termini_noti_A_s(n,k,x,s, A_n_pre,B_n_pre,c_n_pre);
       % set the terms `a_ni0` on the first row, i-th column contains the term which multiplies e^{u} in 'f_n^i
       A_n = [a_ni0s; A_n];
       % calcola i coefficienti b_nij; i,j>0 (pag.7)
       B_n = mtimes(Wb, B_n_pre);
       B_n = [transpose(repelem(0, n-1)),B_n];
       % computation of coefficients b_{n-1,i,j}
       b_ni0s = termini_noti_B_s(n,k,x,s, A_n_pre,B_n_pre,c_n_pre);
       % set the terms `b_ni0` on the first row, i-th column contains the term which multiplies e^{-u} in 'f_n^i
       B_n = [b_ni0s; B_n];
       
       %% Massa in zero 
       % it is equal to the first coefficient a_{n,1,0}
       c_n = A_n(1,1)*s; % 
       
    end
end
       