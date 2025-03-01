%% TODO
% write unittest

function valore = gammaMia(a,b,p,s,k)
    %{
     Calcola l'integrale di ´exp(-x/s)*(x+k)^p´ tra ´a´ e ´b´.
     Con ´s=1´ e ´k=0´ si calcola la gamma incompleta
    %}
    estremo_superiore_integrale = evaluation_above(b,p,s,k);
    estremo_inferiore_integrale = evaluation_below(a,p,s,k);
    valore = estremo_superiore_integrale - estremo_inferiore_integrale;
end


function estremo_superiore_integrale = evaluation_above(b,p,s,k)
    v = 0;
    for d=0:p
        v = v + factorial(p)/factorial(d)*s^(p-d+1)*( - exp(-b/s)*(b+k)^(d) );
    end
    estremo_superiore_integrale = v;
end


function estremo_inferiore_integrale = evaluation_below(a,p,s,k)
   v = 0;
    for d=0:p
        v = v + factorial(p)/factorial(d)*s^(p-d+1)*( - exp(-a/s)*(a+k)^(d) );
    end
    estremo_inferiore_integrale = v;
end

