format long 
% add to path the folders up 2 levels
addpath(fileparts(pwd))
addpath(fileparts((fileparts(pwd))))
% add to path the folders below 1 levels
addpath(genpath(pwd)) % useless


%% Parameters
series_length=40;
num_samples=1000000;
% Mean of incrementes - trend of Linldey process
media=0;
% standard error of increments - it must be positive
s=1;
% initial condition of the process
posizione_init=0;
% barrier - it must be positive
h=2;
% greatest integer for which we want the theorical probability of stopping
n_max=40;


%% Empirical distribution
FETs = repelem(inf,num_samples);
for sim = 1:num_samples
    Wn = posizione_init;
    for i = 1:series_length
        Zn = laprnd(1, 1, media, s*sqrt(2));
        Wn = max(Wn,0) + Zn;
        if Wn>h
            FETs(sim) = i;
            break;
        end
    end
end

edges = 0:series_length;
l2=histogram(FETs,edges, 'Normalization','probability'); hold off;
valoriEmpririci = l2.Values;
bar(valoriEmpririci(2:end)); hold on


%% Theorical distribution
ps_teo = zeros(1, n_max);

if media > 0 && media < h   
    fprintf("0<mu<h \n");
    for i = 1:n_max
            [ps_teo(i),A_n,B_n,C_n] = ProbN(h, i, media, s, posizione_init);
    end
elseif media > 0 && media >= h
    fprintf("h=<mu \n");
    for n = 1:n_max
        [ps_teo(n), b_n] = ProbN_muBiggerThanH(h, n, media, s, posizione_init);
    end
elseif media < 0 && media > -h
    fprintf("-h<mu<0 \n");
    for n = 1:n_max
        ps_teo(n) = ProbN_muSmallerThan0(h, n, media, s, posizione_init);
    end
elseif media < 0 && media <= -h
    fprintf("mu<=-h \n");
    for n = 1:n_max
        ps_teo(n) = ProbN_muSmallerThanMinusH(h, n, media, s, posizione_init);
    end
elseif media == 0
    fprintf("mu=0 \n");
    for n = 1:n_max
        ps_teo(n) = ProbN_muEquals0(h, n, s, posizione_init);
    end
end



%% Plotting
for i=1:n_max
    plot(i,ps_teo(i),'.r','MarkerSize',30); hold on;
end

legend('empirical', 'exact')
xlabel('n')
ylabel('P[N=n]')
title(['h:',num2str(h),', \mu:',num2str(media),', \sigma:', num2str(s), ', x:', num2str(posizione_init)])
