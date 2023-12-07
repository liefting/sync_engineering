% Bonnie Liefting, 5/2023
% VERSION FOR 3D FOURRIER
% Script to call function producing prc, and doing numerical simulations of
% the brusselator with up to four oscillators and feedback

% ------------------  fixed stuff ------------------ %
npartau = 19; % first npartau parameters are taus
npark = 19; % last npark parameters are ks
which_oscillator = "Brusselator"; 
%which_oscillator = "FitzHughNagumo";


%%
% ------------------  loading Brusselator/FHN stuff ------------------ %
% getting Brusselators 'isochrones' from Kyles file to translate between phase and state space
[Hint] = load_isochrones(which_oscillator);
%%
[xn,pn,a0,al,meanz,zl,omega,timevec,yvec,t_prc, prc,sol_old,period] = find_and_plot_waveform_and_prc( ...
  which_oscillator); 

%% Target interaction function
% ------------------  defining H(phi)  ------------------ %

% gains
xi1 = 1;
xi2 = 0.5;
xi3 = 1;  % * 
xi4 = 0;
xi5 = 0;

% * Pairwise model from Thesis: xi3 = 0 
%   Nonpairwise model from Thesis: xi3 = 1

% phase shifts
chi1 = -pi/2; 
chi2 = -pi/2; 
chi3 = 0;
chi4 = 0;  
chi5 = 0;

if xi3 == 0
  pairwise = 1;
else
  pairwise = 0;
end

xis = [xi1 xi2 xi3 xi4 xi5];
chis = [chi1 chi2 chi3 chi4 chi5];

fprintf("-------------------------\n")
printstuff(xis,chis);

H100 = 0.5 * xi1 * exp(-1i * chi1);
H110 = 0.5 * xi3 * exp(-1i * chi3);
H200 = 0.5 * xi2 * exp(-1i * chi2); 

H = [H100 H110 H200];

%% Running optimisation to find feedback parameters matching H

% ------------------  running optimisation to find feedback parameters ------------------ %
fprintf("Finding feedback parameters using optimisation...")

tic
[store_x_nonz] = optimisation(al,zl,omega,H,npartau,npark,period,pairwise);
toc

[nsets1, ~] = size(store_x_nonz) ;

fprintf("  done.\n")

fb_pars1 = round(store_x_nonz,4);  % 4 for pairwise brus. 6 for nonpairwise brus? 4 for pairwise FHN
[~,ia,~] = unique(fb_pars1,'rows'); % find indices of those unique up to 4 or 6 dec.
fb_pars = store_x_nonz(ia,:);
[nsets,~] = size(fb_pars);

fprintf("Found "+num2str(nsets)+" unique feedback parameter sets (originally found "+num2str(nsets1)+")\n.");

fprintf("The feedback parameter sets:\n");
%
format short;
[sumsk] = display_taus_ks(fb_pars,50);
format long;

%% ------------------  initial conditions ------------------ %
startat = "sync";
%startat = "splay";
%startat = "antiphase";

M = 4 ; % number of oscillators. Can be 2, 3, 4

% ------------------  perturbing away from initial conditions before starting simulations ------------------ %
rng(1);
if startat == "sync"
    phases1 = zeros(1,M);                % in-phase
elseif startat == "splay"
    if M == 2
      phases1 = [0,pi];    % splay M=2
    elseif M == 3
      phases1 = [0, 2*pi/3, 4*pi/3];    % splay M=3
    else
      phases1 = [0, pi/2, pi, 3*pi/2];   % splay M=4
    end
elseif startat == "antiphase"
    phases1 = [0, 0, pi, pi];          % antiphase M=4
end
    
%phases = phases1(1:M) + 0.001*(rand(no_in_cond,M)-0.5);  % close to -phases1- (fx full sync 0,0,0)
phases = phases1(1:M) + 0.01*(rand(1,M)-0.5);  % close to -phases1- (fx full sync 0,0,0)


sol_old_tend = sol_old.x(end);
sol_old_tzero = sol_old_tend - period;

zxx = deval(sol_old,sol_old_tzero + phases/(2*pi)*period);
zxx2 = deval(sol_old,sol_old_tzero + phases1/(2*pi)*period);
disp(zxx)
zxx = reshape(zxx,[],1);
zxx2 = reshape(zxx2,[],1);
disp(zxx)
fprintf("Starting at "+startat+": \n" + ...
    "(%3.3f  %3.3f  %3.3f  %3.3f) (phase plane)\n",phases);


%%
figure
plot(sol_old.y(1,:),sol_old.y(2,:))
hold on;
plot(yvec(1,:),yvec(2,:),'o')

for ind = 1:M
  plot(zxx(ind*2-1),zxx(ind*2),'ko','Markersize',10,'Markerfacecolor','k'); hold on;
  plot(zxx2(ind*2-1),zxx2(ind*2),'mo','Markersize',10,'Markerfacecolor','m');
end

if pairwise
  str_wise = "pairwise";
else
  str_wise = "nonpairwise";
end

for par_set = 2
    fprintf("Running simulations for N=" + num2str(M) + " " +which_oscillator+...
      " oscillators for feedback parameter set "+num2str(par_set)+" for the "+str_wise+" target model.\n");
    xp = fb_pars(par_set,:); 

    Dtaus = linspace(0.1,2.5,5);
    
    fprintf("Using "+length(Dtaus)+ " different overall lags Dtau between "+min(Dtaus)+" and "+max(Dtaus)+".\n")

    % ------------------  Coupling strength and simulation time ------------------ %
    tend = 100; %1000 ; % time to run simulations 
    tspan = [0, tend];  
    K = 1e-3; % Overall coupling strength
  
    fprintf("Coupling strength is "+ num2str(K) +" and endtime is "+num2str(tend)+".\n")

    strK = num2str(K);
    name = datestr(datetime("today")) +"_" + which_oscillator +"_"+str_wise +"_parset" + ...
      par_set+"_"+startat +"_K"+strK(3:end)+"_alpha0125_N"+num2str(M);

    fprintf("Coupling strength K=%f\n",K);
    % ------------------  run initial simulation ------------------ %
    %defaults: 'RelTol',1e-3, 'AbsTol', 1e-6
    reltol = 1e-5;  % updated 8 aug
    abstol = 1e-8;  % updated 8 aug
    fprintf("reltol = %f, abstol = %f\n",reltol,abstol);
    tic
    z0 = zxx; 

    [solutions,pars,tspan_old]=simulations(omega,which_oscillator,...
      reltol,abstol,xp,a0,Dtaus,z0,M,K,tspan,npartau,npark);
    fprintf("Done initial run. Filename: "+name)
    toc   
    clear fig    
    if startat == "splay"  % degree of order parameters to store
      deg = M;
    else
      deg = 1;  % in-phase sync and undefined
    end    

    if ~isfolder("data")
       mkdir("data")
    end
    save("data/" + name + ".mat")

    fprintf("writing Rk to file(s)")
    nopoints = 1000;
    for deg1= 1:deg
      writematrix([0,Dtaus],"data/"+name+"_R"+deg1+".txt")
      [R_ordersols,tvecR] = find_order_parameter(Hint,solutions,M,deg1,...
        [0 tspan_old(end)],nopoints); % deg = 1 (in-phase) or deg = 4 (splay)
      writematrix([K*tvecR, R_ordersols],"data/"+name+"_R"+deg1+".txt",'WriteMode','append')
    end
    toc
end

%% Local functions

function printstuff(xi,chi)
    xi1 = xi(1); xi2 = xi(2); xi3 = xi(3); xi4 = xi(4); xi5 = xi(5);
    chi1 = chi(1); chi2 = chi(2); chi3 = chi(3); chi4 = chi(4); chi5 = chi(5);
    fprintf("xi  = (%3.3f  %3.3f  %3.3f  %3.3f  %3.3f)\nchi = (%3.3f  %3.3f  %3.3f  %3.3f  %3.3f)\n\n",xi,chi)
    
    eq16 = xi1*sin(chi1) + 2*xi2*sin(chi2) + 2*xi3*sin(chi3) + xi4*sin(chi4) + xi5*sin(chi5);
    fprintf("In phase bifurcates where lhs eq. (16) from BAR2016 is zero, lhs = %f\n",eq16)
    
    eq19 = xi2*sin(chi2);
    eq20 = xi1*sin(chi1);
    fprintf("Splay bifurcates (steady) where lhs eq. (19) from BAR2016 is zero, lhs = %f\n",eq19)
    fprintf("Splay bifurcates (hopfbi) where lhs eq. (20) from BAR2016 is zero, lhs = %f\n",eq20)
    fprintf("Hopf only holds if cos(chi1) should not be zero, it is: %f\n", cos(chi1))
    
    eq21 = xi1*sin(chi1) - 2*xi2*sin(chi2) - xi4*sin(chi4);
    fprintf("Splay antiphase bifurcation where lhs eq. (21) from BAR2016 is zero, lhs = %f\n\n",eq21)
    
    if eq16 < 0
        fprintf("In phase (sync) stable\n")
        %startat = "sync";
    elseif eq16 > 0
        fprintf("In phase (sync) unstable\n")
    else
        fprintf("In phase (sync) neutral\n")
    end
    
    
    if eq19 > 0
        fprintf("Splay stable, checking HB...\n")
        %startat = "splay";
        if cos(chi1) ~= 0
            if eq20 > 0
                fprintf("   ... limit cycle stable (after Hopf of splay)\n")
            elseif eq20 < 0
                fprintf("   ... limit cycle unstable\n")
            else
                fprintf("   ... limit cycle neutral\n")
            end
        else
            fprintf("...cos(chi1) is zero! (regarding hopf bif)\n")
        end
    elseif eq19 < 0
        fprintf("Splay unstable\n")
    else
        fprintf("Splay neutral\n")
    end
    
    if eq21 > 0
        fprintf("Antiphase stable\n")
        %startat = "antiphase";
    elseif eq21 < 0
        fprintf("Antiphase unstable\n")
    else
        fprintf("Antiphase neutral\n")
    end
    fprintf("\n")

end

function [xn,pn,a0,al,z0,zl,omega,t_waveform,waveform,t_prc, prc,sol_old,...
  period] = find_and_plot_waveform_and_prc(which_oscillator)

    % ------------------  getting waveform and prc ------------------ %
    pulse_x = 0.001;
    pulse =  [pulse_x ; 0 ];
    
    [waveform,t_waveform,prc,t_prc,~,sol_old]=find_prc(pulse,which_oscillator); 
    prc = prc(1:end) / pulse(1);
  
    figure; plot(t_waveform,waveform(1,:),'k','linewidth',2); 
    xlabel('$\theta$','Interpreter', 'latex'); ylabel('$x(\phi$)','Interpreter', 'latex');
    
    nFCs = 5;
    xn = zeros(nFCs, 2);
    pn = zeros(nFCs, 2);
    
    period = t_waveform(end);
    omega = 2*pi/period;

    for n = 1:nFCs
      xn(n,1) = integral(@(t) cos(n*omega*t).*interp1(t_waveform', waveform(1,:)', t), 0, t_waveform(end));
      xn(n,2) = integral(@(t) sin(n*omega*t).*interp1(t_waveform', waveform(1,:)', t), 0, t_waveform(end));
      pn(n,1) = integral(@(t) cos(n*omega*t).*interp1(t_prc', prc', t, 'pchip', 'extrap'), 0, t_prc(end));
      pn(n,2) = integral(@(t) sin(n*omega*t).*interp1(t_prc', prc, t, 'pchip', 'extrap'), 0, t_prc(end));
    end
      
    xn = xn/t_waveform(end);
    pn = pn/t_waveform(end);
    
    xf = zeros(size(t_waveform));
    pf = zeros(size(t_prc));
    
    for n = 1:nFCs
      xf = xf + 2*xn(n,1)*cos(n*omega*t_waveform) + 2*xn(n,2)*sin(n*omega*t_waveform);
      pf = pf + 2*pn(n,1)*cos(n*omega*t_prc) + 2*pn(n,2)*sin(n*omega*t_prc);
    end
    
    xf = xf + mean(waveform(1,:));
    pf = pf + mean(prc);
    
    al = (xn(:,1) + 1i*xn(:,2));
    zl = (pn(:,1) + 1i*pn(:,2));
    
    % Plot result
    plot(t_waveform, waveform(1,:), t_waveform, xf, 'o')
    xlabel('t');
    ylabel('x');
    
    yyaxis right
    plot(t_prc, prc, t_prc, pf, 'o')
    xlabel('t');
    ylabel('PRC');
    
    a0 = mean(waveform(1,:));
    z0 = mean(prc);
end

function [sumsk] = display_taus_ks(store_x_nonz,howmany) 
    % function to display feedback parameters found by optimization

    [hig,~]=size(store_x_nonz);
    disppar = min(howmany,hig);

    disp('     parset    tau10     tau01     tau20     tau02     tau30     tau03     tau12     tau21     tau11')
    disp([(1:disppar)',store_x_nonz(1:disppar,[1,2,4,5,7,8,10,11,17])])
    disp('     parset     k10       k01       k20       k02       k30       k03       k12       k21       k11')
    disp([(1:disppar)',store_x_nonz(1:disppar,19+[1,2,4,5,7,8,10,11,17])])

    % sum over k's to estimate size of perturbation
    sumsk = sum(store_x_nonz(:,20:38),2);     
end