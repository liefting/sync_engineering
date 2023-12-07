function [solutions,pars,tspan] = simulations(omega,which_oscillator,reltol,abstol,x_par,a0,Dtaus1,z0_range,M,K,tspan,npartau,npark)
    % Bonnie, 2-2020
    % Solves the system with phase lag (in the feedback) with dde23 solver.

    [~,wa] = size(Dtaus1);

    solutions(wa).solver = [];solutions(wa).history = [];
    solutions(wa).discont = [];solutions(wa).x = [];
    solutions(wa).y = [];
    solutions(wa).stats = [];
    solutions(wa).yp = [];

    options = ddeset('RelTol',reltol, 'AbsTol', abstol); 
    pars.K = K;
    pars.a0 = a0;
    pars.M = M;
    indx = 1;    
    
    taus = x_par(indx,1:npartau);
    ks = x_par(indx,npartau+1:npartau+npark);
    
    z0 = z0_range(:,1); % use only one initial condition
    
    len_alphas = length(Dtaus1);

    for inda = 1:len_alphas
      Dtau = Dtaus1(inda);
      alpha = Dtau*omega;   
      lags = taus + alpha;
      
      % create history (without feedback) if running new simulation
        K0=0;
        sol_history=dde23(@(t,z,L) oscillators(t,z,L,which_oscillator, ...
        K0,ks,a0,M) ,lags,z0,[-100,0],options);  % K=0
      
      % apply feedback
        solutions(inda)=dde23(@(t,z,L) oscillators(t,z,L,which_oscillator, ...
        K,ks,a0,M) ,lags, sol_history,tspan,options);
    end   
end

