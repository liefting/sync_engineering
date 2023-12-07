function [waveform, t_waveform, PRC_t, t_old_phis, T_old, sol_old] = find_prc(pulse,which_oscillator)

    

    rng(1)
    y0=rand(1,2);
    %y0 = [0.5797    0.3801    0.5509    0.7453]
   

    % Solving the system without adding pulses
    % Also calculating maxima to determinge T_old
    if which_oscillator == "Brusselator"
      % finding the PRC by direct method
      options = ddeset('Events',@findMax_brus,'RelTol',1e-5);

      tend=7*10;        % period: 6.43, 
                      % gives pulses after the
                      % second maximum. and calculates T_new
                      % after the 5th. So should be at least 4*T

      sol_old = ode45(@brus,[0, tend],y0,options);
      
    elseif which_oscillator == "FitzHughNagumo"
      options = ddeset('Events',@findMax_fitshughnagumo,'RelTol',1e-5);

      tend = 35*10;

      sol_old = ode45(@fitzhughnagumo,[0, tend],y0,options);
      
    elseif which_oscillator == "StuartLandau"
      options = ddeset('Events',@findMax_stuartlandau,'RelTol',1e-5);

      tend = 7*10;

      sol_old = ode45(@stuartlandau,[0, tend],y0,options);

    else
      error("Need to pick oscillator type. Define which_oscillator as Brusselator or FitzHughNagumo or StuartLandau")
    end

    t_old=sol_old.x; 
    y_old=sol_old.y; 
    te_old=sol_old.xe; 
    ye_old=sol_old.ye; 
    ie_old=sol_old.ie;

    figure
    plot(t_old,y_old)
    
      
    n1_old = find(ie_old == 1); % indices for maxima of z_1
    n2_old = find(ie_old == 2); % indices for maxima of z_2
 
    period_est = te_old(n1_old(end)) - te_old(n1_old(end-1));

    t_4_old = te_old(n1_old(end)); % time of last maxima of z_1
    t_2_old = te_old(n1_old(end-1)); % time of second last maxima of z_1
    
    hold on;
    plot(te_old(n1_old), ye_old(1,n1_old),'kx'  )
    plot(te_old(n2_old), ye_old(2,n2_old),'ko'  )


    T_old = t_4_old - t_2_old;

    % Times of phases at which to calculate PRC
    %number_of_phases = 500;
    %stepsize = T_old / number_of_phases;
   
    number_of_phases = 1000;
    stepsize = T_old/number_of_phases;

    % stepsize =  0.1;
    % number_of_phases = floor(T_old/stepsize);
     
    
    
    % Peak that is used to compare the difference between original and
    % perturbed trajectory. 
    % *NOTE: The perturbed trajectory runs only from
    % some time between the 2nd and 3rd peak of the original data. So the
    % indexing of the maxima is shifted by "+2"
    peak_calc = 6; % look at shift of ..th peak found (e.g. the first one after the 
                    % interval at which we perterb

    t_old_phis = te_old(n1_old(peak_calc-1)):stepsize:te_old(n1_old(peak_calc)); % pulses between 2nd and 3rd peak
    y_old_phis = deval(sol_old,t_old_phis);
    
    PRC_t = zeros(numel(t_old_phis),1); % initialize PRC_t

    % Run over set of phases at which pulses are added
    for idx = 1:number_of_phases
        % Solve starting at some phase, taking the 
        % t_old at that phase as a starting point
        if which_oscillator == "Brusselator"
          sol = ode45(@brus,[t_old_phis(idx),tend],y_old_phis(:,idx)+pulse,options);
        elseif which_oscillator == "FitzHughNagumo"
          sol = ode45(@fitzhughnagumo,[t_old_phis(idx),tend],y_old_phis(:,idx)+pulse,options);
        elseif which_oscillator == "StuartLandau"
          sol = ode45(@stuartlandau,[t_old_phis(idx),tend],y_old_phis(:,idx)+pulse,options);
        end



        t=sol.x; y=sol.y; te=sol.xe; ye=sol.ye; ie=sol.ie;

        % Also calculate maxima to find T_new
        n1 = find(ie == 1); % indices for maxima of z_1 (=x1)

        tpeak_new = te(n1(peak_calc));
        tpeak_old = te_old(n1_old(peak_calc+2)); % see *NOTE above.

        % Determining PRC by calculating difference between a peak of new
        % trajectory and old trajectory
        %PRC_t(idx) = 2*pi*(tpeak_old -  tpeak_new) / T_old;
        PRC_t(idx) = mod( T_old / 2 + (tpeak_old -  tpeak_new) , T_old) - T_old / 2;

    end
    
    %mod(T/2+peak0-t(peaks(end)),T)-T/2
    
    % plotting 1 waveform
    t_waveform = t_old_phis; % linspace(t_2_num,t_4_num,100);
    waveform = deval(sol_old,t_waveform);
    
    % converting from t to phase by substracting t_0
    t_old_phis = (t_old_phis - t_old_phis(1) + stepsize);
    
    t_waveform = (t_waveform - t_waveform(1));

    %pulsesize = max(pulse); % assuming applying pulse to either x or y
    %PRC_t = PRC_t' * pulsesize;  % scaling by perturbation size

end


function [value,isterminal,direction] = findMax_brus(x,y)
    % finding maxima of solution of ode / dde
    value = brus(x,y);
    isterminal = zeros(2,1);
    direction = -ones(2,1);
end

function [value,isterminal,direction] =findMax_fitshughnagumo(x,y)
    % finding maxima of solution of ode / dde
    value = fitzhughnagumo(x,y);
    isterminal = zeros(2,1);
    direction = -ones(2,1);
end

function [value,isterminal,direction] =findMax_stuartlandau(x,y)
    % finding maxima of solution of ode / dde
    value = stuartlandau(x,y);
    isterminal = zeros(2,1);
    direction = -ones(2,1);
end


