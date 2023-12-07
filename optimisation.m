function [store_x_nonz, h_functions] = optimisation(...
          al,zl,omega,H,ntaus,npark,period_est,pairwise) 
    % Function to find the optimum values for taus, kappas for given 
    % fourier coefficients for the waveform and PRC
    n_par = ntaus + npark; 
    N=100; % number of initial conditions to try
    store_x = zeros(N, n_par);
    store_x0 = zeros(N, n_par);
        
    lb = zeros(1,n_par); ub =ones(1,n_par)*100;   
    h_functions = get_h_functions(al, omega);
    options = optimoptions(@fmincon,'Display','off','Algorithm','sqp');
    f_cost =  @(x) cost(x,ntaus);

    tic
    for n = 1:N
        rng(n)
        x0 = rand(1,n_par).*period_est;
        
        % running the optimisation starting at x0:
        [x,~,exitflag,~] = fmincon(f_cost,x0,[],[],[],[], lb, ub ,@(x) constraints(x, zl,H,period_est,h_functions,ntaus,pairwise), options);
        
        % if there is convergence to a feasible point, store the solution:
        if exitflag > 0
            store_x(n,:) = x;
        end
        store_x0(n,:) = x0;
    end
    toc
    
    % cleaning up array, removing duplicates:
    store_x_nonz = zeros(400,n_par); k = 0;
    for n=1:N
       if store_x(n,1) > 0
           k=k+1;
           store_x_nonz(k,:) = store_x(n,:);
       end
    end
    store_x_nonz = store_x_nonz(1:k,:);
    [n_sols,n_par]=size(store_x_nonz);
    
    for i=1:n_sols  % if k_pq == 0, then set tau_pq to 2
        for j = ntaus+1 : n_par   % running over k_pq
            if store_x_nonz(i,j) < 1e-8 % if k_pq (almost) zero
                store_x_nonz(i,j - ntaus) = 0; % set corresponding t_pq to 0
            end
        end     
    end
    
    for i=1:n_sols
       check = 1; % assume parameter set not found previously
       for j=1:k+1
          if  abs(store_x_nonz(i,:) - store_x_nonz(j,:)) < ones(1,n_par)*0.001
              check = 0; % (approximate) parameter set already found
              break;
          end
       end
       if check
            k=k+1;
            store_x_nonz(k,:)= store_x_nonz(i,:);
       end
    end
    store_x_nonz = store_x_nonz(1:k,:);
end

%% Function to minimize:
function f = cost(x,ntaus)
    % minimize all gains k_{} to get weak feedback
    f1 = sum(x(ntaus+1:end)); 
    f = f1 ; %+ f2;
end

function [c, ceq] = constraints(x, zl,H,period_est,h_functions,ntaus,pairwise)
    fb = vector_to_taus_and_ks(x);
   
    z_1=conj(zl(1)); z_2=conj(zl(2)); z_3=conj(zl(3)); 

    H100 = H(1);  % h10
    H110 = H(2);  % h11
    H200 = H(3);  % h20
    H010 = 0;   % h01
    H020 = 0;   % h02
    H300 = 0;   % h30
    H030 = 0;   % h03
    H210 = 0;   % h21
    H120 = 0;   % h12

    epsilon = 1e-2;
    
    c_pair = [  -x(1:ntaus)';  % all lags must be nonnegative  "c <= 0" --> "-c >= 0" want: c>0 
                x(1:ntaus)' - period_est;   % all lags must me less than or equal to one period    
                
    %% ----------------------  pairwise constraints  ---------------------- %%
                abs( h_functions.h100( fb )*z_1 - H100 )- epsilon; %h100   
                abs( h_functions.h200( fb )*z_2 - H200 )- epsilon; %h200      
    
              ];
    
    
    %% ----------------------  additional nonpairwise constraints  ---------------------- %%
    if ~pairwise   
      c_add_nonp = [  abs( h_functions.h300( fb )*z_3 - H300 ) - epsilon;  %h300
                      abs( h_functions.h110( fb )*z_2 - H110 ) - epsilon; %h110  
                      abs( h_functions.h010( fb )*z_1 - H010 ) - epsilon;  %h010 
                      abs( h_functions.h020( fb )*z_2 - H020 ) - epsilon;  %h020 
                      abs( h_functions.h030( fb )*z_3 - H030 ) - epsilon;  %h030   
                      abs( h_functions.h210( fb )*z_3 - H210 ) - epsilon;  %h210
                      abs( h_functions.h120( fb )*z_3 - H120 ) - epsilon;  %h120
                    ];
      c = [c_pair ; c_add_nonp];
    else
      c = c_pair;
    end
            
    % corresponding to fourway terms:
    ceq_n = [ fb.k001 ;
              fb.k002 ;
              fb.k003 ;
              fb.k111 ;
              fb.k011 ;
              fb.k101 ;
              fb.k012 ;
              fb.k021 ;
              fb.k102 ;
              fb.k201 ;
              fb.tau001 ;
              fb.tau002 ;
              fb.tau003 ;
              fb.tau111 ;
              fb.tau011 ;
              fb.tau101 ;
              fb.tau012 ;
              fb.tau021 ;
              fb.tau102 ;
              fb.tau201 ];

    ceq_add_pair = [  fb.k120; 
                      fb.k210;
                      fb.k110;
                      fb.k300;
                      fb.tau120;
                      fb.tau210;
                      fb.tau110;
                      fb.tau300 ];

    if pairwise
      ceq = [ceq_n; ceq_add_pair]; % if pairwise, then require nonpairwise ones to be zero
    else
      ceq = ceq_n;  % if nonpairwise, allow for nonpairwise feedback parameters to be nonzero
    end     
end
