function dzdt = fitzhughnagumo(t,z)
    % One oscillator model as on http://www.scholarpedia.org/article/FitzHugh-Nagumo_model
    v = z(1);
    w = z(2);
  
   
    A = 1;
    I_app = 0.6;
    epsilon = 0.08;
    gamma = 1;

    dvdt = -w - v*( v - 1 ) * ( v - A ) + I_app;
    dwdt = epsilon*( v - gamma*w);
    
    dzdt = [dvdt; dwdt];

end