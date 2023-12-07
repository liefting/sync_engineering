function [dzdt] = oscillators(t,z,L_vec,which_oscillator,K,k_vec,a0,M)
 
    L_vec = L_vec - a0;

    xlag100 = L_vec(1:2:2*M,1); xlag010 = L_vec(1:2:2*M,2); xlag001 = L_vec(1:2:2*M,3); 
    xlag200 = L_vec(1:2:2*M,4); xlag020 = L_vec(1:2:2*M,5); xlag002 = L_vec(1:2:2*M,6);
    xlag300 = L_vec(1:2:2*M,7); xlag030 = L_vec(1:2:2*M,8); xlag003 = L_vec(1:2:2*M,9);
    xlag120 = L_vec(1:2:2*M,10); xlag210 = L_vec(1:2:2*M,11); xlag012 = L_vec(1:2:2*M,12); 
    xlag021 = L_vec(1:2:2*M,13); xlag102 = L_vec(1:2:2*M,14); xlag201 = L_vec(1:2:2*M,15);
    xlag111 = L_vec(1:2:2*M,16);
    xlag110 = L_vec(1:2:2*M,17); 
    xlag011 = L_vec(1:2:2*M,18); 
    xlag101 = L_vec(1:2:2*M,19);
    k100 = k_vec(1); k010 = k_vec(2); k001 = k_vec(3); k200 = k_vec(4); k020 = k_vec(5); k002 = k_vec(6);
    k300 = k_vec(7); k030 = k_vec(8); k003 = k_vec(9);
    k120 = k_vec(10); k210 = k_vec(11); k012 = k_vec(12); k021 = k_vec(13); k102 = k_vec(14); k201 = k_vec(15);
    k111 = k_vec(16); k110 = k_vec(17); k011 = k_vec(18); k101 = k_vec(19);

    % g,j as in hxgj should match the index of the xlagpq's. 
    % p,q as in xlagpq should match the powers in the products
    feedback_term = 0;
    for i = 1:M
        for j = 1:M
            for k = 1:M
                feedback_term = feedback_term + k100*xlag100(i) + k010*xlag010(j) + k001*xlag001(k) +...
                    k200*xlag200(i)^2 + k020*xlag020(j)^2 + k002*xlag002(k)^2 + ...
                    k110*xlag110(i)*xlag110(j) + k011*xlag011(j)*xlag011(k) + k101*xlag101(i)*xlag101(k) + ...
                    k111*xlag111(i)*xlag111(j)*xlag111(k) + ...
                    k300*xlag300(i)^3 + k030*xlag030(j)^3 + k003*xlag003(k)^3 + ...
                    k021*xlag021(j)^2*xlag021(k) + k012*xlag012(j)*xlag012(k)^2 + ...
                    k201*xlag201(i)^2*xlag201(k) + k102*xlag102(i)*xlag102(k)^2 + ...
                    k210*xlag210(i)^2*xlag210(j) + k120*xlag120(i)*xlag120(j)^2;
            end
        end
    end
    
    feedback_term = K/(M*M*M)*feedback_term; 
    dzdt = zeros(2*M,1);

    if which_oscillator == "Brusselator"
      for i = 1:M
         x_i = z(2*i-1);
         y_i = z(2*i);
         dzdt(2*i-1:2*i,:) =  brus(t,[x_i;y_i]) + [feedback_term ; 0 ];
      end
    elseif which_oscillator == "FitzHughNagumo"
      for i = 1:M
         x_i = z(2*i-1);
         y_i = z(2*i);
         dzdt(2*i-1:2*i,:) = fitzhughnagumo(t,[x_i;y_i]) + [feedback_term ; 0 ];
      end
    else
      error("what oscillator do we simulate??")
    end
end