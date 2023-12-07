function [R_ordersols,timevec] = find_order_parameter(Hint,solutions,M,deg,tspan,nopoints)
     % evaluate orderparameter at nopoints timepoints
    if nargin == 4
      timevec = linspace(0+eps,solutions(1).x(end)-eps,nopoints);
    else  
      timevec = linspace(tspan(1)+eps,tspan(2)-eps,nopoints);
    end
    R_ordersols = zeros(length(solutions),nopoints);
    for soli = 1:length(solutions)
        z = deval(solutions(soli),timevec);
        [zw, zh] = size(z);
        phissoli = zeros(round(zw/2),zh);
          %phis2 = zeros(round(zw/2),zh);
        for i = 1:round(zw/2)
            phissoli(i,:) = Hint(z(2*i-1,:),z(2*i,:));
        end
        R_ordersols(soli,:) = abs( 1/M * sum(exp(1i*deg*phissoli),1));
    end
    R_ordersols = R_ordersols';
    timevec = timevec';
end