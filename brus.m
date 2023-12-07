function dzdt = brus(t,z)
    % One oscillator Brusselator model as in [Rusin, 2010], but without 
    % feedback (and the error fixed)
    % uses function f (see below)
    B=2.3; A=1.0;

    dzdt= [(B-1)*z(1)+A^2*z(2)+f(z(1),z(2),A,B);...  % x
           -B*z(1)-A^2*z(2)-f(z(1),z(2),A,B) ];      % y
end

function [out]=f(x,y,A,B)
    % function used by @brus(t,z)
    out = B/A*x^2 + 2*A*x*y + x^2*y;
end