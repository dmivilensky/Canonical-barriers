function [v, v_prime, v_prime_prime] = asymptotic_near_infty_refined_with_derivatives(t)
    %{
    Description:
      Asymptotics for t -> +\infty of the solution of Painleve III D7
      equation satysfying the condition that \int_{-\infty}^0 \exp(\chi/2) dx diverges, 
      in a form it appears in our problem.
    %}

    s = 3;

    A = log(sqrt(32));
    B = (2^(2/3) * 3^(3/2));
    C = 2^(1/3);
    D = (2*pi)^(-1/2) * 3^(-1/4) * s * sqrt(2/pi);

    x = 2 .* (log(t) - A);
    tau = B .* exp(x./3);
    
    u = C .* exp(2.*x/3) + D .* sqrt(tau) .* exp(tau) .* besselk(0,tau) .* exp(x./2) .* exp(-B.*exp(x./3));
    v = 8 * u ./ t; 

    v_prime = ((32*C.*exp((4.*log(t))./3 - (4*A)/3))./(3.*t) + (8*D.*exp(log(t) - A).*(B.*exp((2.*log(t))./3 - (2*A)/3)).^(1/2).*besselk(0, B.*exp((2.*log(t))./3 - (2*A)/3)))./t + (8*B*D.*exp(log(t) - A).*exp((2.*log(t))./3 - (2*A)/3).*besselk(0, B.*exp((2.*log(t))./3 - (2*A)/3)))./(3.*t.*(B.*exp((2.*log(t))./3 - (2*A)/3)).^(1/2)) - (16*B*D.*exp(log(t) - A).*exp((2.*log(t))./3 - (2*A)/3).*(B.*exp((2.*log(t))./3 - (2*A)/3)).^(1/2).*besselk(1, B.*exp((2.*log(t))./3 - (2*A)/3)))./(3.*t))./t - (8*C.*exp((4.*log(t))./3 - (4*A)/3) + 8*D.*exp(log(t) - A).*(B.*exp((2.*log(t))./3 - (2*A)/3)).^(1/2).*besselk(0, B.*exp((2.*log(t))./3 - (2*A)/3)))./t.^2;
    
    v_prime_prime = (2*(8*C.*exp((4.*log(t))./3 - (4*A)/3) + 8*D.*exp(log(t) - A).*(B.*exp((2.*log(t))./3 - (2*A)/3)).^(1/2).*besselk(0, B.*exp((2.*log(t))./3 - (2*A)/3))))./t.^3 + ((32*C.*exp((4.*log(t))./3 - (4*A)/3))./(9.*t.^2) + (16*B*D.*exp(log(t) - A).*exp((2.*log(t))./3 - (2*A)/3).*(B.*exp((2.*log(t))./3 - (2*A)/3)).^(1/2).*((2.*besselk(1, B.*exp((2.*log(t))./3 - (2*A)/3)))./(3.*t) + (2*B.*exp((2.*log(t))./3 - (2*A)/3).*besselk(0, B.*exp((2.*log(t))./3 - (2*A)/3)))./(3.*t)))./(3.*t) - (32*B^2*D.*exp(log(t) - A).*exp((4.*log(t))./3 - (4*A)/3).*besselk(1, B.*exp((2.*log(t))./3 - (2*A)/3)))./(9.*t.^2.*(B.*exp((2.*log(t))./3 - (2*A)/3)).^(1/2)) - (8*B^2*D.*exp(log(t) - A).*exp((4.*log(t))./3 - (4*A)/3).*besselk(0, B.*exp((2.*log(t))./3 - (2*A)/3)))./(9.*t.^2.*(B.*exp((2.*log(t))./3 - (2*A)/3)).^(3/2)) + (40*B*D.*exp(log(t) - A).*exp((2.*log(t))./3 - (2*A)/3).*besselk(0, B.*exp((2.*log(t))./3 - (2*A)/3)))./(9.*t.^2.*(B.*exp((2.*log(t))./3 - (2*A)/3)).^(1/2)) - (80*B*D.*exp(log(t) - A).*exp((2.*log(t))./3 - (2*A)/3).*(B.*exp((2.*log(t))./3 - (2*A)/3)).^(1/2).*besselk(1, B.*exp((2.*log(t))/3 - (2*A)/3)))./(9.*t.^2))./t - (2.*((32*C.*exp((4.*log(t))./3 - (4*A)/3))./(3.*t) + (8*D.*exp(log(t) - A).*(B.*exp((2.*log(t))./3 - (2*A)/3)).^(1/2).*besselk(0, B.*exp((2.*log(t))./3 - (2*A)/3)))./t + (8*B*D.*exp(log(t) - A).*exp((2.*log(t))./3 - (2*A)/3).*besselk(0, B.*exp((2.*log(t))./3 - (2*A)/3)))./(3.*t.*(B.*exp((2.*log(t))./3 - (2*A)/3)).^(1/2)) - (16*B*D.*exp(log(t) - A).*exp((2.*log(t))./3 - (2*A)/3).*(B.*exp((2.*log(t))./3 - (2*A)/3)).^(1/2).*besselk(1, B.*exp((2.*log(t))./3 - (2*A)/3)))./(3.*t)))./t.^2;
end
