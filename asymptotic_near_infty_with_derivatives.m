function [v, v_prime, v_prime_prime] = asymptotic_near_infty_with_derivatives(t)
    %{
    Description:
      Asymptotics for t -> +\infty of the solution of Painleve III D7
      equation satysfying the condition that \int_{-\infty}^0 \exp(\chi/2) dx diverges, 
      in a form it appears in our problem.
    %}

    c = (3^(3/4))*(pi^(-1/2));
    v = (t.^(1/3)) .* (1 + c.*(t.^(-1/3)).*exp(-(3/2)*sqrt(3).*t.^(2/3)));
    v_prime = ((c.*exp(-(3*3^(1/2).*t.^(2/3))./2))./t.^(1/3) + 1)./(3.*t.^(2/3)) - t.^(1/3).*((c.*exp(-(3*3^(1/2).*t.^(2/3))./2))/(3.*t.^(4/3)) + (3^(1/2)*c.*exp(-(3*3^(1/2).*t.^(2/3))./2))./t.^(2/3));
    v_prime_prime = t.^(1/3).*((3*c.*exp(-(3*3^(1/2).*t.^(2/3))./2))./t + (4*c.*exp(-(3*3^(1/2).*t.^(2/3))./2))./(9.*t.^(7/3)) + (3^(1/2)*c.*exp(-(3*3^(1/2).*t.^(2/3))./2))./t.^(5/3)) - (2*((c.*exp(-(3*3^(1/2).*t.^(2/3))./2))./t.^(1/3) + 1))./(9.*t.^(5/3)) - (2*((c.*exp(-(3*3^(1/2).*t.^(2/3))./2))./(3.*t.^(4/3)) + (3^(1/2)*c.*exp(-(3*3^(1/2).*t.^(2/3))./2))./t.^(2/3)))./(3.*t.^(2/3));
end
