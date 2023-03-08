function [v, v_prime, v_prime_prime] = asymptotic_near_zero_with_derivatives(t)
    %{
    Description:
      Asymptotics for t -> 0+ of the solution of Painleve III D7 equation satysfying
      the condition that \int_{-\infty}^0 \exp(\chi/2) dx diverges, in a form it
      appears in our problem. 
    %}

    v = (2 ./ t) .* (log(t) - 3*log(2) + 3*eulergamma/2).^(-2);
    v_prime = -2 ./ ((t.^2).*((3*eulergamma)/2 + log(t) - log(8)).^2) - 4./((t.^2).*((3*eulergamma)/2 + log(t) - log(8)).^3);
    v_prime_prime = 4 ./ ((t.^3).*((3*eulergamma)/2 + log(t) - log(8)).^2) + 12./((t.^3).*((3*eulergamma)/2 + log(t) - log(8)).^3) + 12./((t.^3).*((3*eulergamma)/2 + log(t) - log(8)).^4);
end