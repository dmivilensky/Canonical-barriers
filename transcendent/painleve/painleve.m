function [v_prime_prime] = painleve(t, v, v_prime)
    %{
    Description:
      Painleve III D7 equation.
    %}
    v_prime_prime = v_prime^2 / v - v_prime / t + v^2 / t - 1 / v;
end
