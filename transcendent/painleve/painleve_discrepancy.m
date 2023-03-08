function [error] = painleve_discrepancy(t, v, v_prime, v_prime_prime)
    %{
    Description:
      Error in given Painleve III D7 transcendent approximation.
    %}
    error = abs(-t .* v .* v_prime_prime + t .* v_prime.^2 - v .* v_prime + v.^3 - t);
end
