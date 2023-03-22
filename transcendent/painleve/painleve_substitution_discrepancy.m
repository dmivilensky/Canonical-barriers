function [error] = painleve_substitution_discrepancy(s, u, u_prime, u_prime_prime)
    %{
    Description:
      Error in given Painleve III D7 transcendent approximation.
    %}
    error = abs(s.^4 .* u_prime.^2 - s.^3 .* u .* u_prime + s .* u.^3 - 1 - s.^4 .* u .* u_prime_prime);
end
