function [error] = painleve_log_substitution_discrepancy(l, w, w_prime, w_prime_prime)
    %{
    Description:
      Error in given Painleve III D7 transcendent approximation.
    %}
    error = abs(w_prime.^2 + exp(l) .* w.^3 - exp(2*l) - w .* w_prime_prime);
end
