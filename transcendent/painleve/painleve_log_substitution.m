function [w_prime_prime] = painleve_log_substitution(l, w, w_prime)
    %{
    Description:
      Painleve III D7 equation after substitution l = log t. In this case,
      dv/dt = e^-l * dw/dl, d^2/dt^2(v) = e^-2l * (d^2/dl^2(w) - dw/dl).
    %}
    w_prime_prime = (w_prime^2) / w - w^2 * exp(l) - exp(2*l) / w;
end