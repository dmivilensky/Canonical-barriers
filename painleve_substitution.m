function [u_prime_prime] = painleve_substitution(s, u, u_prime)
    %{
    Description:
      Painleve III D7 equation after substitution s = 1/t. In this case,
      dv/dt = -s^2 * du/ds, d^2/dt^2(v) = s^4 * d^2/ds^2(u) + 2s^3 * du/ds.
    %}
    u_prime_prime = (u_prime^2) / u - u_prime / s + (u^2) / s^3 - 1 / (u * s^4);
end