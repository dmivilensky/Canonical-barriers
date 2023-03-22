function [v] = asymptotic_near_infty_refined(t)
    %{
    Description:
      Asymptotics for t -> +\infty of the solution of Painleve III D7
      equation satysfying the condition that \int_{-\infty}^0 \exp(\chi/2) dx diverges, 
      in a form it appears in our problem.
    %}

    s = 3;
    x = 2 .* log(t ./ sqrt(32));
    tau = (2^(2/3) * 3^(3/2)) .* exp(x/3);

    v = 2^(1/3) .* exp(2.*x./3) + (2*pi)^(-1/2) * 3^(-1/4) * s * sqrt(2/pi) .* tau .* exp(tau) .* besselk(0,tau) .* exp(x/2) .* exp(-2^(2/3)*3^(3/2).*exp(x/3));
end