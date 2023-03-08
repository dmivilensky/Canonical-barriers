function [v] = asymptotic_near_infty(t)
    %{
    Description:
      Asymptotics for t -> +\infty of the solution of Painleve III D7
      equation satysfying the condition that \int_{-\infty}^0 \exp(\chi/2) dx diverges, 
      in a form it appears in our problem.
    %}

    v = (t.^(1/3)) .* (1 + (3^(3/4))*(pi^(-1/2)).*(t.^(-1/3)).*exp(-(3/2)*sqrt(3).*t.^(2/3)));
end