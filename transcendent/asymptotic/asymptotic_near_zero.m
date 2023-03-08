function [v] = asymptotic_near_zero(t)
    %{
    Description:
      Asymptotics for t -> 0+ of the solution of Painleve III D7 equation satysfying
      the condition that \int_{-\infty}^0 \exp(\chi/2) dx diverges, in a form it
      appears in our problem. 
    %}

    v = (2 ./ t) .* (log(t) - 3*log(2) + 3*eulergamma/2).^(-2);
end