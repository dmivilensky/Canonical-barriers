function [v, v_prime, v_prime_prime] = pade_evaluate(a, b, h)
    %{
    Description:
      Evaluate Pade approximant at point h

    Arguments:
      a --- coefficients of the numerator polynomial
      b --- coefficients of the denominator polynomial

    Returns:
      v, v_prime, v_prime_prime --- values of approximant, its first and second
      derivatives
    %}
    
    numerator = 0;
    numerator_prime = 0;
    numerator_prime_prime = 0;
    for i=1:size(a,2)
        numerator = numerator + a(i) * h^(i-1);
    end
    for i=2:size(a,2)
        numerator_prime = numerator_prime + (i - 1) * a(i) * h^(i-2);
    end
    for i=3:size(a,2)
        numerator_prime_prime = numerator_prime_prime + (i - 2) * (i - 1) * a(i) * h^(i-3);
    end
    denominator = 0;
    denominator_prime = 0;
    denominator_prime_prime = 0;
    for i=1:size(a,2)
        denominator = denominator + b(i) * h^(i-1);
    end
    for i=2:size(a,2)
        denominator_prime = denominator_prime + (i - 1) * b(i) * h^(i-2);
    end
    for i=3:size(a,2)
        denominator_prime_prime = denominator_prime_prime + (i - 2) * (i - 1) * b(i) * h^(i-3);
    end
    v = numerator / denominator;
    v_prime = (denominator * numerator_prime - numerator * denominator_prime) / denominator^2;
    v_prime_prime = (denominator_prime * numerator_prime + denominator * numerator_prime_prime - numerator_prime * denominator_prime - numerator * denominator_prime_prime) / denominator^4;
end