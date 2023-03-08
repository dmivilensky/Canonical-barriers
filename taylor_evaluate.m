function [v, v_prime, v_prime_prime] = taylor_evaluate(coefs, h)
    %{
    Description:
      Evaluate Taylor's series at point h

    Arguments:
      coefs --- coefficients of the series

    Returns:
      v, v_prime, v_prime_prime --- values of series, its first and second
      derivatives
    %}
    
    v = 0;
    for i=1:size(coefs,2)
        v = v + coefs(i) * h^(i-1);
    end
    v_prime = 0;
    for i=2:size(coefs,2)
        v_prime = v_prime + (i - 1) * coefs(i) * h^(i-2);
    end
    v_prime_prime = 0;
    for i=3:size(coefs,2)
        v_prime_prime = v_prime_prime + (i - 2) * (i - 1) * coefs(i) * h^(i-3);
    end
end