function [coefs] = taylor_painleve_log_substitution(l_0, w_0, w_1, max_degree, factorials, truncate_exponent)
    %{
    Description:
      Coefficients in a Taylor series of the solution to Painleve III D7
      equation, after substitution l = log t, at time l_0 with initial point 
      w(l_0) = w_0, dw/dl(l_0) = w_1 (= t * v'(t_0)). 
      One can test a procedure with l_0 = 0, w_0 = 1, w_1 = 1/3, the result
      will correspond to Taylor expansion of w = e^(l/3) at l = 0.

    Arguments:
      l_0 --- point at which Taylor expansion is centered,
      w_0, w_1 --- initial point,
      max_degree --- number of coefficients to calculate,
      factorials --- true if coefficients should be multiplied by k!
      (when we want to calculate just derivatives)
      truncate_exponent --- if greater than zero, sum of w_k w_l w_m divided
      by n! will be calculated for only n <= truncate_exponent (assuming that
      residual is neglectible)

    Returns:
      coefs --- coefficients calculated from (x - l_0)^0 to (x - l_0)^max_degree
    %}

    coefs = zeros(1, max_degree+1);
    coefs(0+1) = w_0; coefs(1+1) = w_1;
    coefs(2+1) = (w_1^2 + exp(l_0) * w_0^3 - exp(2*l_0)) / (2 * w_0);
    power_of_2_div_factorial = 1;
    for d=1:max_degree
        power_of_2_div_factorial = power_of_2_div_factorial * 2 / d;

        sum_0 = 0;
        for k=0:d
            l = d - k;
            sum_0 = sum_0 + (k + 1) * (l + 1) * coefs(k + 1 + 1) * coefs(l + 1 + 1);
        end

        sum_2 = 0;
        for k=0:d-1
            l = d - 1 - k;
            sum_2 = sum_2 + (l + 1) * (l + 2) * coefs(k + 1 + 1) * coefs(l + 2 + 1);
        end

        sum_1 = 0;
        factorial_n = 1;
        if truncate_exponent <= 0
            max_n = d;
        else
            max_n = min(d, truncate_exponent);
        end
        for n=0:max_n
            subsum = 0;
            for k=0:(d-n)
                for l=0:(d-n-k)
                    m = d - n - k - l;
                    subsum = subsum + coefs(k + 1) * coefs(l + 1) * coefs(m + 1);
                end
            end
            sum_1 = sum_1 + subsum / factorial_n;
            factorial_n = factorial_n * (n + 1);
        end
        
        coefs(d + 2 + 1) = (sum_0 - sum_2 + exp(l_0) * sum_1 - power_of_2_div_factorial * exp(2*l_0)) / ((d + 1) * (d + 2) * w_0);
    end
    if factorials
        factorial = 1;
        for i=1:max_degree+1
            coefs(i) = coefs(i) * factorial;
            factorial = factorial * i;
        end
    end
end
