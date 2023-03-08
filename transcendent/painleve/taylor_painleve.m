function [coefs] = taylor_painleve(t_0, u_0, u_1, max_degree, factorials)
    %{
    Description:
      Coefficients in a Taylor series of the solution to Painleve III D7
      equation at time t_0 with initial point u(t_0) = u_0, u'(t_0) = u_1.
      One can test a procedure with t_0 = 1, u_0 = 1, u_1 = 1/3, the result
      will correspond to Taylor expansion of u = t^(1/3) at t = 1.

    Arguments:
      t0 --- point at which Taylor expansion is centered,
      u_0, u_1 --- initial point,
      max_degree --- number of coefficients to calculate,
      factorials --- true if coefficients should be multiplied by k!
      (when we want to calculate just derivatives)

    Returns:
      coefs --- coefficients calculated from (x - t_0)^0 to (x - t_0)^max_degree
    %}

    coefs = zeros(1, max_degree+1);
    coefs(0+1) = u_0; coefs(1+1) = u_1;
    coefs(2+1) = (t_0 * u_1^2 - u_0 * u_1 + u_0^3 - t_0) / (2 * t_0 * u_0);
    coefs(3+1) = (-4 * u_0 * coefs(2+1) + 2 * t_0 * u_1 * coefs(2+1) + 3 * u_0^2 * u_1 - 1) / (6 * t_0 * u_0);
    for d=4:max_degree
        q = d - 2;
        sum_2 = 0;
        sum_4 = 0;
        for k = 0:q
            l = q - k;
            sum_2 = sum_2 + (k + 1) * (l + 1) * coefs(k + 1 + 1) * coefs(l + 1 + 1);
            sum_4 = sum_4 + (k + 1) * coefs(k + 1 + 1) * coefs(l + 1);
        end
        sum_1 = 0;
        sum_3 = 0;
        sum_6 = 0;
        for k = 0:q-1
            l = q - 1 - k;
            sum_1 = sum_1 + (k + 1) * (k + 2) * coefs(k + 2 + 1) * coefs(l + 1);
            sum_3 = sum_3 + (k + 1) * (l + 1) * coefs(k + 1 + 1) * coefs(l + 1 + 1);
            sum_6 = sum_6 + (k + 1) * (k + 2) * coefs(k + 2 + 1) * coefs(l + 1 + 1);
        end
        sum_5 = 0;
        for k = 0:q
            for l = 0:(q-k)
                m = q - k - l;
                sum_5 = sum_5 + coefs(k + 1) * coefs(l + 1) * coefs(m + 1);
            end
        end
        coefs(d+1) = (-sum_1 + t_0 * sum_2 + sum_3 - sum_4 + sum_5 - t_0 * sum_6) / ((q + 1) * (q + 2) * t_0 * u_0);
    end
    if factorials
        factorial = 1;
        for i=1:max_degree+1
            coefs(i) = coefs(i) * factorial;
            factorial = factorial * i;
        end
    end
end
