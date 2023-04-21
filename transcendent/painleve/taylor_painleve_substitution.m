function [coefs] = taylor_painleve_substitution(s_0, u_0, u_1, max_degree, factorials)
    %{
    Description:
      Coefficients in a Taylor series of the solution to Painleve III D7
      equation, after substitution s = 1/t, at time s_0 with initial point 
      u(s_0) = u_0, du/ds(s_0) = u_1 (= -t0^2 * v'(t_0)). 
      One can test a procedure with s_0 = 1, u_0 = 1, u_1 = -1/3, the result
      will correspond to Taylor expansion of u = s^(-1/3) at s = 1.

    Arguments:
      s_0 --- point at which Taylor expansion is centered,
      u_0, u_1 --- initial point,
      max_degree --- number of coefficients to calculate,
      factorials --- true if coefficients should be multiplied by k!
      (when we want to calculate just derivatives)

    Returns:
      coefs --- coefficients calculated from (x - s_0)^0 to (x - s_0)^max_degree
    %}

    coefs = zeros(1, max_degree+1);
    coefs(0+1) = u_0; coefs(1+1) = u_1;
    coefs(2+1) = (s_0^4 * u_1^2 - s_0^3 * u_0 * u_1 + s_0 * u_0^3 - 1) / (2 * s_0^4 * u_0);
    coefs(3+1) = (-2 * s_0^4 * u_1 * coefs(2+1) - 10 * s_0^3 * u_0 * coefs(2+1) + 4 * s_0^4 * u_1 * coefs(2+1) + 3 * s_0^3 * u_1^2 - 3 * s_0^2 * u_0 * u_1 + 3 * s_0 * u_0^2 * u_1 + u_0^3) / (6 * s_0^4 * u_0);
    for d=4:max_degree
        q = d - 2;
        sum_6 = 0; sum_14 = 0;
        for k=0:q
            l = q - k;
            sum_6 = sum_6 + (k+1) * (l+1) * coefs(k + 1 + 1) * coefs(l + 1 + 1);
            sum_14 = sum_14 + (k+1) * coefs(l + 1) * coefs(k + 1 + 1);
        end
        sum_1 = 0; sum_2 = 0; sum_7 = 0; sum_13 = 0;
        for k=0:(q-1)
            l = q - 1 - k;
            sum_1 = sum_1 + (k+1) * (k+2) * coefs(k + 2 + 1) * coefs(l + 1 + 1);
            sum_2 = sum_2 + (k+1) * (k+2) * coefs(k + 2 + 1) * coefs(l + 1);
            sum_7 = sum_7 + (k+1) * (l+1) * coefs(k + 1 + 1) * coefs(l + 1 + 1);
            sum_13 = sum_13 + (k+1) * coefs(k + 1 + 1) * coefs(l + 1);
        end
        sum_3 = 0; sum_8 = 0; sum_12 = 0;
        for k=0:(q-2)
            l = q - 2 - k;
            sum_3 = sum_3 + (k+1) * (k+2) * coefs(k + 2 + 1) * coefs(l + 1);
            sum_8 = sum_8 + (k+1) * (l+1) * coefs(k + 1 + 1) * coefs(l + 1 + 1);
            sum_12 = sum_12 + (k+1) * coefs(k + 1 + 1) * coefs(l + 1);
        end
        sum_4 = 0; sum_9 = 0; sum_11 = 0;
        for k=0:(q-3)
            l = q - 3 - k;
            sum_4 = sum_4 + (k+1) * (k+2) * coefs(k + 2 + 1) * coefs(l + 1);
            sum_9 = sum_9 + (k+1) * (l+1) * coefs(k + 1 + 1) * coefs(l + 1 + 1);
            sum_11 = sum_11 + (k+1) * coefs(k + 1 + 1) * coefs(l + 1); 
        end
        sum_5 = 0; sum_10 = 0;
        for k=0:(q-4)
            l = q - 4 - k;
            sum_5 = sum_5 + (k+1) * (k+2) * coefs(k + 2 + 1) * coefs(l + 1);
            sum_10 = sum_10 + (k+1) * (l+1) * coefs(k + 1 + 1) * coefs(l + 1 + 1);
        end
        sum_16 = 0;
        for k=0:(q-1)
            for l=0:(q-1-k)
                m = q - 1 - k - l;
                sum_16 = sum_16 + coefs(k + 1) * coefs(l + 1) * coefs(m + 1);
            end
        end
        sum_15 = 0;
        for k=0:q
            for l=0:(q-k)
                m = q - k - l;
                sum_15 = sum_15 + coefs(k + 1) * coefs(l + 1) * coefs(m + 1);
            end
        end
        coefs(d+1) = ( ...
            -s_0^4 * sum_1 - 4 * s_0^3 * sum_2 - 6 * s_0^2 * sum_3 - ...
            4 * s_0 * sum_4 - sum_5 + s_0^4 * sum_6 + 4 * s_0^3 * sum_7 + ...
            6 * s_0^2 * sum_8 + 4 * s_0 * sum_9 + sum_10 - sum_11 - ...
            3 * s_0 * sum_12 - 3 * s_0^2 * sum_13 - s_0^3 * sum_14 + ...
            s_0 * sum_15 + sum_16) / ((q+1) * (q+2) * s_0^4 * u_0);
    end
    if factorials
        factorial = 1;
        for i=1:max_degree+1
            coefs(i) = coefs(i) * factorial;
            factorial = factorial * i;
        end
    end
end
