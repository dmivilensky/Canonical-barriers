function [root] = brent(f, a, b, error)
    %{
    Description:
      Brent root-finding method.

    Arguments:
      f --- function which root the method is finding,
      a, b --- left and right bounds of segment to find root on, f(a) and
      f(b) must be of the different signs,
      error --- positive threshold which bounds the absolute value of
      function at point we consider to be the approximate root.

    Returns:
      root --- point from segment, absolute function at which is less than given
      threshold, or center of the segment of candidate points if its length
      is less than computer zero (10^(-16) in our case).
    %}

    x1 = a;
    f1 = f(x1);
    if abs(f1) < error
        root = x1;
        return;
    end
    x2 = b;
    f2 = f(x2);
    if abs(f2) < error
        root = x2;
        return;
    end
    assert(f1 * f2 < 0);
    x3 = (a + b) / 2;

    while true
        f3 = f(x3);
        % clc;
        % fprintf('%0.2e, %0.16f\n', b-a, f3);
        if abs(f3) < error
            root = x3;
            return;
        end
        if f1 * f3 < 0
            b = x3;
        else
            a = x3;
        end
        if b - a < 1e-16
            root = (a + b) / 2;
            return;
        end
        
        denominator = (f2 - f1) * (f3 - f1) * (f2 - f3);
        if denominator == 0
            dx = b - a;
        else
            numerator = x3 * (f1 - f2) * (f2 - f3 + f1) + f2 * x1 * (f2 - f3) + f1 * x2 * (f3 - f1);
            dx = f3 * numerator / denominator;
        end
        x = x3 + dx;
        
        if (b - x) * (x - a) < 0
            dx = (b - a) / 2;
            x = a + dx;
        end
        if x < x3
            x2 = x3;
            f2 = f3;
        else 
            x1 = x3;
            f1 = f3;
        end
        x3 = x;
    end
end