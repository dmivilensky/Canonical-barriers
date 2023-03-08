function [is_in_cone, x] = canonize(x_)
    %{
    Description:
      Find the symmetry M * T * S \in GL(3, R), mapping the point x_ \in R^3 
      from the cone into the fundamental region (x = (1, x_1, x_2), 
      where 0 <= x_1 <= 1/2, x_1 < x_2), and the image x 
      of the corresponding mapping.

    Arguments:
      x_ --- point from the cone in R^3 stretched on integer points of the
      form (1, n, n^2).

    Returns:
      is_in_cone --- true if x_ lies in the given cone,
      x --- image of x_ lying in the fundamental region.

    Usage example:
      x = [2; 5; 13];
      [is_in_cone, x] = canonize(x);
    %}

    A = 2 * floor(x_(2) / x_(1)) + 1;
    B = -1;
    C = -floor(x_(2) / x_(1)) * (floor(x_(2) / x_(1)) + 1);

    is_in_cone = (x_(1) > 0) & ((A * (x_(2) / x_(1)) + B * (x_(3) / x_(1)) + C) / B >= 0);

    if is_in_cone
        S = eye(3) / x_(1);
        n = -round(x_(2) / x_(1));
        T = [1 0 0; n 1 0; n^2 2*n 1];
        M = [1 0 0; 0 sign(x_(2) / x_(1) + n) 0; 0 0 1];

        x = (M * T * S) * x_;
    else
        x = nan;
    end
end