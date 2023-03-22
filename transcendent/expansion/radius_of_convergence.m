function [R] = radius_of_convergence(coefs)
    %{
    Description:
      Mercer--Roberts procedure to calculate radius of convergence
      of the series with coefficients `coefs`

    Arguments:
      coefs --- coefficients of the series

    Returns:
      R --- radius of convergence of the series
    %}
    
    b_sqr = abs((coefs(4:end) .* coefs(2:end-2) - coefs(3:end-1).^2) ./ (coefs(3:end-1) .* coefs(1:end-3) - coefs(2:end-2).^2));
    reciprocal_R = interp1(1 ./ (3:(size(b_sqr, 2)+2)), sqrt(b_sqr), 0, 'linear', 'extrap');
    hold on
    plot(1 ./ (3:(size(sqrt(b_sqr), 2)+2)), sqrt(b_sqr));
    plot(0, reciprocal_R, '*r')
    hold off
    R = 1 / reciprocal_R;
%     if reciprocal_R > 0
%         R = 1 / (reciprocal_R + 1e-16);
%     elseif reciprocal_R > -1e-16
%         R = 1 / 1e-16;
%     else
%         R = -1;
%     end
end