function [coefs_numerator, coefs_denominator] = taylor_to_pade(coefs, L, M)
    %{
    Description:
      Calculation of Pade approximant numerator and denominator
      coefficients by definition (using linear system based on coefficients
      of Taylor's expansion). See Baker, Graves-Morris §1.1

    Arguments:
      coefs --- coefficients of Taylor's expansion
      L --- maximal degree of numerator polynomial
      M --- maximal degree of denominator polynomial

    Returns:
      coefs_numerator --- coefficients of numerator from h^0 to h^L
      coefs_denominator --- coefficients of denominator from h^0 to h^M
    %}
    
    system_denominator_matrix = zeros(M, M);
    system_denominator_rhs = zeros(M, 1);
    for i=1:M
        system_denominator_rhs(i, 1) = -coefs(L + i + 1);
    end
    for i=1:M
        for j=1:M
            if L - M + i + j < 1
                system_denominator_matrix(i, j) = 0;
            else
                system_denominator_matrix(i, j) = coefs(L - M + i + j);
            end
        end
    end
    coefs_denominator = linsolve(system_denominator_matrix, system_denominator_rhs);
    coefs_denominator = coefs_denominator.';
    coefs_denominator = [1, flip(coefs_denominator, 1)];
    coefs_numerator = zeros(1, L+1);
    for i=1:(L+1)
        convolution = 0;
        for j=1:min(i-1, M)
            convolution = convolution + coefs_denominator(j+1) .* coefs(L-j+1);
        end
        coefs_numerator(i) = coefs(i) + convolution;
    end
end
