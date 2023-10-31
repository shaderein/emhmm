function [hNull X2 pValue] = chi2ind(o, alpha)
    %#  CHISQUARETEST  Pearson's Chi-Square test of independence
    %#
    %#    @param o          The Contignecy Table of the joint frequencies
    %#                      of the two events (attributes)
    %#    @param alpha      Significance level for the test
    %#    @return hNull     hNull = 1: null hypothesis accepted (independent)
    %#                      hNull = 0: null hypothesis rejected (dependent)
    %#    @return pValue    The p-value of the test (the prob of obtaining
    %#                      the observed frequencies under hNull)
    %#    @return X2        The value for the chi square statistic
    %#

    %# o:     observed frequency
    %# e:     expected frequency
    %# dof:   degree of freedom

    [r c] = size(o);
    dof = (r-1)*(c-1);

    %# e = (count(A=ai)*count(B=bi)) / N
    e = sum(o,2)*sum(o,1) / sum(o(:));

    %# [ sum_r [ sum_c ((o_ij-e_ij)^2/e_ij) ] ]
    X2 = sum(sum( (o-e).^2 ./ e ));

    %# p-value needed to reject hNull at the significance level with dof
    pValue = 1 - chi2cdf(X2, dof);
    hNull = (pValue > alpha);

    %# X2 value needed to reject hNull at the significance level with dof
    %#X2table = chi2inv(1-alpha, dof);
    %#hNull = (X2table > X2);

end