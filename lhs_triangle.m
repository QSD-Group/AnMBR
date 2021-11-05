function s=lhs_triangle(min,probable,max,nsample)
% LHS from triangular distribution
% Input (see figure below):
%   min:        minimum value
%   probable:   most probable value
%   max:        maximum value
%   nsample:    number of samples
%
% Output:
%   s:          vector of random variables using LHS
%
%            *
%          * |      *
%        *   |              *
%      *     |                     *    
%    *       |                            *
%  *---------|-----------------------------------*
%  |         |                                   |
% MIN     PROBABLE                              MAX

% If questions, email Diana Byrne, byrne5@illinois.edu

x = lhsdesign(nsample,2);
pd = makedist('triangular','a',min,'b',probable,'c',max);
x(:,1) = icdf(pd,x(:,1));
s = x(:,1);

end

