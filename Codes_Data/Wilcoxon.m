function [p_values] = Wilcoxon (array) 

%NOTE: This is specifically for the work done in B. Shoener's AnMBR
%comparisons. The function will need to be adapted to other research.
%Data used in this analysis was copied from Excel by creating a variable 
%[e.g., x = []).

%Contact: shoener2@illinois.edu

% Input Parameters: 
% array: Array of values to be compared

p_values = [1,46];
for i = 2:24
        a = 1;
        p = ranksum((array(:,a)),(array(:,i)));
        p_values(i-1) = p;
end

for i = 26:48
        a = 25;
        p = ranksum((array(:,a)),(array(:,i)));
        p_values(i-1) = p;
en(
                       
end