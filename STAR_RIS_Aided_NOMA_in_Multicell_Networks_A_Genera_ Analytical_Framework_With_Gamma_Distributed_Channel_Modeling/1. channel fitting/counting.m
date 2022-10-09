function [x, value_pdf, value_cdf] = counting(data, range_min, range_max, N)
%% Function£ºcalculation of PDF and CDF
% Ziyi Xie
% Harbin Institute of Technology
%Input:
%   data                    Data to be counted
%   range_min               Lower bound 
%   range_max               Upper bound
%   N                       Number of scale value
%Output
%   x                       Scale values at x-axis
%   value_pdf               Probability distribution function 
%   value_cdf               Cumulative distribution function
%--------------------------------------------------------------------------
scale = (range_max - range_min)/N;     % scale interval

x =  range_min:scale:range_max;
num = hist(data, x);

[value_pdf, value_cdf] = deal(zeros(size(num)));

value_pdf(1) = num(1)/sum(num)/(scale/2);
value_pdf(end) = num(end)/sum(num)/(scale/2);
value_pdf(2:end-1) = num(2:end-1)/sum(num)/scale;

temp_num = cumsum(num);

value_cdf(1) = temp_num(1)/sum(num);
value_cdf(end) = temp_num(end)/sum(num);
value_cdf(2:end-1) = temp_num(2:end-1)/sum(num);


end