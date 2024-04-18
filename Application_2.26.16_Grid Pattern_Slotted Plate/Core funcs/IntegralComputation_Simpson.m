function [IntegralValue] = IntegralComputation_Simpson(FuncVec, LowValue, UpperValue, StepLength)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IntegralComputation_Simpson: compute the function's integral using
%                              Simpson' rule
%
% syntax: IntegralComputation_Simpson(...)
%
% Input arguments:
%   FuncVec: the values of the function to be integral(a vector)
%   LowValue: the lower value(scale)
%   UpperValue: the upper value(scale)
%   StepLength: the step length of the independant variable
% Output arguments
%   IntegralValue: the compuation integral result
% 
% Zemin Cai, 2008.6.27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
if nargin < 4
    error('The functions parameters does not match!');
end

PointsNum = max(size(FuncVec));
if PointsNum < 2
    error('The integral function should contains more than one values!');
end

RegionsNum = PointsNum - 1;

sum = 0;
if mod(RegionsNum, 2)~=0   
    if RegionsNum == 1           % use the Trapezoid sums directly
        sum = sum + 0.5*(FuncVec(1) + FuncVec(2))*StepLength;
    else
        for k = 1:2:(RegionsNum-2)
            sum = sum + StepLength/3*(FuncVec(k) + 4*FuncVec(k+1) + FuncVec(k+2));  % Simpson sum
        end
        sum = sum + 0.5*(FuncVec(end-1) + FuncVec(end))*StepLength;   % Trapezoid sums to compute the last region
    end
else
    for k = 1:2:RegionsNum
        sum = sum + StepLength/3*(FuncVec(k) + 4*FuncVec(k+1) + FuncVec(k+2));
    end
end
IntegralValue = sum;