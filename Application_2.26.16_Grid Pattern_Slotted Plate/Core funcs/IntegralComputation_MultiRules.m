function [IntegralValue] = IntegralComputation_MultiRules(FuncVec, LowValue, UpperValue, StepLength)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IntegralComputation_MultiRuless: compute the function's integral using
%                                  multi rules: Left Riemanm Sum, Trapezoid
%                                  Sum and Simpson's rule
%
% syntax: IntegralComputation_MultiRuless(...)
%
% Input arguments:
%   FuncVec: the values of the function to be integral(a vector)
%   LowValue: the lower value(scale)
%   UpperValue: the upper value(scale)
%   StepLength: the step length of the independant variable
% Output arguments
%   IntegralValue: the compuation integral result
% 
% Zemin Cai, 2008.5.27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
if nargin < 4
    error('The functions parameters does not match!');
end

PointsNum = max(size(FuncVec));
Step = (UpperValue-LowValue)/(PointsNum-1);
if PointsNum < 2
    error('The integral function should contains more than one values!');
end

RegionsNum = PointsNum - 1;

sum = 0;
if mod(RegionsNum,2)~=0
    if RegionsNum==1
        sum = sum + FuncVec(1)*StepLength;             % Left Riemann sum
    else 
        for k = 1:2:(RegionsNum-2)
            sum = sum + StepLength/3*(FuncVec(k) + 4*FuncVec(k+1) + FuncVec(k+2));  % Simpson sum
        end
        sum = sum + FuncVec(end-1)*StepLength;
    end
else
    if RegionsNum==2
        sum = sum + 0.5*(FuncVec(1)+FuncVec(2))*StepLength;     % Trapezoid sums
        sum = sum + FuncVec(2)*StepLength;
    else
        for k = 1:2:(RegionsNum-2)
            sum = sum + StepLength/3*(FuncVec(k) + 4*FuncVec(k+1) + FuncVec(k+2));
        end
        sum = sum + 0.5*(FuncVec(RegionsNum-1)+FuncVec(RegionsNum))*StepLength;
        sum = sum + FuncVec(end-1)*StepLength;
    end
end
IntegralValue = sum;