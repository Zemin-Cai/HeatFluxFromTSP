function [FuncValue, Legal] = IntegralComputation(FuncStr,LowValue,UpperValue,FuncVariate)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IntegralComputation: compute the function's integral while a string of the
%                  function is given(with the lower and upper integral value)
%
% syntax: IntegralComputation(...)
%
% Input arguments:
%   FuncStr: the function (a string)
%   LowValue: the lower value(scale)
%   UpperValue: the upper value(scale)
%   FuncVariate: the variate in the function (string)
% Output arguments:
%   FuncValue: the compuation integral result
%   Legal: to decide whether the result is legal, 0: not legal; 1: legal
% 
% Zemin Cai, 2008.5.16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Function = inline(FuncStr, FuncVariate);
FuncValue = quadl(Function,LowValue,UpperValue);        % 
%FuncValue = quad(Function,LowValue,UpperValue);
Legal = 1;
if(imag(FuncValue) ~= 0)
    Legal = 0;
end