function [FuncValue] = IntegralComputation_Infinite_CloseZero(FuncStr,FuncVariate, torrence)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IntegralComputation_Infinite_CloseZero: compute the function's integral while a string of the
%                  function is given(the low is close to zero but not zero because zero is a sigular point)
%
% syntax: IntegralComputation_Infinite_CloseZero(...)
%
% Input arguments:
%   FuncStr: the function (a string)
%   FuncVariate: the variate in the function (string)
%   torrence: error torrence
% Output arguments:
%   FuncValue: the compuation integral result
% 
% Zemin Cai, 2008.6.25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Function = inline(FuncStr, FuncVariate);
r = 2;        % step length;
i = 1;
%epsilon = 1.0e-32;
epsilon = torrence;
CloseZero = 10e-8;

r1 = r^i;
r2 = r^(i+1);
%disp('before');
FuncIncre = abs(real(quadl(Function, r1, r2)));
%disp('after...')
%FuncValueTemp = quadl(Function, 0, r2);
FuncValueTemp = quadl(Function, CloseZero, r2);
while(FuncIncre > epsilon)
%while(i<10)
    %FuncValueTemp = quadl(Function, 0, r2);
    FuncValueTemp = quadl(Function, CloseZero, r2);
    i = i + 1;
    r1 = r^i;
    r2 = r^(i+1);
    %FuncIncre = quadl(Function, r1, r2);
    FuncIncre = abs(real(quadl(Function, r1, r2)));
end
% r2 = r^torrence;
% FuncValueTemp = quadl(Function, 0, r2);

FuncValue = FuncValueTemp;

