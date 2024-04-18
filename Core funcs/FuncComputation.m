function [FuncValue] = FuncComputation(FuncStr, FuncVariate, VariateVec)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FuncComputation: compute the function value while a string of the
%                  function is given(one dimension)
%
% syntax: FuncValue(...)
%
% Input arguments:
%   FuncStr: the function (a string)
%   FuncVariate: the variate in the function (string)
%   VariateVec: the variate value vector
% Output arguments:
%   FuncValues: the compuation function value vector
% 
% Zemin Cai, 2008.5.6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
Function = inline(FuncStr, FuncVariate);
Function_Value = feval(Function, VariateVec);

FuncValue = Function_Value;