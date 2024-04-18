function [AbsolutePath] = GetAbsolutePath(RelativePath)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GetAbsolutePath: get the absolute path of a file while the relative 
%                  path is given
% 
% Input argument:
%  RelativePath: the relative path of the file
%
% Ouput argument:
%  AbsolutePath: the absolute path in which the file is saved
%
% 2008.07.21, Zemin Cai
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[stat,mess] = fileattrib(RelativePath);
AbsolutePath = mess.Name;