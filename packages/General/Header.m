function [header,myline1,myline2] = Header
%Author: Feng-Yi Liao & Yang Zheng
%        SOC Lab @UC San Diego
% Set stuff
myline1 = [repmat('=',1,78),'\n'];
myline2 = [repmat('-',1,78),'\n'];
header  = [' iter|   Obj   | CostDrop |  PFeasi  |  DFeasi  |   Gap   |    rho  | time (s) \n'];
