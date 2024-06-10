function [restineqx,resteqx]=rest1_x(d1)
 global Delta1
 global k
 restineqx =  d1'*d1 - (Delta1(k))^2;
 resteqx=[];
end %end function

