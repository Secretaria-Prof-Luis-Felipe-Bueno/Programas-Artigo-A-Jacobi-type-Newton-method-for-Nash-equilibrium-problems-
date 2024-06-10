function [restineqy,resteqy]=rest1_y(d1)
 global Delta2
 global k
 restineqy =  d1'*d1 - (Delta2(k))^2;
 resteqy=[];
end %end function
