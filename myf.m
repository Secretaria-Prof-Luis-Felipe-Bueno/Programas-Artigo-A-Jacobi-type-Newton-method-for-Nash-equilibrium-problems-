function [z] = myf(x,y,p)

z = 0;
for i=1:length(p)
   z = z + (x-p(i))^2/((x-p(i))^2+(y(i)-p(i))^2);
end

end

