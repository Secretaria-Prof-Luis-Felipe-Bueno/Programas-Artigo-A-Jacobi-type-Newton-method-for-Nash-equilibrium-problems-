function [x,y,lambda] = resolvesistema(A,rho)
 x=zeros(2,1); y=zeros(2,1); lambda=zeros(2,1);
 a=A(1,1); b=A(1,2); c=A(2,1); d=A(2,2);
 
    function [g] = testefun(s)
        g1 = (a-b)*s(2) + b + s(3) -rho*max([0,-s(1)]);
        g2 = (c-d)*s(2) + d + s(3) -rho*max([0,s(1)-1]);
        g3 = (a-b)*s(1) + b + s(4) -rho*max([0,-s(2)]);
        g4 = (c-d)*s(1) + d + s(4) -rho*max([0,s(2)-1]);
        g = [g1;g2;g3;g4];
    end

   x0=[x;y];
   fun = @testefun;
   newx = fsolve(fun,x0);
   
   x=[newx(1);1-newx(1)];
   y=[newx(2);1-newx(2)];
   lambda = [newx(3);newx(4)];

end

