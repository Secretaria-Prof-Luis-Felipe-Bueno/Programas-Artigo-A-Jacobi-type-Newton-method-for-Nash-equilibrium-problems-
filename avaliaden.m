function [coefs] = avaliaden(x,y)

n=length(x);

if n==2

%f1, z=(1,0)
z=[1;0];
den = norm(x-z)^2 + norm(y-z)^2;
f(1,:)=(-2/den^2)*(norm(y-z)^2)*(x-z);
z=[0;1];
den = norm(x-z)^2 + norm(y-z)^2;
f(2,:)=(-2/den^2)*(norm(y-z)^2)*(x-z);
z=[-1;0];
den = norm(x-z)^2 + norm(y-z)^2;
f(3,:)=(-2/den^2)*(norm(y-z)^2)*(x-z);
z=[0;-1];
den = norm(x-z)^2 + norm(y-z)^2;
f(4,:)=(-2/den^2)*(norm(y-z)^2)*(x-z);
  
%f1,
  
z=[1;0];
den = norm(x-z)^2 + norm(y-z)^2;
f(5,:)=(-2/den^2)*(norm(x-z)^2)*(y-z);
z=[0;1];
den = norm(x-z)^2 + norm(y-z)^2;
f(6,:)=(-2/den^2)*(norm(x-z)^2)*(y-z);
z=[-1;0];
den = norm(x-z)^2 + norm(y-z)^2;
f(7,:)=(-2/den^2)*(norm(x-z)^2)*(y-z);
z=[0;-1];
den = norm(x-z)^2 + norm(y-z)^2;
f(8,:)=(-2/den^2)*(norm(x-z)^2)*(y-z);
  
f=f';

f1=f(:,1:4);
f2=f(:,5:8);




%acha os coefs agora
coefs1 = f1(:,2:4)\(-f1(:,1));
coefs2 = f2(:,2:4)\(-f2(:,1));
coefs1=coefs1';
coefs2=coefs2';

coefs=[coefs1;coefs2];

%provareal1 = f1(:,1) + coefs1(1)*f1(:,2) + coefs1(2)*f1(:,3) + coefs1(3)*f1(:,4)
%provareal2 = f2(:,1) + coefs2(1)*f2(:,2) + coefs2(2)*f2(:,3) + coefs2(3)*f2(:,4)

else

    %f1, z=(1,0)
z=[1];
den = norm(x-z)^2 + norm(y-z)^2;
f(1,:)=(-2/den^2)*(norm(y-z)^2)*(x-z);
z=[-1];
den = norm(x-z)^2 + norm(y-z)^2;
f(2,:)=(-2/den^2)*(norm(y-z)^2)*(x-z);
z=[2];
den = norm(x-z)^2 + norm(y-z)^2;
f(3,:)=(-2/den^2)*(norm(y-z)^2)*(x-z);
%z=[0;-1];
%den = norm(x-z)^2 + norm(y-z)^2;
%f(4,:)=(-2/den^2)*(norm(y-z)^2)*(x-z);
  
%f1,
  
  
f=f';

f1=f(:,1:3);




%acha os coefs agora
coefs1 = f1(:,2:3)\(-f1(:,1));
%coefs2 = f2(:,2:4)\(-f2(:,1));
coefs=coefs1';
%coefs2=coefs2';

%coefs=[coefs1;coefs2];

%provareal1 = f1(:,1) + coefs1(1)*f1(:,2) + coefs1(2)*f1(:,3) + coefs1(3)*f1(:,4)
%provareal2 = f2(:,1) + coefs2(1)*f2(:,2) + coefs2(2)*f2(:,3) + coefs2(3)*f2(:,4)

end

