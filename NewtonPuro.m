function [x,y,f1k,f2k,g1k,g2k,H1k,H2k]=NewtonPuro(x,y,teste,itmax,tol)

 if teste ==1

f1 = @(x,y) -x^3*y^2/3+x^2/2;
f2 = @(x,y) -y^3*x^2/3+y^2/2;

g1 = @(x,y) -x^2*y^2+x;
g2 = @(x,y) -x^2*y^2+y;

H1 = @(x,y) -2*x*y^2+1;
H2 = @(x,y) -2*x^2*y+1;

H1xy = @(x,y) -2*x^2*y;
H2xy = @(x,y) -2*x*y^2;

end


if teste ==2

z1=2;
z2=1;
z3=-1;
h=10^-5;
h2=2*h;
hq=h^2;

f1 = @(x,y) norm(x-z1)^2/(norm(x-z1)^2+norm(y-z1)^2)+norm(x-z2)^2/(norm(x-z2)^2+norm(y-z2)^2)+norm(x-z3)^2/(norm(x-z3)^2+norm(y-z3)^2);
f2 = @(x,y)  norm(y-z1)^2/(norm(x-z1)^2+norm(y-z1)^2)+norm(y-z2)^2/(norm(x-z2)^2+norm(y-z2)^2)+norm(y-z3)^2/(norm(x-z3)^2+norm(y-z3)^2);

g1 = @(x,y) (f1(x+h,y)-f1(x-h,y))/h2;
g2 = @(x,y) (f2(x,y+h)-f2(x,y-h))/h2;

H1 = @(x,y) (f1(x+h,y)-2*f1(x,y)+f1(x-h,y))/hq;
H2 = @(x,y) (f2(x,y+h)-2*f2(x,y)+f2(x,y-h))/hq;

H1xy = @(x,y) (g1(x,y+h)-g1(x,y-h))/h2;
H2xy = @(x,y) (g2(x+h,y)-g2(x+h,y))/h2;


end

    g1k=g1(x,y);
    g2k=g2(x,y);
    it=1;

while it<itmax && g1k^2+g2k^2 > tol
    it=it+1;
    H1k=H1(x,y);
    H2k=H2(x,y);


    % descomentar para fazer Newton mas com modelos convexos.
    %if H1k<0.1
    %    H1k=1;
    %end

    %if H2k<0.1
    %     H2k=1;
    % end




    H1xyk=H1xy(x,y);
    H2xyk=H2xy(x,y);

    g1k=g1(x,y);
    g2k=g2(x,y);

    b=-[g1k;g2k];
    M=[H1k H1xyk; H2xyk H2k];


    d=M\b;

    x=x+d(1);
    y=y+d(2);

end

f1k=f1(x,y);
f2k=f2(x,y);

end
