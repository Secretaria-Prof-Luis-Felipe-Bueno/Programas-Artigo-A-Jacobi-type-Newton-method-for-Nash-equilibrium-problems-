function [x,k,erro] = nikaidoisoda(x0,tol,mystr)
%manda ver os x,k,t,erro d nikaidoisoda

k=1; psimin=1; psiint=1;
x=x0;
while (psiint>tol)&&(k<1000)
    
%computa ymax e psimax
if strcmp(mystr,'exartigo1_mbom')
    %fazer a funcao psi
    x1=x(1); x2=x(2);
    psyfun = @(y)  (x1^2 + x1*x2 - 5*x1 + (3/2)*x2^2 - x2*x1 - x2) - (y(1)^2 + y(1)*x2 - 5*y(1) + (3/2)*y(2)^2 - y(2)*x1 - y(2)); 
    [ymax,psimax] = fminunc(psyfun,[5;1]);

end
if strcmp(mystr,'exnovo6')
    x1=x(1); x2=x(2);
    psyfun = @(y) ((1/2)*(x1+x2)^2 + sin(x1)+sin(x2) ) - ((1/4)*(y(1)+x2)^2 + sin(y(1)) + (1/4)*(x1+y(2))^2 + sin(y(2)) ); 
    [ymax,psimax] = fminunc(psyfun,[5;1]);
    
end

if abs(psimax)<psimin
    psimin = abs(psimax);
    %x = ymax;
end
psiint = psimax;
x = ymax;
k = k+1;


end %end while
erro = psimin;


end

