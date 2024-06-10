function [x0,k,t] = jacobi_quad(A,b,x0,tol)
%resolve Ax=b pelo jacobi classico
n = length(b);
k=1;
t=0;
tic
diagmax = zeros(n,n); lmax = diagmax; umax=diagmax;
erro=norm(A*x0-b); kmax=1000; nb=zeros(n,1);

variab=1;

while (erro>tol)&&(k<kmax)
    
    for i=1:n
        diagmax(i,i) = 1/A(i,i);
    for j=1:n
         if i>j
          lmax(i,j) = A(i,j);
         elseif i<j
             umax(i,j) = A(i,j);
         end
     end
    end
        
newd = diagmax*(b-(lmax+umax)*x0);  
x0 = newd
pause
erro = norm(A*x0-b);
k = k+1;
t=t+toc;

if variab==1
             nb(1) = 2*x0(1) + x0(2) - 5;
             nb(2) = 3*x0(2) - x0(1) - 1; 
             b = -nb;
end

end %do while