function [x0,erro] = Elvis(n,m,gama,x0,lamb0)
%esqueleto do laco externo do elvis (IALEM)
%efinir gama, sigma, x0, lambda0
kmaxfora = 1000;
kcont=0; condicao=0;

   alfa = 1e-2; %alfa em (0,1)
   tol = 1e-4; kmax=1000;
   minvaltol = 1e-6;
   fixvaltol=1;
   maxvaltol=1e10;
   zerodograd=0;
   lambdai=0;
   mystr = 'elvis';
  [x0,A1,A2,b1,b2,c1,c2] = montadados(mystr,n,m);


while (k<kmaxfora)&&(condicao==0)

 xvelho = x0;
 %resolva o negocio sem o e
 x0 = rodar geral...
     
 %atualiza lambda
 for j=1:n+m
 lamb0(j) = lamb0(j)*exp((-x0(j))/gama);
 end
 
 erro = norm(x0-xvelho);
 if erro<1e-4
    condicao=1; 
 end
 
 kcont=kcont+1;
end
 
 
  
  


end

