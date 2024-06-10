% script para rodar o geral, sem argumentos. qqr argumento mudamos aqui msm
% clear all

%flaghess = 'hessiana';
%flaghess = 'bfgs';
%flaghess = 'identidade'; multhess = 1;
flaghess = 'ifdois'; multhess = 1;
%flaghess = 'tfixo1';
%flaghess = 'newtonpuro';
%flaghess = 'jacobiintern';
%flaghess = 'jacobidois';
%flaghess = 'quadrado'; multhess=1;

%mystr = 'quadratica';
%mystr = 'testepaper';
%mystr = 'problemayuan';
%mystr = 'quadmisto';
%mystr = 'exponencial';
%mystr = 'separavel';
%mystr = 'cubico';
%mystr = 'diferentes';
%mystr = 'quadratica_yuan';
%mystr = 'quartico';
%mystr = 'quadtestes1';
%mystr = 'quadtestes2';
%mystr = 'quadtestes3';
%mystr = 'quadtestes4';
%mystr = 'quadtestes5';
%mystr = 'exp1';
%mystr = 'sin2';
%mystr = 'arcsin';
%mystr = 'sincos';
%mystr = 'testapol';
%mystr = 'testacirc';
%mystr = 'gradlimmin';
%mystr = 'vacina';
%mystr = 'cubicasem';
%mystr = 'exbueno';
%mystr = 'exbueno2';
%mystr='exacho';
%mystr = 'aplicacao';
%mystr = 'hartung';
%mystr = 'hartunglagrangeano';
%mystr = 'penalcubic';
%mystr = 'formation';
%mystr = 'formationzero';
%mystr = 'aplbueno';
%mystr = 'aplbuenovdd';
%mystr = 'aplbuenosem';
%mystr = 'aplbuenoy';
%mystr = 'aplbuenox';
%mystr = 'aplbuenoyoutro';

%mystr = 'exartigo1_mbom';
%mystr = 'exartigo2_mruim';
%mystr = 'exartigo3_indef';
%mystr = 'exartigo4_cubic';
%mystr = 'exartigo5_futebol';
%mystr = 'futebolzero';
%mystr = 'exartigo6_vacina';
%mystr = 'exnovo6';
mystr = 'exrevisao';
mystr = 'exrevisao2';
%mystr = 'maxes';

if (strcmp(mystr,'aplicacao')==0)
    r1 = 1; r2=1;
    ro1 = 1; ro2=1;
    %vaimontadados = 1;
end
if(strcmp(mystr,'hartung')==1)||( (strcmp(mystr,'hartunglagrangeano')==1)||(strcmp(mystr,'penalcubic'))  )
    vaimontadados = 0;
end

treset=1;
maint=0;
fazcomlambdas=1;
mytic = cputime;

fazcholmod=0;

vaidisplay = 0;
vaigrafico=1;
%vaimontadados = 0;

if vaimontadados==1
    alfa = 1e-6; %alfa em (0,1)
    %alfa=0;
    tol = 1e-5; kmax=1000;
    minvaltol = 1e-6;
    fixvaltol=1;
    maxvaltol=1e10;
    %zerodograd=1e-4;
    zerodograd=0;
    lambdai=0;
    n = 1;
    m = 1;
    [A1,A2,b1,b2,c1,c2,x0] = montadados(mystr,n,m);
end

%calcula lambdamin e lambdamax
betta = 100000;
tetta = 0.01;
meugama=0.00001;

%tetta=0;
%meugama=0;

%if (n==1)&&(m==1)&&(vaimontadados==1)
xvec=[x0(1)];
yvec=[x0(2)];
if (n>1)
    xvec2=[x0(3)];
    yvec2=[x0(4)];
end
%end

if strcmp(flaghess,'tfixo1')
    xvecnewton=[x0(1)];
    yvecnewton=[x0(2)];
end

grad1 = Feval(x0,A1,b1,c1,1,mystr,1,r1,ro1);
grad2 = Feval(x0,A2,b2,c2,1,mystr,2,r2,ro2);
erro = norm(grad1)+norm(grad2);

if vaidisplay==1
    disp('x0: '); disp(x0);
    disp('grad1: '); disp(grad1);
    disp('grad2: '); disp(grad2);
    disp('erro: '); disp(erro);
    pause
end

srstring='sr';

k=0; t=1; tau=0.99;
maint=cputime-mytic;
x0jacobi = x0;
while((erro>tol)&&(k<kmax))
    
    mytic = cputime;
    
    if strcmp(mystr,'formation3')
        b1 = x0(1:end/2);
        b2 = x0(end/2+1:end);
    end
    
    %matrizes pos defas sao hessianas mesmo
    if strcmp(flaghess,'hessiana')||strcmp(flaghess,'ifdois')||strcmp(flaghess,'tfixo1')||strcmp(flaghess,'jacobiintern')||strcmp(flaghess,'jacobidois')||strcmp(flaghess,'quadrado')||strcmp(flaghess,'newtonpuro')
        H1 = Feval(x0,A1,b1,c1,2,mystr,1,r1,ro1);
        H2 = Feval(x0,A2,b2,c2,2,mystr,2,r2,ro2);
    end
    
    %bfgs
    if strcmp(flaghess,'bfgs')
        if k==0
            H1 = Feval(x0,A1,b1,c1,2,mystr,1,r1,ro1);
            H2 = Feval(x0,A2,b2,c2,2,mystr,2,r2,ro2);
            %H1 = 1*eye(n,n);
            %H2 = 1*eye(m,m);
        else
            sk1 = t*newd(1:n);
            sk2 = t*newd(n+1:n+m);
            yk1 = grad1-gradvelho1;
            yk2 = grad2-gradvelho2;
            vk1 = H1*sk1;
            vk2 = H2*sk2;
            H1 = H1 + (yk1*yk1')/(yk1'*sk1) - (vk1*vk1')/(sk1'*vk1);
            H2 = H2 + (yk2*yk2')/(yk2'*sk2) - (vk2*vk2')/(sk2'*vk2);
        end %end if do bfgs
    end  %end do flaghess do bfgs
    
    
    %testando identidade
    if strcmp(flaghess,'identidade')
        H1 = multhess*eye(n,n);
        H2 = multhess*eye(m,m);
    end
    
    if vaidisplay==1
        disp('H1: '); disp(H1);
        disp('H2: '); disp(H2);
        pause
    end
    
    mista1 = Feval(x0,A1,b1,c1,3,mystr,1,r1,ro1);
    mista2 = Feval(x0,A2,b2,c2,3,mystr,2,r2,ro2);
    bigmatrix = [H1 t*mista1; t*mista2 H2];
    biggrad = [grad1;grad2];
    
    if vaidisplay==1
        disp('mista1: ');disp(mista1);
        disp('mista2: ');disp(mista2);
        disp('bigmatrix: ');disp(bigmatrix);
        disp('biggrad: ');disp(biggrad);
        pause
    end
    
    if (strcmp(flaghess,'ifdois')||strcmp(flaghess,'quadrado') )%||strcmp(flaghess,'tfixo1')
        
        if fazcholmod==0
            
            %condnumber = cond(bigmatrix);
            if (strcmp(mystr,'problemayuan')==0)&&(strcmp(mystr,'vacina')==0)&&(strcmp(mystr,'sin2')==0)
                minval = eigs(bigmatrix,1,srstring);
            else
                minval = min(eig(bigmatrix));
            end
            
            if vaidisplay==1
                disp('cond: '); disp(cond(bigmatrix));
                disp('minval: '); disp(minval);
            end
            
            increment=1e-4;
            while(real(minval)<minvaltol)&&(norm(grad1)>zerodograd)&&(norm(grad2)>zerodograd)
                
                bigmatrix=bigmatrix+(lambdai*abs(minval)+increment*fixvaltol)*eye(n+m,n+m);
                %bigmatrix = (abs(minval)+minvaltol)*eye(n+m,n+m);
                
                if (strcmp(mystr,'problemayuan')==0)&&(strcmp(mystr,'vacina')==0)&&(strcmp(mystr,'sin2')==0)
                    minval = eigs(bigmatrix,1,srstring);
                else
                    minval = min(eig(bigmatrix));
                end
                %minval = min(eig(bigmatrix));
                increment=2*increment;
                
                if vaidisplay==1
                    disp('cond: '); disp(cond(bigmatrix));
                    disp('minvaltol: '); disp(minvaltol);
                    disp('maxvaltol: '); disp(maxvaltol);
                end
            end
            
        else %ai cholmod nao eh zero
            bigmatrix = cholmod(bigmatrix,1e-6);
            minval = min(eig(bigmatrix));
        end %end cholmod
        
    end
    
    if vaidisplay==1
        disp('mista1: ');disp(mista1);
        disp('mista2: ');disp(mista2);
        disp('bigmatrix: ');disp(bigmatrix);
        disp('biggrad: ');disp(biggrad);
        pause
    end
    
    if (strcmp(flaghess,'jacobiintern')==0 )&&(strcmp(flaghess,'jacobidois')==0 )%&&(strcmp(flaghess,'tfixo1')==0)
        
        if (norm(grad1)>zerodograd)&&(norm(grad2)>zerodograd)
            newd = -bigmatrix\biggrad;
        end
        
        if norm(grad1)<=zerodograd
            newd(1:n)=zeros(n,1);
            disp('H2: '); disp(H2);
            minlam2 = eigs(H2,1,srstring);
            while minlam2<=0.5
                H2=H2+eye(m,m);
                minlam2 = eigs(H2,1,srstring);
            end
            newd(n+1:n+m) = -H2\grad2;
        end
        if norm(grad2)<=zerodograd
            newd(n+1:n+m)=zeros(m,1);
            minlam1 = eigs(H1,1,srstring);
            while minlam1<=0.5
                H1=H1+eye(m,m);
                minlam1 = eigs(H1,1,srstring);
            end
            newd(1:n) = -H1\grad1;
        end
        
    elseif strcmp(flaghess,'jacobiintern')
        vessel=zeros(n+m,n+m);
        diagmax = vessel; lmax = vessel; umax=vessel;
        for i=1:n+m
            diagmax(i,i) = 1/bigmatrix(i,i);
            for j=1:n+m
                if i>j
                    lmax(i,j) = bigmatrix(i,j);
                elseif i<j
                    umax(i,j) = bigmatrix(i,j);
                end
            end
        end
        
        newd = diagmax*(-biggrad-(lmax+umax)*x0);
        
        x0 = newd;
        
        %pause
        grad1 = Feval(x0,A1,b1,c1,1,mystr,1,r1,ro1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2,r2,ro2);
        erro = norm(grad1)+norm(grad2);
        k = k+1;
        
        maint=maint+(cputime-mytic);
        
    else  %ai strcmp(flaghess,'jacobidois')
        newd = x0;
        
        newd(1:n) = H1\(-biggrad(1:n)-mista1*x0(n+1:n+m));
        newd(n+1:n+m) = H2\(-biggrad(n+1:n+m)-mista2*x0(1:n));
        
        x0 = newd;
        
        grad1 = Feval(x0,A1,b1,c1,1,mystr,1,r1,ro1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2,r2,ro2);
        erro = norm(grad1)+norm(grad2);
        k = k+1;
        maint=maint+(cputime-mytic);
    end %end os 3 strcmps
    
    if vaidisplay==1
        disp('newd: '); disp(newd);
        pause
    end
    
    gradval1 = Feval([x0(1:n);x0(n+1:n+m)+t*newd(n+1:n+m)],A1,b1,c1,1,mystr,1,r1,ro1);
    gradval2 = Feval([x0(1:n)+t*newd(1:n);x0(n+1:n+m)],A2,b2,c2,1,mystr,2,r2,ro2);
    produto1 = (gradval1)'*newd(1:n);
    produto2 = (gradval2)'*newd(n+1:n+m);
    direito1 = -tetta*norm(gradval1)*norm(newd(1:n));
    direito2 = -tetta*norm(gradval2)*norm(newd(n+1:n+m));
    numero1 = Feval(x0+t*newd,A1,b1,c1,0,mystr,1,r1,ro1)-Feval([x0(1:n);x0(n+1:n+m)+t*newd(n+1:n+m)],A1,b1,c1,0,mystr,1,r1,ro1)-alfa*t*produto1;
    numero2 = Feval(x0+t*newd,A2,b2,c2,0,mystr,2,r2,ro2)-Feval([x0(1:n)+t*newd(1:n);x0(n+1:n+m)],A2,b2,c2,0,mystr,2,r2,ro2)-alfa*t*produto2;
    if strcmp(flaghess,'quadrado')
        numero1 = Feval(x0+t*newd,A1,b1,c1,0,mystr,1,r1,ro1)-Feval([x0(1:n);x0(n+1:n+m)+t*newd(n+1:n+m)],A1,b1,c1,0,mystr,1,r1,ro1)+alfa*t*(newd(1:n))'*newd(1:n);
        numero2 = Feval(x0+t*newd,A2,b2,c2,0,mystr,2,r2,ro2)-Feval([x0(1:n)+t*newd(1:n);x0(n+1:n+m)],A2,b2,c2,0,mystr,2,r2,ro2)+alfa*t*(newd(n+1:n+m))'*newd(n+1:n+m);
    end
    
    %coloquei esse artificial pra caso a norma dos grad for zero
    if (norm(grad1)<=zerodograd)&&(t<tau)
        numero1=0;
    end
    if (norm(grad2)<=zerodograd)&&(t<tau)
        numero2=0;
    end
    
    if (strcmp(flaghess,'tfixo1')==0)&&(strcmp(flaghess,'newtonpuro')==0)&&( strcmp(flaghess,'jacobiintern')==0 )&&( strcmp(flaghess,'jacobidois')==0 )
        passa=1;
    else
        passa=0;
    end
    
    if (strcmp(flaghess,'tfixo1')==0)&&(strcmp(flaghess,'newtonpuro')==0)&&( strcmp(flaghess,'jacobiintern')==0 )&&( strcmp(flaghess,'jacobidois')==0 )
        
        if (  (  (real(minval)<minvaltol)&&(norm(grad1)>zerodograd)&&(norm(grad2)>zerodograd)   )||(...
                (numero1>0)&&(abs(numero2)>zerodograd)   )||(  (numero2>0)&&(abs(numero1)>zerodograd) )||(...
                (produto1>0)&&(abs(numero2)>zerodograd)&&(fazcomlambdas==0)   )||(  (produto2>0)&&(abs(numero1)>zerodograd)&&(fazcomlambdas==0) )||(...
                (produto1>direito1)&&(abs(numero2)>zerodograd)&&(fazcomlambdas~=0)   )||(  (produto2>direito2)&&(...
                abs(numero1)>zerodograd)&&(fazcomlambdas~=0) )||(  (produto1>0)&&(abs(numero2)>zerodograd)&&(fazcomlambdas==0)   )||(...
                ((norm(newd(1:n))>betta*norm(gradval1)  )||(norm(newd(1:n))<meugama*norm(gradval1))  )&&(fazcomlambdas~=0) )||(  ((norm(newd(n+1:n+m))<meugama*norm(gradval2)  )||(norm(newd(n+1:n+m))>betta*norm(gradval2))  )&&(fazcomlambdas~=0) )     )
            
            if vaidisplay==1
                disp('gradval1: '); disp(gradval1);
                disp('gradval2: '); disp(gradval2);
                disp('numero1: '); disp(numero1);
                disp('numero2: '); disp(numero2);
                pause
            end
            
            %disp('passei no if');
            %pause
            
            
            %aqui nao aceita o passo, apenas diminui o t
            if treset == 0
                t=t/2;
                if vaidisplay==1
                    disp('t: '); disp(t);
                    disp('recusou o passo');
                end
                %pause
            else %if treset=1
                
                %nesse teste comentado, vc reseta o t
                %passa = 0;
                while (  (  (real(minval)<minvaltol)&&(norm(grad1)>zerodograd)&&(norm(grad2)>zerodograd)   )||(...
                (numero1>0)&&(abs(numero2)>zerodograd)   )||(  (numero2>0)&&(abs(numero1)>zerodograd) )||(...
                (produto1>0)&&(abs(numero2)>zerodograd)&&(fazcomlambdas==0)   )||(  (produto2>0)&&(abs(numero1)>zerodograd)&&(fazcomlambdas==0) )||(...
                (produto1>direito1)&&(abs(numero2)>zerodograd)&&(fazcomlambdas~=0)   )||(  (produto2>direito2)&&(...
                abs(numero1)>zerodograd)&&(fazcomlambdas~=0) )||(  (produto1>0)&&(abs(numero2)>zerodograd)&&(fazcomlambdas==0)   )||(...
                ((norm(newd(1:n))>betta*norm(gradval1)  )||(norm(newd(1:n))<meugama*norm(gradval1))  )&&(fazcomlambdas~=0) )||(  ((norm(newd(n+1:n+m))<meugama*norm(gradval2)  )||(norm(newd(n+1:n+m))>betta*norm(gradval2))  )&&(fazcomlambdas~=0) )     )
                    %(produto1>0)||(produto2>0)||(numero1>0)||(numero2>0)||(  (real(minval)<minvaltol)&&(norm(grad1)>zerodograd)&&(norm(grad2)>zerodograd)   )
                    
                    
                    t=t/2;
                    %disp('recusou o passo');
                    bigmatrix = [H1 t*mista1; t*mista2 H2];
                    if (strcmp(mystr,'problemayuan')==0)&&(strcmp(mystr,'vacina')==0)&&(strcmp(mystr,'sin2')==0)
                        minval = eigs(bigmatrix,1,srstring);
                    else
                        minval = min(eig(bigmatrix));
                    end
                    increment=1e-4;
                    
                    if (real(minval)<minvaltol)
                        if fazcholmod==0
                            
                            while (real(minval) < minvaltol)
                                bigmatrix=bigmatrix+(lambdai*abs(minval)+(increment*fixvaltol))*eye(n+m,n+m);
                                if (strcmp(mystr,'problemayuan')==0)&&(strcmp(mystr,'vacina')==0)&&(strcmp(mystr,'sin2')==0)
                                    minval = eigs(bigmatrix,1,srstring);
                                else
                                    minval = min(eig(bigmatrix));
                                end
                                increment=2*increment;
                            end
                            
                        else
                            bigmatrix = cholmod(bigmatrix,1e-5);
                        end
                        
                    end
                    newd = -bigmatrix\biggrad;
                    
                    
                    gradval1 = Feval([x0(1:n);x0(n+1:n+m)+t*newd(n+1:n+m)],A1,b1,c1,1,mystr,1,r1,ro1);
                    gradval2 = Feval([x0(1:n)+t*newd(1:n);x0(n+1:n+m)],A2,b2,c2,1,mystr,2,r2,ro2);
                    produto1 = (gradval1)'*newd(1:n);
                    produto2 = (gradval2)'*newd(n+1:n+m);
                    direito1 = -tetta*norm(gradval1)*norm(newd(1:n));
                    direito2 = -tetta*norm(gradval2)*norm(newd(n+1:n+m));
                    numero1 = Feval(x0+t*newd,A1,b1,c1,0,mystr,1,r1,ro1)-Feval([x0(1:n);x0(n+1:n+m)+t*newd(n+1:n+m)],A1,b1,c1,0,mystr,1,r1,ro1)-alfa*t*produto1;
                    numero2 = Feval(x0+t*newd,A2,b2,c2,0,mystr,2,r2,ro2)-Feval([x0(1:n)+t*newd(1:n);x0(n+1:n+m)],A2,b2,c2,0,mystr,2,r2,ro2)-alfa*t*produto2;
                    if strcmp(flaghess,'quadrado')
                        numero1 = Feval(x0+t*newd,A1,b1,c1,0,mystr,1,r1,ro1)-Feval([x0(1:n);x0(n+1:n+m)+t*newd(n+1:n+m)],A1,b1,c1,0,mystr,1,r1,ro1)+alfa*t*(newd(1:n))'*newd(1:n);
                        numero2 = Feval(x0+t*newd,A2,b2,c2,0,mystr,2,r2,ro2)-Feval([x0(1:n)+t*newd(1:n);x0(n+1:n+m)],A2,b2,c2,0,mystr,2,r2,ro2)+alfa*t*(newd(n+1:n+m))'*newd(n+1:n+m);
                    end
                    %coloquei esse artificial pra caso a norma dos grad for zero
                    
                    %disp(produto1>direito1);
                    %disp(produto2>direito2);
                    %disp(real(minval)<minvaltol);
                    %disp(   norm(newd(n+1:n+m))<meugama*norm(gradval2)   );
                    %disp(norm(newd(1:n))>betta*norm(gradval1));
                    %disp(norm(newd(1:n))<meugama*norm(gradval1) );
                    %disp(norm(newd(n+1:n+m))<meugama*norm(gradval2)  );
                    %disp(norm(newd(n+1:n+m))>betta*norm(gradval2));
                    %pause
                    
                end
                

                
            end %treset
            maint = maint+(cputime-mytic);
            
        else %nesse caso aceita o passo
            if passa == 1
                
                if (norm(grad1)<=zerodograd)&&(t<tau)
                    numero1=0;
                end
                if (norm(grad2)<=zerodograd)&&(t<tau)
                    numero2=0;
                end
                
                if norm(grad1)<=zerodograd
                    newd(1:n)=zeros(n,1);
                    minlam2 = eigs(H2,1,srstring);
                    while minlam2<=0
                        H2=H2+eye(m,m);
                        minlam2 = eigs(H2,1,srstring);
                    end
                    newd(n+1:n+m) = -H2\grad2;
                end
                if norm(grad2)<=zerodograd
                    newd(n+1:n+m)=zeros(m,1);
                    minlam1 = eigs(H1,1,srstring);
                    while minlam1<=0
                        H1=H1+eye(m,m);
                        minlam1 = eigs(H1,1,srstring);
                    end
                    newd(1:n) = -H1\grad1;
                end
                
                x0 = x0+t*newd
                grad1
                grad2
                t
                erro
                %teste = Feval([x0(1:n)+newd(1:n);x0(n+1:n+m)],A1,b1,c1,0,mystr,1,r1,ro1)-Feval(x0,A1,b1,c1,0,mystr,1,r1,ro1)-Feval(x0,A1,b1,c1,1,mystr,1,r1,ro1)'*(newd(1:n)) - (newd(1:n))'*(Feval(x0,A1,b1,c1,2,mystr,1,r1,ro1)*newd(1:n))
                pause

                
                if strcmp(mystr,'aplicacao')
                    x0(2)=1-x0(1);
                    x0(4)=1-x0(3);
                    x0
                end
                
                if vaidisplay==1
                    disp('------------------------');
                    disp('f1: '); disp(Feval(x0,A1,b1,c1,0,mystr,1,r1,ro1));
                    disp('f2: '); disp(Feval(x0,A2,b2,c2,0,mystr,2,r2,ro2));
                    disp('grad1: ' ); disp(grad1);
                    disp('grad2: ' ); disp(grad2);
                    disp('gradval1: '); disp(gradval1'*newd(1:n));
                    disp('gradval2: '); disp(gradval2'*newd(n+1:n+m));
                    %disp('newd: '); disp(newd);
                    disp('erro: '); disp(erro);
                    %disp('x0: '); disp(x0');
                    disp('t: '); disp(t);
                    pause
                end
                
                %esses sao pra se tiver bfgs
                if strcmp(flaghess,'bfgs')
                    gradvelho1 = grad1;
                    gradvelho2 = grad2;
                end
                
                grad1 = Feval(x0,A1,b1,c1,1,mystr,1,r1,ro1);
                grad2 = Feval(x0,A2,b2,c2,1,mystr,2,r2,ro2);
                erro = norm(grad1)+norm(grad2);
                
                k = k+1;
                                
                if vaidisplay==1
                    disp('grad1: '); disp(grad1);
                    disp('grad2: '); disp(grad2);
                    disp('erro: '); disp(erro);
                    pause
                end
                
                maint=maint+(cputime-mytic);
                
                if ((n==1)&&(m==1)&&(vaimontadados==1)&&(vaigrafico==1))
                    xvec = [xvec; x0(1)];
                    yvec = [yvec; x0(2)];
                    
                    scatter(xvec,yvec);
                    hold on
                    switch mystr
                        case 'quadratica'
                            fun1 = @(x) x'*(A1*x)+(b1')*x + c1;
                            fun2 = @(x) x'*(A2*x)+(b2')*x + c2;
                        case 'testepaper'
                            fun1 = @(x,y) (x'*(A1*x))/2+((b1*y - c1)')*x;
                            fun2 = @(x,y) (x'*(A2*x))/2+((b2*y - c2)')*x;
                        case 'problemayuan'
                            fun1 = @(x,y) ((b1*y - c1)')*x;
                            fun2 = @(x,y) ((b2*y - c2)')*x;
                        case 'quadmisto'
                            fun1 = @(x,y) (x'*(A1*x) - y'*(A2*y))/2 + x'*(b1*y) + c1'*x + c2'*y;
                            fun2 = @(x,y) (-x'*(A1*x) + y'*(A2*y))/2 + x'*(b2*y) + c1'*x + c2'*y;
                        case 'exponencial'
                            fun1 = @(x,y) sum(A1.*exp(b1.*(x.*y)));
                            fun2 = @(x,y) sum(A2.*exp(b2.*(x.*y)));
                        case 'separavel'
                            fun1 = @(x,y) sum(exp(x).*(c1+x-log(exp(x)+exp(y))));
                            fun2 = @(x,y) sum(exp(x).*(c2+x-log(exp(x)+exp(y))));
                        case 'cubico'
                            fun1 = @(x,y) (x.^3)'*(y.^2)/3 + (norm(x))^2/2;
                            fun2 = @(x,y) (y.^3)'*(x.^2)/3 + (norm(y))^2/2;
                    end
                    %pause
                end
                
                if ((n==2)&&(m==2)&&((strcmp(mystr,'aplbueno'))||strcmp(mystr,'aplbuenovdd')))
                    %xvec = [xvec; x0(1)];
                    %yvec = [yvec; x0(2)];
                    
                    %hold on
                    %xvec2 = [xvec2; x0(3)];
                    %yvec2 = [yvec2; x0(4)];
                    %pause
                    
                    xvec = [xvec; x0(1)];
                    yvec = [yvec; x0(3)];
                    
                    hold on
                    
                end
                
                
            end %if passa
            
            if treset==1
                t = 1;
            end
            
        end % end if
        
    end % end if das strings
    
    
    if strcmp(flaghess,'tfixo1')
        %nesse caso atualiza simples
        x0 = x0+t*newd
        %pause
        
        
        xvecnewton = [xvecnewton; x0(1)];
        yvecnewton = [yvecnewton; x0(2)];
        %x0(2)=1-x0(1);
        %x0(4)=1-x0(3);
        %x0;
        
        if (n==1)&&(m==1)
        xvec = [xvec; x0(1)];
        yvec = [yvec; x0(2)];
        end
         if (n==2)&&(m==2)
        xvec = [xvec; x0(1)];
        yvec = [yvec; x0(3)];
        end
        %scatter(xvec,yvec,'red');
        if n>1
            xvec2 = [xvec2; x0(3)];
            yvec2 = [yvec2; x0(4)];
        end
        hold on
        
        
        grad1 = Feval(x0,A1,b1,c1,1,mystr,1,r1,ro1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2,r2,ro2);
        erro = norm(grad1)+norm(grad2);
        k = k+1;
        
    end
    
    
    %disp('k, erro, norm(d), t');
    %disp([k;erro;norm(newd);t]');
    %pause
    
    
    
end %end while

if (strcmp(flaghess,'ifdois'))&&(n>1)
scatter([xvec;xvec2],[yvec;yvec2],'green');
%if n>1
%    scatter(xvec2,yvec2,'cyan');
%end
end

if (n==11)
    hold on

    [gridx,gridy]=meshgrid(-5:0.05:5,-5:0.05:5);
    [s1,s2]=size(gridx);
    z=zeros(s1,s2); w=zeros(s1,s2);
    %s1=100; s2=100;
    %z=zeros(100,1); w=zeros(100,1); t=linspace(-20,20);
    
    for i=1:s1
        for j=1:s2
            z(i,j) = Feval([gridy(i),gridy(j)],A1,b1,c1,0,mystr,1,r1,ro1);
            %z(i,j) = Feval([gridx(i)/sqrt(2)-gridy(i)/sqrt(2),-gridx(j)/sqrt(2)+gridy(j)/sqrt(2)],A1,b1,c1,0,mystr,1,r1,ro1);
            w(i,j) = Feval([gridy(i),gridy(j)],A2,b2,c2,0,mystr,2,r2,ro2);
            %w(i,j) = Feval([gridx(i)/sqrt(2)-gridy(i)/sqrt(2),-gridx(j)/sqrt(2)+gridy(j)/sqrt(2)],A2,b2,c2,0,mystr,2,r2,ro2);
            %tes(i,j) = gridy(i)^2 + gridy(j)^2;
            %z(i) =  Feval([t(i);1],A1,b1,c1,0,'aplbuenovdd',1);
            %w(i) =  Feval([t(i);1.2],A2,b2,c2,0,'aplbuenovdd',2);
        end
    end
    surf(gridx,gridy,z)
    %surf(gridx,gridy,w)
    %gridx=linspace(-20,20);
    %gridy=linspace(-20,20);
    zvec=xvec;
    for i=1:length(xvec)
      zvec(i) = Feval([xvec(i);1],A1,b1,c1,0,'aplbuenovdd',1);
    end
    %scatter(xvec,zvec,'green');

    %plot(t,z)
    %plot(t,w)
end

if (n==22)
    hold on

    [gridx,gridy]=meshgrid(-5:0.05:5,-5:0.05:5);
    [s1,s2]=size(gridy);
    z=zeros(s1,s2); w=zeros(s1,s2);
    %s1=100; s2=100;
    %z=zeros(100,1); w=zeros(100,1); t=linspace(-20,20);
    
    for i=1:s1
        for j=1:s2
            z(i,j) = Feval([gridy(i);gridy(j);0;0],A1,b1,c1,0,mystr,1,r1,ro1);
            %z(i,j) = Feval([gridx(i)/sqrt(2)-gridy(i)/sqrt(2),-gridx(j)/sqrt(2)+gridy(j)/sqrt(2)],A1,b1,c1,0,mystr,1,r1,ro1);
            w(i,j) = Feval([0;0;gridy(i);gridy(j)],A2,b2,c2,0,mystr,2,r2,ro2);
            %w(i,j) = Feval([gridx(i)/sqrt(2)-gridy(i)/sqrt(2),-gridx(j)/sqrt(2)+gridy(j)/sqrt(2)],A2,b2,c2,0,mystr,2,r2,ro2);
            %tes(i,j) = gridy(i)^2 + gridy(j)^2;
            %z(i) =  Feval([t(i);1],A1,b1,c1,0,'aplbuenovdd',1);
            %w(i) =  Feval([t(i);1.2],A2,b2,c2,0,'aplbuenovdd',2);
        end
    end
    %surf(gridx,gridy,z)
    %surf(gridx,gridy,w)
    %gridx=linspace(-20,20);
    %gridy=linspace(-20,20);
    %zvec=xvec;
    %for i=1:length(xvec)
    %  zvec(i) = Feval([xvec(i);1],A1,b1,c1,0,'aplbuenovdd',1);
    %end
    %scatter(xvec,zvec,'green');

    %plot(t,z)
    %plot(t,w)
    
    teste1 = Feval(x0,A1,b1,c1,2, mystr,1,r1,ro1);
                teste2 = Feval(x0,A2,b2,c2,2, mystr,2,r2,ro2);
                teste=[min(eigs(teste1)),min(eigs(teste2))];
                disp(teste);
    
end



%end

%fazendo compara jacobi pra ver como vai
if strcmp(mystr,'quadratica')||strcmp(mystr,'testepaper')||strcmp(mystr,'problemayuan')||strcmp(mystr,'quadtestes1')||...
        strcmp(mystr,'cubico')||strcmp(mystr,'quadtestes5')||strcmp(mystr,'quadratica_yuan')||strcmp(mystr,'vacina')||...
        strcmp(mystr,'exartigo1_mbom')||strcmp(mystr,'exartigo2_mruim')||strcmp(mystr,'exartigo3_indef')||...
        strcmp(mystr,'exartigo4_cubic')||strcmp(mystr,'formation')||(strcmp(mystr,'aplbuenovdd'))%||(strcmp(mystr,'exnovo6'))
    [xjacobi,kjacobi,tjacobi,ejacobi]=ComparaJacobi(x0jacobi,A1,b1,c1,A2,b2,c2,tol,mystr);
end
if strcmp(mystr,'exartigo1_mbom')||(strcmp(mystr,'exnovo6'))
    [xnika,knika,erronika] = nikaidoisoda(x0jacobi,tol,mystr);
end
%fimplicit(fun1,[-100 100]);
%fimplicit(fun2,[-100 100]);
%disp('meut');

%scatter
%vecbom1 = [0.60385,0.7971;-0.88771,-0.4604];
%vecbom2 = [-0.7971,-0.60385;0.4604,0.88771];
%vecbom1 = [0.60385; -0.88771; -0.7971; 0.4604];
%vecbom2 = [0.7971; -0.4604; -0.60385; 0.88771];
%scatter(vecbom1,vecbom2,'filled');

%vecruim1 = [-0.24504,-0.50284,-0.104,-0.48123,0.087088,-0.12219,-0.0657]';
%vecruim2 = [0.24504,0.41887,0.12047,0.53458,0.51903,0.10589,-3.9025]';
%scatter(vecruim1,vecruim2,'s');

if (n==1)&&(m==1)&&(vaigrafico==1)
    %contour %(meio caro computacionalmente...)
    %vecgridx = linspace(x0(1)-10,x0(1)+10);
    %vecgridy = linspace(x0(2)-10,x0(2)+10);
    %matrixval1 = zeros(100,100);
    %matrixval2 = zeros(100,100);
    %for i=1:100
    %    for j=1:100
    %        matrixval1(i,j) = Feval([vecgridx(i);vecgridy(j)],A1,b1,c1,1,mystr,1);
    %matrixval2(i,j) = Feval([vecgridx(i);vecgridy(j)],A2,b2,c2,1,mystr,2);
    %    end
    %end
    %contour(vecgridx,vecgridy,matrixval1);
    %contour(vecgridx,vecgridy,matrixval2);
    %hold off
end

if (n==1)&&(m==1)&&(strcmp(mystr,'aplbuenovdd'))
    t = linspace(-3,3); ty=-0; tx=0;
    y=t;
    %z=t;
    scatvecx=xvec
    scatvecy=yvec
    
    if strcmp(flaghess,'tfixo1')
        x0newton=x0;
    end
    
    for i=1:100
        y(i) =  Feval([t(i)-tx;x0(2)-tx],A1,b1,c1,0,'aplbuenovdd',1);
        %z(i) =  Feval([t(i)-tx;x0newton(2)-tx],A1,b1,c1,0,'aplbuenovdd',1);
        z(i) =  Feval([x0(1)-ty;t(i)-tx],A2,b2,c2,0,'aplbuenovdd',2);
    end
    for i=1:length(xvec)
        scatvecx(i) = Feval([xvec(i)-tx;x0(2)-tx],A1,b1,c1,0,'aplbuenovdd',1);
        scatvecy(i) = Feval([x0(1)-ty;yvec(i)-ty],A2,b2,c2,0,'aplbuenovdd',2);
    end
    meut
    %if strcmp(flaghess,'tfixo1')
    %    for i=1:length(xvec)
    %    scatvecy(i) = Feval([ty-x0(1);ty-yvec(i)],A1,b1,c1,0,'aplbuenovdd',1);
    %    end
    %    for i=1:100
    %    z(i) =  Feval([x0(1)-ty;t(i)-tx],A2,b2,c2,0,'aplbuenovdd',2);
    %    end
    %end
    %for i=1:length(xvecnewton)
    %scatvecy(i) = Feval([xvecnewton(i)-tx;x0newton(2)-tx],A1,b1,c1,0,'aplbuenovdd',1);
    %end
    
    plot(t,y)
    plot(t,z)
    hold on
    %plot(t,z)
    %size(xvec)
    %size(xvecnewton)
    
    scatter(xvec,scatvecx);
    scatter(yvec,scatvecy);
    
    teste1 = Feval(x0,A1,b1,c1,2, mystr,1,r1,ro1);
                teste2 = Feval(x0,A2,b2,c2,2, mystr,2,r2,ro2);
                teste=[eigs(teste1),eigs(teste2)];
                disp(teste);
    
end

%pra fazer o corte
if (n==2)&&(m==2)&&(strcmp(mystr,'aplbuenovdd'))
    t = linspace(-3,3); ty=-0; tx=0;
    y=t;
    %z=t;
    scatvecx=xvec
    scatvecy=yvec
    
    if strcmp(flaghess,'tfixo1')
        x0newton=x0;
    end
    
    for i=1:100
        y(i) =  Feval([t(i)-tx;0;x0(3)-tx;0],A1,b1,c1,0,'aplbuenovdd',1);
        %z(i) =  Feval([t(i)-tx;x0newton(2)-tx],A1,b1,c1,0,'aplbuenovdd',1);
        z(i) =  Feval([x0(1)-ty;0;t(i)-tx;0],A2,b2,c2,0,'aplbuenovdd',2);
    end
    for i=1:length(xvec)
        scatvecx(i) = Feval([xvec(i)-tx;0;x0(3)-tx;0],A1,b1,c1,0,'aplbuenovdd',1);
        scatvecy(i) = Feval([x0(1)-ty;0;yvec(i)-ty;0],A2,b2,c2,0,'aplbuenovdd',2);
    end
    %for i=1:length(xvecnewton)
    %scatvecy(i) = Feval([xvecnewton(i)-tx;x0newton(2)-tx],A1,b1,c1,0,'aplbuenovdd',1);
    %end
    
    plot(t,y)
    plot(t,z)
    hold on
    %plot(t,-z)
    %size(xvec)
    %size(xvecnewton)
    
    scatter(xvec,scatvecx);
    scatter(yvec,scatvecy);
    
    teste1 = Feval(x0,A1,b1,c1,2, mystr,1,r1,ro1);
                teste2 = Feval(x0,A2,b2,c2,2, mystr,2,r2,ro2);
                teste=[eigs(teste1),eigs(teste2)];
                disp(teste);
    
end

%descomente esse pra fazer o grafico 2 por 2 depois
if (strcmp(mystr,'aplbuenovdd'))&&(m==11)
    scatponto=[x0(1),x0(3);x0(2),x0(4)];
    %scatter(x0(1),x0(2),'o','MarkerFaceColor','red');
    %scatter(x0(3),x0(4),'o','MarkerFaceColor','red');
    scatter(scatponto(1,:),scatponto(2,:),'o','MarkerFaceColor','red');
    scatcruz=[1,0,0,-1;0,-1,1,0];
    scatter(scatcruz(1,:),scatcruz(2,:),'+','magenta');
    %scatter(0,-1,'+','magenta');
    %scatter(0,1,'+','magenta');
    %scatter(-1,0,'+','magenta');
end


if (n==1)&&(m==1)&&(vaimontadados==1)&&(vaigrafico==1)
    %   hold off
end
