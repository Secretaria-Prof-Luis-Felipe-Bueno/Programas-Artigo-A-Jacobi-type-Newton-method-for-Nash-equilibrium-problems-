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
mystr = 'quadfacil';

%mystr = 'exartigo1_mbom';
%mystr = 'exartigo2_mruim';
%mystr = 'exartigo3_indef';
%mystr = 'exartigo4_cubic';
%mystr = 'exartigo5_futebol';
%mystr = 'futebolzero';
%mystr = 'exartigo6_vacina';

 r1 = 1; r2=1;
 ro1 = 1; ro2=1;

treset=1;
maint=0;
fazcomlambdas=1;
tic

fazcholmod=0;

vaidisplay = 0;
vaigrafico=1;
vaimontadados = 1;

if vaimontadados==1 
   alfa = 1e-2; %alfa em (0,1)
   tol = 1e-4; kmax=100;
   minvaltol = 1e-6;
   fixvaltol=1;
   maxvaltol=1e10;
   zerodograd=1e-4;
   %zerodograd=0;
   lambdai=0;
   n = 1;
   m = 1;
  [A1,A2,b1,b2,c1,c2,x0,rest1,rest2,bdarest1,bdarest2] = montadados(mystr,n,m);
end

 gamma = 1;
 gammamax = gamma;

grad1 = Feval(x0,A1,b1,c1,1,mystr,1,r1,ro1);
grad2 = Feval(x0,A2,b2,c2,1,mystr,2,r2,ro2);
gradlam1 = bdarest1 - rest1*x0(1:n);
gradlam2 = bdarest2 - rest2*x0(n+1:n+m);
%erro = norm(grad1)+norm(grad2);

srstring='sr';


k=0; t=1; tau=0.99;
maint=toc;
aux=1;
%x0jacobi = x0;
while((aux==1)&&(k<kmax))
    
    tic
       
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
    
    [sizer1n,sizer1m] = size(rest1');
    [sizer2n,sizer2m] = size(rest2');
      
    
    mista1 = Feval(x0,A1,b1,c1,3,mystr,1,r1,ro1);
    mista2 = Feval(x0,A2,b2,c2,3,mystr,2,r2,ro2);
    bigmatrix = [H1 rest1 t*mista1  zeros(n,sizer2n); rest1' zeros(sizer1n,sizer1n) zeros(sizer1n,m) zeros(sizer1n,sizer2n);  t*mista2 zeros(m,sizer1n) H2 rest2; zeros(sizer2n,n) zeros(sizer2n,sizer1n) rest2' zeros(sizer2n,sizer2n) ];
    
    biggrad = [grad1;gradlam1';grad2;gradlam2'];
    
        
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
            
            [sizeeyen,sizeeyem]=size(bigmatrix);
            bigmatrix=bigmatrix+(lambdai*abs(minval)+increment*fixvaltol)*eye(sizeeyen,sizeeyem);
            %bigmatrix = (abs(minval)+minvaltol)*eye(n+m,n+m);
            
            if (strcmp(mystr,'problemayuan')==0)&&(strcmp(mystr,'vacina')==0)&&(strcmp(mystr,'sin2')==0)
               minval = eigs(bigmatrix,1,srstring);
            else
               minval = min(eig(bigmatrix));
            end
            %minval = min(eig(bigmatrix));
            increment=2*increment;
            
        end
        
        else %ai cholmod nao eh zero
            bigmatrix = cholmod(bigmatrix,1e-6);
            minval = min(eig(bigmatrix));
        end %end cholmod
        
    end
        
    if (strcmp(flaghess,'jacobiintern')==0 )&&(strcmp(flaghess,'jacobidois')==0 )
       
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
            newd(n+1+sizer1n:n+m+sizer1n)=zeros(m,1);
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
            maint=maint+toc;
        
    else  %ai strcmp(flaghess,'jacobidois')
            newd = x0;
        
            newd(1:n) = H1\(-biggrad(1:n)-mista1*x0(n+1:n+m));
            newd(n+1+sizer1n:n+m+sizer1n) = H2\(-biggrad(n+1+sizer1n:n+m+sizer1n)-mista2*x0(1:n));

            x0 = newd;
            grad1 = Feval(x0,A1,b1,c1,1,mystr,1,r1,ro1);
            grad2 = Feval(x0,A2,b2,c2,1,mystr,2,r2,ro2);
            erro = norm(grad1)+norm(grad2);
            k = k+1;
            maint=maint+toc;
    end %end os 3 strcmps
    
    gradval1 = Feval([x0(1:n);x0(n+1:n+m)+t*newd(n+1+sizer1n:n+m+sizer1n)],A1,b1,c1,0,mystr,1,r1,ro1);
    gradval2 = Feval([x0(1:n)+t*newd(1:n);x0(n+1:n+m)],A2,b2,c2,0,mystr,2,r2,ro2);
    produto1 = -(newd(1:n))'*H1*newd(1:n);
    produto2 = -(newd(n+1+sizer1n:n+m+sizer1n))'*H2*newd(n+1+sizer1n:n+m+sizer1n);
    %direito1 = -tetta*norm(gradval1)*norm(newd(1:n));
    %direito2 = -tetta*norm(gradval2)*norm(newd(n+1:n+m));
    newds = [newd(1:n);newd(n+1+sizer1n:n+m+sizer1n)];
    numero1 = Feval(x0+t*newds,A1,b1,c1,0,mystr,1,r1,ro1)-Feval([x0(1:n);x0(n+1:n+m)+t*newd(n+1+sizer1n:n+m+sizer1n)],A1,b1,c1,0,mystr,1,r1,ro1)-alfa*t*produto1;
    numero2 = Feval(x0+t*newds,A2,b2,c2,0,mystr,2,r2,ro2)-Feval([x0(1:n)+t*newd(1:n);x0(n+1:n+m)],A2,b2,c2,0,mystr,2,r2,ro2)-alfa*t*produto2;
    %if strcmp(flaghess,'quadrado')
    %    numero1 = Feval(x0+t*newd,A1,b1,c1,0,mystr,1,r1,ro1)-Feval([x0(1:n);x0(n+1:n+m)+t*newd(n+1:n+m)],A1,b1,c1,0,mystr,1,r1,ro1)+alfa*t*(newd(1:n))'*newd(1:n);
    %    numero2 = Feval(x0+t*newd,A2,b2,c2,0,mystr,2,r2,ro2)-Feval([x0(1:n)+t*newd(1:n);x0(n+1:n+m)],A2,b2,c2,0,mystr,2,r2,ro2)+alfa*t*(newd(n+1:n+m))'*newd(n+1:n+m);
    %end
    
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
                
       if (  (  (real(minval)<minvaltol)||(...
               (numero1>0)&&(numero2>zerodograd)   )||(  (numero2>0)&&(numero1>zerodograd) )||(...
               (produto1>0)&&(numero2>zerodograd)&&(fazcomlambdas==0)   )||(  (produto2>0)&&(numero1>zerodograd)&&(fazcomlambdas==0) )||(...
               (numero2>zerodograd)&&(fazcomlambdas~=0)   )||(...
               numero1>zerodograd)&&(fazcomlambdas~=0) )||(  (produto1>0)&&(numero2>zerodograd)&&(fazcomlambdas==0)   )    )
            

       
        
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
            while (  (  (real(minval)<minvaltol)||(...
               (numero1>0)&&(numero2>zerodograd)   )||(  (numero2>0)&&(numero1>zerodograd) )||(...
               (produto1>0)&&(numero2>zerodograd)&&(fazcomlambdas==0)   )||(  (produto2>0)&&(numero1>zerodograd)&&(fazcomlambdas==0) )||(...
               (numero2>zerodograd)&&(fazcomlambdas~=0)   )||(...
               numero1>zerodograd)&&(fazcomlambdas~=0) )||(  (produto1>0)&&(numero2>zerodograd)&&(fazcomlambdas==0)   )  )
               
               

                t=t/2;
                %disp('recusou o passo');
                bigmatrix = [H1 rest1 t*mista1  zeros(n,sizer2n); rest1' zeros(sizer1n,sizer1n) zeros(sizer1n,m) zeros(sizer1n,sizer2n);  t*mista2 zeros(m,sizer1n) H2 rest2; zeros(sizer2n,n) zeros(sizer2n,sizer1n) rest2' zeros(sizer2n,sizer2n) ];
                if (strcmp(mystr,'problemayuan')==0)&&(strcmp(mystr,'vacina')==0)&&(strcmp(mystr,'sin2')==0)
                 minval = eigs(bigmatrix,1,srstring);
                else
                 minval = min(eig(bigmatrix));
                end
                increment=1e-4;
                
                if (real(minval)<minvaltol)
                    if fazcholmod==0
                    
                while (real(minval) < minvaltol)
                    [sizeeyen,sizeeyem] = size(bigmatrix);
                    bigmatrix=bigmatrix+(lambdai*abs(minval)+(increment*fixvaltol))*eye(sizeeyen,sizeeyem);
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
                lamb1 = newd(n+1:n+sizer1n);
                lamb2 = newd(n+m+sizer1n+1:end);
                
                gradval1 = Feval([x0(1:n);x0(n+1:n+m)+t*newd(n+1+sizer1n:n+m+sizer1n)],A1,b1,c1,0,mystr,1,r1,ro1);
                gradval2 = Feval([x0(1:n)+t*newd(1:n);x0(n+1:n+m)],A2,b2,c2,0,mystr,2,r2,ro2);
                produto1 = -(newd(1:n))'*H1*newd(1:n);
                produto2 = -(newd(n+1:n+m))'*H2*newd(n+1+sizer1n:n+m+sizer1n);
                newds = [newd(1:n);newd(n+1+sizer1n:n+m+sizer1n)];
                numero1 = Feval(x0+t*newds,A1,b1,c1,0,mystr,1,r1,ro1)-Feval([x0(1:n);x0(n+1:n+m)+t*newd(n+1+sizer1n:n+m+sizer1n)],A1,b1,c1,0,mystr,1,r1,ro1)-alfa*t*produto1;
                numero2 = Feval(x0+t*newds,A2,b2,c2,0,mystr,2,r2,ro2)-Feval([x0(1:n)+t*newd(1:n);x0(n+1:n+m)],A2,b2,c2,0,mystr,2,r2,ro2)-alfa*t*produto2;
                %if strcmp(flaghess,'quadrado')
                %    numero1 = Feval(x0+t*newd,A1,b1,c1,0,mystr,1,r1,ro1)-Feval([x0(1:n);x0(n+1:n+m)+t*newd(n+1:n+m)],A1,b1,c1,0,mystr,1,r1,ro1)+alfa*t*(newd(1:n))'*newd(1:n);
                %    numero2 = Feval(x0+t*newd,A2,b2,c2,0,mystr,2,r2,ro2)-Feval([x0(1:n)+t*newd(1:n);x0(n+1:n+m)],A2,b2,c2,0,mystr,2,r2,ro2)+alfa*t*(newd(n+1:n+m))'*newd(n+1:n+m);
                %end
                %coloquei esse artificial pra caso a norma dos grad for zero
            end
            
        end %treset
        maint = maint+toc;
        
    else %nesse caso aceita o passo
        if passa == 1
            
               %if (norm(grad1)<=zerodograd)&&(t<tau)
               %    numero1=0;
               %end
               %if (norm(grad2)<=zerodograd)&&(t<tau)
               %    numero2=0;
               %end
               
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
                 newd(n+1+sizer1n:n+m+sizer1n)=zeros(m,1);
                 minlam1 = eigs(H1,1,srstring);
                 while minlam1<=0
                     H1=H1+eye(m,m);
                     minlam1 = eigs(H1,1,srstring);
                 end
                 newd(1:n) = -H1\grad1;
               end
            
               newds = [newd(1:n);newd(n+1+sizer1n:n+m+sizer1n)];
            x0 = x0+t*newds;

    
            grad1 = Feval(x0,A1,b1,c1,1,mystr,1,r1,ro1);
            grad2 = Feval(x0,A2,b2,c2,1,mystr,2,r2,ro2);
            erro = norm(grad1+rest1*lamb1 + H1*newd(1:n))+norm(grad2+rest2*lamb2 + H2*newd(n+1:n+m));
            %disp('t: '); disp(t);
            %disp('[x,y]: '); disp(x0);
            %disp('d: '); disp(newd);
            %pause
            k = k+1;
            
            if erro<=tol
                aux=0; 
            end
            
            %update params
            lab = max(norm(lamb1),norm(lamb2));
            if gamma>=1.1*(lab+gammamax)
               gamma = (gamma + gammamax + lab)/2; 
            else
                if gamma < lab + gammamax
                    gamma = max(1.5*gamma,lab+gammamax);
                end
                
            end
          
                       
            maint=maint+toc;
            
            if ((n==1)&&(m==1)&&(vaimontadados==1)&&(vaigrafico==1))&&(strcmp(mystr,'quadfacil')==0)
                xvec = [xvec; x0(1)];
                yvec = [yvec; x0(2)];

                %scatter(xvec,yvec);
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
                xvec = [xvec; x0(1)];
                yvec = [yvec; x0(2)];
                
                %hold on
                xvec2 = [xvec2; x0(3)];
                yvec2 = [yvec2; x0(4)];
                %pause
                
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
        x0 = x0+t*newd;
        
        xvecnewton = [xvecnewton; x0(1)];
        yvecnewton = [yvecnewton; x0(2)];
        %x0(2)=1-x0(1);
        %x0(4)=1-x0(3);
        %x0;
        
        
        xvec = [xvec; x0(1)];
        yvec = [yvec; x0(2)];
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
   
    
    
end %end while


