
flaghess = 'ifdois'; multhess = 1;
%flaghess = 'tfixo1';

mystr = 'quadfacil'; %4.1do paper
%mystr = 'quadmlau';
%mystr = 'quadfacilmista'; %ex2 do paper
mystr = 'tenta1'; %4.3 do paper
%mystr='quadhessmista';
%mystr = 'dozero';
mystr = 'exmisto'; %4.4 do paper
%mystr = 'naonosso';
%mystr = 'aplbuenovdd';
%mystr = 'aplbuenovddiq';
%mystr = 'quadtranspo';
mystr = 'apalfa1';
%mystr='quadfacillau'

treset=1;
fazcomlambdas=0;
fazzerograd=1;
paramlamb = 100000;
tic

fazcholmod=0;

vaidisplay = 0;
vaigrafico=1;

alfa = 1e-6; %alfa em (0,1)
tol = 1e-4; kmax=1000;
minvaltol = 1e-10;
fixvaltol=1e-4;
maxvaltol=1e10;
cholmodtol=1e-4;
passaassim=0;
%zerodograd=1e-10;
zerodograd=0;
lambdai=0;
n = 2;
m = 2;

%tem q definir vaimontadados e os ris fora
if vaimontadados==1
    r1 = 1; r2=1;
    ro1 = 1; ro2=1;
    [A1,A2,b1,b2,c1,c2,x0,rest1,rest2,bdarest1,bdarest2] = montadados(mystr,n,m);
end

[y1,y2] = gradrest(x0,mystr,1,bdarest1,bdarest2,rest1,rest2);
rest1 = y1; rest2=y2;

gamma = 0.1;
gammamax = 10000;

grad1 = Feval(x0,A1,b1,c1,1,mystr,1,r1,ro1);
grad2 = Feval(x0,A2,b2,c2,1,mystr,2,r2,ro2);
[y1,y2] = gradrest(x0,mystr,0,bdarest1,bdarest2,rest1,rest2);
gradlam1 = y1; gradlam2=y2;

srstring='sr';
%if fazcholmod~=0
%   srstring='sr';
%end

k=0; t=1; tau=1;
maint=toc;
aux=1;

[sizer1n,sizer1m] = size(rest1);
[sizer2n,sizer2m] = size(rest2);

lamb1 = zeros(sizer1n,1);
lamb2 = zeros(sizer2n,1);

lamb1velho = lamb1;
lamb2velho = lamb2;

xvec=[]; yvec=[];
xvec2=[]; yvec2=[];


while((aux==1)&&(k<kmax))
    
    tic
    
    H1 = Feval(x0,A1,b1,c1,2,mystr,1,r1,ro1);
    H2 = Feval(x0,A2,b2,c2,2,mystr,2,r2,ro2);
    [y1,y2] = gradrest(x0,mystr,2,bdarest1,bdarest2,rest1,rest2,H1,H2,lamb1,lamb2);
    H1 = y1; H2=y2;
    
    mista1 = Feval(x0,A1,b1,c1,3,mystr,1,r1,ro1);
    mista2 = Feval(x0,A2,b2,c2,3,mystr,2,r2,ro2);
    if strcmp(mystr,'exmisto')
        mista1=mista1+lamb1*gradrest(x0,mystr,3,bdarest1,bdarest2,rest1,rest2,H1,H2,lamb1,lamb2);
        mista2=mista2+lamb2*gradrest(x0,mystr,3,bdarest1,bdarest2,rest1,rest2,H1,H2,lamb1,lamb2);
    end
    if strcmp(mystr,'apalfa1')==0
      bigmatrix = [H1 rest1' t*mista1  zeros(n,sizer2n); rest1 zeros(sizer1n,sizer1n) zeros(sizer1n,m) zeros(sizer1n,sizer2n);  t*mista2 zeros(m,sizer1n) H2 rest2'; zeros(sizer2n,n) zeros(sizer2n,sizer1n) rest2 zeros(sizer2n,sizer2n) ];
    end
    if strcmp(mystr,'apalfa1')
        a1 = Feval(x0,A1,b1,c1,3,mystr,1,r1,ro1);
        a1=a1(3,:);
        a2 = Feval(x0,A2,b2,c2,3,mystr,2,r2,ro2);
        a2=a2(3,:);
        bigmatrix = [H1 rest1' t*mista1  zeros(n+1,sizer2n); rest1 zeros(sizer1n,sizer1n) a1 zeros(sizer1n,sizer2n);  t*mista2 zeros(m+1,sizer1n) H2 rest2'; a2 zeros(sizer2n,sizer1n) rest2 zeros(sizer2n,sizer2n) ];
    end
    
    biggrad = [grad1+rest1'*lamb1;gradlam1;grad2+rest2'*lamb2;gradlam2];
    
    if (strcmp(flaghess,'ifdois')||strcmp(flaghess,'quadrado') )%||strcmp(flaghess,'tfixo1')
        
        minval1 = eigs(H1,1,srstring);
        minval2 = eigs(H2,1,srstring);
        minval = min(minval1,minval2);
        
        increment=1;
        while(real(minval)<minvaltol)
            
            [sizeeyen,sizeeyem]=size(bigmatrix);
            if fazcholmod==0
                H1 = H1+(lambdai*abs(minval)+increment*fixvaltol)*eye(size(H1));
                H2 = H2+(lambdai*abs(minval)+increment*fixvaltol)*eye(size(H2));
            else
                H1 = cholmod(H1,cholmodtol);
                H2 = cholmod(H2,cholmodtol);
            end
            if strcmp(mystr,'apalfa1')==0
                bigmatrix=[H1 rest1' t*mista1  zeros(n,sizer2n); rest1 zeros(sizer1n,sizer1n) zeros(sizer1n,m) zeros(sizer1n,sizer2n);  t*mista2 zeros(m,sizer1n) H2 rest2'; zeros(sizer2n,n) zeros(sizer2n,sizer1n) rest2 zeros(sizer2n,sizer2n) ];
            end
            if strcmp(mystr,'apalfa1')
                a1 = Feval(x0,A1,b1,c1,3,mystr,1,r1,ro1);
                a1=a1(3,:);
                a2 = Feval(x0,A2,b2,c2,3,mystr,2,r2,ro2);
                a2=a2(3,:);
                bigmatrix = [H1 rest1' t*mista1  zeros(n+1,sizer2n); rest1 zeros(sizer1n,sizer1n) a1 zeros(sizer1n,sizer2n);  t*mista2 zeros(m+1,sizer1n) H2 rest2'; a2 zeros(sizer2n,sizer1n) rest2 zeros(sizer2n,sizer2n) ];
            end
            %bigmatrix = (abs(minval)+minvaltol)*eye(n+m,n+m);
            
            minval1 = eigs(H1,1,srstring);
            minval2 = eigs(H2,1,srstring);
            minval = min(minval1,minval2);
            %minval = min(eig(bigmatrix));
            increment=2*increment;
            
        end
        
        if (norm([grad1+rest1'*lamb1;gradlam1])<=zerodograd)&&(fazzerograd==1)&&(strcmp(mystr,'ifdois'))
            if strcmp(mystr,'apalfa1')==0
                bigmatrix = [H1 rest1' zeros(size(mista1))  zeros(n,sizer2n); rest1 zeros(sizer1n,sizer1n) zeros(sizer1n,m) zeros(sizer1n,sizer2n);  t*mista2 zeros(m,sizer1n) H2 rest2'; zeros(sizer2n,n) zeros(sizer2n,sizer1n) rest2 zeros(sizer2n,sizer2n) ];
            end
            if strcmp(mystr,'apalfa1')
                a1 = Feval(x0,A1,b1,c1,3,mystr,1,r1,ro1);
                a1=a1(3,:);
                a2 = Feval(x0,A2,b2,c2,3,mystr,2,r2,ro2);
                a2=a2(3,:);
                bigmatrix = [H1 rest1' t*mista1  zeros(n+1,sizer2n); rest1 zeros(sizer1n,sizer1n) a1 zeros(sizer1n,sizer2n);  t*mista2 zeros(m+1,sizer1n) H2 rest2'; a2 zeros(sizer2n,sizer1n) rest2 zeros(sizer2n,sizer2n) ];
            end
        end
        if (norm([grad2+rest2'*lamb2;gradlam2])<=zerodograd)&&(fazzerograd==1)&&(strcmp(mystr,'ifdois'))
            if strcmp(mystr,'apalfa1')==0
                bigmatrix = [H1 rest1' t*mista1  zeros(n,sizer2n); rest1 zeros(sizer1n,sizer1n) zeros(sizer1n,m) zeros(sizer1n,sizer2n);  zeros(size(mista2)) zeros(m,sizer1n) H2 rest2'; zeros(sizer2n,n) zeros(sizer2n,sizer1n) rest2 zeros(sizer2n,sizer2n) ];
            end
            if strcmp(mystr,'apalfa1')
                a1 = Feval(x0,A1,b1,c1,3,mystr,1,r1,ro1);
                a1=a1(3,:);
                a2 = Feval(x0,A2,b2,c2,3,mystr,2,r2,ro2);
                a2=a2(3,:);
                bigmatrix = [H1 rest1' t*mista1  zeros(n+1,sizer2n); rest1 zeros(sizer1n,sizer1n) a1 zeros(sizer1n,sizer2n);  t*mista2 zeros(m+1,sizer1n) H2 rest2'; a2 zeros(sizer2n,sizer1n) rest2 zeros(sizer2n,sizer2n) ];
            end
        end
        
    end %%ifdois ou quadrado
    
    newd = -bigmatrix\biggrad;
    
    if strcmp(mystr,'apalfa1') n=3; m=3; end
    
    gradval1 = Feval([x0(1:n);x0(n+1:n+m)+t*newd(n+1+sizer1n:n+m+sizer1n)],A1,b1,c1,0,mystr,1,r1,ro1);
    gradval2 = Feval([x0(1:n)+t*newd(1:n);x0(n+1:n+m)],A2,b2,c2,0,mystr,2,r2,ro2);
    
    %esses produtos sao velhos, vamos pegar um multiplo da norma do
    %gradiente
    %produto1 = -(newd(1:n))'*H1*newd(1:n);
    %produto2 = -(newd(n+1+sizer1n:n+m+sizer1n))'*H2*newd(n+1+sizer1n:n+m+sizer1n);
    produto1 = (1/paramlamb)*norm([grad1+rest1'*lamb1;gradlam1]);
    produto2 = (1/paramlamb)*norm([grad2+rest2'*lamb2;gradlam2]);
    
    newds = [newd(1:n);newd(n+1+sizer1n:n+m+sizer1n)];
    
    x0agr = x0+t*newds;
    [y1,y2] = gradrest(x0agr,mystr,0,bdarest1,bdarest2,rest1,rest2,H1,H2,lamb1,lamb2);
    hnovo1 = y1; hnovo2 = y2;
    
    numero1 = Feval(x0+t*newds,A1,b1,c1,0,mystr,1,r1,ro1)-Feval([x0(1:n);x0(n+1:n+m)+t*newd(n+1+sizer1n:n+m+sizer1n)],A1,b1,c1,0,mystr,1,r1,ro1)+gamma*(norm(hnovo1) - norm(gradlam1))-alfa*t*produto1;
    numero2 = Feval(x0+t*newds,A2,b2,c2,0,mystr,2,r2,ro2)-Feval([x0(1:n)+t*newd(1:n);x0(n+1:n+m)],A2,b2,c2,0,mystr,2,r2,ro2)+gamma*(norm(hnovo2) - norm(gradlam2))-alfa*t*produto2;
    
    if strcmp(mystr,'apalfa1') n=2; m=2; end
    
    if (strcmp(flaghess,'tfixo1')==0)&&(strcmp(flaghess,'newtonpuro')==0)&&( strcmp(flaghess,'jacobiintern')==0 )&&( strcmp(flaghess,'jacobidois')==0 )
        
        if (  (det(bigmatrix)<minvaltol)||(real(minval)<minvaltol)||(  (( numero1>zerodograd )||( numero2>zerodograd ))&&(passaassim==0)  )  )
            %aqui nao aceita o passo, apenas diminui o t
            
            while  (  (det(bigmatrix)<minvaltol)||(real(minval)<minvaltol)||(  (( numero1>zerodograd )||( numero2>zerodograd ))&&(passaassim==0)  )  )
                
                %disp('recusou o passo');
                
                minval1 = eigs(H1,1,srstring);
                minval2 = eigs(H2,1,srstring);
                minval = min(minval1,minval2);
                
                increment=1;
                while(real(minval)<minvaltol)
                    %disp('bugou o minval');
                    [sizeeyen,sizeeyem]=size(bigmatrix);
                    if fazcholmod==0
                        H1 = H1+(lambdai*abs(minval)+increment*fixvaltol)*eye(size(H1));
                        H2 = H2+(lambdai*abs(minval)+increment*fixvaltol)*eye(size(H2));
                    else
                        H1 = cholmod(H1,cholmodtol);
                        H2 = cholmod(H2,cholmodtol);
                    end
                    if strcmp(mystr,'apalfa1')==0
                        bigmatrix=[H1 rest1' t*mista1  zeros(n,sizer2n); rest1 zeros(sizer1n,sizer1n) zeros(sizer1n,m) zeros(sizer1n,sizer2n);  t*mista2 zeros(m,sizer1n) H2 rest2'; zeros(sizer2n,n) zeros(sizer2n,sizer1n) rest2 zeros(sizer2n,sizer2n) ];
                    end
                    if strcmp(mystr,'apalfa1')
                        a1 = Feval(x0,A1,b1,c1,3,mystr,1,r1,ro1);
                        a1=a1(3,:);
                        a2 = Feval(x0,A2,b2,c2,3,mystr,2,r2,ro2);
                        a2=a2(3,:);
                        bigmatrix = [H1 rest1' t*mista1  zeros(n+1,sizer2n); rest1 zeros(sizer1n,sizer1n) a1 zeros(sizer1n,sizer2n);  t*mista2 zeros(m+1,sizer1n) H2 rest2'; a2 zeros(sizer2n,sizer1n) rest2 zeros(sizer2n,sizer2n) ];
                    end
                    %bigmatrix = (abs(minval)+minvaltol)*eye(n+m,n+m);
                    
                    minval1 = eigs(H1,1,srstring);
                    minval2 = eigs(H2,1,srstring);
                    minval = min(minval1,minval2);
                    %minval = min(eig(bigmatrix));
                    increment=2*increment;
                    
                end
                
                if (det(bigmatrix)<minvaltol)
                    while  (det(bigmatrix)<minvaltol)
                        t=t/2;
                        if strcmp(mystr,'apalfa1')==0
                            bigmatrix = [H1 rest1' t*mista1  zeros(n,sizer2n); rest1 zeros(sizer1n,sizer1n) zeros(sizer1n,m) zeros(sizer1n,sizer2n);  t*mista2 zeros(m,sizer1n) H2 rest2'; zeros(sizer2n,n) zeros(sizer2n,sizer1n) rest2 zeros(sizer2n,sizer2n) ];
                        end
                        if strcmp(mystr,'apalfa1')
                            a1 = Feval(x0,A1,b1,c1,3,mystr,1,r1,ro1);
                            a1=a1(3,:);
                            a2 = Feval(x0,A2,b2,c2,3,mystr,2,r2,ro2);
                            a2=a2(3,:);
                            bigmatrix = [H1 rest1' t*mista1  zeros(n+1,sizer2n); rest1 zeros(sizer1n,sizer1n) a1 zeros(sizer1n,sizer2n);  t*mista2 zeros(m+1,sizer1n) H2 rest2'; a2 zeros(sizer2n,sizer1n) rest2 zeros(sizer2n,sizer2n) ];
                        end
                    end
                end
                
                if (norm([grad1+rest1'*lamb1;gradlam1])<=zerodograd)&&(fazzerograd==1)&&(strcmp(mystr,'ifdois'))
                    if strcmp(mystr,'apalfa1')==0
                        bigmatrix = [H1 rest1' zeros(size(mista1))  zeros(n,sizer2n); rest1 zeros(sizer1n,sizer1n) zeros(sizer1n,m) zeros(sizer1n,sizer2n);  t*mista2 zeros(m,sizer1n) H2 rest2'; zeros(sizer2n,n) zeros(sizer2n,sizer1n) rest2 zeros(sizer2n,sizer2n) ];
                    end
                    if strcmp(mystr,'apalfa1')
                        a1 = Feval(x0,A1,b1,c1,3,mystr,1,r1,ro1);
                        a1=a1(3,:);
                        a2 = Feval(x0,A2,b2,c2,3,mystr,2,r2,ro2);
                        a2=a2(3,:);
                        bigmatrix = [H1 rest1' t*mista1  zeros(n+1,sizer2n); rest1 zeros(sizer1n,sizer1n) a1 zeros(sizer1n,sizer2n);  t*mista2 zeros(m+1,sizer1n) H2 rest2'; a2 zeros(sizer2n,sizer1n) rest2 zeros(sizer2n,sizer2n) ];
                    end
                end
                if (norm([grad2+rest2'*lamb2;gradlam2])<=zerodograd)&&(fazzerograd==1)&&(strcmp(mystr,'ifdois'))
                    if strcmp(mystr,'apalfa1')==0
                        bigmatrix = [H1 rest1' t*mista1  zeros(n,sizer2n); rest1 zeros(sizer1n,sizer1n) zeros(sizer1n,m) zeros(sizer1n,sizer2n);  zeros(size(mista2)) zeros(m,sizer1n) H2 rest2'; zeros(sizer2n,n) zeros(sizer2n,sizer1n) rest2 zeros(sizer2n,sizer2n) ];
                    end
                    if strcmp(mystr,'apalfa1')
                        a1 = Feval(x0,A1,b1,c1,3,mystr,1,r1,ro1);
                        a1=a1(3,:);
                        a2 = Feval(x0,A2,b2,c2,3,mystr,2,r2,ro2);
                        a2=a2(3,:);
                        bigmatrix = [H1 rest1' t*mista1  zeros(n+1,sizer2n); rest1 zeros(sizer1n,sizer1n) a1 zeros(sizer1n,sizer2n);  t*mista2 zeros(m+1,sizer1n) H2 rest2'; a2 zeros(sizer2n,sizer1n) rest2 zeros(sizer2n,sizer2n) ];
                    end
                end
                
                if strcmp(mystr,'apalfa1') n=3; m=3; end
                
                if (min([numero1-zerodograd,numero2-zerodograd])<1e-8)||(passaassim==1)
                    numero1=zerodograd;numero2=zerodograd;
                    passaassim=1;
                    disp('fez o treco');
                    newd = -bigmatrix\biggrad;
                    lamb1 = lamb1velho+newd(n+1:n+sizer1n);
                    lamb2 = lamb2velho+newd(n+m+sizer1n+1:end);
                    lab = max(norm(lamb1)+2*norm(newd(n+1:n+sizer1n)),norm(lamb2)+2*norm(newd(n+m+sizer1n+1:end)));
                    if gamma>=1.1*(lab+gammamax)
                        gamma = (gamma + gammamax + lab)/2;
                    else
                        if gamma < lab + gammamax
                            gamma = max(1.5*gamma,lab+gammamax);
                        end
                        
                    end
                    
                    gradval1 = Feval([x0(1:n);x0(n+1:n+m)+t*newd(n+1+sizer1n:n+m+sizer1n)],A1,b1,c1,0,mystr,1,r1,ro1);
                    gradval2 = Feval([x0(1:n)+t*newd(1:n);x0(n+1:n+m)],A2,b2,c2,0,mystr,2,r2,ro2);
                    produto1 = -(newd(1:n))'*H1*newd(1:n);
                    produto2 = -(newd(n+1:n+m))'*H2*newd(n+1+sizer1n:n+m+sizer1n);
                    newds = [newd(1:n);newd(n+1+sizer1n:n+m+sizer1n)];
                    x0agr = x0+t*newds
                    [y1,y2] = gradrest(x0agr,mystr,0,bdarest1,bdarest2,rest1,rest2,H1,H2,lamb1,lamb2);
                    hnovo1 = y1;  hnovo2 = y2;
                    
                    if strcmp(mystr,'apalfa1') n=2; m=2; end
                    
                end
                
                while (( numero1>zerodograd )||( numero2>zerodograd ))&&(passaassim==0)
                    
                    
                    t=t/2;
                    disp('bugou o numero');
                    if ( numero1>zerodograd )
                        disp('n-z: ' ); disp(numero1-zerodograd);
                    end
                    if ( numero2>zerodograd )
                        disp('n-z: ' ); disp(numero2-zerodograd);
                    end
                    %pause
                    
                    if strcmp(mystr,'apalfa1') n=3; m=3; end
                    
                    newd = -bigmatrix\biggrad;
                    lamb1 = lamb1velho+newd(n+1:n+sizer1n);
                    lamb2 = lamb2velho+newd(n+m+sizer1n+1:end);
                    
                    %update params
                    lab = max(norm(lamb1)+2*norm(newd(n+1:n+sizer1n)),norm(lamb2)+2*norm(newd(n+m+sizer1n+1:end)));
                    if gamma>=1.1*(lab+gammamax)
                        gamma = (gamma + gammamax + lab)/2;
                    else
                        if gamma < lab + gammamax
                            gamma = max(1.5*gamma,lab+gammamax);
                        end
                        
                    end
                    
                    gradval1 = Feval([x0(1:n);x0(n+1:n+m)+t*newd(n+1+sizer1n:n+m+sizer1n)],A1,b1,c1,0,mystr,1,r1,ro1);
                    gradval2 = Feval([x0(1:n)+t*newd(1:n);x0(n+1:n+m)],A2,b2,c2,0,mystr,2,r2,ro2);
                    produto1 = -(newd(1:n))'*H1*newd(1:n);
                    produto2 = -(newd(n+1:n+m))'*H2*newd(n+1+sizer1n:n+m+sizer1n);
                    newds = [newd(1:n);newd(n+1+sizer1n:n+m+sizer1n)];
                    x0agr = x0+t*newds
                    [y1,y2] = gradrest(x0agr,mystr,0,bdarest1,bdarest2,rest1,rest2,H1,H2,lamb1,lamb2);
                    hnovo1 = y1;  hnovo2 = y2;
                    
                    numero1 = Feval(x0+t*newds,A1,b1,c1,0,mystr,1,r1,ro1)-Feval([x0(1:n);x0(n+1:n+m)+t*newd(n+1+sizer1n:n+m+sizer1n)],A1,b1,c1,0,mystr,1,r1,ro1)+gamma*(norm(hnovo1) - norm(gradlam1))-alfa*t*produto1;
                    numero2 = Feval(x0+t*newds,A2,b2,c2,0,mystr,2,r2,ro2)-Feval([x0(1:n)+t*newd(1:n);x0(n+1:n+m)],A2,b2,c2,0,mystr,2,r2,ro2)+gamma*(norm(hnovo2) - norm(gradlam2))-alfa*t*produto2;
                    
                    if strcmp(mystr,'apalfa1') n=2; m=2; end
                    
                end
                
            end
            
            maint = maint+toc;
            
        else
            
            if strcmp(mystr,'apalfa1') n=3; m=3; end
            
            newds = [newd(1:n);newd(n+1+sizer1n:n+m+sizer1n)];
            x0 = x0+t*newds
            erro
            
            lamb1 = lamb1velho+newd(n+1:n+sizer1n);
            lamb2 = lamb2velho+newd(n+m+sizer1n+1:end);
            
            grad1 = Feval(x0,A1,b1,c1,1,mystr,1,r1,ro1);
            grad2 = Feval(x0,A2,b2,c2,1,mystr,2,r2,ro2);
            [y1,y2] = gradrest(x0,mystr,1,bdarest1,bdarest2,rest1,rest2,H1,H2,lamb1,lamb2);
            rest1 = y1; rest2 = y2;
            [y1,y2] = gradrest(x0,mystr,0,bdarest1,bdarest2,rest1,rest2,H1,H2,lamb1,lamb2);
            gradlam1 = y1; gradlam2 = y2;
            erro = norm(biggrad)
            
            if strcmp(mystr,'apalfa1') n=2; m=2; end
            
            k = k+1
            
            lamb1velho = lamb1;
            lamb2velho = lamb2;
            
            
            if ((n==2)&&(m==2)&&((strcmp(mystr,'aplbueno'))||strcmp(mystr,'aplbuenovdd')))
                xvec = [xvec; x0(1)];
                yvec = [yvec; x0(2)];
                
                xvec2 = [xvec2; x0(3)];
                yvec2 = [yvec2; x0(4)];
                
                scatter([xvec;xvec2],[yvec;yvec2],'green');
                hold on
                
            end
            
            if erro<=tol
                aux=0;
            end
            
            maint=maint+toc;
            
            if treset==1
                t = 1;
            end
            
        end
        
    end % end if das strings
    
    if (strcmp(flaghess,'tfixo1')==1)
        
        if strcmp(mystr,'apalfa1') n=3; m=3; end
        newd = -bigmatrix\biggrad;
        lamb1 = lamb1velho+newd(n+1:n+sizer1n);
        lamb2 = lamb2velho+newd(n+m+sizer1n+1:end);
        
        newds = [newd(1:n);newd(n+1+sizer1n:n+m+sizer1n)];
        x0 = x0+t*newds
        %pause
        
        if strcmp(mystr,'apalfa1') n=2; m=2; end
        
        if ((n==2)&&(m==2)&&((strcmp(mystr,'aplbueno'))||strcmp(mystr,'aplbuenovdd')))
            xvec = [xvec; x0(1)];
            yvec = [yvec; x0(2)];
            
            xvec2 = [xvec2; x0(3)];
            yvec2 = [yvec2; x0(4)];
            
            scatter([xvec;xvec2],[yvec;yvec2],'green');
            hold on
            
        end
        
        grad1 = Feval(x0,A1,b1,c1,1,mystr,1,r1,ro1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2,r2,ro2);
        [y1,y2] = gradrest(x0,mystr,1,bdarest1,bdarest2,rest1,rest2,H1,H2,lamb1,lamb2);
        rest1 = y1; rest2 = y2;
        [y1,y2] = gradrest(x0,mystr,0,bdarest1,bdarest2,rest1,rest2,H1,H2,lamb1,lamb2);
        gradlam1 = y1; gradlam2 = y2;
        erro = norm(biggrad);
        
        if erro<=tol
            aux=0;
        end
        
        k = k+1
        
        lamb1velho = lamb1;
        lamb2velho = lamb2;
        maint = maint+toc;
    end
    
    
    
end %end while

if (strcmp(mystr,'aplbuenovdd'))%&&(m==11)
    scatponto=[x0(1),x0(3);x0(2),x0(4)];
    scatter(scatponto(1,:),scatponto(2,:),'o','MarkerFaceColor','red');
    scatcruz=[1,0,0,-1;0,-1,1,0];
    scatter(scatcruz(1,:),scatcruz(2,:),'+','magenta');
    
    hold on
    
    t=linspace(-2,2);
    plot(t,-sqrt((1.1727-t).^2-1));
    plot(t,sqrt((1.3994-t).^2-1));
    
    %plot(t,sqrt((0-t).^2-1));
    %plot(t,-sqrt((0-t).^2-1));
    
    %x = linspace(-2,2); y=x;
    %myn=length(x);
    %for myi = 1:myn
    %    for myj=1:myn
    %        z(myi,myj) = -x(myi) + sqrt(1+y(myj)^2);
    %    end
    %end
    %z = x + sqrt(1+y.^2);
    %contourf1 = @(x,y) x + sqrt(1+y^2) - 1.1727;
    %contourf2 = @(x,y) x + sqrt(1+y^2) - 1.3994;
    %contour(x,y,z,[1.1727;1.3994])
    
    hold off
end



