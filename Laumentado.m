function [x0,errofora,contfora,xyuan,eyuan,kforayuan] = Laumentado(flaghess,n,m,ro1,ro2,r1,r2)
tau = 0.9; gama=2; mu1max = 10000*ones(n,1); mu2max = 10000*ones(m,1);
contfora=0; kmaxfora=100; condc=0; kforayuan=0; contador=0;
flagyuan=0; flagesse=0; %1 significa q parou

if nargin<5
    r1 = 0; r2 = 0;
end

alfa = 1e-6; %alfa em (0,1)
tolfora=1e-4;
tol = 1e-3; kmax=1000;
minvaltol = 1e-6;
fixvaltol=1;
maxvaltol=1e10;
zerodograd=1e-6;
cholmodtol=1e-4;
lambdai=0;

%mystr = 'quadfacil';
mystr = 'quadmlau';
%mystr = 'quadfacilmista';
%mystr = 'tenta1';
%mystr='quadhessmista';
%mystr = 'dozero';
%mystr = 'exmisto';
%mystr = 'naonosso';
%mystr = 'aplbuenovdd';
mystr = 'aplbuenovddiq';

penalstr = 'penal';
%penalstr = 'lagrang';


[A1,A2,b1,b2,c1,c2,x0,rest1,rest2,bdarest1,bdarest2] = montadados(mystr,n,m);
x0=[0.9;0.1;0.9;0.1];
x0=[0.9;0.1;0;0.9;0.1;0];
if strcmp(mystr,'aplbuenovddiq')
    x0=[5;3;0;-5;2;0];
    x0=[0.1540; 0.5427;0; -0.7960; -0.1559;0];
    x0=[1; 0;0; -1; -0;0];
    %x0=[1; 4;0; 1; -49;0];
    %x0=[0.1540; 0.5427;0; -0.7960; -0.1559;0];
end

%parametros yuan
xyuan = x0;
ry1 = r1; ry2 = r2; roy1 = ro1; roy2 = ro2;
by1 = b1; by2 = b2; cy1 = c1; cy2 = c2;

lambvelho1=b1;
lambvelho2=b2;
lambvelhoy1 = by1;
lambvelhoy2 = by2;
if strcmp(flaghess,'ifdois')
            x0=[0.1540; 0.5427;0; -0.7960; -0.1559;0];
end

while (contfora<kmaxfora)&&(condc==0)
    xvelho = x0;
    xvelhoyuan = xyuan;
    %lambvelho1=b1;
    %lambvelho2=b2;
    %lambvelhoy1 = by1;
    %lambvelhoy2 = by2;
    vaimontadados=0;
    
    if flagesse==0
        %Geral
        Geral_rest_simpl
        
        lamb1
        lamb2
        pause
        
        contador = contador+k
    end
    
    %parar yuan
    if (nargout>3)&&(flagyuan==0)
        kforayuan = kforayuan+1
        [xyuan,~,~,eyuan] = Yuan_Varios_Fun(xvelhoyuan,n,m,A1,b1,c1,A2,b2,c2,mystr,r1,r2,ro1,ro2);
        if eyuan<tolfora
            flagyuan=1;
        end
    else
        flagyuan=1;
    end
    
    %atualiza parametros para esse
    if flagesse==0
        
        if strcmp(mystr,'penalcubic')==1
            b1 = (1/ro1)*(b1.^2);
            b2 = (1/ro2)*(b2.^2);
        end
        
        %c1 = c1 + r1*(sum(x0(1:n))-1);
        %c2 = c2 + r2*(sum(x0(n+1:n+m))-1);
        
        if strcmp(penalstr,'penal')
            
            if strcmp(mystr,'aplbuenovddiq')
                menor1=0; menor2=0;
                if x0(3)<0
                    menor1=x0(3);
                end
                if x0(6)<0
                    menor1=x0(6);
                end
            else
                menor1=0; menor2=0;
                for i=1:n
                    if x0(i)<0
                        menor1 = menor1 + x0(i);
                    end
                end
                
                for i=n+1:n+m
                    if x0(i)<0
                        menor2 = menor2 + x0(i);
                    end
                end
                
            end
            
            if (norm(menor1) > tolfora)
                ro1 = gama*ro1;
            end
            if (norm(menor2) > tolfora)
                ro2 = gama*ro2;
            end
            
        end
        
        if strcmp(penalstr,'lagrang')
            
            if strcmp(mystr,'aplbuenovddiq')
                menor1=0; menor2=0;
                if x0(3)<0
                    menor1=x0(3);
                end
                if x0(6)<0
                    menor1=x0(6);
                end
            else
                menor1=0; menor2=0;
                for i=1:n
                    if x0(i)<0
                        menor1 = menor1 + x0(i);
                    end
                end
                
                for i=n+1:n+m
                    if x0(i)<0
                        menor2 = menor2 + x0(i);
                    end
                end
                
            end
            
            if strcmp(mystr,'aplbuenovddiq')
                lambfora1 = max(0,c1 - ro1*x0(3));
                lambfora2 = max(0,c2 - ro2*x0(6));
                
                if (norm(max(-x0(3),-lambfora1))) >= tau*(norm(max(-xvelho(3),-lambvelho1)))
                    ro1 = ro1*gama;
                end
                if (norm(max(-x0(6),-lambfora2))) >= tau*(norm(max(-xvelho(6),-lambvelho2)))
                    ro2 = ro2*gama;
                end
                lambvelho1 = lambfora1;
                lambvelho2 = lambfora2;
                c1 = min(1e39,lambfora1);
                c2 = min(1e39,lambfora2);
            else
                lambfora1 = max(zeros(n,1),b1-ro1*x0(1:n));
                lambfora2 = max(zeros(m,1),b2-ro2*x0(n+1:n+m));
                
                if (norm(max(-x0(1:n),-lambfora1))) >= tau*(norm(max(-xvelho(1:n),-lambvelho1)))
                    ro1 = ro1*gama;
                end
                if (norm(max(-x0(n+1:n+m),-lambfora2))) >= tau*(norm(max(-xvelho(n+1:n+m),-lambvelho2)))
                    ro2 = ro2*gama;
                end
                lambvelho1 = lambfora1;
                lambvelho2 = lambfora2;
                b1 = min(1e39*ones(n,1),lambfora1);
                b2 = min(1e39*ones(m,1),lambfora2);
                
            end
            
        end
        
    end%flagesse
    
    %atualiza parametros para yuan
    if flagyuan==0
        for j=1:n
            by1(j) = max([0;by1(j)-roy1*xyuan(j)]);
        end
        for j=1:m
            by2(j) = max([0;by2(j)-roy2*xyuan(n+j)]);
        end
        if strcmp(mystr,'penalcubic')==1
            by1 = (1/roy1)*(by1.^2);
            by2 = (1/roy2)*(by2.^2);
        end
        cy1 = cy1 + ry1*(sum(xyuan(1:n))-1);
        cy2 = cy2 + ry2*(sum(xyuan(n+1:n+m))-1);
        
        if (norm(min(by1,xyuan(1:n))) > tau*norm(min(lambvelhoy1,xvelhoyuan(1:n))))
            roy1 = gama*roy1;
        end
        if (norm(min(by2,xyuan(n+1:n+m))) > tau*norm(min(lambvelhoy2,xvelhoyuan(n+1:n+m))))
            roy2 = gama*roy2;
        end
        if (  norm(1-sum(xyuan(1:n))) > tau*norm(1-sum(xvelhoyuan(1:n)))   )
            ry1 = gama*ry1;
            if nargin<=4
                cy1 = gama*cy1;
            end
        end
        if (   norm(1-sum(xyuan(n+1:n+m))) > tau*norm(1-sum(xvelhoyuan(n+1:n+m)))    )
            ry2 = gama*ry2;
            if nargin<=4
                cy2 = gama*cy2;
            end
        end
        by1 = max(zeros(n,1),min(by1,mu1max));
        by2 = max(zeros(m,1),min(by2,mu2max));
        
    end%flagyuan
    
    
    contfora=contfora+1
    
    erro11 = norm(menor1);
    erro12 = norm(menor2);
    erro21 = norm(erro);
    erro22 = norm(gradlam1+gradlam2);
    if (erro11<tolfora)&&(erro12<tolfora)&&(erro21<tolfora)&&(erro22<tolfora)
        flagesse=1;
    end
    errofora = [erro11;erro12;erro21;erro22];
    
    if (flagesse==1)&&(flagyuan==1)
        condc=1;
    end
    
    %fazer o teste da hessiana da lag
    if strcmp(mystr,'aplbuenovddiq')
        teste1 = Feval(x0,A1,b1,c1,2, mystr,1,r1,ro1);
    teste1 = teste1(1:2,1:2);
    teste1 = teste1+ [0,0;0,lamb1/((1+x0(2)^2)^(3/2))];
    teste2 = Feval(x0,A2,b2,c2,2, mystr,2,r2,ro2);
    teste2 = teste2(1:2,1:2);
    teste2 = teste2+ [0,0;0,lamb2/((1+x0(5)^2)^(3/2))];
    end
    teste=[eigs(teste1),eigs(teste2)]
    %pause
    
end

%fazer o teste da hessiana da lag
if strcmp(mystr,'aplbuenovddiq')
    teste1 = Feval(x0,A1,b1,c1,2, mystr,1,r1,ro1);
    teste1 = teste1(1:2,1:2);
    teste1 = teste1+ [0,0;0,lamb1/((1+x0(2)^2)^(3/2))];
    teste2 = Feval(x0,A2,b2,c2,2, mystr,2,r2,ro2);
    teste2 = teste2(1:2,1:2);
    teste2 = teste2+ [0,0;0,lamb2/((1+x0(5)^2)^(3/2))];
end
teste=[eigs(teste1),eigs(teste2)]
pause

