function [ retx,retk,t,erro ] = Yuan_Varios_Fun(x0,n,m,A1,b1,c1,A2,b2,c2,mystr,r1,r2,ro1,ro2)

vaiplot=0;

tic

%pkg load optim
global k
k = 2; % primeira casa do vetor ?1 =>  primeira iteração k = 2;
cont = 1; % contador de iterações

%CONSTANTES
beta1 = 1/2; %entre 0 e 1
beta2 = 1/2;

delt1 = 0.01; %positivo
delt2 = 0.01; %previous 0.09

tau1 = 1; %positivo
tau2 = 1;

%PONTO INICIAL
x = x0(1:n);
y = x0(n+1:end);

s=[x;y];
s1 = s(1:n);
s2 = s(n+1:end);


%PRIMEIRA PARTE

%Problema Original:

%escolhe string
%mystr = 'vacina';
%mystr = 'problemayuan';
%mystr = 'cubico';
%mystr = 'quadratica_yuan';
%mystr = 'quadtestes1';
%mystr = 'quartico';
%mystr = 'exartigo1_mbom';
%mystr = 'exartigo2_mruim';
%mystr = 'exartigo3_indef';
%mystr = 'exartigo4_cubic';
%mystr = 'exartigo5_futebol';
mystr = 'exnovo6';
%mystr = 'formation';
%mystr = 'formationzero';
%mystr = 'aplicacao';
%mystr='aplbueno';

%com restricao
%mystr = 'quadfacil';
%mystr = 'quadmlau';
%mystr = 'quadfacilmista';
%mystr = 'tenta1';
%mystr='quadhessmista';
%mystr = 'dozero';

switch mystr
    case 'vacinaold'
        
        a = 0.45;
        b1 = 0.3;
        b2 = 0.2;
        
        u1 = @(s)  -(s(1:n)')*(a*s(n+1:end) - b1*ones(n,1));
        u2 = @(s)  (s(n+1:end)')*(a*s(1:n) - b2*ones(m,1));
        
        %Derivadas
        gl_x = @(s)  -(a*s(n+1:end) - b1*ones(n,1)); %primeira derivada de u1 em relação a x
        gll_x = @(s) 0; %segunda derivada de u1 em relação a x
        
        gl_y = @(s) (a*s(1:n) - b2*ones(m,1)); %primeira derivada de u2 em relação a y
        gll_y = @(s) 0; %segunda derivada de u2 em relação a y
        
    case 'problemayuan'
        
        
        u1 = @(s)  (s(1:n)')*(b1*s(n+1:end) - c1);
        u2 = @(s)  (s(n+1:end)')*(b2*s(1:n) - c2);
        
        %Derivadas
        gl_x = @(s)  (b1*s(n+1:end) - c1); %primeira derivada de u1 em relação a x
        gll_x = @(s) zeros(n,n); %segunda derivada de u1 em relação a x
        
        gl_y = @(s) (b2*s(1:n) - c2); %primeira derivada de u2 em relação a y
        gll_y = @(s) zeros(m,m); %segunda derivada de u2 em relação a y
        
        
    case 'vacina'
        
        u1 = @(s)  -(s(1:n)')*(b1*s(n+1:end) - c1);
        u2 = @(s)  (s(n+1:end)')*(b2*s(1:n) - c2);
        
        %Derivadas
        gl_x = @(s)  (c1 - b1*s(n+1:end)); %primeira derivada de u1 em relação a x
        gll_x = @(s) zeros(n,n); %segunda derivada de u1 em relação a x
        
        gl_y = @(s) (b2*s(1:n) - c2); %primeira derivada de u2 em relação a y
        gll_y = @(s) zeros(m,m); %segunda derivada de u2 em relação a y
        
        
    case 'testepaper'
        
        u1 = @(s)  (s(1:n)'*(A1*s(1:n)))/2+((b1*s(n+1:end) - c1)')*s(1:n);
        u2 = @(s)  (s(n+1:end)'*(A2*s(n+1:end)))/2+((b2*s(1:n) - c2)')*s(n+1:end);
        
        %Derivadas
        gl_x = @(s)  A1*s(1:n) + b1*s(n+1:end) - c1; %primeira derivada de u1 em relação a x
        gll_x = @(s) A1; %segunda derivada de u1 em relação a x
        
        gl_y = @(s) A2*s(n+1:end) + b2*s(1:n) - c2; %primeira derivada de u2 em relação a y
        gll_y = @(s) b2; %segunda derivada de u2 em relação a y
        
        
    case 'cubico'
        n = length(s1);
        
        u1 = @(s)  (s(1:n).^3)'*(s(n+1:end).^2)/3 + (norm(s(1:n)))^2/2;
        u2 = @(s)  (s(n+1:end).^3)'*(s(1:n).^2)/3 + (norm(s(n+1:end)))^2/2;
        
        %Derivadas
        gl_x = @(s)  (s(1:n).^2).*(s(n+1:end).^2) + s(1:n); %primeira derivada de u1 em relação a x
        gll_x = @(s) diag(  2*(s(1:n)).*(s(n+1:end).^2) + ones(n,1)   ); %segunda derivada de u1 em relação a x
        
        gl_y = @(s) (s(1:n).^2).*(s(n+1:end).^2) + s(n+1:end); %primeira derivada de u2 em relação a y
        gll_y = @(s) diag(  2*(s(n+1:end)).*(s(1:n).^2) + ones(n,1)   ); %segunda derivada de u2 em relação a y
        
    case 'quadtestes1'
        
        u1 = @(s)  (s(1:n)-A1*ones(n,1))'*(s(1:n)-A1*ones(n,1)) + (s(n+1:end)-b1*ones(m,1))'*(s(n+1:end)-b1*ones(m,1));
        u2 = @(s)  (s(1:n)-A2*ones(n,1))'*(s(1:n)-A2*ones(n,1)) - (s(n+1:end)-b2*ones(m,1))'*(s(n+1:end)-b2*ones(m,1));
        
        %Derivadas
        gl_x = @(s)  2*(s(1:n)-A1*ones(n,1)); %primeira derivada de u1 em relação a x
        gll_x = @(s) 2*eye(n); %segunda derivada de u1 em relação a x
        
        gl_y = @(s) -2*(s(n+1:end)-b2*ones(m,1)); %primeira derivada de u2 em relação a y
        gll_y = @(s) 2*eye(m); %segunda derivada de u2 em relação a y
        
    case 'quadratica_yuan'
        
        u1 = @(s)  (  (s(1:n)-c1)'*(s(1:n)-c1) + (s(n+1:end)-c1)'*(s(n+1:end)-c1)   )/2;
        u2 = @(s)  (  (s(1:n)-c2)'*(s(1:n)-c2) + (s(n+1:end)-c2)'*(s(n+1:end)-c2)   )/2;
        
        %Derivadas
        gl_x = @(s)  s(1:n)-c1; %primeira derivada de u1 em relação a x
        gll_x = @(s) eye(n); %segunda derivada de u1 em relação a x
        
        gl_y = @(s) s(n+1:end)-c2; %primeira derivada de u2 em relação a y
        gll_y = @(s) eye(m); %segunda derivada de u2 em relação a y
        
        
    case 'quartico'
        
        u1 = @(s)  (  ((s(1:n)-b1).^2)'*((s(1:n)-b1).^2)  )/4; %- ( (s1)'*(s2) );
        u2 = @(s)  (  ((s(n+1:end)-c1).^2)'*((s(n+1:end)-c1).^2)  )/4; %- ( (s1)'*(s2) );
        
        %Derivadas
        gl_x = @(s)  (s(1:n)-b1).^3; %- s2; %primeira derivada de u1 em relação a x
        gll_x = @(s)  3*diag((s(1:n)-b1).^2); %segunda derivada de u1 em relação a x
        
        gl_y = @(s) (s(n+1:end)-c1).^3; %- s1; %primeira derivada de u2 em relação a y
        gll_y = @(s) 3*diag((s(n+1:end)-c1).^2); %segunda derivada de u2 em relação a y
        
        
    case 'exartigo1_mbom'
        
        u1 = @(s)  s(1)^2 + s(1)*s(2) - 5*s(1);
        u2 = @(s)  3*(s(2)^2/2) - s(1)*s(2) - s(2);
        
        %Derivadas
        gl_x = @(s)  2*s(1) + s(2) - 5; %primeira derivada de u1 em relação a x
        gll_x = @(s)  2; %segunda derivada de u1 em relação a x
        
        gl_y = @(s) 3*s(2) - s(1) - 1; %primeira derivada de u2 em relação a y
        gll_y = @(s) 3; %segunda derivada de u2 em relação a y
        
        
    case 'exartigo2_mruim'
        
        u1 = @(s)  (s(1)^2)/4 + s(1)*s(2) - 5*s(1);
        u2 = @(s)  (s(2)^2)/6 - s(1)*s(2) - s(2);
        
        %Derivadas
        gl_x = @(s)  s(1)/2 + s(2) - 5; %primeira derivada de u1 em relação a x
        gll_x = @(s)  1/2; %segunda derivada de u1 em relação a x
        
        gl_y = @(s) s(2)/3 - s(1) - 1; %primeira derivada de u2 em relação a y
        gll_y = @(s) 1/3; %segunda derivada de u2 em relação a y
        
    case 'exartigo3_indef'
        
        u1 = @(s)  s(1)^2 + s(1)*s(2) - 5*s(1);
        u2 = @(s)  -3*(s(2)^2/2) - s(1)*s(2) - s(2);
        
        %Derivadas
        gl_x = @(s)  2*s(1) + s(2) - 5; %primeira derivada de u1 em relação a x
        gll_x = @(s)  2; %segunda derivada de u1 em relação a x
        
        gl_y = @(s) -3*s(2) - s(1) - 1; %primeira derivada de u2 em relação a y
        gll_y = @(s) -3; %segunda derivada de u2 em relação a y
        
        
    case 'exartigo4_cubic'
        
        u1 = @(s)  ((s(1)^3)*(s(2)^2))/3 + (s(1)^2)/2;
        u2 = @(s)  ((s(1)^2)*(s(2)^3))/3 + (s(2)^2)/2;
        
        %Derivadas
        gl_x = @(s)  s(1)^2*s(2)^2 + s(1); %primeira derivada de u1 em relação a x
        gll_x = @(s) 2*s(1)*s(2)^2 + 1; %segunda derivada de u1 em relação a x
        
        gl_y = @(s)  s(1)^2*s(2)^2 + s(2); %primeira derivada de u2 em relação a y
        gll_y = @(s) 2*s(2)*s(1)^2 + 1; %segunda derivada de u2 em relação a y
        
    case 'exartigo5_futebol'
        
        u1 = @(s)  -s(1)*(0.6-s(2));
        u2 = @(s)  s(2)*(0.7-s(1));
        
        %Derivadas
        gl_x = @(s)  s(2) - 0.6; %primeira derivada de u1 em relação a x
        gll_x = @(s)  0; %segunda derivada de u1 em relação a x
        
        gl_y = @(s) 0.7 - s(1); %primeira derivada de u2 em relação a y
        gll_y = @(s) 0; %segunda derivada de u2 em relação a y
        
    case 'exartigo6_vacina'
        
        u1 = @(s)  s(1)*(0.45*s(2)-0.3);
        u2 = @(s)  -s(2)*(0.45*s(1)-0.2);
        
        %Derivadas
        gl_x = @(s)  0.45*s(2) - 0.3; %primeira derivada de u1 em relação a x
        gll_x = @(s)  0; %segunda derivada de u1 em relação a x
        
        gl_y = @(s) 0.2 - 0.45*s(1); %primeira derivada de u2 em relação a y
        gll_y = @(s) 0; %segunda derivada de u2 em relação a y
        
         case 'exartigo1_mbom'
        
        u1 = @(s)  s(1)^2 + s(1)*s(2) - 5*s(1);
        u2 = @(s)  3*(s(2)^2/2) - s(1)*s(2) - s(2);
        
        %Derivadas
        gl_x = @(s)  2*s(1) + s(2) - 5; %primeira derivada de u1 em relação a x
        gll_x = @(s)  2; %segunda derivada de u1 em relação a x
        
        gl_y = @(s) 3*s(2) - s(1) - 1; %primeira derivada de u2 em relação a y
        gll_y = @(s) 3; %segunda derivada de u2 em relação a y
        
        
    case 'exnovo6'
        
        u1 = @(s) (1/4)*(s(1)+s(2))^2 + cos(s(1));
        u2 = @(s) (1/4)*(s(1)+s(2))^2 + cos(s(2));
        
        %Derivadas
        
        gl_x = @(s) (1/2)*(s(1)+s(2)) - sin(s(1)); %primeira derivada de u1 em relação a x
        gll_x = @(s) (1/2) - cos(s(1)); %segunda derivada de u1 em relação a x
        
        gl_y = @(s) (1/2)*(s(1)+s(2)) - sin(s(2)); %primeira derivada de u2 em relação a y
        gll_y = @(s) (1/2) - cos(s(2)); %segunda derivada de u2 em relação a y
        
    case 'formation'
        
        u1 = @(s)  (1+c1)*(s(1:n)')*s(1:n) - 2*(s(1:n)')*(b1 + c1*s(n+1:end)) + (b1')*b1 + c1*(s(n+1:end)')*s(n+1:end);
        u2 = @(s)  (1+c2)*(s(n+1:end)')*s(n+1:end) - 2*(s(n+1:end)')*(b2 + c2*s(1:n)) + (b2')*b2 + c2*(s(1:n)')*s(1:n);
        
        %Derivadas
        gl_x = @(s)  2*(1+c1)*s(1:n) - 2*(b1 + c1*s(n+1:end)); %primeira derivada de u1 em relação a x
        gll_x = @(s) 2*(1+c1)*eye(n,n); %segunda derivada de u1 em relação a x
        
        gl_y = @(s) 2*(1+c2)*s(n+1:end) - 2*(b2 + c2*s(1:n)); %primeira derivada de u2 em relação a y
        gll_y = @(s) 2*(1+c2)*eye(n,n); %segunda derivada de u2 em relação a y
        
    case 'formationzero'
        
        u1 = @(s)  (c1)*(s(1:n)')*s(1:n) - 2*(s(1:n)')*( c1*s(n+1:end)) + c1*(s(n+1:end)')*s(n+1:end);
        u2 = @(s)  (c2)*(s(n+1:end)')*s(n+1:end) - 2*(s(n+1:end)')*( c2*s(1:n)) + c2*(s(1:n)')*s(1:n);
        
        %Derivadas
        gl_x = @(s)  2*(c1)*s(1:n) - 2*( c1*s(n+1:end)); %primeira derivada de u1 em relação a x
        gll_x = @(s) 2*(c1)*eye(n,n); %segunda derivada de u1 em relação a x
        
        gl_y = @(s) 2*(c2)*s(n+1:end) - 2*( c2*s(1:n)); %primeira derivada de u2 em relação a y
        gll_y = @(s) 2*(c2)*eye(n,n); %segunda derivada de u2 em relação a y
        
        
        
    case 'aplicacao'
        
        u1 = @(s) (s(1:n)')*(A1*s(n+1:end)) + c1*(sum(s(1:n))-1) + (ro1/2)*sum((max(zeros(n,1),b1/ro1 - s(1:n))).^2);
        u2 = @(s) (s(n+1:end)')*(A2*s(1:n)) + c2*(sum(s(n+1:end))-1) + (ro2/2)*sum((max(zeros(n,1),b2/ro2 - s(n+1:end))).^2);
        
        gl_x = @(s) A1*s(n+1:end) + c1*ones(n,1) - max(zeros(n,1),b1 - ro1*s(1:n));
        gll_x = @(s)  ro1*eye(n,n);
        
        gl_y = @(s) (A2')*s(1:n) + c2*ones(n,1) - max(zeros(n,1),b2 - ro2*s(n+1:end));
        gll_y = @(s)  ro2*eye(n,n);
        
        
    case 'quadfacil'
        
        
        bdarest1 = 1;
        bdarest2 = 1;
        rest1 = ones(1,2);
        rest2 = ones(1,2);
        
        u1 = @(s)  (s(1:n))'*s(1:n) + [-2, 2]*s(1:n) + 2;
        u2 = @(s)  (s(n+1:end))'*s(n+1:end) + [-2, 2]*s(n+1:end) + 2;
        
        %Derivadas
        gl_x = @(s)  2*s(1:n)+[-2;2]; %primeira derivada de u1 em relação a x
        gll_x = @(s) 2*eye(2,2); %segunda derivada de u1 em relação a x
        
        gl_y = @(s) 2*s(n+1:end)+[-2;2]; %primeira derivada de u2 em relação a y
        gll_y = @(s) 2*eye(2,2); %segunda derivada de u2 em relação a y
        
        case 'dozero'
        
        
        bdarest1 = 1;
        bdarest2 = 1;
        rest1 = ones(1,2);
        rest2 = ones(1,2);
        
        u1 = @(s)  (s(1:n))'*s(1:n) ;
        u2 = @(s)  (s(n+1:end))'*s(n+1:end) ;
        
        %Derivadas
        gl_x = @(s)  2*s(1:n); %primeira derivada de u1 em relação a x
        gll_x = @(s) 2*eye(2,2); %segunda derivada de u1 em relação a x
        
        gl_y = @(s) 2*s(n+1:end); %primeira derivada de u2 em relação a y
        gll_y = @(s) 2*eye(2,2); %segunda derivada de u2 em relação a y
        
         
        
end %switch mystr



%DELTA INICIAL e t inicial

t1(2) = 1;
t2(2) = 1;

if (strcmp(mystr,'aplbueno'))||(strcmp(mystr,'formation'))||(strcmp(mystr,'aplicacao'))
    u1s = Feval(s,A1,b1,c1,0,mystr,1,r1,ro1);
    u2s = Feval(s,A2,b2,c2,0,mystr,2,r2,ro2);
    gl_xs = Feval(s,A1,b1,c1,1,mystr,1,r1,ro1);
    gll_xs = Feval(s,A1,b1,c1,2,mystr,1,r1,ro1);
    gl_ys = Feval(s,A2,b2,c2,1,mystr,2,r2,ro2);
    gll_ys = Feval(s,A2,b2,c2,2,mystr,2,r2,ro2);
end


global Delta1
global Delta2
if (strcmp(mystr,'aplbueno')==0)&&(strcmp(mystr,'formation')==0)&&(strcmp(mystr,'aplicacao')==0)
    Delta1(2) = (1/(tau1 + t1(2)))*norm(gl_x(s))
    Delta2(2) = (1/(tau2 + t2(2)))*norm(gl_y(s))
else
    Delta1(2) = (1/(tau1 + t1(2)))*norm(gl_xs)
    Delta2(2) = (1/(tau2 + t2(2)))*norm(gl_ys)
end


saidowhile=1;
while (cont < 1000)&&(saidowhile>1e-4)
      
    
    s=[x;y];
    s1 = s(1:n);
    s2 = s(n+1:end);
    
    if (n==1)&&(m==1)&&(vaiplot==1)
        p(k) = x;
        q(k) = y;
    end
    
    %PROBLEMA AUXILIAR, JOGADOR 1:
    
    if (strcmp(mystr,'aplbueno')==0)&&(strcmp(mystr,'formation')==0)&&(strcmp(mystr,'aplicacao')==0)
        phi1 = @(d1) u1(s) + d1'*gl_x(s) + (d1'*(gll_x(s)*d1))/2; % função objetivo
    else
        phi1 = @(d1) u1s + d1'*gl_xs + (d1'*(gll_xs*d1))/2;
    end
    d1k(:,1) = ones(n,1);
    d1 = d1k(:,1); %ponto inicial exigido pela função nonlin_min
    %Otimização de phi1
    %comentado era do otccve:
    %[d1k(k), objf1(k), cvg, outp] = nonlin_min (phi1, d1, optimset ("inequc", {rest1_x}));
    if (strcmp(mystr,'aplbueno')==0)&&(strcmp(mystr,'formation')==0)&&(strcmp(mystr,'aplicacao')==0)
        if (n==1)&&(m==1)
            gll_x(s)
            gl_x(s)
            pause
            [d1k(:,k),objf1(k)] = quadprog(gll_x(s),gl_x(s),[eye(1,1);-eye(1,1)],[Delta1(k);Delta1(k)]);
            objf1(k) = objf1(k) + u1(s);
        else
            if (strcmp(mystr,'quadfacil')==0)||(strcmp(mystr,'dozero')==0)
                [d1k(:,k),objf1(k)] = quadprog(gll_x(s),gl_x(s),[eye(n,n);-eye(n,n)],Delta1(k)*ones(2*n,1));
                objf1(k) = objf1(k) + u1(s);
            else
                vec1=bdarest1-rest1*s(1:n);
                [d1k(:,k),objf1(k)] = quadprog(gll_x(s),gl_x(s),[eye(n,n);-eye(n,n)],Delta1(k)*ones(2*n,1),rest1,vec1);
                objf1(k) = objf1(k) + u1(s);
            end
        end
    else
        if (n==1)&&(m==1)
            [d1k(:,k),objf1(k)] = quadprog(gll_xs,gl_xs,[1;-1],[Delta1(k),Delta1(k)]);
            objf1(k) = objf1(k) + u1s;
        else
            %[d1k(:,k), objf1(k), cvg, outp] = fmincon(phi1, d1,[],[],[],[],[],[],@rest1_x);
            [d1k(:,k),objf1(k)] = quadprog(gll_xs,gl_xs,[eye(n,n);-eye(n,n)],Delta1(k)*ones(2*n,1));
            objf1(k) = objf1(k) + u1s;
        end
        
    end
    
    %Variáveis para calcular r1(k)
    Pred1(k) = phi1(0*ones(n,1)) - phi1(d1k(:,k));
    s_novo_x = [(x + d1k(:,k));y];
    
    if (strcmp(mystr,'aplbueno'))||(strcmp(mystr,'formation'))||(strcmp(mystr,'aplicacao'))
        u1snovo = Feval(s_novo_x,A1,b1,c1,0,mystr,1,r1,ro1);
        gl_xsnovo = Feval(s_novo_x,A1,b1,c1,1,mystr,1,r1,ro1);
    end
    
    if (strcmp(mystr,'aplbueno')==0)&&(strcmp(mystr,'formation')==0)&&(strcmp(mystr,'aplicacao')==0)
        Ared1(k) = u1(s) - u1(s_novo_x);
    else
        Ared1(k) = u1s - u1snovo;
    end
    
    
    %PROBLEMA AUXILIAR, JOGADOR 2:
    if (strcmp(mystr,'aplbueno')==0)&&(strcmp(mystr,'formation')==0)&&(strcmp(mystr,'aplicacao')==0)
        phi2 = @(d2) u2(s) + d2'*gl_y(s) + (d2'*(gll_y(s)*d2))/2; % função objetivo
    else
        phi2 = @(d2) u2s + d2'*gl_ys + (d2'*(gll_ys*d2))/2; % função objetivo
    end
    d2k(:,1) = ones(m,1);
    d2 = d2k(:,1);

    
    %Otimização de phi2
    %comentado era do octave:
    %[d2k(k), objf2(k), cvg, outp] = nonlin_min (phi2, d2, optimset ("inequc", {rest1_y}));
    if (strcmp(mystr,'aplbueno')==0)&&(strcmp(mystr,'formation')==0)&&(strcmp(mystr,'aplicacao')==0)
        if  (n==1)&&(m==1)
            [d2k(:,k),objf2(k)] = quadprog(gll_y(s),gl_y(s),[1;-1],[Delta2(k),Delta2(k)]);
            objf2(k) = objf2(k) + u2(s);
        else
            if (strcmp(mystr,'quadfacil')==0)||(strcmp(mystr,'dozero')==0)
                [d2k(:,k),objf2(k)] = quadprog(gll_y(s),gl_y(s),[eye(n,n);-eye(n,n)],Delta1(k)*ones(2*n,1));
                objf2(k) = objf2(k) + u2(s);
            else
                vec2=bdarest2-rest2*s(n+1:end);
                [d2k(:,k),objf2(k)] = quadprog(gll_y(s),gl_y(s),[eye(n,n);-eye(n,n)],Delta1(k)*ones(2*n,1),rest2,vec2);
                objf2(k) = objf2(k) + u2(s);
            end
        end
    else
        if  (n==1)&&(m==1)
            [d2k(:,k),objf2(k)] = quadprog(gll_ys,gl_ys,[1;-1],[Delta2(k),Delta2(k)]);
            objf2(k) = objf2(k) + u2s;
        else
            %[d2k(:,k), objf2(k), cvg, outp] = fmincon(phi2, d2,[],[],[],[],[],[],@rest2_x);
            [d2k(:,k),objf2(k)] = quadprog(gll_ys,gl_ys,[eye(n,n);-eye(n,n)],Delta1(k)*ones(2*n,1));
            objf2(k) = objf2(k) + u2s;
        end
    end
    
    
    %Variáveis para calcular r2(k)
    Pred2(k) = phi2(0*ones(m,1)) - phi2(d2k(:,k));
    s_novo_y = [x; (y + d2k(:,k))];
    
    if (strcmp(mystr,'aplbueno'))||(strcmp(mystr,'formation'))||(strcmp(mystr,'aplicacao'))
        u2snovo = Feval(s_novo_y,A2,b2,c2,0,mystr,2,r2,ro2);
        gl_ysnovo = Feval(s_novo_y,A2,b2,c2,1,mystr,2,r2,ro2);
    end
    
    if (strcmp(mystr,'aplbueno')==0)&&(strcmp(mystr,'formation')==0)&&(strcmp(mystr,'aplicacao')==0)
        Ared2(k) = u2(s) - u2(s_novo_y);
    else
        Ared2(k) = u2s - u2snovo;
    end
    
    
    % VETOR F(X)
    % variáveis equivalentes a y na função pi
    if (strcmp(mystr,'aplbueno')==0)&&(strcmp(mystr,'formation')==0)&&(strcmp(mystr,'aplicacao')==0)
        zx = s1 - gl_x(s);
        zy = s2 - gl_y(s);
    else
        zx = s1 - gl_xs;
        zy = s2 - gl_ys;
    end
         %definindo o vetor F
    if (strcmp(mystr,'quadfacil'))||(strcmp(mystr,'dozero'))
        [zx] = quadprog(2*eye(2,2),-2*(s1-gl_x(s)),[],[],rest1,bdarest1);
        %myobj1 = myobj1+norm(gl_y(s))^2;
        [zy] = quadprog(2*eye(2,2),-2*(s2-gl_y(s)),[],[],rest2,bdarest2);
        %myobj2 = myobj2+norm(gl_y(s))^2;
    end

    
    %Eu acho q n precisa nada daquilo ali de cima, da p fazer simplsmente:
    arg1 = zx;
    arg2 = zy;

    
    Fx = [(arg1 - s1); (arg2 - s2)];
    saidowhile = norm(Fx);
    
    %Calculo do Eta_k
    if (strcmp(mystr,'aplbueno')==0)&&(strcmp(mystr,'formation')==0)&&(strcmp(mystr,'aplicacao')==0)
        psi = @(s) sum((gl_x(s)).^2) + sum((gl_y(s)).^2); % somatório das derivads ao quadrado pois não tem restrições.
    else
        psi = @(s) sum((gl_xs).^2) + sum((gl_ys).^2);
    end
    
    %Psi(1) = 0;
    Psi(k) = psi(s);
    if (strcmp(mystr,'quadfacil'))||(strcmp(mystr,'dozero'))
        psis =  norm(arg1 - s1)^2 + norm(arg2 - s2)^2;
        Psi(k) = psis;
    end
    [eta(k), position] = min(Psi);
    
    
    % Atualização do ponto de iteração, rvk>0
    
    if( Ared1(k)*Pred1(k) > 0) % supondo que Pred>0
        x = x + d1k(k);
    end%if
    
    if( Ared2(k)*Pred2(k) > 0)
        y = y + d2k(k);
    end%if
    
    s_novo = [x;y];
    Pred(k) = Pred1(k) + Pred2(k);
    ro(k) = (eta(k) - psi(s_novo))/(Pred(k)); %supondo  sempre Pred>0
    
    % se valhe 2.16
    if ( ro(k) >= beta1)
        % v=1
        %RC aumenta
        if (Ared1(k) > beta2 * Pred1(k))
            [t1(k+1), position] = max([t1(k) - delt1,0]);
        end%if
        %RC constante
        if (0 < Ared1(k)*Pred1(k) && Ared1(k) <= beta2 * Pred1(k))
            t1(k+1) = t1(k);
        end%if
        %RC diminui
        if (Ared1(k)*Pred1(k) <= 0 )
            t1(k+1) = t1(k) + delt1;
        end%if
        % v = 2
        if (Ared2(k) > beta2 * Pred2(k))
            [t2(k+1), position] = max([t2(k) - delt2,0]);
        end%if
        
        if (0 < Ared2(k)*Pred2(k) && Ared2(k) <= beta2 * Pred2(k))
            t2(k+1) = t2(k);
        end%if
        
        if (Ared2(k)*Pred2(k) <= 0 )
            t2(k+1) = t2(k) + delt2;
        end%if
        %se não valhe 2.16
    else
        t1(k+1) = t1(k) + delt1;
        t2(k+1) = t2(k) + delt2;
    end%if
    
    if (strcmp(mystr,'aplbueno')==0)&&(strcmp(mystr,'formation')==0)&&(strcmp(mystr,'aplicacao')==0)
        Delta1(k+1) = (1/(tau1 + t1(k+1))).*norm(gl_x(s_novo));
        Delta2(k+1) = (1/(tau2 + t2(k+1))).*norm(gl_y(s_novo));
    else
        if (strcmp(mystr,'quadfacil')==0)||(strcmp(mystr,'dozero')==0)
            Delta1(k+1) = (1/(tau1 + t1(k+1))).*norm(gl_xsnovo);
            Delta2(k+1) = (1/(tau2 + t2(k+1))).*norm(gl_ysnovo);
        else
            glxprojetado= s1 - quadprog(2*eye(2,2),-2*(s1-gl_x(s)),[],[],rest1,bdarest1);
            glyprojetado= s2 - quadprog(2*eye(2,2),-2*(s2-gl_y(s)),[],[],rest2,bdarest2);
            Delta1(k+1) = (1/(tau1 + t1(k+1))).*norm(glxprojetado);
            Delta2(k+1) = (1/(tau2 + t2(k+1))).*norm(glyprojetado);
        end
    end

    
    %pause
    k = k + 1;
    cont = cont + 1;
    
    if strcmp(mystr,'formation3')
        b1 = s1;
        b2 = s2;
    end
    
end%while


if (n==1)&&(m==1)&&(vaiplot==1)
    p(1)= 0;
    q(1)= 0;
    
    plot (p,q,'.');
end



t = toc;
retx = [x;y];
retk = cont;
erro = saidowhile;

end

