clear

%pkg load optim
global k
k = 2; % primeira casa do vetor � 1 =>  primeira itera��o k = 2;
cont = 1; % contador de itera��es

%CONSTANTES 
beta1 = 1/2; %entre 0 e 1
beta2 = 1/2;

delt1 = 0.01; %positivo
delt2 = 0.09;

tau1 = 1; %positivo
tau2 = 1;

%PONTO INICIAL
n=1;
m=1;
x = ones(n,1); 
y = ones(m,1); 

s = [x;y];
s1 = s(1:n);
s2 = s(n+1:end);


%PRIMEIRA PARTE

%Problema Original:

%escolhe string
%mystr = 'vacina';
mystr = 'problemayuan';
%mystr = 'cubico';
%mystr = 'quadratica_yuan';
%mystr = 'quadtestes1';
%mystr = 'quartico';

switch mystr
    case 'vacina'
        
     a = 0.45;
     b1 = 0.3;
     b2 = 0.2;

     u1 = @(s)  -(s1')*(a*s2 - b1*ones(n,1)); 
     u2 = @(s)  (s2')*(a*s1 - b2*ones(m,1)); 

     %Derivadas
     gl_x = @(s)  -(a*s2 - b1*ones(n,1)); %primeira derivada de u1 em rela��o a x
     gll_x = @(s) 0; %segunda derivada de u1 em rela��o a x

     gl_y = @(s) (a*s1 - b2*ones(m,1)); %primeira derivada de u2 em rela��o a y
     gll_y = @(s) 0; %segunda derivada de u2 em rela��o a y
     
    case 'problemayuan'
        
     b1 = 0.6;
     b2 = 0.7;
        
     u1 = @(s)  (s1')*(s2 - b1*ones(n,1)); 
     u2 = @(s)  (s2')*(s1 + b2*ones(m,1)); 

     %Derivadas     
     gl_x = @(s)  (s2 - b1*ones(n,1)); %primeira derivada de u1 em rela��o a x
     gll_x = @(s) zeros(n,n); %segunda derivada de u1 em rela��o a x

     gl_y = @(s) (s1 + b2*ones(m,1)); %primeira derivada de u2 em rela��o a y
     gll_y = @(s) zeros(m,m); %segunda derivada de u2 em rela��o a y  
     
    case 'cubico'
        
     u1 = @(s)  (s1.^3)'*(s2.^2)/3 + (s1'*s1)/2; 
     u2 = @(s)  (s2.^3)'*(s1.^2)/3 + (s2'*s2)/2;

     %Derivadas
     gl_x = @(s)  (s1.^2).*(s2.^2)+s1; %primeira derivada de u1 em rela��o a x
     gll_x = @(s) eye(n)+2*s1*(s2.^2)'; %segunda derivada de u1 em rela��o a x

     gl_y = @(s) (s1.^2).*(s2.^2)+s2; %primeira derivada de u2 em rela��o a y
     gll_y = @(s) eye(m)+2*s2*(s1.^2)'; %segunda derivada de u2 em rela��o a y     
     
     case 'quadtestes1'
         
     a1 = 1*ones(n,1);
     a2 = 2*ones(n,1);     
     b1 = 1*ones(m,1);
     b2 = 2*ones(m,1);    
        
     u1 = @(s)  (s1-a1)'*(s1-a1) + (s2-b1)'*(s2-b1); 
     u2 = @(s)  (s1-a2)'*(s1-a2) - (s2-b2)'*(s2-b2);

     %Derivadas
     gl_x = @(s)  2*(s1-a1); %primeira derivada de u1 em rela��o a x
     gll_x = @(s) 2*eye(n); %segunda derivada de u1 em rela��o a x

     gl_y = @(s) -2*(s2-b2); %primeira derivada de u2 em rela��o a y
     gll_y = @(s) 2*eye(m); %segunda derivada de u2 em rela��o a y
     
     case 'quadratica_yuan'
         
     b1 = 0*ones(n,1);
     b2 = 1*ones(m,1);

     u1 = @(s)  (  (s1-b1)'*(s1-b1) + (s2-b1)'*(s2-b1)   )/2; 
     u2 = @(s)  (  (s1-b2)'*(s1-b2) + (s2-b2)'*(s2-b2)   )/2;

     %Derivadas
     gl_x = @(s)  s1-b1; %primeira derivada de u1 em rela��o a x
     gll_x = @(s) eye(n); %segunda derivada de u1 em rela��o a x

     gl_y = @(s) s2-b2; %primeira derivada de u2 em rela��o a y
     gll_y = @(s) eye(m); %segunda derivada de u2 em rela��o a y
     
     
    case 'quartico'
        
     b1 = 0*ones(n,1);
     b2 = 1*ones(m,1);
        
     u1 = @(s)  (  ((s1-b1).^2)'*((s1-b1).^2)  )/4 - ( (s1.^2)'*(s2^2) )/2; 
     u2 = @(s)  (  ((s2-b2).^2)'*((s2-b2).^2)  )/4 - ( (s1.^2)'*(s2^2) )/2; 

     %Derivadas
     gl_x = @(s)  (s1-b1).^3 - s1.*(s2.^2); %primeira derivada de u1 em rela��o a x
     gll_x = @(s)  3*diag((s1-b1).^2) - diag(s2.^2); %segunda derivada de u1 em rela��o a x

     gl_y = @(s) (s2-b2).^3 - s2.*(s1.^2); %primeira derivada de u2 em rela��o a y
     gll_y = @(s) 3*diag((s2-b2).^2) - diag(s1.^2); %segunda derivada de u2 em rela��o a y
        
        
     
     
end %switch mystr



%DELTA INICIAL e t inicial

t1(2) = 1;
t2(2) = 1;

global Delta1
global Delta2
Delta1(2) = (1/(tau1 + t1(2)))*norm(gl_x(s));
Delta2(2) = (1/(tau2 + t2(2)))*norm(gl_y(s));


while (cont < 100)

fprintf ("\n\n K = %d \n",k)


s = [x;y]

if (n==1)&&(m==1)
 p(k) = x;
 q(k) = y;
end

%PROBLEMA AUXILIAR, JOGADOR 1: 

phi1 = @(d1) u1(s) + d1'*gl_x(s) + (d1'*(gll_x(s)*d1))/2; % fun��o objetivo
d1k(:,1) = 0*ones(n,1);
d1 = d1k(:,1); %ponto inicial exigido pela fun��o nonlin_min

%(Restri��o1: delta da regi�o de confian�a) 
%rest1_x = @(d1)  norm(d1) - delt # <= 0; � ACEITA pelo inequc

%anterior
%restineqx = @(d1)  - norm(d1) + Delta1(k); % >= 0; ACEITA pelo inequc
%resteqx=[];
%rest1_x = [restineqx,resteqx];

%antes
%function [restineqx,resteqx]=rest1_x(d1)
% global Delta1
% restineqx = - norm(d1) + Delta1(k);
% resteqx=[];
%end %end function


% (Restri��o 2 : o valor com a mudan�a deve pertencer ao conjunto) 
%h1_x = @(d1) - s(1) - d1 + 0.5# >= 0 # solu��o � menor que 3
%h2_x = @(d1) + s(1) + d1 - 0.4 # >= 0 # solu��o � maior que -3
%Neste exemplo o conjunto de possibilidades s�o valores entre -3 e 3. 

%Otimiza��o de phi1
%comentado era do otccve:
%[d1k(k), objf1(k), cvg, outp] = nonlin_min (phi1, d1, optimset ("inequc", {rest1_x}));
[d1k(:,k), objf1(k), cvg, outp] = fmincon(phi1, d1,[],[],[],[],[],[],@rest1_x);


%Vari�veis para calcular r1(k)
Pred1(k) = phi1(0*ones(n,1)) - phi1(d1k(:,k));
s_novo_x = [(x + d1k(:,k));y];
Ared1(k) = u1(s) - u1(s_novo_x);


%PROBLEMA AUXILIAR, JOGADOR 2: 
phi2 = @(d2) u2(s) + d2'*gl_y(s) + (d2'*(gll_y(s)*d2))/2; % fun��o objetivo
d2k(:,1) = 0*ones(m,1);
d2 = d2k(:,1);

% (Restri��o1) 
%estava aqui antes
%restineqy = @(d2)  - norm(d2) + Delta2(k); % >= 0; 
%resdteqy=[];
%rest1_y = [restineqy,resteqy];

% (Restri��o 2) 
%Restri��o 2 � sobre a soma da coordenada com a mudan�a ainda permanecer no conjunto. Neste caso a soma certamente estara em R1.
%h1_y = @(d2) - s(2) - d2 + 3 # >= 0 # solu��o � menor que 3
%h2_y = @(d2) + s(2) + d2 + 3 # >= 0 # solu��o � maior que -3
%Exemplo o conjunto de possibilidades ser� valores reais entre -3 e 3. 

%Otimiza��o de phi2
%comentado era do octave:
%[d2k(k), objf2(k), cvg, outp] = nonlin_min (phi2, d2, optimset ("inequc", {rest1_y}));
[d2k(:,k), objf2(k), cvg, outp] = fmincon(phi2, d2,[],[],[],[],[],[], @rest1_y);


%Vari�veis para calcular r2(k)
Pred2(k) = phi2(0*ones(m,1)) - phi2(d2k(:,k));
s_novo_y = [x, (y + d2k(:,k))];
Ared2(k) = u2(s) - u2(s_novo_y);


% VETOR F(X)
% vari�veis equivalentes a y na fun��o pi
zx = s1 - gl_x(s);
zy = s2 - gl_y(s);


%%%%%%%%% aqui estava, mas acho q n precisa ne, no caso q X eh Rn
%fun��o pi
%arg_pi1 = @(x) norm(x - zx);
%arg_pi2 = @(y) norm(y - zy);

%pin = 0;%ponto inicial exigido pela fun��o nonlin_min

%argN s�o pontos do dom�nio em que arg_piN � m�nimo 
%comentadas eram do octave
%[arg1, objf_pi1, cvg, outp] = nonlin_min (arg_pi1, pin);
%[arg2, objf_pi2, cvg, outp] = nonlin_min (arg_pi2, pin);
%[arg1, objf_pi1, cvg, outp] = fmincon(arg_pi1, pin);
%[arg2, objf_pi2, cvg, outp] = fmincon(arg_pi2, pin);

%Eu acho q n precisa nada daquilo ali de cima, da p fazer simplsmente:
arg1 = zx; arg2 = zy;


%disp('testando a teoria1: '); disp(norm(arg1-zx)); pause
%disp('testando a teoria2: '); disp(norm(arg2-zy)); pause


Fx = [(arg1 - s1); (arg2 - s2)];

%Calculo do Eta_k
psi = @(s) sum((gl_x(s)).^2) + sum((gl_y(s)).^2); % somat�rio das derivads ao quadrado pois n�o tem restri��es.
%Psi(1) = 0;
Psi(k) = psi(s);
[eta(k), position] = min (Psi);


% Atualiza��o do ponto de itera��o, rvk>0
 
if( Ared1(k) > 0) % supondo que Pred>0
    x = x + d1k(k);
end%if
 
if( Ared2(k) > 0)
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
    if (0 <= Ared1(k) && Ared1(k) <= beta2 * Pred1(k))
        t1(k+1) = t1(k); 
    end%if
    %RC diminui
    if (Ared1(k) <= 0 )
        t1(k+1) = t1(k) + delt1; 
    end%if
% v = 2
    if (Ared2(k) > beta2 * Pred2(k))
        [t2(k+1), position] = max([t2(k) - delt2,0]); 
    end%if
    
    if (0 <= Ared2(k) && Ared2(k) <= beta2 * Pred2(k))
        t2(k+1) = t2(k); 
    end%if
    
    if (Ared2(k) <= 0 )
        t2(k+1) = t2(k) + delt2; 
    end%if
%se n�o valhe 2.16
else
  t1(k+1) = t1(k) + delt1;
  t2(k+1) = t2(k) + delt2;   
end%if

  Delta1(k+1) = (1/(tau1 + t1(k+1))).*norm(gl_x(s_novo));
  Delta2(k+1) = (1/(tau2 + t2(k+1))).*norm(gl_y(s_novo));


fprintf ("\t Jog 1: \n");
fprintf ("Pred1(%d) = %f \n",k,Pred1(k));
fprintf ("Ared1(%d) = %f \n",k,Ared1(k));
fprintf ("Delta1(%d) = %f \n",k,Delta1(k));
fprintf ("d1k(%d) = %f \n",k,d1k(k));

fprintf ("\t Jog 2: \n");
fprintf ("Pred2(%d) = %f \n",k,Pred2(k));
fprintf ("Ared2(%d) = %f \n",k,Ared2(k));
fprintf ("Delta2(%d) = %f \n",k,Delta2(k));
fprintf ("d2k(%d) = %f \n",k,d2k(k));


k = k + 1;
cont = cont + 1;
    end%while

if (n==1)&&(m==1)&&(vaiplot=1)
 p(1)= 0;
 q(1)= 0;

 plot (p,q,'.');
end