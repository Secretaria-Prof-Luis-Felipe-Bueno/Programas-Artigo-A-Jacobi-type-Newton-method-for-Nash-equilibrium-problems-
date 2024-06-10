clear
%pkg load optim
k = 2; % primeira casa do vetor é 1 =>  primeira iteração k = 2;
cont = 1; % contador de iterações

%CONSTANTES 
beta1 = 1/2; %entre 0 e 1
beta2 = 1/2;

delt1 = 0.0001; %positivo
delt2 = 0.0001;

tau1 = 1; %positivo
tau2 = 1;

%PONTO INICIAL
x = 3; 
y = 0; 

s = [x,y];


%PRIMEIRA PARTE

%Problema Original:

u1 = @(s) ((s(1).^2) + s(1).*s(2) - 5.*s(1)); %Jog 1 : min u1(s) = x² + xy - 5x  x = s(1) in R1
u2 = @(s) ((3.*(s(2).^2))./2 - (s(1).*s(2)) - s(2)); %Jog 2 : min u2(s) = 3y²/2 - xy  - y  y = s(2) in R1


%Derivadas
gl_x = @(s) 2.*s(1) + s(2) - 5; %primeira derivada de u1 em relação a x
gll_x = @(s) 2; %segunda derivada de u1 em relação a x

gl_y = @(s) -s(1) + 3.*s(2) -1; %primeira derivada de u2 em relação a y
gll_y = @(s) 3; %segunda derivada de u2 em relação a y



%DELTA INICIAL e t inicial

t1(2) = 1;
t2(2) = 1;

Delta1(2) = (1/(tau1 + t1(2))).*norm(gl_x(s));
Delta2(2) = (1/(tau2 + t2(2))).*norm(gl_y(s));


while (cont < 100)

fprintf ("\n\n K = %d \n",k)


s = [x;y]
p(k) = x;
q(k) = y;

%PROBLEMA AUXILIAR, JOGADOR 1:   
phi1 = @(d1) u1(s) + d1.*gl_x(s) + (d1.*gll_x(s).*d1)./2; % função objetivo
d1k(1) = 0;
d1 = d1k(1); %ponto inicial exigido pela função nonlin_min

%(Restrição1: delta da região de confiança) 
%antes
%rest1_x = @(d1)  norm(d1) - delt % <= 0; Ñ ACEITA pelo inequc
%rest1_x = @(d1)  - norm(d1) + Delta1(k); % >= 0; ACEITA pelo inequc

% (Restrição 2 : o valor com a mudança deve pertencer ao conjunto) 
%h1_x = @(d1) - s(1) - d1 + 0.5% >= 0 # solução é menor que 3
%h2_x = @(d1) + s(1) + d1 - 0.4 % >= 0 # solução é maior que -3
%Neste exemplo o conjunto de possibilidades são valores entre -3 e 3. 

%Otimização de phi1
%estava aqui do octave
%[d1k(k), objf1(k), cvg, outp] = nonlin_min (phi1, d1, optimset ("inequc", {rest1_x}));
[d1k(k), objf1(k), cvg, outp] = fmincon(phi1, d1, [],[],[],[],[],[],@rest1_x);


%Variáveis para calcular r1(k)
Pred1(k) = phi1(0) - phi1(d1k(k));
s_novo_x = [(x + d1k(k)),y];
Ared1(k) = u1(s) - u1(s_novo_x);


%PROBLEMA AUXILIAR, JOGADOR 2: 
phi2 = @(d2) u2(s) + d2.*gl_y(s) + (d2.*gll_y(s).*d2)./2; % função objetivo
d2k(1) = 0;
d2 = d2k(1);

% (Restrição1) 
%estava aqui antes
%rest1_y = @(d2)  - norm(d2) + Delta2(k); % >= 0; 

% (Restrição 2) 
%Restrição 2 é sobre a soma da coordenada com a mudança ainda permanecer no conjunto. Neste caso a soma certamente estara em R1.
%h1_y = @(d2) - s(2) - d2 + 3 % >= 0 % solução é menor que 3
%h2_y = @(d2) + s(2) + d2 + 3 % >= 0 % solução é maior que -3
%Exemplo o conjunto de possibilidades será valores reais entre -3 e 3. 

%Otimização de phi2
%estava aqui do octava:
%[d2k(k), objf2(k), cvg, outp] = nonlin_min (phi2, d2, optimset ("inequc", {rest1_y}));
[d2k(k), objf2(k), cvg, outp] = fmincon(phi2, d2, [],[],[],[],[],[],@rest1_y);


%Variáveis para calcular r2(k)
Pred2(k) = phi2(0) - phi2(d2k(k));
s_novo_y = [x, (y + d2k(k))];
Ared2(k) = u2(s) - u2(s_novo_y);


% VETOR F(X)
% variáveis equivalentes a y na função pi
zx = s(1) - gl_x(s);
zy = s(2) - gl_y(s);

%função pi
arg_pi1 = @(x) norm(x - zx);
arg_pi2 = @(y) norm(y - zy);

pin = 0; %ponto inicial exigido pela função nonlin_min

%argN são pontos do domínio em que arg_piN é mínimo 
%comentados eram do octave:
%[arg1, objf_pi1, cvg, outp] = nonlin_min (arg_pi1, pin);
%[arg2, objf_pi2, cvg, outp] = nonlin_min (arg_pi2, pin);
[arg1, objf_pi1, cvg, outp] = fmincon(arg_pi1, pin);
[arg2, objf_pi2, cvg, outp] = fmincon(arg_pi2, pin);

Fx = [(arg1 - s(1)); (arg2 - s(2))];

%Calculo do Eta_k
psi = @(s) (gl_x(s)).^2 + (gl_y(s)).^2; % somatório das derivads ao quadrado pois não tem restrições.
%Psi(1) = NA;
Psi(k) = psi(s);
[eta(k), position] = min (Psi);


% Atualização do ponto de iteração, rvk>0
 
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
%se não valhe 2.16
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

p(1)= 0;
q(1)= 0;

plot (p,q);
