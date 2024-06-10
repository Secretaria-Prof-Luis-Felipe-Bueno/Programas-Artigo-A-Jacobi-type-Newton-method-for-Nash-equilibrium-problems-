function [A1,A2,b1,b2,c1,c2,x0,rest1,rest2,bdarest1,bdarest2] = montadados(mystr,n,m)

switch mystr
    case 'quadratica'
        %quad separavel
        % f1 = x'A1x + b1'x + c1, f2 = x'A2x + b2'x + c2
        
                A1 = rand(n,n);
                A1=(A1')*A1;
                A2 = rand(m,m);
                A2=(A2')*A2;
                b1=ones(n,1);
                b2 = ones(m,1);
                c1=0;
                c2=0;
                %descomente para salvar os dados
                %mystr =
                %strcat('.\Salvar\matrizesquadratica',num2str(n),num2str(m));
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                %quadfile = matfile(mystr);
                %A1 = quadfile.A1;
                %A2 = quadfile.A2;
                %b1 = quadfile.b1;
                %b2 = quadfile.b2;
                %c1 = quadfile.c1;
                %c2 = quadfile.c2;
                
         %end switch case quadratica

         
         case 'naonosso'
        %quad separavel
        % f1 = x'A1x + b1'x + c1, f2 = x'A2x + b2'x + c2
        
                
                A1 = eye(2,2);
                A2 = -eye(2,2);
                b1= [-2;2];
                b2 = [2;-2];
                c1=2;
                c2=-2;

         %end switch case naonosso
         
         case 'quadfacil'
        %quad separavel
        % f1 = x'A1x + b1'x + c1, f2 = x'A2x + b2'x + c2
        
                
                A1 = eye(2,2);
                A2 = eye(2,2);
                b1= [-2;2];
                b2 = [-2;2];
                c1=2;
                c2=-2;

         %end switch case quadfacil
         
         case 'quadtranspo'
        %quad separavel
        % f1 = x'A1x + b1'x + c1, f2 = x'A2x + b2'x + c2
        
                
                A1 = [1,3,2,1;3,1,1,2;2,1,1,3;1,2,3,1];
                A2 = A1;
                b1= [0;0];
                b2 = [0;0];
                c1=0;
                c2=0;

         %end switch case quadtranspo
         
         case 'locationeq'
        %quad separavel
        % f1 = x'A1x + b1'x + c1, f2 = x'A2x + b2'x + c2
        
                %A1 seria o p dos clientes
                A1 = [-1,0,0,1;0,-1,1,0];
                A2 = A1;
                b1= [1;1;1;1]; %demanda dos clientes
                b2 = [1;1;1;1];
                c1=zeros(6,1); %multiplicador, no caso
                c2=zeros(6,1);

         %end switch case locationeq
         
          case 'dozero'
        %quad separavel      
                
                A1 = eye(2,2);
                A2 = eye(2,2);
                b1= [0;0];
                b2 = [0;0];
                c1=0;
                c2=0;

                
         %end switch case dozero
         
         case 'quadhessmista'
        %quad separavel
        % f1 = x'A1x + b1'x + c1, f2 = x'A2x + b2'x + c2
 
                A1 = 0.5*eye(2,2);
                A2 = 0.5*eye(2,2);
                b1= [0;0];
                b2 = [0;0];
                c1=0;
                c2=0;
                
                
         %end switch case quadhessmista
         
         
          case 'tenta1'
        %quad separavel
                
                
                A1 = [2,0;0,1];
                A2 = [2,0;0,1];
                b1= [0;0];
                b2 = [0;0];
                c1=-1;
                c2=-1;

                
         %end switch case tenta1
         
          case 'exmisto'
        %quad separavel
                
                
                A1 = [2,0;0,1];
                A2 = [2,0;0,1];
                b1= [0;0];
                b2 = [0;0];
                c1=-1;
                c2=-1;

                
         %end switch case exmisto
         
          case 'quadfacillau'
        %quad separavel
        % f1 = x'A1x + b1'x + c1, f2 = x'A2x + b2'x + c2
        
                A1 = eye(1,1);
                A2 = eye(1,1);
                b1=ones(1,1);
                b2 = ones(1,1);
                c1=0;
                c2=0;
                
                
                A1 = eye(2,2);
                A2 = eye(2,2);
                b1= [-2;2];
                b2 = [-2;2];
                c1=1;
                c2=1;
                %descomente para salvar os dados
                %mystr =
                %strcat('.\Salvar\matrizesquadratica',num2str(n),num2str(m));
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                %quadfile = matfile(mystr);
                %A1 = quadfile.A1;
                %A2 = quadfile.A2;
                %b1 = quadfile.b1;
                %b2 = quadfile.b2;
                %c1 = quadfile.c1;
                %c2 = quadfile.c2;
                
         %end switch case quadfacillau
         
          case 'quadmlau'
               
                A1 = [1,0;0,1];
                A2 = [1,0;0,1];
                b1= [0;0];
                b2 = [0;0];
                c1=0;
                c2=0;
                
                A1 = [1,-1,0;-1,1,1;0,1,0];
                A2 = [1,-1,0;-1,1,1;0,1,0];
                b1= [0;0;0];
                b2 = [0;0;0];
                c1=0;
                c2=0;
                
                %descomente para salvar os dados
                
                %A1 = [4,3;3,1];
                %A2 = [4,3;3,1];
                %b1= [0;0];
                %b2 = [0;0];
                %c1=0;
                %c2=0;
               
                
         %end switch case quadmlau
         
          case 'wanet'
        %quad separavel
               
        %nesse caso A1 eh a matriz dos Ls
                A1 = [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0;
                    0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0;
                    0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0;
                    0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0;
                    0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0;
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0;
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0;
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0;
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1;
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1;
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1;
                    0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0;
                    1,1,0,0,0,0,0,0,0,0,0,1,1,0,0,0;
                    0,0,1,1,1,0,0,0,0,0,0,0,1,0,1,0;
                    0,0,0,0,0,1,1,1,0,0,0,0,0,1,0,1;
                    0,0,0,0,0,0,0,0,1,1,1,0,1,0,1,0;
                    ];
                A2 = A1;
                b1= zeros(16,1);
                b2 = zeros(16,1);
                c1=0;
                c2=0;
                
         %end switch case wanet
         
          case 'quadfacilmista'
        %quad separavel
        % f1 = x'A1y + b1'x + c1, f2 = x'A2y + b2'y + c2
        
               
                A1 = [1,0;0,1];
                A2 = [1,0;0,1];
                b1= [0;0];
                b2 = [0;0];
                c1=0;
                c2=0;
                
                
         %end switch case quadfacilmista
                
        case 'testepaper'
            %quad mais geral do paper
            %f_i = xi'Axi/2 +(bixj-ci)'xi

                A1 = rand(n,n);
                A1=(A1')*A1;
                A2 = rand(m,m);
                A2=(A2')*A2;
                b1 = rand(n,m);
                b2 = rand(m,n);
                c1=zeros(n,1);
                c2=zeros(m,1);
                
                %testando valores fixos
                %A1 = 1;
                %A2 = 2;
                %b1 = 2;
                %b2 = 3;
                %c1 = 1;
                %c2 = 0;
                
                %descomente para salvar os dados
                %mystr =
                %strcat('.\Salvar\matrizestestepaper','num2str(n)','num2str(m)');
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                %quadfile = matfile(mystr);
                %A1 = quadfile.A1;
                %A2 = quadfile.A2;
                %b1 = quadfile.b1;
                %b2 = quadfile.b2;
                %c1 = quadfile.c1;
                %c2 = quadfile.c2;
                
       %end switch case testepaper
        
        case 'problemayuan'
            %problema q o yuan n consegue
            %min_x x(y-0.6), min_y y(0.7-x)
            
                A1 = zeros(n,n);
                A2 = zeros(m,m);
                b1 = eye(n,m);
                b2 = -eye(m,n);
                c1 = 0.6*ones(n,1);
                c2 = -0.7*ones(m,1);
                
                %descomente para salvar os dados
                %mystr =
                %strcat('.\Salvar\matrizesproblemayuan','num2str(n)','num2str(m)');
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                %quadfile = matfile(mystr);
                %A1 = quadfile.A1;
                %A2 = quadfile.A2;
                %b1 = quadfile.b1;
                %b2 = quadfile.b2;
                %c1 = quadfile.c1;
                %c2 = quadfile.c2;
       
        %end switch case problemayuan
        
        case 'quadmisto'
            %quadratica mais geral
            %f(x,y) = x'A1x/2+ y'A2y/2+x'By +c1x + c2x 
            
                A1 = rand(n,n);
                A1=(A1')*A1;
                A2 = rand(m,m);
                A2=(A2')*A2;
                A = zeros(n+m,n+m);
                A(1:n,1:n) = A1;
                A(n+1:end,n+1:end) = A2;
                b = rand(n,m);
                c1 = rand(n,1);
                c2 = rand(m,1);
                c = zeros(n+m,1);
                c(1:n) = c1;
                c(n+1:end) = c2;
                
                
                %so para ficar ok no principal, vc coloca td no msm a e c
                %nesse caso
                A1 = A; A2 = A;
                c1 = c; c2 = c;
                b1 = b; b2 = b;
                
                %descomente para salvar os dados
                %mystr =
                %strcat('.\Salvar\matrizesquadmisto','num2str(n)','num2str(m)');
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                %quadfile = matfile(mystr);
                %A1 = quadfile.A1;
                %A2 = quadfile.A2;
                %b1 = quadfile.b1;
                %b2 = quadfile.b2;
                %c1 = quadfile.c1;
                %c2 = quadfile.c2;
                
        %end switch number case quadmisto
        
        case 'exponencial'
            %exponencial sem minimizador
            %f(x) = sum(ai*exp(bi xi yi))

                A1 = rand(n,1);
                A2 = rand(m,1);
                b1 = rand(n,1);
                b2 = rand(m,1);
                c1 = zeros(n,1);
                c2 = zeros(m,1);
                
                %descomente para salvar os dados
                %mystr =
                %strcat('.\Salvar\matrizesexponencial','num2str(n)','num2str(m)');
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                %quadfile = matfile(mystr);
                %A1 = quadfile.A1;
                %A2 = quadfile.A2;
                %b1 = quadfile.b1;
                %b2 = quadfile.b2;
                %c1 = quadfile.c1;
                %c2 = quadfile.c2;
                
       %end switch number case exponencial
       
        case 'separavel'
            %problema nao convexo
            %fi(x) = sum(e^(xi)*(ci+xi-log(e^xi+e^yi)))
            
                A1 = [];
                A2 = [];
                b1 = [];
                b2 = [];
                c1 = rand(n,1);
                c2 = rand(n,1);
                
                %descomente para salvar os dados
                %mystr =
                %strcat('.\Salvar\matrizesseparavel','num2str(n)','num2str(m)');
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                %quadfile = matfile(mystr);
                %A1 = quadfile.A1;
                %A2 = quadfile.A2;
                %b1 = quadfile.b1;
                %b2 = quadfile.b2;
                %c1 = quadfile.c1;
                %c2 = quadfile.c2;
       
        %end switch case separavel
        
        
        case 'formation'
            %x^3y^2/3 + x^2/2
               n=3;
               m=n;
               c=0;
            
                A1 = [];
                A2 = [];
                b1 = ones(n,1);
                b2 = -ones(n,1);
                c1 = c;
                c2 = c;
       
        %end switch case formation
        
        case 'formationzero'
            %x^3y^2/3 + x^2/2
               n=2;
               m=n;
               c=10;
            
                A1 = [];
                A2 = [];
                b1 = ones(n,1);
                b2 = -ones(n,1);
                c1 = c;
                c2 = c;
       
        %end switch case formationzero
              
        
        case 'quartico'
            %(  ((x-b).^2)'*((x-b).^2)  )/4 - ( (x)'*(y) )/2;
            %(  ((y-c).^2)'*((y-c).^2)  )/4 - ( (x)'*(y) )/2;
            
                A1 = [];
                A2 = [];
                b1 = [0];
                b2 = [0];
                c1 = [1];
                c2 = [1];
                
                %descomente para salvar os dados
                %mystr =
                %strcat('.\Salvar\matrizesquartico','num2str(n)','num2str(m)');
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                %quadfile = matfile(mystr);
                %A1 = quadfile.A1;
                %A2 = quadfile.A2;
                %b1 = quadfile.b1;
                %b2 = quadfile.b2;
                %c1 = quadfile.c1;
                %c2 = quadfile.c2;
       
        %end switch case quartico
        
        case 'quadratica_yuan'
            %f1 = (  (x-b)'*(x-b) + (y-b)'*(y-b)   )/2; 
            %f2 = (  (y-c)'*(y-c) + (x-c)'*(x-c)   )/2;
            
                A1 = [];
                A2 = [];
                b1 = [0];
                b2 = [0];
                c1 = [1];
                c2 = [1];
                
                %descomente para salvar os dados
                %mystr =
                %strcat('.\Salvar\matrizesquadratica_yuan','num2str(n)','num2str(m)');
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                %quadfile = matfile(mystr);
                %A1 = quadfile.A1;
                %A2 = quadfile.A2;
                %b1 = quadfile.b1;
                %b2 = quadfile.b2;
                %c1 = quadfile.c1;
                %c2 = quadfile.c2;
       
        %end switch case quadratica_yuan
        
        case 'diferentes'
            %f1: (x-y)^2
            %f2: (x-(y-a))^2
            
                A1 = [1];
                A2 = [10];
                b1 = [];
                b2 = [];
                c1 = [];
                c2 = [];
                
                %descomente para salvar os dados
                %mystr =
                %strcat('.\Salvar\matrizesdiferentes','num2str(n)','num2str(m)');
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                %quadfile = matfile(mystr);
                %A1 = quadfile.A1;
                %A2 = quadfile.A2;
                %b1 = quadfile.b1;
                %b2 = quadfile.b2;
                %c1 = quadfile.c1;
                %c2 = quadfile.c2;
       
        %end switch case diferentes


       case 'quadtestes1'
            %f1: A1(x-c1)^2+b1(y-c1)^2 (minimo c1,c1)
            %f2: A2(x-c2)^2 + b2(y-c2)^2 (minimo c2,c2)
            %testar se b<0 (uma funcao n tem minimizador)
            
                A1 = [1];
                A2 = [1];
                b1 = [1];
                b2 = [-1];
                c1 = [0];
                c2 = [0];
                
                %descomente para salvar os dados
                %mystr =
                %strcat('.\Salvar\matrizesquadtestes1','num2str(n)','num2str(m)');
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                %quadfile = matfile(mystr);
                %A1 = quadfile.A1;
                %A2 = quadfile.A2;
                %b1 = quadfile.b1;
                %b2 = quadfile.b2;
                %c1 = quadfile.c1;
                %c2 = quadfile.c2;
       
          %end switch case quadtestes1
          
          case 'quadtestes2'
            %f1: A1x^3*y 
            %f2: A2x*y^3 (se A1=-A2 um e menos infinito e o outro mais)
            
            
                A1 = [1];
                A2 = [1];
                b1 = [];
                b2 = [];
                c1 = [];
                c2 = [];
                
                %descomente para salvar os dados
                %mystr =
                %strcat('.\Salvar\matrizesquadtestes2','num2str(n)','num2str(m)');
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                %quadfile = matfile(mystr);
                %A1 = quadfile.A1;
                %A2 = quadfile.A2;
                %b1 = quadfile.b1;
                %b2 = quadfile.b2;
                %c1 = quadfile.c1;
                %c2 = quadfile.c2;
                
       
          %end switch case quadtestes2
          
          case 'quadtestes3'
            %f1: A1(x-b1)
            %f2: A2(y-b2)
            %aqui o grad ?cte, logo n pode zerar
            
                A1 = [1];
                A2 = [-1];
                b1 = [1];
                b2 = [1];
                c1 = [];
                c2 = [];
                
                %descomente para salvar os dados
                %mystr =
                %strcat('.\Salvar\matrizesquadtestes3','num2str(n)','num2str(m)');
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                %quadfile = matfile(mystr);
                %A1 = quadfile.A1;
                %A2 = quadfile.A2;
                %b1 = quadfile.b1;
                %b2 = quadfile.b2;
                %c1 = quadfile.c1;
                %c2 = quadfile.c2;
       
          %end switch case quadtestes3
          
          case 'quadtestes4'
            %f1: A1x^3/3 + b1x
            %f2: A2y^3/3 + b2y
            %grad nao zera mas n eh cte
            
                A1 = [1];
                A2 = [1];
                b1 = [1];
                b2 = [1];
                c1 = [];
                c2 = [];
                
                %descomente para salvar os dados
                %mystr =
                %strcat('.\Salvar\matrizesquadtestes4','num2str(n)','num2str(m)');
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                %quadfile = matfile(mystr);
                %A1 = quadfile.A1;
                %A2 = quadfile.A2;
                %b1 = quadfile.b1;
                %b2 = quadfile.b2;
                %c1 = quadfile.c1;
                %c2 = quadfile.c2;
       
          %end switch case quadtestes4
          
          case 'quadtestes5'
            %f1: sin(b1x) + A1x
            %f2: cos(b2y) + A2y
            %grad n zera mas n eh cte
            
                A1 = [10];
                A2 = [10];
                b1 = [5];
                b2 = [5];
                c1 = [];
                c2 = [];
                
                %descomente para salvar os dados
                %mystr =
                %strcat('.\Salvar\matrizesquadtestes5','num2str(n)','num2str(m)');
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                %quadfile = matfile(mystr);
                %A1 = quadfile.A1;
                %A2 = quadfile.A2;
                %b1 = quadfile.b1;
                %b2 = quadfile.b2;
                %c1 = quadfile.c1;
                %c2 = quadfile.c2;
       
          %end switch case quadtestes5
                    
          
          case 'exp1'
            %f1: -exp(-(x^2+y^2-c^2)^2) quando x=y, minimo 
            %f2: -exp(-(x^2+(y-a)^2)^2-c^2) quando x=y-a, minimo
            
                A1 = [0];%2 1.1
                A2 = [1];%2
                b1 = [];
                b2 = [];
                c1 = [1];%1 0.5
                c2 = [3];%1
                
                %descomente para salvar os dados
                %mystr =
                %strcat('.\Salvar\matrizesexp1','num2str(n)','num2str(m)');
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                %quadfile = matfile(mystr);
                %A1 = quadfile.A1;
                %A2 = quadfile.A2;
                %b1 = quadfile.b1;
                %b2 = quadfile.b2;
                %c1 = quadfile.c1;
                %c2 = quadfile.c2;
       
          %end switch case exp1
                    
          case 'exp2'
            %f1: -exp(-(x^2+y^2-A)^2) quando x^2+y^2=1
            %f2: -exp(-(x^2+y^2-c)^2) quando x^2+y^2=6
            
                A1 = [1];
                A2 = [1];
                b1 = [];
                b2 = [];
                c1 = [6];
                c2 = [6];
                
                %descomente para salvar os dados
                %mystr =
                %strcat('.\Salvar\matrizesexp1','num2str(n)','num2str(m)');
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                %quadfile = matfile(mystr);
                %A1 = quadfile.A1;
                %A2 = quadfile.A2;
                %b1 = quadfile.b1;
                %b2 = quadfile.b2;
                %c1 = quadfile.c1;
                %c2 = quadfile.c2;
       
          %end switch case exp2
          
          case 'arcsin'
            %f1: asin(x^3/3 + x*(y^2-c^2))
            %f2: asin((y-A)^3/3 + y*(x^2-c^2))
            
                A1 = [1];
                A2 = [1];
                b1 = [];
                b2 = [];
                c1 = [1];
                c2 = [1];
                
                %descomente para salvar os dados
                %mystr =
                %strcat('.\Salvar\matrizesarcsin','num2str(n)','num2str(m)');
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                %quadfile = matfile(mystr);
                %A1 = quadfile.A1;
                %A2 = quadfile.A2;
                %b1 = quadfile.b1;
                %b2 = quadfile.b2;
                %c1 = quadfile.c1;
                %c2 = quadfile.c2;
       
          %end switch case arctg
          
          case 'sincos'
            %f1 = sin(valx^3/3 + (valy^2-A)*valx);
            %f2 = cos(valy^3/3 + (valx^2-A)*valy);
            %grad n zera mas n eh cte
            
                A1 = [1];
                A2 = [1];
                b1 = [];
                b2 = [];
                c1 = [];
                c2 = [];
       
          %end switch case sincos
          
           case 'sin2'
            %f1:  sin(valx^3/3 + (valy^2-A)*valx);
            %f2:  sin((valy)^3/3 + (valx^2-c)*valy);
            %grad n zera mas n eh cte
            
                A1 = [1];
                A2 = [1];
                b1 = [0];
                b2 = [0];
                c1 = [1];
                c2 = [1];
       
          %end switch case sin2
          
          case 'testapol'
            %f1: (x-y)^2(x+y)
            %f2: (x-(y-a))^2(x+(y-a))
            %gradiente no caso fica
            
            % grad_x f1 = 2(x-y)(x+y) + (x-y)^2  =  (x-y)[2(x+y-1)]
            % grad_y f2 = 2(x-(y-1))(1-y-x) + (x-(y-1))^2 = (x-(y-1))[2(1-y-x)+1]
            
                A1 = [1];
                A2 = [1];
                b1 = [];
                b2 = [];
                c1 = [];
                c2 = [];
       
          %end switch case testapol
        
          
          case 'testacirc'
            %f1: x^3/3 -y^3/3 + xy^2
            %f2: y^3/3 -x^3/3 + (x-a)^2y
            %gradiente no caso fica
            
            % grad_x f1 = x^2 + y^2
            % grad_y f2 = (x-a)^2 + y^2
            
                A1 = [1];
                A2 = [1];
                b1 = [];
                b2 = [];
                c1 = [];
                c2 = [];
       
          %end switch case testacirc
          
           case 'gradlimmin'
            %f1: ((x-c11)^2 + (y-c12)^2 - A1)^2
            %f2: ((x-c21)^2 + (y-c22)^2 - A2)^2
            
                A1 = [1];
                A2 = [3];
                b1 = [];
                b2 = [];
                c1 = [0;0];
                c2 = [3;4];
       
          %end switch case gradlimmin
          
          case 'exbueno'
            %f1: xy
            %f2: -y
            
                A1 = [];
                A2 = [];
                b1 = [];
                b2 = [];
                c1 = [];
                c2 = [];
       
          %end switch case exbueno
          
          case 'aplicacao'
            %f1: a fucao da aplicacao
            %f1: x'A1y + lamb(sum(xi)-1) + (r/2)sq(sum(xi)-1)^2 + (ro/2)max
            
            %bi vai fazer papel do mui
            %ci vai fazer papel do lambdai
            %gama=2, ro e r comecam com 1 n?
                n = 1; m = 1; strmat='1';
                x0 = [1;1];
                %descomente para salvar os dados
                myloadstr = strcat('.\Salvar\matmistas',num2str(n),num2str(m),strmat);
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                quadfile = matfile(myloadstr);
                A1 = quadfile.A1;
                A2 = quadfile.A2;
                b1 = quadfile.b1;
                b2 = quadfile.b2;
                c1 = quadfile.c1;
                c2 = quadfile.c2;
                
                x0 = [0.9;0.1;0.9;0.1];
                A1 = [1,2;2,77777];
                A2 = [1,2;2,77777];
                b1 = [0;0];
                b2 = [0;0];
                c1=0;
                c2=0;
                
                %x0 = [5;-1];
                %A1 = [1];
                %A2 = [-1];
                %b1 = [0];
                %b2 = [0];
                %c1=0;
                %c2=0;
                
       
          %end switch case aplicacao
          
           case 'apalfa1'
            %f1: a fucao da aplicacao
            %f1: x'A1y + lamb(sum(xi)-1) + (r/2)sq(sum(xi)-1)^2 + (ro/2)max
            
            %bi vai fazer papel do mui
            %ci vai fazer papel do lambdai
            %gama=2, ro e r comecam com 1 n?
                n = 2; m = 2; %strmat='1';
                x0 = [1;1;1;1;1;1];
                %descomente para salvar os dados
                %myloadstr = strcat('.\Salvar\matmistas',num2str(n),num2str(m),strmat);
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                %quadfile = matfile(myloadstr);
                %A1 = quadfile.A1;
                %A2 = quadfile.A2;
                %b1 = quadfile.b1;
                %b2 = quadfile.b2;
                %c1 = quadfile.c1;
                %c2 = quadfile.c2;
                
                %b eh os ws, A eh as posicoes dos clientes
                %x0=(x,q1,y,q2)
                x0 = [0.9;0.1;1;0.9;0.1;1];
                A1 = [1,-1,0,0;0,0,1,-1];
                A2 = [1,-1,0,0;0,0,1,-1];
                b1 = [1;2;1;1];
                b2 = [1;2;1;1];
                c1=1; %ci eh o alfa
                c2=1;
                
                %x0 = [5;-1];
                %A1 = [1];
                %A2 = [-1];
                %b1 = [0];
                %b2 = [0];
                %c1=0;
                %c2=0;
                
       
          %end switch case apalfa1
          
          case 'penalcubic'
            %f1: x'A1y + lamb(sum(xi)-1) + (r/2)sq(sum(xi)-1)^2 + (ro/2)max
            
            %bi vai fazer papel do mui
            %ci vai fazer papel do lambdai
            %gama=2, ro e r comecam com 1 n?
                n = 1; m = 1; strmat='1';
                x0 = [1;1];
                %descomente para salvar os dados
                myloadstr = strcat('.\Salvar\matmistas',num2str(n),num2str(m),strmat);
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                quadfile = matfile(myloadstr);
                A1 = quadfile.A1;
                A2 = quadfile.A2;
                b1 = quadfile.b1;
                b2 = quadfile.b2;
                c1 = quadfile.c1;
                c2 = quadfile.c2;
                
                x0 = [1;0;0;1];
                A1 = [1,2;2,77777];
                A2 = [1,2;2,77777];
                b1 = [1;1];
                b2 = [1;1];
                
       
          %end switch case penalcubic
          
          case 'aplbueno'
            %f1: funcao do bueno
            
            %bi vai fazer papel do mui
            %ci vai fazer papel do lambdai
            %gama=2, ro e r comecam com 1 n?
                n = 2; m = 2; strmat='2';
                x0 = [1;1;-1;-1];
                %descomente para salvar os dados
                %myloadstr = strcat('.\Salvar\matmistas',num2str(n),num2str(m),strmat);
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                
                A1 = [1,0,-1,0;0,1,0,-1];
                A2 = [1,0,-1,0;0,1,0,-1];
                b1 = [1;2;1;1];
                b2 = [1;2;2;3];
                c1=[];
                c2=[];
                
                %unidimensional
                A1 = [1,-1];
                A2 = [1,-1];
                b1 = [1;1];
                b2 = [1;1];
                
       
          %end switch case aplbueno
          
           case 'aplbuenoy'
            %f1: funcao do bueno
            
            %bi vai fazer papel do mui
            %ci vai fazer papel do lambdai
            %gama=2, ro e r comecam com 1 n?
                n = 2; m = 2; strmat='2';
                x0 = [1;1;-1;-1];
                %descomente para salvar os dados
                %myloadstr = strcat('.\Salvar\matmistas',num2str(n),num2str(m),strmat);
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                
                A1 = [1,0,-1,0;0,1,0,-1];
                A2 = [1,0,-1,0;0,1,0,-1];
                b1 = [200;0;1;0];
                b2 = [1;0;200;0];
                c1=[];
                c2=[];
                
                %unidimensional
                A1 = [1,-1];
                A2 = [1,-1];
                b1 = [1;1];
                b2 = [1;1];
                
       
          %end switch case aplbuenoy
          
          case 'aplbuenoyoutro'
            %f1: funcao do bueno
            
            %bi vai fazer papel do mui
            %ci vai fazer papel do lambdai
            %gama=2, ro e r comecam com 1 n?
                n = 2; m = 2; strmat='2';
                x0 = [1;1;-1;-1];
                %descomente para salvar os dados
                %myloadstr = strcat('.\Salvar\matmistas',num2str(n),num2str(m),strmat);
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                
                A1 = [1,0,-1,0;0,1,0,-1];
                A2 = [1,0,-1,0;0,1,0,-1];
                b1 = [200;0;1;0];
                b2 = [1;0;200;0];
                c1=[];
                c2=[];
                
                %unidimensional
                A1 = [1,-1];
                A2 = [1,-1];
                b1 = [1;1];
                b2 = [1;1];
                
       
          %end switch case aplbuenoyoutro
          
           case 'aplbuenox'
            %f1: funcao do bueno
            
            %bi vai fazer papel do mui
            %ci vai fazer papel do lambdai
            %gama=2, ro e r comecam com 1 n?
                n = 2; m = 2; strmat='2';
                x0 = [1;1;-1;-1];
                %descomente para salvar os dados
                %myloadstr = strcat('.\Salvar\matmistas',num2str(n),num2str(m),strmat);
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                
                A1 = [1,0,-1,0;0,1,0,-1];
                A2 = [1,0,-1,0;0,1,0,-1];
                b1 = [200;0;1;0];
                b2 = [1;0;200;0];
                c1=[];
                c2=[];
                
                %unidimensional
                A1 = [1,-1];
                A2 = [1,-1];
                b1 = [1;1];
                b2 = [1;1];
                
       
          %end switch case aplbuenox
          
          case 'aplbuenosem'
            %f1: funcao do bueno
            
            %bi vai fazer papel do mui
            %ci vai fazer papel do lambdai
            %gama=2, ro e r comecam com 1 n?
                n = 2; m = 2; strmat='2';
                x0 = [1;1;-1;-1];
                %descomente para salvar os dados
                %myloadstr = strcat('.\Salvar\matmistas',num2str(n),num2str(m),strmat);
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                
                A1 = [1,0,-1,0;0,1,0,-1];
                A2 = [1,0,-1,0;0,1,0,-1];
                b1 = [200;0;1;0];
                b2 = [1;0;200;0];
                c1=[];
                c2=[];
                
                %unidimensional
                A1 = [1,-1];
                A2 = [1,-1];
                b1 = [1;-1];
                b2 = [1;-1];
                
       
          %end switch case aplbuenosem
          
           case 'aplbuenovdd'
            %f1: funcao do bueno
            
            %bi vai fazer papel do mui
            %ci vai fazer papel do lambdai
            %gama=2, ro e r comecam com 1 n?
                n = 2; m = 2; strmat='2';
                %x0 = [1;1;-1;-1];
                %descomente para salvar os dados
                %myloadstr = strcat('.\Salvar\matmistas',num2str(n),num2str(m),strmat);
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                
                A1 = [1,0,-1,0;0,1,0,-1];
                A2 = [1,0,-1,0;0,1,0,-1];
                b1 = [1;0;0.1429;0];
                b2 = [1;0;0.1429;0];
                b1 = [1;2;1;1];
                b2 = [1;2;2;3];
                
                A1 = [1,-3,3;0,0,0];
                A2 = [1,-3,2;0,0,0];
                %A1 = [1,-2,4;0,0,0];
                %A2 = [1,-2,3;0,0,0];
                b1 = [1;1;1];
                b2 = [1;1;1];
                c1=[];
                c2=[];
                
                %unidimensional
                %A1 = [1,-3,3];%aql meio estranho, e inflexao p 1,-1
                %A2 = [1,-3,2];
                %A1 = [1,-2,4]; %tb estranho
                %A2 = [1,-2,3];
                %paper
                A1 = [1,-1,2];
                A2 = [1,-1,2];
                b1 = [1;1;1];
                b2 = [1;1;1];
                
       
          %end switch case aplbuenovdd
          
          case 'aplbuenovddiq'
            %f1: funcao do bueno
            
            %bi vai fazer papel do mui
            %ci vai fazer papel do lambdai
            %gama=2, ro e r comecam com 1 n?
                n = 2; m = 2; strmat='2';
                x0 = [1;1;-1;-1];
                %descomente para salvar os dados
                %myloadstr = strcat('.\Salvar\matmistas',num2str(n),num2str(m),strmat);
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                
                A1 = [1,0,-1,0;0,1,0,-1];
                A2 = [1,0,-1,0;0,1,0,-1];
                b1 = [1;0;0.1429;0];
                b2 = [1;0;0.1429;0];
                b1 = [1;2;1;1];
                b2 = [1;2;2;3];
                c1=[0]; %multiplicadores (so para o s)
                c2=[0];
                
                %unidimensional
                %A1 = [1,-1,2];
                %A2 = [1,-1,2];
                %b1 = [1;1;1];
                %b2 = [1;1;1];
                
       
          %end switch case aplbuenovddiq
          
          case 'hartunglagrangeano'
            %f1: x'A1y + lamb(sum(xi)-1) + (r/2)sq(sum(xi)-1)^2 + (ro/2)exp
            
            %bi vai fazer papel do mui
            %ci vai fazer papel do lambdai
            %gama=2, ro e r comecam com 1 n?
                n = 1; m = 1; strmat='1';
                x0 = [1;1];
                %descomente para salvar os dados
                myloadstr = strcat('.\Salvar\matmistas',num2str(n),num2str(m),strmat);
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                quadfile = matfile(myloadstr);
                A1 = quadfile.A1;
                A2 = quadfile.A2;
                b1 = quadfile.b1;
                b2 = quadfile.b2;
                c1 = quadfile.c1;
                c2 = quadfile.c2;
                
                x0 = [1;0;0;1];
                A1 = [1,2;2,1];
                A2 = [1,2;2,1];
                b1 = [1;1];
                b2 = [1;1];
                
       
          %end switch case hartunglagrangeano
          
          case 'hartung'
            %f1: x'A1y + lamb(sum(xi)-1) + (r/2)sq(sum(xi)-1)^2 + (ro/2)max
            
            %bi vai fazer papel do mui
            %ci vai fazer papel do lambdai
            %gama=2, ro e r comecam com 1 n?
                n = 1; m = 1; strmat='1';
                %x0 = [1;1];
                %descomente para salvar os dados
                myloadstr = strcat('.\Salvar\matmistas',num2str(n),num2str(m),strmat);
                %save(mystr,'A1','A2','b1','b2','c1','c2')
                
                %para dar load em vez de fazer randomico
                quadfile = matfile(myloadstr);
                A1 = quadfile.A1;
                A2 = quadfile.A2;
                b1 = quadfile.b1;
                b2 = quadfile.b2;
                c1 = quadfile.c1;
                c2 = quadfile.c2;
                
                x0 = [0;1;0.5;0.5];
                n=2;
                m=2;
                A1 = [1,2;2,1];
                A2 = [1,2;2,1];
                c1 = 1;
                c2 = 1;
                %teria q mudar se fosse n diferente de m
                b1 = (1e-4)/(length(x0)/2);
                b2 = (1e-4)/(length(x0)/2);
                
       
          %end switch case hartung
          
         otherwise
                disp('otherwise');
                A1 = [];
                A2 = [];
                b1 = [];
                b2 = [];
                c1 = [];
                c2 = [];
        
end  %end switch mystr

%x0 inicial eh randomico
%x0 = rand(n+m,1);
%if strcmp(mystr,'separavel')||strcmp(mystr,'quimica')||strcmp(mystr,'cubico')
%   x0=abs(x0); 
%end
%x0 = [4;2.9]; %[1 3]


%if (strcmp(mystr,'aplicacao')==0)&&(strcmp(mystr,'hartung')==0)&&(strcmp(mystr,'penalcubic')==0)&&(strcmp(mystr,'hartunglagrangeano')==0)
% x0=[1;2;3;-2;1;1];
% x0=[1;2;4;-3];
%end

%atualiza rest e bda
x0=[0.9;0.1;0.9;0.1];

if strcmp(mystr,'apalfa1')
   x0 = [0.9;0.1;1;0.9;0.1;1]; 
end

%x0=[-3.14;3.14];
%x0 = [-1e-4;sqrt(1e-4)];
%x0=[2;1];
if (strcmp(mystr,'tenta1'))||(strcmp(mystr,'exmisto'))
    x0=[0.3713; -0.5934; 0.3713; -0.5934];
    x0=[2; sqrt(3); 2; sqrt(3)];
end
if strcmp(mystr,'aplbuenovdd')
%x0=[2;3;-2;2];
%x0=[5;3;-5;2];
%x0=[4.5972; 14.2695; -0.1660; -0.1117];
%x0=[-164.5890;  -368.7072;   -14.3585;    37.1194];
x0 = [-1e-4,sqrt(1e-4)];
%x0=[0;0]; %[2,1] outro;  [1,-1]: pto de inflexao
%x0=[2;0;1;0];
%x0=[1;0;-1;0];
end
%x0=[1;1];
bdarest1 = 1; 
bdarest2 = 1;
if strcmp(mystr,'dozero')
    bdarest1=0;
    bdarest2=0;
    x0=[-1; 0; -1; 0];
end
%n = length(x0);
rest1 = ones(1,2);
rest2 = ones(1,2);

if strcmp(mystr,'quadmlau')
    x0=[2;1;3;-2;1;2];
    bdarest1=1;
    bdarest2=1;
    rest1 = ones(1,3);
    rest2 = ones(1,3);
end

if strcmp(mystr,'quadtranspo')
    bdarest1=1;
    bdarest2=1;
    rest1 = ones(1,2);
    rest2 = ones(1,2);
end

if strcmp(mystr,'locationeq')
 x0=[2;3;1;-2;2;1];
end

if strcmp(mystr,'maxes')
 x0=[3;0];
end
if strcmp(mystr,'apalfa1')
 x0=[0.9;0.1;1;0.9;0.1;1];
 bdarest1=1;
    bdarest2=1;
    rest1 = 1;
    rest2 = 1;
end

if strcmp(mystr,'aplbuenovddiq')
    x0=[1;-1;0;-1;1;0]; %os x3 e x6 sao os s
    x0=[2;1;0;-2;1;0];
    x0=[-0.1540; 0.5427; 0.1889; -0.7960; -0.1559; 1.1833];
    x0=[4.5972; 14.2695; 0; -0.1660; -0.1117; 0];
    bdarest1=1;
    bdarest2=1;
    rest1 = ones(1,3);
    rest2 = ones(1,3);
end


%rest1=rest1';
%rest2=rest2';

%x0=[0;2.7]; %exp1
%x0 = [5;-1]; %amanda

%x0 = ones(n+m,1);

end %endfunction