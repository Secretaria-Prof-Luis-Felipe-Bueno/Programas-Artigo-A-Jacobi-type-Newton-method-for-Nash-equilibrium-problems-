function [x,k,t,erro]=ComparaJacobi(x0,A1,b1,c1,A2,b2,c2,tol,mystr)
% faz o metodo do jacobi


tic
kmax=10000; k=1;

switch mystr
    
    case 'quadratica'
        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);
        lildim=length(b1);
        valx=x0(1:lildim);
        valy=x0(end-lildim+1:end);

       while(k<kmax)&&(erro>tol)
        valx = (A1\(-b1/2));
        valy = (A2\(-b2/2));
        x0=[valx;valy];

        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);   
        k=k+1; 
       end

    %fim do case quadratica
    
    case 'aplbuenovdd'
        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);
        lildim=length(A1(:,1));
        valx=x0(1:lildim);
        valy=x0(end-lildim+1:end);


       while(k<kmax)&&(erro>tol)
       %montar as fcs manualmente pq n da pra fazer geral...  
       f1 = @(x) b1(1)*( (norm(x-A1(:,1))^2)/(norm(x-A1(:,1))^2 + norm(valy-A1(:,1))^2) ) + b1(2)*( (norm(x-A1(:,2))^2)/(norm(x-A1(:,2))^2+ norm(valy-A1(:,2))^2) )...
           + b1(3)*( (norm(x-A1(:,3))^2)/(norm(x-A1(:,3))^2 + norm(valy-A1(:,3))^2) ) + b1(4)*( (norm(x-A1(:,4))^2)/(norm(x-A1(:,4))^2 + norm(valy-A1(:,4))^2) );
       
       f2 = @(x) b1(1)*( (norm(x-A1(:,1))^2)/(norm(x-A1(:,1))^2 + norm(valx-A1(:,1))^2) ) + b1(2)*( (norm(x-A1(:,2))^2)/(norm(x-A1(:,2))^2+ norm(valx-A1(:,2))^2) )...
           + b1(3)*( (norm(x-A1(:,3))^2)/(norm(x-A1(:,3))^2 + norm(valx-A1(:,3))^2) ) + b1(4)*( (norm(x-A1(:,4))^2)/(norm(x-A1(:,4))^2 + norm(valx-A1(:,4))^2) );
           
        nvalx = fminunc(f1,valx);
        nvaly = fminunc(f2,valy);
        x0=[nvalx;nvaly];

        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);   
        k=k+1; 
       end

    %fim do case aplbuenovdd
    
    case 'quadratica_yuan'
        
     n = length(x0)/2;   
     valx = x0(1:n);
     valy = x0(n+1:end);
     grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
     grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
     erro = norm(grad1)+norm(grad2);
     
     while(k<kmax)&&(erro>tol)
        valx = b1;
        valy = c2;
        x0=[valx;valy];

        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);   
        k=k+1; 
     end
     
     %end case quadratica yuan

    case 'testepaper'
    
    %nesse caso bi eh matriz, as funcoes sao f_i = xi'Axi/2 +
     %(bixj-ci)'xi

        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);
        [n,m]= size(b1);
        valx = x0(1:n);
        valy = x0(n+1:end);

     while(k<kmax)&&(erro>tol)
        valxnovo = (A1\(-b1*valy+c1));
        valynovo = (A2\(-b2*valx+c2));
        x0=[valxnovo;valynovo];
        valx=valxnovo;
        valy=valynovo; 

        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);  
        k=k+1;  
       end

    %fim do case testepaper
 
        
   case 'problemayuan'
        %nesse caso bi eh matriz, as funcoes sao f_i = (bixj-ci)'xi 
        %com bi=+-id, e ci=0.6 ou -0.7 pra segunda
        
        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);
        [m,n] = size(b1);
        valx = x0(1:n);
        valy = x0(n+1:end);


     while(k<kmax)&&(erro>tol)
        valxnovo = (b1\c1);
        valynovo = (b2\c2);
        x0=[valxnovo;valynovo];
        valx=valxnovo;
        valy=valynovo; 

        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);  
        k=k+1;  
     end

   % end case problemayuan
   
   
   case 'vacina'
        %nesse caso bi eh matriz, as funcoes sao f_i = +-(bixj-ci)'xi 
        %com bi=0.3, 0.2, e Ai = 0.45 
        
        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);
        [m,n] = size(b1);
        valx = x0(1:n);
        valy = x0(n+1:end);


     while(k<kmax)&&(erro>tol)
        valxnovo = (b1\c1);
        valynovo = (b2\c2);
        x0=[valxnovo;valynovo];
        valx=valxnovo;
        valy=valynovo; 

        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);  
        k=k+1;  
     end

   % end case vacina
   
    case 'exnovo6'
        %nesse caso bi eh matriz, as funcoes sao f_i = +-(bixj-ci)'xi 
        %com bi=0.3, 0.2, e Ai = 0.45 
        
        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);
        [m,n] = size(b1);
        valx = x0(1:n);
        valy = x0(n+1:end);


     while(k<kmax)&&(erro>tol)
        valxnovo = fsolve(@(x) novo6(x,valy),[0]);
        valynovo = fsolve(@(x) novo6(x,valx),[0]);
        x0=[valxnovo;valynovo]
        valx=valxnovo;
        valy=valynovo; 

        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2)  
        k=k+1;
        pause
     end

   % end case exnovo6
  
  
 case 'quadmisto'
     % nesse caso as duas fc sao quadraticas gerais tipo
     % f(x,y) = (xA1x)/2 + (yA2y)/2 + xB1y + c1'x + c2'y.

        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);
     
      m = length(b(1,:));
      n = length(b(:,1));
      %aqui passamos os dados tds, a matrix A eh na vdd A=[A1 A2]
            A1 = A(1:n,1:n);
            A2 = A(n+1:end,n+1:end);
            c1 = c(1:n);
            c2 = c(n+1:end);
            
            valx = x0(1:n);
            valy = x0(n+1:end);   
        
    while(k<kmax)&&(erro>tol)
            valnovox = A1\(-b'*valy-c1); 
            valnovoy = A2\(-b'*valx - c2);
            x0=[valxnovo;valynovo];
            valx=valxnovo;
            valy=valynovo; 
            
            grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
            grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
            erro = norm(grad1)+norm(grad2);  
            k=k+1; 
        end
    % end case quadmisto
    
    case 'cubico'
        lildim=length(x0)/2;
            valx = x0(1:lildim);
            valy = x0(lildim+1:end);

        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);

     while(k<kmax)&&(erro>tol)
        funx = @(x) (x.^2).*(valy.^2) + x;
        funy = @(y) (valx.^2).*(y.^2) + y;
         
        valxnovo = fzero(funx,0);
        valynovo = fzero(funy,0);
        x0=[valxnovo;valynovo];
        valx=valxnovo;
        valy=valynovo; 

        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);  
        k=k+1;  
       end
       
    % end case cubico
    
    
    case 'quartico'
        n=length(x0)/2;
        valx = x0(1:n);
        valy = x0(n+1:end);
            
        %b1 = 0*ones(n,1);
        %b2 = 1*ones(n,1);

        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);

     while(k<kmax)&&(erro>tol)
        
        funx = @(x) (x-b1).^3 - valy;
        funy = @(y) (y-b2).^3 - valx;
         
        valxnovo = fzero(funx,0);
        valynovo = fzero(funy,0);
        x0=[valxnovo;valynovo];
        valx=valxnovo;
        valy=valynovo; 

        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);  
        k=k+1;  
     end
       
    % end case quartico


case 'quadtestes1'
        %f1: A1(x-c1)^2+b1(y-c1)^2 (minimo c1,c1)
        %f2: A2(x-c2)^2 + b2(y-c2)^2 (minimo c2,c2)

        valx=x0(1);valy=x0(2);

       grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);

      while(k<kmax)&&(erro>tol)
        valxnovo = c1/A1;
        valynovo = c2/b2;
        x0=[valxnovo;valynovo];
        valx=valxnovo;
        valy=valynovo; 

        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);  
        k=k+1;  
       end
       
    % end case quadtestes1
    
    
    case 'quadtestes5'
        %%f1: sin(b1x) + A1x
        %f2: cos(b2y) + A2y
        
        valx=x0(1);valy=x0(2);

       grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);

      while(k<kmax)&&(erro>tol)
        funx = @(x) b1*cos(b1*x) + A1;
        funy = @(y) -b2*sin(b2*y) + A2;
         
        valxnovo = fzero(funx,0);
        valynovo = fzero(funy,0);
        x0=[valxnovo;valynovo];
        valx=valxnovo;
        valy=valynovo; 

        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);  
        k=k+1;  
       end
       
    % end case quadtestes5
    
    
    case 'exartigo1_mbom'
        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);
        valx=x0(1);
        valy=x0(2);

       while(k<kmax)&&(erro>tol)
        valx = (5 - valy)/2;
        valy = (1+valx)/3;
        x0=[valx;valy];

        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);   
        k=k+1; 
       end

    %fim do case exartigo1_mbom
    
    case 'exartigo2_mruim'
        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);
        valx=x0(1);
        valy=x0(2);

       while(k<kmax)&&(erro>tol)
        valx = (5 - valy)*2;
        valy = (1+valx)*3;
        x0=[valx;valy];

        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);   
        k=k+1; 
       end

    %fim do case exartigo2_mruim
    
    
    case 'gradlimmin'
        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);
        valx=x0(1);
        valy=x0(2);

       while(k<kmax)&&(erro>tol)
        valx = c1(1) + sqrt((valy - c1(2))^2 - A1);
        valy = c2(2) + sqrt((valx - c2(1))^2 - A2);
        x0=[valx;valy];

        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);   
        k=k+1; 
       end
       
     %fim do case  gradlimmin
    
    case 'exartigo3_indef'
        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);
        valx=x0(1);
        valy=x0(2);

       while(k<kmax)&&(erro>tol)
        valx = (5 - valy)/2;
        valy = -((1+valx)/3);
        x0=[valx;valy];

        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);   
        k=k+1; 
       end

    %fim do case exartigo3_indef
    
    case 'exartigo4_cubic'
        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);
        valx=x0(1);
        valy=x0(2);

       while(k<kmax)&&(erro>tol)
        valx = -1/(valy)^2;
        valy = -1/(valx)^2;
        x0=[valx;valy];

        grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
        grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);   
        k=k+1; 
       end

    %fim do case exartigo4_cubic
    
    
    
    case 'formation'
        
     n = length(x0)/2;   
     valx = x0(1:n);
     valy = x0(n+1:end);
     grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
     grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
     erro = norm(grad1)+norm(grad2);
     
     while(k<kmax)&&(erro>tol)
        valx = (x0(1:n)+c1*x0(n+1:end))/(1+c1);
        valy = (x0(n+1:end)+c2*x0(1:n))/(1+c2);
        x0=[valx;valy];

        grad1 = Feval(x0,A1,x0(1:n),c1,1,mystr,1);
        grad2 = Feval(x0,A2,x0(n+1:end),c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);   
        k=k+1; 
     end
     
     %end case formation
     
     case 'formationzero'
        
     n = length(x0)/2;   
     valx = x0(1:n);
     valy = x0(n+1:end);
     grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
     grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
     erro = norm(grad1)+norm(grad2);
     
     while(k<kmax)&&(erro>tol)
        valx = x0(n+1:end);
        valy = x0(1:n);
        x0=[valx;valy];

        grad1 = Feval(x0,A1,x0(1:n),c1,1,mystr,1);
        grad2 = Feval(x0,A2,x0(n+1:end),c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);   
        k=k+1; 
     end
     
     %end case formationzero
     
     case 'aplbueno'
        
     n = length(x0)/2;   
     valx = x0(1:n);
     valy = x0(n+1:end);
     grad1 = Feval(x0,A1,b1,c1,1,mystr,1);
     grad2 = Feval(x0,A2,b2,c2,1,mystr,2);
     erro = norm(grad1)+norm(grad2);
     
     while(k<kmax)&&(erro>tol)
        valx = x0(n+1:end);
        valy = x0(1:n);
        x0=[valx;valy];

        grad1 = Feval(x0,A1,x0(1:n),c1,1,mystr,1);
        grad2 = Feval(x0,A2,x0(n+1:end),c2,1,mystr,2);
        erro = norm(grad1)+norm(grad2);   
        k=k+1; 
     end
     
     %end case aplbueno
     
     
     y = zeros(lildim,1);
           for j=1:cols
            den = norm(valx - A(:,j))^2 + norm(valx-valy)^2;
            y = y + (2*den*(valx-A(:,j)) - 2*norm(valx-A(:,j))^2*(2*valx-(valy+A(:,j))))/den^2;
           end
    


end %end switch

x = x0;
t = toc ;


end
        