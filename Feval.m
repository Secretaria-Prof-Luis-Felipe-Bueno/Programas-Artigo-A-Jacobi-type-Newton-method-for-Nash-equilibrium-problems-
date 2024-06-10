function y=Feval(x0,A,b,c,ordem, mystr,qual,r,ro)
% avalia a funcao quadratica f(x) = x'Ax + b'x + c ou a derivada ou a
% hessiana
%talvez eu faca outras fevais de outras funcoes, ai ponho de argumento
%qual indica se eh o primeiro ou segundo vetor

switch mystr
    
    case 'quadratica'
        lildim=length(b);
        if qual==1
            valx=x0(1:lildim);
        else %ai qual=2
            valx=x0(end-lildim+1:end);
        end

     switch ordem
       case 0 %avalia f
            y = valx'*(A*valx)+(b')*valx + c;
       case 1 %avalia grad
            y = 2*(A*valx)+b;
       case 2 %avalia hess
            y = 2*A;
       case 3 %avalia hess mista
            y = zeros(lildim,length(x0)-lildim);
    end % end switch ordem quadratica
    
    case 'quadfacil'
        lildim=length(b);
        if qual==1
            valx=x0(1:lildim);
        else %ai qual=2
            valx=x0(end-lildim+1:end);
        end

     switch ordem
       case 0 %avalia f
            y = valx'*(A*valx)+(b')*valx + c;
       case 1 %avalia grad
            y = 2*(A*valx)+b;
       case 2 %avalia hess
            y = 2*A;
       case 3 %avalia hess mista
            y = zeros(lildim,length(x0)-lildim);
    end % end switch ordem quadfacil
    
    case 'quadtranspo'
        lildim=length(b);
        if qual==1
            valx=x0(1:lildim);
            valy=x0(end-lildim+1:end);
        else %ai qual=2
            valx=x0(end-lildim+1:end);
            valy=x0(1:lildim);
        end

     switch ordem
       case 0 %avalia f
            y = ([valx; valy]')*(A*[valx; valy]);
       case 1 %avalia grad
           if qual==1
              y = 2*(A(1:2,1:2)*valx + A(1:2,3:4)*valy); 
           else
              y = 2*(A(3:4,3:4)*valy + A(3:4,1:2)*valx); 
           end
       case 2 %avalia hess
           if qual==1
               y = 2*A(1:2,1:2);
           else
               y = 2*A(3:4,3:4);
           end
       case 3 %avalia hess mista
            if qual==1
               y = 2*A(1:2,3:4);
           else
               y = 2*A(3:4,1:2);
           end
    end % end switch ordem quadtranspo
    
    case 'apalfa1'
        %lembre 1 n=m=2, e o tanto de z=4 e tem os qs a mais (entao t \'e
        %3x3)
        if qual==1
        valx = x0(1:2);
         q1 = x0(3);
         valy = x0(4:5);
         q2 = x0(6);
        else
            valx = x0(4:5);
         q1 = x0(6);
         valy = x0(1:2);
         q2 = x0(3);
        end

     switch ordem
         
       case 0 %avalia f
           y = 0;
           for i=1:4
               Den = (q2*norm(valx-A(:,i))^2 + q1*norm(valy-A(:,i))^2);
               y = y+ (1/Den)*b(i)*q1*(norm(valy-A(:,i))^2);
           end
           y = y - c*q1;
           
       case 1 %avalia grad
            gradx = zeros(2,1);
            gradq = 0;
            for i=1:4
                Den = (q2*norm(valx-A(:,i))^2 + q1*norm(valy-A(:,i))^2);
                gradx = gradx-2*(1/Den^2)*b(i)*q1*q2*(norm(valy-A(:,i))^2)*(valx-A(:,i));
                gradq = gradq + (1/Den^2)*b(i)*q2*(norm(valy-A(:,i))^2)*(norm(valx-A(:,i))^2);
            end
            gradq = gradq - c;
            y = [gradx;gradq];
            
       case 2 %avalia hess
            hess1 = zeros(2,2);
            vechess = zeros(2,1);
            restim = 0;
            for i=1:4
                Den = (q2*norm(valx-A(:,i))^2 + q1*norm(valy-A(:,i))^2);
                hess1 = hess1 -2*((b(i)*q1*q2*(norm(valy-A(:,i)))^2)/Den^3)*(Den*eye(2,2) - 4*q2*(valx-A(:,i))*(valx-A(:,i))');
                vechess = vechess -2*((b(i)*q2*(norm(valy-A(:,i)))^2)/Den^3)*(Den - 2*q1*norm(valy-A(:,i))^2)*(valx - A(:,i));
                restim = restim -2*(b(i)*q1/Den^3)*(norm(valx-A(:,i)))^4*(norm(valy-A(:,i)))^2;
            end
            restim = restim + 0;
            y = [hess1, vechess; vechess', restim];
            
       case 3 %avalia hess mista
            hess1 = zeros(2,2);
            vechess = zeros(2,1);
            vet = zeros(2,1);
            restim = 0;
            for i=1:4
                Den = (q2*norm(valx-A(:,i))^2 + q1*norm(valy-A(:,i))^2);
                hess1 = hess1 +4*((b(i)*q1*q2*(Den-2*q1*norm(valy-A(:,i)))^2)/Den^3)*(valx-A(:,i))*(valy-A(:,i))';
                vechess = vechess -2*((b(i)*q1*(norm(valy-A(:,i)))^2)/Den^3)*(Den - 2*q2*norm(valx-A(:,i))^2)*(valx - A(:,i));
                vet = vet -2*((b(i)*q2*(norm(valx-A(:,i)))^2)/Den^3)*(Den - 2*q1*norm(valy-A(:,i))^2)*(valy - A(:,i));
                restim = restim +(b(i)/Den^3)*(Den-2*q2*norm(valx-A(:,i))^2)*(norm(valx-A(:,i)))^2*(norm(valy-A(:,i)))^2;
            end
            restim = restim + 0;
            y = [hess1, vechess; vet', restim];
    end % end switch ordem apalfa1
    
    
    case 'locationeq'
            valx=x0(1:2);
            valy=x0(4:5);
     switch ordem
       case 0 %avalia f
           if qual==1
               y = x0(3);%q1
               for j=1:4%4 clientes
                  y = y - b(j)*((x0(3)*norm(valy-A(:,j))^2)/(x0(3)*norm(valy-A(:,j))^2 + x0(6)*norm(valx-A(:,j))^2)); 
               end
           else
               y = x0(6);%q1
               for j=1:4%4 clientes
                  y = y + b(j)*((x0(3)*norm(valy-A(:,j))^2)/(x0(3)*norm(valy-A(:,j))^2 + x0(6)*norm(valx-A(:,j))^2) - 1); 
               end
           end
       case 1 %avalia grad
           y = zeros(3,1);
            if qual == 1
                for j=1:4
                   Den = (x0(3)*norm(valy-A(:,j))^2 + x0(6)*norm(valx-A(:,j))^2);
                   y = y-(1/Den^2)*[-2*x0(3)*x0(6)*norm(valy-A(:,j))^2*(valx-A(:,j)) ; (norm(valx-A(:,j))^2)*(norm(valy-A(:,j))^2)*x0(6)];
                end
            else
                for j=1:4
                   Den = (x0(3)*norm(valy-A(:,j))^2 + x0(6)*norm(valx-A(:,j))^2);
                   y = y+ (1/Den^2)*[2*x0(3)*x0(6)*norm(valx-A(:,j))^2*(valy-A(:,j)) ; -norm(valx-A(:,j))^2*(norm(valy-A(:,j))^2)*x0(3)]; 
                end
            end
       case 2 %avalia hess
            y = zeros(3,3);
            if qual==1
               for j=1:4
                  Den = (x0(3)*norm(valy-A(:,j))^2 + x0(6)*norm(valx-A(:,j))^2);
                  y = y -  (1/Den^4)*[-2*x0(3)*x0(6)*Den^2*((norm(valy-A(:,j)))^2)*eye(2,2) + 8*x0(3)*Den*(x0(6)^2)*(norm(valy-A(:,j))^2)*((valx-A(:,j))*(valx-A(:,j))'), (-2*Den^2*x0(6)*(norm(valy-A(:,j))^2) + 4*x0(3)*Den*x0(6)*(norm(valy-A(:,j))^4)  )*(valx-A(:,j)); -2*x0(6)*(norm(valy-A(:,j))^2)*(Den^2)*(valx-A(:,j))' + 4*Den*x0(3)*x0(6)*(norm(valy-A(:,j))^4)*(valx-A(:,j))' , -2*Den*x0(6)*(norm(valy-A(:,j))^4)*(norm(valx-A(:,j))^2)];
               end
            else
                for j=1:4
                   Den = (x0(3)*norm(valy-A(:,j))^2 + x0(6)*norm(valx-A(:,j))^2);
                   y = y + (1/Den^4)*[-2*x0(3)*x0(6)*Den^2*((norm(valx-A(:,j)))^2)*eye(2,2) + 8*x0(6)*Den*(x0(3)^2)*(norm(valx-A(:,j))^2)*((valy-A(:,j))*(valy-A(:,j))'), (-2*Den^2*x0(3)*(norm(valx-A(:,j))^2) + 4*x0(6)*Den*x0(3)*(norm(valx-A(:,j))^4)  )*(valy-A(:,j)); -2*x0(3)*(norm(valx-A(:,j))^2)*(Den^2)*(valy-A(:,j))' + 4*Den*x0(3)*x0(6)*(norm(valx-A(:,j))^4)*(valy-A(:,j))' , -2*Den*x0(3)*(norm(valx-A(:,j))^4)*(norm(valy-A(:,j))^2)];
               end
            end
       case 3 %avalia hess mista
            y = zeros(3,3);
            if qual==1
               for j=1:4
                  Den = (x0(3)*norm(valy-A(:,j))^2 + x0(6)*norm(valx-A(:,j))^2);
                  y = y - (1/Den^4)*[-4*x0(3)*x0(6)*Den^2*(valx-A(:,j))*(valy-A(:,j))' + 8*(x0(3)^2)*x0(6)*Den*(norm(valy-A(:,j))^2)*(valx-A(:,j))*(valy-A(:,j))',  -2*(norm(valy-A(:,j))^2)*(Den^2)*x0(3)*(valx-A(:,j)) + 4*Den*x0(3)*x0(6)*(norm(valy-A(:,j))^2)*(norm(valx-A(:,j))^2)*(valx-A(:,j));  2*(Den^2)*x0(6)*(norm(valx-A(:,j))^2)*(valy-A(:,j)) - 4*Den*x0(3)*x0(6)*(norm(valy-A(:,j))^2)*(norm(valx-A(:,j))^2)*(valy-A(:,j)),   (Den^2)*(norm(valy-A(:,j))^2)*(norm(valx-A(:,j))^2) - 2*Den*x0(6)*(norm(valy-A(:,j))^2)*(norm(valx-A(:,j))^4)  ];
               end
            else
                for j=1:4
                   Den = (x0(3)*norm(valy-A(:,j))^2 + x0(6)*norm(valx-A(:,j))^2);
                    y = y - (1/Den^4)*[4*x0(3)*x0(6)*Den^2*(valy-A(:,j))*(valx-A(:,j))' - 8*(x0(3)^2)*x0(6)*Den*(norm(valx-A(:,j))^2)*(valy-A(:,j))*(valx-A(:,j))',  2*(Den^2)*x0(6)*(norm(valx-A(:,j))^2)*(valy-A(:,j)) - 4*Den*x0(3)*x0(6)*(norm(valy-A(:,j))^2)*(norm(valx-A(:,j))^2)*(valy-A(:,j));  -2*(norm(valy-A(:,j))^2)*(Den^2)*x0(3)*(valx-A(:,j)) + 4*Den*x0(3)*x0(6)*(norm(valy-A(:,j))^2)*(norm(valx-A(:,j))^2)*(valx-A(:,j)),   -(Den^2)*(norm(valy-A(:,j))^2)*(norm(valx-A(:,j))^2) + 2*Den*x0(3)*(norm(valy-A(:,j))^4)*(norm(valx-A(:,j))^2)  ];
               end
            end
    end % end switch ordem locationeq
    
    
    case 'naonosso'
        lildim=length(b);
        if qual==1
            valx=x0(1:lildim);
        else %ai qual=2
            valx=x0(end-lildim+1:end);
        end

     switch ordem
       case 0 %avalia f
            y = valx'*(A*valx)+(b')*valx + c;
       case 1 %avalia grad
            y = 2*(A*valx)+b;
       case 2 %avalia hess
            y = 2*A;
       case 3 %avalia hess mista
            y = zeros(lildim,length(x0)-lildim);
    end % end switch ordem naonosso
    
    case 'exmisto'
        lildim=length(b);
        if qual==1
            valx=x0(1:lildim);
            valy=x0(end-lildim+1:end);
        else %ai qual=2
            valx=x0(end-lildim+1:end);
            valy=x0(1:lildim);
        end

     switch ordem
       case 0 %avalia f
            y = valx(1)^2 - valx(2)^2 + valx(1)*valy(1);
       case 1 %avalia grad
            y = [2*valx(1) + valy(1); -2*valx(2)];
       case 2 %avalia hess
            y = [2, 0; 0, -2];
       case 3 %avalia hess mista
            y = [1, 0; 0, 0];
    end % end switch ordem exmisto
    
    case 'dozero'
        lildim=length(b);
        if qual==1
            valx=x0(1:lildim);
        else %ai qual=2
            valx=x0(end-lildim+1:end);
        end

     switch ordem
       case 0 %avalia f
            y = valx'*(A*valx)+(b')*valx + c;
       case 1 %avalia grad
            y = 2*(A*valx)+b;
       case 2 %avalia hess
            y = 2*A;
       case 3 %avalia hess mista
            y = zeros(lildim,length(x0)-lildim);
    end % end switch ordem dozero
    
    case 'quadhessmista'
        lildim=length(b);
        if qual==1
            valx=x0(1:lildim);
            valy=x0(end-lildim+1:end);
        else %ai qual=2
            valx=x0(end-lildim+1:end);
            valy=x0(1:lildim);
        end

     switch ordem
       case 0 %avalia f
            y = valx(1)*valy(1) + valx(2)*valy(2);
       case 1 %avalia grad
            y = valy;
       case 2 %avalia hess
            y = zeros(2,2);
       case 3 %avalia hess mista
            y = eye(2,2);
    end % end switch ordem quadhessmista
    
    case 'tenta1'
        lildim=length(b);
        if qual==1
            valx=x0(1:lildim);
        else %ai qual=2
            valx=x0(end-lildim+1:end);
        end

     switch ordem
       case 0 %avalia f
            y = valx'*(A*valx)+(b')*valx + c;
       case 1 %avalia grad
            y = 2*(A*valx)+b;
       case 2 %avalia hess
            y = 2*A;
       case 3 %avalia hess mista
            y = zeros(lildim,length(x0)-lildim);
    end % end switch ordem tenta1
    
     case 'quadfacilmista'
        lildim=length(b);
        if qual==1
            valx=x0(1:lildim);
            valy=x0(end-lildim+1:end);
        else %ai qual=2
            valx=x0(end-lildim+1:end);
            valy=x0(1:lildim);
        end

     switch ordem
       case 0 %avalia f
            y = valx'*(A*valy)+(b')*valx + c;
       case 1 %avalia grad
            y = 2*(A*valx)+b;
       case 2 %avalia hess
            y = 2*A;
       case 3 %avalia hess mista
            y = zeros(lildim,length(x0)-lildim);
    end % end switch ordem quadfacilmista
    
    
    case 'quadfacillau'
        lildim=length(b);
        if qual==1
            valx=x0(1:lildim);
        else %ai qual=2
            valx=x0(end-lildim+1:end);
        end

     switch ordem
       case 0 %avalia f
            y = valx'*(A*valx)+(b')*valx + c + (ro/2)*(valx'*valx);
       case 1 %avalia grad
            y = 2*(A*valx)+b + ro*valx;
       case 2 %avalia hess
           [tn,tm]=size(A);
            y = 2*A+ro*eye(tn,tm);
       case 3 %avalia hess mista
            y = zeros(lildim,length(x0)-lildim);
    end % end switch ordem quadfacillau
    
    case 'quadmlau'
        lildim=length(b);
        if qual==1
            valx=x0(1:lildim);
            valy = x0(end-lildim+1:end);
        else %ai qual=2
            valx=x0(end-lildim+1:end);
            valy = x0(1:lildim);
        end

     switch ordem
       case 0 %avalia f
           penalvec = max(zeros(lildim,1),valx+(1/ro)*b); 
            y = valx'*(A*valy)  + (ro/2)*(penalvec'*penalvec);
       case 1 %avalia grad
           penalvec = max(zeros(lildim,1),valx+(1/ro)*b) ;
            y = 2*(A*valy)- ro*penalvec;
       case 2 %avalia hess
           [tn,tm]=size(A);
            y = ro*eye(tn,tm);
            for i=1:lildim
                if valx(i)<0
                   y(i)=0; 
                end
            end
       case 3 %avalia hess mista
            y = 2*A;
    end % end switch ordem quadmlau
    
    case 'aplbuenox'
        lildim=length(A(:,1));
        cols = length(A(1,:));
        if qual==1
            valx = x0(1:lildim);
            valy = x0(lildim+1:end);
        else %ai qual=2
            valx = x0(end-lildim+1:end);
            valy = x0(1:end-lildim);
        end

     switch ordem
       case 0 %avalia f
           y = 0;
           for j=1:cols
            y = y + b(j)*((norm(valx-A(:,j)))^2/(norm(valx - A(:,j))^2 + norm(valx-valy)^2));
           end
       case 1 %avalia grad
           y = zeros(lildim,1);
           for j=1:cols
            den = norm(valx - A(:,j))^2 + norm(valx-valy)^2;
            y = y + b(j)*(  (2*den*(valx-A(:,j)) - 2*norm(valx-A(:,j))^2*((1+1)*valx-(valy+A(:,j))))/den^2  );
           end
       case 2 %avalia hess
            y = zeros(lildim,lildim);
           for j=1:cols
            den = norm(valx - A(:,j))^2 + norm(valx-valy)^2;
            denl = 2*((1+1)*valx-valy-A(:,j));
            num = 2*den*(valx-A(:,j)) - (norm(valx-A(:,j))^2)*2*((1+1)*valx-(valy+A(:,j)));
            y = y + b(j)*((den^2*(   2*(den*eye(lildim,lildim)+(valx-A(:,j))*denl') - (2*(valx-A(:,j))*denl' + 2*(1+1)*(norm(valx-A(:,j))^2)*eye(lildim,lildim) )    ) - 2*den*num*denl'  )/(den^4));
           end
       case 3 %avalia hess mista
            y = zeros(lildim,lildim);
            for j=1:cols
               den = norm(valx - A(:,j))^2 + norm(valx-valy)^2;
               denl = 2*((1+1)*valx-valy-A(:,j));
               num = 2*den*(valx-A(:,j)) - (norm(valx-A(:,j))^2)*2*((1+1)*valx-(valy+A(:,j)));
               y=y+ (b(j)/(den^4))*( 2*den^2*(norm(valx-A(:,j))^2)*eye(lildim,lildim) - 4*den*num*(valy-A(:,j))' ); 
            end
    end % end switch ordem aplbuenox
    
    case 'aplbueno'
        lildim=length(A(:,1));
        cols = length(A(1,:));
        if qual==1
            valx = x0(1:lildim);
            valy = x0(lildim+1:end);
        else %ai qual=2
            valx = x0(end-lildim+1:end);
            valy = x0(1:end-lildim);
        end

     switch ordem
       case 0 %avalia f
           y = 0;
           for j=1:cols
            y = y + b(j)*((norm(valx-A(:,j)))^2/(norm(valx - A(:,j))^2 + norm(valx-valy)^2));
           end
       case 1 %avalia grad
           y = zeros(lildim,1);
           for j=1:cols
            den = norm(valx - A(:,j))^2 + norm(valx-valy)^2;
            y = y + b(j)*(  (2*den*(valx-A(:,j)) - 2*norm(valx-A(:,j))^2*((1+1)*valx-(valy+A(:,j))))/den^2  );
           end
       case 2 %avalia hess
            y = zeros(lildim,lildim);
           for j=1:cols
            den = norm(valx - A(:,j))^2 + norm(valx-valy)^2;
            denl = 2*((1+1)*valx-valy-A(:,j));
            num = 2*den*(valx-A(:,j)) - (norm(valx-A(:,j))^2)*2*((1+1)*valx-(valy+A(:,j)));
            y = y + b(j)*((den^2*(   2*(den*eye(lildim,lildim)+(valx-A(:,j))*denl') - (2*(valx-A(:,j))*denl' + 2*(1+1)*(norm(valx-A(:,j))^2)*eye(lildim,lildim) )    ) + 8*num*denl'  )/(den^4));
           end
       case 3 %avalia hess mista
            y = zeros(lildim,lildim);
            for j=1:cols
               den = norm(valx - A(:,j))^2 + norm(valx-valy)^2;
               y=y+ b(j)*( ((2*(norm(valx-A(:,j)))^2)/(den^2)) )*eye(lildim,lildim); 
            end
    end % end switch ordem aplbueno
    
    case 'aplbuenoy'
        lildim=length(A(:,1));
        cols = length(A(1,:));
        if qual==1
            valx = x0(1:lildim);
            valy = x0(lildim+1:end);
        else %ai qual=2
            valx = x0(end-lildim+1:end);
            valy = x0(1:end-lildim);
        end

     switch ordem
       case 0 %avalia f
           y = 0;
           for j=1:cols
            y = y + b(j)*((norm(valx-A(:,j)))^2/(norm(valy - A(:,j))^2 + norm(valx-valy)^2));
           end
       case 1 %avalia grad
           y = zeros(lildim,1);
           for j=1:cols
            den = norm(valy - A(:,j))^2 + norm(valx-valy)^2;
            y = y + (b(j)/den^2)*(  2*den*(valx-A(:,j)) - 2*norm(valx-A(:,j))^2*(valx-valy)  );
           end
       case 2 %avalia hess
            y = zeros(lildim,lildim);
           for j=1:cols
            den = norm(valy - A(:,j))^2 + norm(valx-valy)^2;
            denl = 2*(valx-valy);
            num = 2*den*(valx-A(:,j)) - 2*(norm(valx-A(:,j))^2)*(valx-valy);
            y = y + (2*b(j)/(den^4))*(den^2*( 2*(valx-A(:,j))*(valx-valy)' + 2*norm(valy-A(:,j))^2*eye(lildim,lildim)- 2*(valx-valy)*(valx-A(:,j))'    ) - 2*den*num*denl'  );
           end
       case 3 %avalia hess mista
            y = zeros(lildim,lildim);
            for j=1:cols
               den = norm(valy - A(:,j))^2 + norm(valx-valy)^2;
               num = 2*den*(valx-A(:,j)) - 2*(norm(valx-A(:,j))^2)*(valx-valy);
               y=y+ (b(j)/(den^4))*(den^2*norm(valx-A(:,j))^2*eye(lildim,lildim) - 2*den*num*(2*valy -valx - A(:,j))'); 
            end
    end % end switch ordem aplbuenoy
    
    case 'aplbuenoyoutro'
        lildim=length(A(:,1));
        cols = length(A(1,:));
        if qual==1
            valx = x0(1:lildim);
            valy = x0(lildim+1:end);
        else %ai qual=2
            valx = x0(end-lildim+1:end);
            valy = x0(1:end-lildim);
        end

     switch ordem
       case 0 %avalia f
           y = 0;
           for j=1:cols
            y = y + b(j)*((norm(valx-A(:,j)))^2/(norm(valy - A(:,j))^2 + norm(valx-valy)^2));
           end
       case 1 %avalia grad
           y = zeros(lildim,1);
           for j=1:cols
            den = norm(valy - A(:,j))^2 + norm(valx-valy)^2;
            y = y + (b(j)/den^2)*(  2*den*(valx-A(:,j)) - 2*norm(valx-A(:,j))^2*(valx-valy)  );
           end
       case 2 %avalia hess
            y = zeros(lildim,lildim);
           for j=1:cols
            den = norm(valy - A(:,j))^2 + norm(valx-valy)^2;
            denl = 2*(valx-valy);
            num = 2*den*(valx-A(:,j)) - 2*(norm(valx-A(:,j))^2)*(valx-valy);
            y = y + (2*b(j)/(den^4))*(den^2*( 2*(valx-A(:,j))*(valx-valy)' - 2*(valx-valy)*(valx-A(:,j))'    ) + 8*den*num*denl'  );
           end
       case 3 %avalia hess mista
            y = zeros(lildim,lildim);
            for j=1:cols
               den = norm(valy - A(:,j))^2 + norm(valx-valy)^2;
               num = 2*den*(valx-A(:,j)) - 2*(norm(valx-A(:,j))^2)*(valx-valy);
               y=y + (2*b(j)/(den^2))*(2*norm(valx-A(:,j))^2*eye(lildim,lildim) ); 
            end
    end % end switch ordem aplbuenoyoutro
    
    
    case 'aplbuenosem'
        lildim=length(A(:,1));
        cols = length(A(1,:));
        if qual==1
            valx = x0(1:lildim);
            valy = x0(lildim+1:end);
        else %ai qual=2
            valx = x0(end-lildim+1:end);
            valy = x0(1:end-lildim);
        end
        
        const=1;

     switch ordem
       case 0 %avalia f
           y = 0;
           for j=1:cols
            y = y + b(j)*((norm(valx-A(:,j)))^2/(const+norm(valy-valx)^2));
           end
       case 1 %avalia grad
           y = zeros(lildim,1);
           for j=1:cols
            %if valx==valy
                den = const+norm(valy-valx)^2;
            %else
            %    den=1;
            %end
            y = y + ((2*b(j))/den^2)*( den*( valx-A(:,j) ) - (norm(valx-A(:,j)))^2*( valx-valy )  );
           end
       case 2 %avalia hess
            y = zeros(lildim,lildim);
           for j=1:cols
            %if valx==valy
                den = (const+norm(valy-valx)^2)^2;
            %else
            %    den=1;
            %end
            num = den*( valx-A(:,j) ) - (norm(valx-A(:,j)))^2*( valx-valy );
            y = y + ((2*b(j))/(den^2))*( (den^2)*(-2*(valx-A(:,j))*(valy-valx)' + (const+(norm(valy-valx)^2)-(norm(valx-valy)^2) )*eye(lildim,lildim) - 2*(valx-valy)*(valx-A(:,j))' ) - 4*(const+norm(valx-valy))^3*num*(valx-valy)');
           end
       case 3 %avalia hess mista
            y = zeros(lildim,lildim);
            for j=1:cols
               %if valx==valy
                den = (const+norm(valy-valx)^2)^2;
               %else
               % den=1;
               %end
               num = den^2*( valx-A(:,j) ) - (norm(valx-A(:,j)))^2*( valx-valy );
               y=y +  (2*b(j)/(den^2))*( (den^2)*(2*(valx-A(:,j))*(valy-valx)'+(norm(valx-A(:,j))^2)*eye(lildim,lildim)) + 4*(const+norm(valx-valy))^3*num*(valx-valy)' ); 
            end
    end % end switch ordem aplbuenosem
    
    
    case 'aplbuenovdd'
        lildim=length(A(:,1));
        cols = length(A(1,:));
        if qual==1
            valx = x0(1:lildim);
            valy = x0(lildim+1:end);
        else %ai qual=2
            valx = x0(end-lildim+1:end);
            valy = x0(1:end-lildim);
        end

     switch ordem
       case 0 %avalia f
           y = 0;
           for j=1:cols
            y = y + b(j)*((norm(valx-A(:,j)))^2/(norm(valx - A(:,j))^2 + norm(valy-A(:,j))^2));
           end
       case 1 %avalia grad
           y = zeros(lildim,1);
           for j=1:cols
            den = norm(valx - A(:,j))^2 + norm(valy-A(:,j))^2;
            y = y + (2*b(j)*(norm(valy-A(:,j))^2)/den^2)*( valx-A(:,j)  );
           end
       case 2 %avalia hess
            y = zeros(lildim,lildim);
           for j=1:cols
            den = norm(valx - A(:,j))^2 + norm(valy-A(:,j))^2;
            denl = 2*(valx-A(:,j));
            %num = -2*(norm(valy-A(:,j))^2)*(valx-A(:,j));

              y = y + (2*b(j)/(den^4))*( 2*(norm(valy-A(:,j))^2)*den^2*eye(lildim,lildim)  - 2*(norm(valy-A(:,j))^2)*(den)*(valx-A(:,j))*(valx-A(:,j))' );
              %y = y + (2*b(j)/(den^4))*((valy-A(:,j))^2)*(-3*valx^2+6*valx*A(:,j)+valy^2-2*valy*A(:,j)-2*A(:,j)^2);
       
           end
       case 3 %avalia hess mista
            y = zeros(lildim,lildim);
            for j=1:cols
               den = norm(valx - A(:,j))^2 + norm(valy-A(:,j))^2;
               y=y+  (2*b(j)/(den^4))*( 2*den^2*(valx-A(:,j))*(valy-A(:,j))' -  2*(den*norm(valy-A(:,j))^2)*(valx-A(:,j))*(valy-A(:,j))' ); 
               %y=y-  (4*b(j)/(den^4))*( (valx-A(:,j))^3*(valy-A(:,j))); 
            end
    end % end switch ordem aplbuenovdd
    
    case 'aplbuenovddiq'
        lildim=length(A(:,1));
        cols = length(A(1,:));
        if qual==1
            valx = x0(1:2);
            valy = x0(4:5);
            vals=x0(3);
        else %ai qual=2
            valx = x0(4:5);
            valy = x0(1:2);
            vals=x0(6);
        end

     switch ordem
       case 0 %avalia f
           y = 0;
           for j=1:cols
            y = y + b(j)*((norm(valx-A(:,j)))^2/(norm(valx - A(:,j))^2 + norm(valy-A(:,j))^2)) + (ro/2)*(max([vals;0]).^2);
           end
       case 1 %avalia grad
           y = zeros(lildim,1);
           
           for j=1:cols
            den = norm(valx - A(:,j))^2 + norm(valy-A(:,j))^2;
            y = y + (2*b(j)*(norm(valy-A(:,j))^2)/den^2)*( valx-A(:,j)  );
           end
           y = [y;ro*max([vals;0])];
       case 2 %avalia hess
            y = zeros(lildim,lildim);
           for j=1:cols
            den = norm(valx - A(:,j))^2 + norm(valy-A(:,j))^2;
            denl = 2*(valx-A(:,j));
            %num = -2*(norm(valy-A(:,j))^2)*(valx-A(:,j));
              y = y + (2*b(j)/(den^4))*( 2*(norm(valy-A(:,j))^2)*den^2*eye(lildim,lildim)  - 2*(norm(valy-A(:,j))^2)*(den)*(valx-A(:,j))*(valx-A(:,j))' );
              %y = y + (2*b(j)/(den^4))*((valy-A(:,j))^2)*(-3*valx^2+6*valx*A(:,j)+valy^2-2*valy*A(:,j)-2*A(:,j)^2);
       
           end
           y = [y,[0;0];[0,0,0]];
           if vals>0
               y(3,3)=ro;
           end
       case 3 %avalia hess mista
            y = zeros(lildim,lildim);
            for j=1:cols
               den = norm(valx - A(:,j))^2 + norm(valy-A(:,j))^2;
               y=y+  (2*b(j)/(den^4))*( 2*den^2*(valx-A(:,j))*(valy-A(:,j))' -  2*(den*norm(valy-A(:,j))^2)*(valx-A(:,j))*(valy-A(:,j))' ); 
               %y=y-  (4*b(j)/(den^4))*( (valx-A(:,j))^3*(valy-A(:,j))); 
            end
            y = [y,[0;0];[0,0,0]];
    end % end switch ordem aplbuenovddiq
    
    case 'aplicacao'
        lildim=length(b);
        if qual==1
            valx = x0(1:lildim);
            valy = x0(lildim+1:end);
        else %ai qual=2
            valx = x0(end-lildim+1:end);
            valy = x0(1:end-lildim);
        end

     switch ordem
       case 0 %avalia f
            parenteses = (sum(valx)-1);
            soma=0;
            for j=1:lildim
               soma=soma+ (max([0;(b(j)/ro) - valx(j)]))^2;
            end
            if qual == 1
              y = valx'*(A*valy) + c*parenteses + (r/2)*(parenteses^2) + (ro/2)*soma;
            else
              y = valy'*(A*valx) + c*parenteses + (r/2)*(parenteses^2) + (ro/2)*soma;  
            end
       case 1 %avalia grad
           parenteses = (sum(valx)-1);
           base = eye(lildim,lildim); soma=zeros(lildim,1);
           for j=1:lildim
              soma = soma + (max([0;b(j)+ro*valx(j)]))*base(:,j); 
           end
           if qual == 1
               y = A*valy +(c+r*parenteses)*ones(lildim,1) - soma;
           else
               y = A'*valy +(c+r*parenteses)*ones(lildim,1) - soma;
           end
       case 2 %avalia hess
            y = r*ones(lildim,lildim)+ro*eye(lildim,lildim);
       case 3 %avalia hess mista
           if qual==1
               y = A;
           else
               y = A';
           end
    end % end switch ordem aplicacao
    
    case 'penalcubic'
        lildim=length(b);
        if qual==1
            valx = x0(1:lildim);
            valy = x0(lildim+1:end);
        else %ai qual=2
            valx = x0(end-lildim+1:end);
            valy = x0(1:end-lildim);
        end

     switch ordem
       case 0 %avalia f
            parenteses = (sum(valx)-1);
            soma=0;
            for j=1:lildim
               soma=soma+ (max([0;sign(b(j))*sqrt(b(j)) - (valx(j))/(ro)]))^3 - (b(j))^(3/2);
            end
            if qual == 1
              y = valx'*(A*valy) + c*parenteses + (r/2)*(parenteses^2) + (ro/3)*soma;
            else
              y = valy'*(A*valx) + c*parenteses + (r/2)*(parenteses^2) + (ro/2)*soma;  
            end
       case 1 %avalia grad
           parenteses = (sum(valx)-1);
           base = eye(lildim,lildim); soma=zeros(lildim,1);
           for j=1:lildim
              soma = soma + (max([0;sign(b(j))*sqrt(b(j)) - (valx(j))/(ro)])^2)*base(:,j); 
           end
           if qual == 1
               y = A*valy +(c+r*parenteses)*ones(lildim,1) - (soma)/ro;
           else
               y = A'*valy +(c+r*parenteses)*ones(lildim,1) - (soma)/ro;
           end
       case 2 %avalia hess
           vetor=zeros(lildim,1);
           for j=1:lildim
               vetor(j) = max([0;sign(b(j))*sqrt(b(j)) - (valx(j))/(ro)]);
           end
            y = r*ones(lildim,lildim)+(2/ro)*diag(vetor);
       case 3 %avalia hess mista
           if qual==1
               y = A;
           else
               y = A';
           end
    end % end switch ordem penalcubic
    
    case 'hartunglagrangeano'
        %tem q ver como vai fazer o lildim se n for diferente de m
        lildim=length(x0)/2;
        if qual==1
            valx = x0(1:lildim);
            valy = x0(lildim+1:end);
        else %ai qual=2
            valx = x0(end-lildim+1:end);
            valy = x0(1:end-lildim);
        end
        parenteses = (sum(valx)-1);

     switch ordem
       case 0 %avalia f
            soma=0;
            for j=1:lildim
               soma=soma+ exp(-valx(j)*b(j));
            end
            if qual == 1
              y = valx'*(A*valy) + c*parenteses + (r/2)*(parenteses^2) + (ro/2)*soma;
            else
              y = valy'*(A*valx) + c*parenteses + (r/2)*(parenteses^2) + (ro/2)*soma;  
            end
       case 1 %avalia grad
           if qual==1
              y = A*valy + (c + (r/2)*(parenteses))*ones(lildim,1);
              for j=1:lildim
                 y(j) = y(j) - b(j)*(ro/2)*exp(-b(j)*valx(j)); 
              end
           else
               y = (A')*valy + (c + (r/2)*(parenteses))*ones(lildim,1);
              for j=1:lildim
                 y(j) = y(j) - b(j)*(ro/2)*exp(-b(j)*valx(j)); 
              end
           end
       case 2 %avalia hess
            y = (r/2)*ones(lildim,lildim);
            for j=1:lildim
               y(j,j) = y(j,j)+(ro/2)*(b(j)^2)*exp(-b(j)*valx(j)); 
            end
                
       case 3 %avalia hess mista
          if qual==1
              y = A;
          else
              y = A';
          end
    end % end switch ordem hartunglagrangeano 
    
    case 'hartung'
        %tem q ver como vai fazer o lildim se n for diferente de m
        lildim=length(x0)/2;
        if qual==1
            valx = x0(1:lildim);
            valy = x0(lildim+1:end);
        else %ai qual=2
            valx = x0(end-lildim+1:end);
            valy = x0(1:end-lildim);
        end

     switch ordem
       case 0 %avalia f
            parenteses = (1-sum(valx(1:lildim-1)));
            soma=0; somaexp=0;
            for j=1:lildim-1
               soma=soma+ (valx(j))^2;
               somaexp = somaexp + exp(-(c^2)*valx(j));
            end
            if qual == 1
              y = valx'*(A*valy) + (b/sqrt(c))*(soma + parenteses^2) + b*(somaexp - exp(-(c^2)*parenteses));
            else
              y = valy'*(A*valx) + (b/sqrt(c))*(soma + parenteses^2) + b*(somaexp - exp(-(c^2)*parenteses));  
            end
       case 1 %avalia grad
           parenteses = (1-sum(valx(1:lildim-1))); y=zeros(lildim,1); 
           if qual==1
               v = A*valy;
           else
               v = A'*valy;
           end
           for j=1:lildim-1
              y(j) = v(j)-v(lildim) + (2*b/sqrt(c))*(valx(j) - parenteses) + b*(c^2)*(exp(-c^2*parenteses) - exp(-c^2*valx(j))); 
           end
       case 2 %avalia hess
            parenteses = (1-sum(valx(1:lildim-1))); 
            y = zeros(lildim,lildim);
            for i=1:lildim-1
               for j=1:lildim
                   if i==j
                       y(i,j) = (4*b/sqrt(c)) + (b*c^4)*(exp(-c^2*valx(i)) + exp(-c^2*parenteses) );
                   else
                       y(i,j) = (2*b/sqrt(c)) + (b*c^4)*exp(-c^2*parenteses);
                   end
               end
            end
                
       case 3 %avalia hess mista
           y = zeros(lildim,lildim);
            if qual==1
                for i=1:lildim-1
                  for j=1:lildim
                      y(i,j) = A(i,j)-A(i,lildim);
                  end
                end
            else
                for i=1:lildim-1
                  for j=1:lildim
                      y(i,j) = A(j,i)-A(lildim,i);
                  end
                end
            end
    end % end switch ordem hartung
  
    case 'testepaper'
        %nesse caso bi eh matriz, as funcoes sao f_i = xi'Axi/2 +
        %(bixj-ci)'xi
        
        n = length(A(1,:));
        
        if qual==1
            valx = x0(1:n);
            valy = x0(n+1:end);
        else %ai qual=2
            valx = x0(end-n+1:end);
            valy = x0(1:end-n);
        end

  switch ordem
    case 0 %avalia f
        y = (valx'*(A*valx))/2+((b*valy - c)')*valx;
    case 1 %avalia grad
        y = (A*valx)+b*valy - c;
    case 2 %avalia hess
        y = A;
    case 3 %avalia hess mista
        y = b;
  end % end switch ordem teste paper
        
        
  case 'problemayuan'
        %nesse caso bi eh matriz, as funcoes sao f_i = (bixj-ci)'xi 
        %com bi=+-id, e ci=0.6 ou -0.7 pra segunda
        
        if qual==1
            [n,m] = size(b);
            valx = x0(1:n);
            valy = x0(n+1:end);
        else %ai qual=2
            [m,n] = size(b);
            valx = x0(1:n);
            valy = x0(n+1:end);
        end

  switch ordem
    case 0 %avalia f
        if qual == 1
          y = ((b*valy - c)')*valx;
        else
          y = ((b*valx - c)')*valy;
        end
    case 1 %avalia grad
        if qual == 1  
          y = b*valy - c;
        else
          y = b*valx - c;
        end
    case 2 %avalia hess
        if qual == 1
            y = zeros(n,n);
        else
            y = zeros(m,m);
        end
    case 3 %avalia hess mista
        y = b;
  end % end switch ordem problemayuan
  
  
  case 'vacina'
        %nesse caso bi eh matriz, as funcoes sao f_i = +-(bixj-ci)'xi 
        %com bi=0.3, 0.2, e Ai = 0.45 
        
        if qual==1
            [n,m] = size(b);
            valx = x0(1:n);
            valy = x0(n+1:end);
        else %ai qual=2
            [m,n] = size(b);
            valx = x0(1:n);
            valy = x0(n+1:end);
        end

  switch ordem
    case 0 %avalia f
        if qual == 1
          y = -((b*valy - c)')*valx;
        else
          y = ((b*valx - c)')*valy;
        end
    case 1 %avalia grad
        if qual == 1  
          y = -b*valy + c;
        else
          y = b*valx - c;
        end
    case 2 %avalia hess
        if qual == 1
            y = zeros(n,n);
        else
            y = zeros(m,m);
        end
    case 3 %avalia hess mista
        if qual == 1
            y = -b;
        else
            y = b;
        end
  end % end switch ordem vacina
  
  
 case 'quadmisto'
     % nesse caso as duas fc sao quadraticas gerais tipo
     % f(x,y) = (xA1x)/2 + (yA2y)/2 + xB1y + c1'x + c2'y.
     
     m = length(b(1,:));
     n = length(b(:,1));
     %aqui passamos os dados tds, a matrix A eh na vdd A=[A1 A2]
            A1 = A(1:n,1:n);
            A2 = A(n+1:end,n+1:end);
            c1 = c(1:n);
            c2 = c(n+1:end);
            
            valx = x0(1:n);
            valy = x0(n+1:end);   
        
        switch ordem
    case 0 %avalia f
        if qual == 1
            y = (valx'*(A1*valx) - valy'*(A2*valy))/2 + valx'*(b*valy) + c1'*valx + c2'*valy;
        else
            y = (-valx'*(A1*valx) + valy'*(A2*valy))/2 + valx'*(b*valy) + c1'*valx + c2'*valy;
        end
    case 1 %avalia grad
        if qual==1
            %deriva no x
            y = A1*valx + b*valy + c1;
        else
            y = A2*valy + b'*valx + c2;
        end
    case 2 %avalia hess
        if qual==1
            y = A1;
        else
            y = A2;
        end
    case 3 %avalia hess mista
        if qual==1
            y = b;
        else
            y = b';
        end
  end % end switch ordem quadmisto
  
  
  case 'exponencial'
        lildim=length(b);
        %nesse caso deve-se ter n=m
        if qual==1
            valx = x0(1:lildim);
            valy = x0(lildim+1:end);
        else %ai qual=2
            valx = x0(lildim+1:end);
            valy = x0(1:lildim);
        end

     switch ordem
       case 0 %avalia f
            y = sum(A.*exp(b.*(valx.*valy)));
       case 1 %avalia grad
            y = A.*(b.*valy).*exp(b.*(valx.*valy));
       case 2 %avalia hess
            y = diag(A.*((b.^2).*(valy.^2)).*exp(b.*(valx.*valy)));
       case 3 %avalia hess mista
            y = diag(  (A.*b).*(ones(lildim,1)+ b.*(valx.*valy)).*exp(b.*(valx.*valy))  );
    end % end switch ordem exponencial
    
    case 'quadratica_yuan'
         
     n = length(x0)/2;   
     valx = x0(1:n);
     valy = x0(n+1:end);

     switch ordem
         case 0
             if qual==1
              y = (  (valx-b)'*(valx-b) + (valy-b)'*(valy-b)   )/2;
             else
              y = (  (valx-c)'*(valx-c) + (valy-c)'*(valy-c)   )/2; 
             end
             
         case 1
             if qual==1
              y = valx-b;
             else
              y = valy-c; 
             end     
             
         case 2
             y = eye(n,n);
             
         case 3
             y = zeros(n,n);
             
     end
      
    case 'quartico'
        
     n = length(x0)/2;   
     valx = x0(1:n);
     valy = x0(n+1:end);   
     
     switch ordem
         
         case 0
             if qual==1
              y = (  ((valx-b).^2)'*((valx-b).^2)  )/4; %- ( (valx)'*(valy) );
             else
              y = (  ((valy-c).^2)'*((valy-c).^2)  )/4; %- ( (valx)'*(valy) );  
             end
         case 1
             if qual==1
              y = (valx-b).^3; %- valy;
             else
              y = (valy-c).^3; %- valx;
             end
         case 2
             if qual==1
              y = 3*diag((valx-b).^2);
             else
              y = 3*diag((valy-c).^2);
             end
             
         case 3
              y = zeros(n,n); %-eye(n,n);
             
     end    
    
    case 'separavel'
        lildim=length(c);
        if qual==1
            valx = x0(1:lildim);
            valy = x0(lildim+1:end);
        else %ai qual=2
            valx = x0(lildim+1:end);
            valy = x0(1:lildim);
        end

     switch ordem
       case 0 %avalia f
            y = sum(exp(valx).*(c+valx-log(exp(valx)+exp(valy))));
       case 1 %avalia grad
            y = exp(valx).*(c+valx-log(exp(valx)+exp(valy))+ones(lildim,1)-valx./(exp(valx)+exp(valy)) );
       case 2 %avalia hess
            y = exp(valx).*(c+valx-log(exp(valx)+exp(valy))+2*ones(lildim,1)-(2*valx)./(exp(valx)+exp(valy)) - (exp(valx).*exp(valy))./((exp(valx)+exp(valy)).^2) );
            y = diag(y);
       case 3 %avalia hess mista
            y = -diag( (exp(valx).*exp(valy))./(exp(valx)+exp(valy))  + (exp(valx).*exp(valy))./((exp(valx)+exp(valy)).^2)  );
    end % end switch ordem separavel
    
    
    case 'formation'
        lildim=length(x0)/2;
        if qual==1
            valx = x0(1:lildim);
            valy = x0(lildim+1:end);
        else %ai qual=2
            valx = x0(lildim+1:end);
            valy = x0(1:lildim);
        end

     switch ordem
       case 0 %avalia f
            y = (1+c)*(valx'*valx) -2*valx'*(b+c*valy) + b'*b + c*(valy'*valy);
       case 1 %avalia grad
            y = 2*(1+c)*valx -2*(b+c*valy);
       case 2 %avalia hess
            y = 2*(1+c)*eye(lildim,lildim);
       case 3 %avalia hess mista
            y = -2*c*eye(lildim,lildim);
    end % end switch ordem formation
    
    case 'formationzero'
        lildim=length(x0)/2;
        if qual==1
            valx = x0(1:lildim);
            valy = x0(lildim+1:end);
        else %ai qual=2
            valx = x0(lildim+1:end);
            valy = x0(1:lildim);
        end

     switch ordem
       case 0 %avalia f
            y = (c)*(valx'*valx) -2*valx'*(c*valy) + c*(valy'*valy);
       case 1 %avalia grad
            y = 2*(c)*valx -2*(c*valy);
       case 2 %avalia hess
            y = 2*(c)*eye(lildim,lildim);
       case 3 %avalia hess mista
            y = -2*c*eye(lildim,lildim);
    end % end switch ordem formationzero
    
    
     case 'diferentes'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
         if qual == 1
             y = (valx-valy)^2;
         else
             y = (valx-(valy-A))^2;
         end
       case 1 %avalia grad
            if qual == 1
             y = 2*(valx-valy);
            else
             y = -2*(valx-(valy-A));
            end
       case 2 %avalia hess
            if qual == 1
             y = 2;
            else
             y = 2;
            end
       case 3 %avalia hess mista
            if qual == 1
             y = -2;
            else
             y = -2;
            end
    end % end switch ordem diferentes


case 'quadtestes1'
        %f1: A1(x-c1)^2+b1(y-c1)^2 (minimo c1,c1)
        %f2: A2(x-c2)^2 + b2(y-c2)^2 (minimo c2,c2)

        lildim=1;
        valx=x0(1:lildim);valy=x0(lildim+1:end);

  switch ordem
    case 0 %avalia f
        y = A*(valx-c)^2+b*(valy-c)^2;
    case 1 %avalia grad
        if qual == 1  
          y = 2*A*(valx-c);
        else
          y = 2*b*(valy - c);
        end
    case 2 %avalia hess
        if qual == 1
            y = 2*A;
        else
            y = 2*b;
        end
    case 3 %avalia hess mista
        y = 0;
  end % end switch ordem quadtestes1
  
  case 'quadtestes2'
        %f1: A1x^3*y
        %f2: A2y^3*x

        lildim=1;
        valx=x0(1:lildim);
        valy=x0(lildim+1:end);

  switch ordem
    case 0 %avalia f
        if qual==1
           y = A*(valx.^3)'*valy;
        else
           y = A*(valy.^3)'*valx; 
        end
    case 1 %avalia grad
        if qual == 1  
          y = 3*A*(valy)*(valx)^2;
        else
          y = 3*A*(valx)*(valy)^2;
        end
    case 2 %avalia hess
        y = 6*A*valx*valy;
    case 3 %avalia hess mista
        if qual == 1  
          y = 3*A*(valx)^2;
        else
          y = 3*A*(valy)^2;
        end
  end % end switch ordem quadtestes2
  
  case 'quadtestes3'
        %f1: A1(x-b1)
        %f2: A2(y-b2)
        
        if qual==1
          valx=x0(1);
        else
          valx=x0(2);
        end

  switch ordem
      
    case 0 %avalia f
        y = A*(valx-b);
    case 1 %avalia grad
        y = A;
    case 2 %avalia hess
        y = 0;
    case 3 %avalia hess mista
        y = 0;
  end % end switch ordem quadtestes3
  
  case 'quadtestes4'
        %f1: A1(x^3/3) + b1x
        %f2: A2(y^3/3) + b2y
        
        if qual==1
          valx=x0(1);
        else
          valx=x0(2);
        end

  switch ordem
      
    case 0 %avalia f
        y = A*(valx^3/3) + b*valx;
    case 1 %avalia grad
        y = A*(valx^2) + b;
    case 2 %avalia hess
        y = 2*A*valx;
    case 3 %avalia hess mista
        y = 0;
  end % end switch ordem quadtestes4
  
  
  case 'quadtestes5'
        %f1: sin(b1x) + A1x
        %f2: cos(b2y) + A2y
        
        if qual==1
          valx=x0(1);
        else
          valx=x0(2);
        end

  switch ordem
      
    case 0 %avalia f
        if qual==1
            y = sin(b*valx)+A*valx;
        else
            y = cos(b*valx)+A*valx;
        end
    case 1 %avalia grad
        if qual==1
            y = b*cos(b*valx)+A;
        else
            y = -b*sin(b*valx)+A;
        end
    case 2 %avalia hess
        if qual==1
            y = (b^2)*sin(b*valx);
        else
            y = -(b^2)*cos(b*valx);
        end
    case 3 %avalia hess mista
        y = 0;
  end % end switch ordem quadtestes5
  
  
  case 'exp1'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
               y =  -exp(-(valx^2+valy^2-(c^2))^2);
           else
               y =  -exp(-(valx^2+(valy-A)^2-(c^2))^2);
           end
       case 1 %avalia grad
            if qual == 1
             y = 4*valx*(valx^2+valy^2-(c^2))*exp(-(valx^2+valy^2-(c^2))^2);
            else
             y = 4*(valy-A)*(valx^2+(valy-A)^2-(c^2))*exp(-(valx^2+(valy-A)^2-(c^2))^2);
            end
       case 2 %avalia hess
            if qual == 1
             y = (8*valx^2-16*valx^2*(valx^2+valy^2-(c^2))^2 + 4*(valx^2+valy^2-(c^2)) )*exp(-(valx^2+valy^2-(c^2))^2);
            else
             y = (8*(valy-A)^2 -16*((valy-A)^2)*(valx^2+(valy-A)^2-(c^2))^2 + 4*(valx^2+(valy-A)^2-(c^2)) )*exp(-(valx^2+(valy-A)^2-(c^2))^2);
            end
       case 3 %avalia hess mista
            if qual == 1
             y = 8*valx*(valy)*(1-2*(valx^2+valy^2-(c^2))^2)*exp(-(valx^2+valy^2-(c^2))^2);
            else
             y = 8*valx*(valy-A)*(1-2*(valx^2+(valy-A)^2-(c^2))^2)*exp(-(valx^2+(valy-A)^2-(c^2))^2);
            end
    end % end switch ordem exp1
     
    
    case 'exp2'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
               y =  -exp(-(valx^2+valy^2-(A^2))^2);
           else
               y =  -exp(-(valx^2+valy^2-(c^2))^2);
           end
       case 1 %avalia grad
            if qual == 1
             y = 4*valx*(valx^2+valy^2-(A^2))*exp(-(valx^2+valy^2-(A^2))^2);
            else
             y = 4*valy*(valx^2+valy^2-(c^2))*exp(-(valx^2+valy^2-(c^2))^2);
            end
       case 2 %avalia hess
            if qual == 1
             y = (8*valx^2-16*valx^2*(valx^2+valy^2-(A^2))^2 + 4*(valx^2+valy^2-(A^2)) )*exp(-(valx^2+valy^2-(A^2))^2);
            else
             y = (8*valy^2 -16*valy^2*(valx^2+valy^2-(c^2))^2 + 4*(valx^2+valy^2-(c^2)) )*exp(-(valx^2+valy^2-(c^2))^2);
            end
       case 3 %avalia hess mista
            if qual == 1
             y = 8*valx*(valy)*(1-2*(valx^2+valy^2-(A^2))^2)*exp(-(valx^2+valy^2-(A^2))^2);
            else
             y = 8*valx*(valy)*(1-2*(valx^2+valy^2-(c^2))^2)*exp(-(valx^2+valy^2-(c^2))^2);
            end
    end % end switch ordem exp2
    
    
    case 'arcsin'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
               y = asin(valx^3/3 + valx*(valy^2-c^2));
           else
               y = asin((valy-A)^3/3 + valy*(valx^2-c^2));
           end
       case 1 %avalia grad
            if qual == 1
             y = (valx^2 + valy^2 - c^2)/sqrt(1 - (valx^3/3+valx*(valy^2-c^2))^2);
            else
             y = (valx^2 + (valy-A)^2 - c^2)/sqrt(1 - ((valy-A)^3/3+valy*(valx^2-c^2))^2);
            end
       case 2 %avalia hess
            if qual == 1
             y = (2*valx)/sqrt(1 - (valx^3/3+valx*(valy^2-c^2))^2) + ( (valx^3/3+valx*(valy^2-c^2))*(valx^2 + valy^2 - c^2) )/(1 - (valx^3/3+valx*(valy^2-c^2))^(3/2) );
            else
             y = (2*(valy-A))/sqrt(1 - ((valy-A)^3/3+valy*(valx^2-c^2))^2) + ( ((valy-A)^3/3+valy*(valx^2-c^2))*(valx^2 + (valy-A)^2 - c^2) )/(1 - ((valy-A)^3/3+valy*(valx^2-c^2))^(3/2) );
            end
       case 3 %avalia hess mista
            if qual == 1
             y = (2*valy)/sqrt(1 - (valx^3/3+valx*(valy^2-c^2))^2) + (2*valx*valy*(valx^3/3+valx*(valy^2-c^2))*(valx^2 + valy^2 - c^2))/(1 - (valx^3/3+valx*(valy^2-c^2))^(3/2));
            else
             y = (2*valx)/sqrt(1 - ((valy-A)^3/3+valy*(valx^2-c^2))^2) + (2*valx*valy*((valy-A)^3/3+valy*(valx^2-c^2))*(valx^2 + (valy-A)^2 - c^2))/(1 - ((valy-A)^3/3+valy*(valx^2-c^2))^(3/2));
            end
    end % end switch ordem arcsin
    
    
    case 'sin2'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
               y = sin(valx^3/3 + (valy^2-A)*valx);
           else
               y = sin((valy-b)^3/3 + (valx^2-c)*valy);
           end
       case 1 %avalia grad
            if qual == 1
             y = (valx^2+(valy)^2-A)*cos(valx^3/3 + (valy^2-A)*valx);
            else
             y = (valx^2+(valy-b)^2-c)*cos((valy-b)^3/3 + (valx^2-c)*valy);
            end
       case 2 %avalia hess
            if qual == 1
             y = 2*valy*cos(valx^3/3 + (valy^2-A)*valx) - ((valx^2+valy^2-A)^2)*sin(valx^3/3 + (valy^2-A)*valx);
            else
             y = 2*(valy-b)*cos((valy-b)^3/3 + (valx^2-c)*valy) - ((valx^2+(valy-b)^2-c)^2)*sin((valy-b)^3/3 + (valx^2-c)*valy);
            end
       case 3 %avalia hess mista
            if qual == 1
             y = 2*valy*(  cos(valx^3/3 + (valy^2-A)*valx) - valx*(valx^2+valy^2-A)*sin(valx^3/3 + (valy^2-A)*valx)   );
            else
             y = 2*valx*(  cos((valy-b)^3/3 + (valx^2-c)*valy) - valy*(valx^2+(valy-b)^2-c)*sin((valy-b)^3/3 + (valx^2-c)*valy)   );
            end
    end % end switch ordem sincos
    
    case 'sincos'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
               y = sin(valx^3/3 + (valy^2-A)*valx);
           else
               y = cos(valy^3/3 + (valx^2-A)*valy);
           end
       case 1 %avalia grad
            if qual == 1
             y = (valx^2+valy^2^2-A)*cos(valx^3/3 + (valy^2-A)*valx);
            else
             y = -(valx^2+valy^2^2-A)*sin(valy^3/3 + (valx^2-A)*valy);
            end
       case 2 %avalia hess
            if qual == 1
             y = 2*valx*cos(valx^3/3 + (valy^2-A)*valx) - ((valx^2+valy^2-A)^2)*sin(valx^3/3 + (valy^2-A)*valx);
            else
             y = -2*valy*sin(valx^3/3 + (valy^2-A)*valx) - ((valx^2+valy^2-A)^2)*cos(valx^3/3 + (valy^2-A)*valx);
            end
       case 3 %avalia hess mista
            if qual == 1
             y = 2*valy*(  cos(valx^3/3 + (valy^2-A)*valx) - valx*(valx^2+valy^2-A)*sin(valx^3/3 + (valy^2-A)*valx)   );
            else
             y = -2*valx*(  sin(valy^3/3 + (valx^2-A)*valy) + valy*(valx^2+valy^2-A)*cos(valy^3/3 + (valx^2-A)*valy)   ); 
            end
    end % end switch ordem sincos
   
    
    
    case 'testapol'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
               y = ((valx-valy)^2)*(valx+valy);
           else
               y = ((valx-(valy-A))^2)*(valx+(valy-A));
           end
       case 1 %avalia grad
            if qual == 1
             y = (valx-valy)*(3*valx+valy);
            else
             y = -(valx-(valy-A))*(valx+3*valy-3*A);
            end
       case 2 %avalia hess
            if qual == 1
             y = 6*valx-2*valy;
            else
             y = -2*(valx-3*valy+3*A);
            end
       case 3 %avalia hess mista
            if qual == 1
             y = -2*(valx+valy);
            else
             y = -2*(valx+(valy-A));
            end
    end % end switch ordem testapol
    
    
    case 'testacirc'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
               y = valx^3/3 - valy^3/3 + (valx*valy^2);
           else
               y = valy^3/3 - valx^3/3 + ((valx-A)^2*valy);
           end
       case 1 %avalia grad
            if qual == 1
             y = valx^2+ valy^2;
            else
             y = (valx-A)^2 + valy^2;
            end
       case 2 %avalia hess
            if qual == 1
             y = 2*valx;
            else
             y = 2*valy;
            end
       case 3 %avalia hess mista
            if qual == 1
             y = 2*valy;
            else
             y = 2*(valx-A);
            end
    end % end switch ordem testacirc
    
    case 'gradlimmin'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           y = ( (valx-c(1))^2 + (valy-c(2))^2-A   )^2;
       case 1 %avalia grad
            if qual == 1
             y = 4*(valx-c(1))*((valx-c(1))^2 + (valy-c(2))^2-A);
            else
             y = 4*(valy-c(2))*((valx-c(1))^2 + (valy-c(2))^2-A);
            end
       case 2 %avalia hess
            if qual == 1
             y = 4*( ((valx-c(1))^2 + (valy-c(2))^2-A) + 2*(valx-c(1))^2  );
            else
             y = 4*( ((valx-c(1))^2 + (valy-c(2))^2-A) + 2*(valy-c(2))^2  );
            end
       case 3 %avalia hess mista
             y = 8*(valx-c(1))*(valy-c(2));
    end % end switch ordem gradlimmin
    
    
    case 'cubicasem'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
              ra=sqrt(1+valx^2); 
              y = (-valx+(valx*ra + log(valx+ra))/2)^2;
           else
              ra=sqrt(1+valy^2); 
              y = (-valy+(valy*ra + log(valy+ra))/2)^2;
           end
       case 1 %avalia grad
            if qual == 1
                ra= sqrt(1+valx^2);
                y = -(  (1-valx^2+ra)*(valx*(ra-2) + log(valx+ra))   )/(ra);
            else
                ra= sqrt(1+valy^2);
                y = -(  (1-valy^2+ra)*(valy*(ra-2) + log(valy+ra))   )/(ra);
            end
       case 2 %avalia hess
            if qual == 1
                ra= sqrt(1+valx^2);
                y = (1/ra)*(3*valx^2*(ra-2)+4*(ra-1)+valx*(log(ra+valx))   );
            else
                ra= sqrt(1+valy^2);
                y = (1/ra)*(3*valy^2*(ra-2)+4*(ra-1)+valy*(log(ra+valy))   );
            end
       case 3 %avalia hess mista
             y = 0;
    end % end switch ordem cubicasem
    
    case 'exbueno'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
              y = valx*valy;
           else
              y = -valy;
           end
       case 1 %avalia grad
            if qual == 1 
                y = valy;
            else
                y = -1;
            end
       case 2 %avalia hess
            y = 0;
       case 3 %avalia hess mista
            if qual == 1 
                y = 1;
            else
                y = 0;
            end
    end % end switch ordem exbueno
    
    case 'exbueno2'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
              y = (valy-valx)^2/4;
           else
              y = (valy-1)^2;
           end
       case 1 %avalia grad
            if qual == 1 
                y = (valy-valx)/2;
            else
                y = 2*valy;
            end
       case 2 %avalia hess
            if qual == 1 
                y = -1/2;
            else
                y = 2;
            end
       case 3 %avalia hess mista
            if qual == 1 
                y = 1/2;
            else
                y = 0;
            end
    end % end switch ordem exbueno2
    
    case 'maxes'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual == 1
               y = (valx^3)/3 + (valy^2-4)*valx + (valy^2+1)*(  (max([valx-3;0]))^4 + (min([1-valx;0]))^4 );
           else
               y = (valy^3)/3 + ((valx-5)^2)*valy + (3/2)*valy^2 - valy + ((valx-2)^2+1)*( (min([1+valy;0]))^4 );
           end
       case 1 %avalia grad
            if qual == 1
             y = valx^2 + valy^2-4 + 4*(valy^2+1)*(  (max([valx-3;0]))^3 - (min([1-valx;0]))^3 );
            else
             y = valy^2 + (valx-5)^2 + 3*valy - 1 + 4*((valx-2)^2+1)*( (min([1+valy;0]))^3 );
            end
       case 2 %avalia hess
            if qual == 1
             y = 2*valx + 12*(valy^2+1)*(  (max([valx-3;0]))^2 + (min([1-valx;0]))^2 );
            else
             y = 2*valy +3 + 12*((valx-2)^2+1)*( (min([1+valy;0]))^2 );
            end
       case 3 %avalia hess mista
            if qual == 1
             y = 2*valy + 8*valy*(  (max([valx-3;0]))^3 - (min([1-valx;0]))^3 );
            else
             y = 2*(valx-5)+ 8*(valx-2)*( (min([1+valy;0]))^3 );
            end
    end % end switch ordem maxes
    
    
    case 'exacho'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
              y = (1-valx+valy)^2;
           else
              y = (valx^2+(valy-1)^2)/2;
           end
       case 1 %avalia grad
            if qual == 1 
                y = 2*(1-valx+valy);
            else
                y = valy-1;
            end
       case 2 %avalia hess
            if qual == 1 
                y = -2;
            else
                y = 1;
            end
       case 3 %avalia hess mista
            if qual == 1 
                y = 2;
            else
                y = 0;
            end
    end % end switch ordem exacho

    case 'exartigo1_mbom'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
               y = valx^2 + valx*valy - 5*valx;
           else
               y = 3*(valy^2/2) - valx*valy - valy;
           end
       case 1 %avalia grad
            if qual == 1
             y = 2*valx + valy - 5;
            else
             y = 3*valy - valx - 1;
            end
       case 2 %avalia hess
            if qual == 1
             y = 2;
            else
             y = 3;
            end
       case 3 %avalia hess mista
            if qual == 1
             y = 1;
            else
             y = -1;
            end
    end % end switch exartigo1_mbom
    
    
    case 'exnovo6'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
               y = (1/4)*(valx+valy)^2 + sin(valx);
           else
               y = (1/4)*(valx+valy)^2 + sin(valy);
           end
       case 1 %avalia grad
            if qual == 1
             y = (1/2)*(valx + valy) + cos(valx);
            else
             y = (1/2)*(valx + valy) + cos(valy);
            end
       case 2 %avalia hess
            if qual == 1
             y = (1/2)-sin(valx);
            else
             y = (1/2)-sin(valy);
            end
       case 3 %avalia hess mista
            if qual == 1
             y = 1/2;
            else
             y = 1/2;
            end
    end % end switch exnovo6
    
    case 'exrevisao'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
               y = valx^2 - 4*valx*valy^2;
           else
               y = (1/2)*(valy)^2 + valy;
           end
       case 1 %avalia grad
            if qual == 1
             y = 2*valx - 4*valy^2;
            else
             y = valy+1;
            end
       case 2 %avalia hess
            if qual == 1
             y = 2;
            else
             y = 1;
            end
       case 3 %avalia hess mista
            if qual == 1
             y = -8*valy;
            else
             y = 0;
            end
    end % end switch exrevisao
    
     case 'exrevisao2'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
               y = (valx^2)/2 - 2*valx*valy^2;
           else
               y = (1/2)*(valy)^2;
           end
       case 1 %avalia grad
            if qual == 1
             y = valx - 2*valy^2;
            else
             y = valy;
            end
       case 2 %avalia hess
            if qual == 1
             y = 1;
            else
             y = 1;
            end
       case 3 %avalia hess mista
            if qual == 1
             y = -4*valy;
            else
             y = 0;
            end
    end % end switch exrevisao2
 
    
    case 'exartigo2_mruim'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
               y = (valx^2)/4 + valx*valy - 5*valx;
           else
               y = (valy^2/6) - valx*valy - valy;
           end
       case 1 %avalia grad
            if qual == 1
             y = valx/2 + valy - 5;
            else
             y = valy/3 - valx - 1;
            end
       case 2 %avalia hess
            if qual == 1
             y = 1/2;
            else
             y = 1/3;
            end
       case 3 %avalia hess mista
            if qual == 1
             y = 1;
            else
             y = -1;
            end
    end % end switch exartigo2_mruim
    
    case 'exartigo3_indef'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
               y = (valx^2) + valx*valy - 5*valx;
           else
               y = -3*(valy^2/2) - valx*valy - valy;
           end
       case 1 %avalia grad
            if qual == 1
             y = 2*valx + valy - 5;
            else
             y = -3*valy - valx - 1;
            end
       case 2 %avalia hess
            if qual == 1
             y = 2;
            else
             y = -3;
            end
       case 3 %avalia hess mista
            if qual == 1
             y = 1;
            else
             y = -1;
            end
    end % end switch exartigo3_indef
    
    
    case 'exartigo4_cubic'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
               y = ((valx^3)*(valy^2))/3 + (valx^2)/2;
           else
               y = ((valx^2)*(valy^3))/3 + (valy^2)/2;
           end
       case 1 %avalia grad
            if qual == 1
             y = valx^2*valy^2 + valx;
            else
             y = valx^2*valy^2 + valy;
            end
       case 2 %avalia hess
            if qual == 1
             y = 2*valx*valy^2 + 1;
            else
             y = 2*valy*valx^2 + 1;
            end
       case 3 %avalia hess mista
            if qual == 1
             y = 2*valy*valx^2;
            else
             y = 2*valx*valy^2;
            end
    end % end switch exartigo4_cubic
    
    
    case 'exartigo5_futebol'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
               y = -valx*(0.6-valy);
           else
               y = valy*(0.7-valx);
           end
       case 1 %avalia grad
            if qual == 1
             y = valy - 0.6;
            else
             y = 0.7 - valx;
            end
       case 2 %avalia hess
            if qual == 1
             y = 0;
            else
             y = 0;
            end
       case 3 %avalia hess mista
            if qual == 1
             y = 1;
            else
             y = -1;
            end
    end % end switch exartigo5_futebol
    
    
    case 'exartigo6_vacina'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
               y = valx*(0.45*valy-0.3);
           else
               y = -valy*(0.45*valx-0.2);
           end
       case 1 %avalia grad
            if qual == 1
             y = 0.45*valy - 0.3;
            else
             y = 0.2 - 0.45*valx;
            end
       case 2 %avalia hess
            if qual == 1
             y = 0;
            else
             y = 0;
            end
       case 3 %avalia hess mista
            if qual == 1
             y = 0.45;
            else
             y = -0.45;
            end
    end % end switch exartigo6_vacina
    
     case 'futebolzero'      
            valx = x0(1);
            valy = x0(2);

     switch ordem
       case 0 %avalia f
           if qual==1
               y = -valx*(0.0-valy);
           else
               y = valy*(0.0-valx);
           end
       case 1 %avalia grad
            if qual == 1
             y = valy - 0.0;
            else
             y = 0.0 - valx;
            end
       case 2 %avalia hess
            if qual == 1
             y = 0;
            else
             y = 0;
            end
       case 3 %avalia hess mista
            if qual == 1
             y = 1;
            else
             y = -1;
            end
    end % end switch futebolzero
    
  

end %end switch

end
        