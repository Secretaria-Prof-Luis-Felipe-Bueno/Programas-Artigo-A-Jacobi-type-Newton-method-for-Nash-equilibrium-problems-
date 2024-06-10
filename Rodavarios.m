%aplbueno eh 2 por 2
clear x0;


total = 3;

fazuni=1;
if fazuni==0
    
    
    x0=[-1;-1;0;-1;-1;0]; comeco=-2; retr = 1; %x0=[-1;-1;-1;-1];
    %Geral
    %salvastr = strcat('.\Dados_apl1\newton',num2str(-3),num2str(-3),num2str(-3),num2str(-3));
    %save(salvastr,'x0','k','erro');
    
    for indi=1:total
        for indj=1:total
            for indk=1:total
                for indl=1:total
                    
                    clear x0;
                    disp([indi,indj,indk,indl]);
                    %x0=comeco*[1;1;0;1;1;0]+retr*[indi;indj;0;indk;indl;0];
                    x0=comeco*[1;1;1;1]+retr*[indi;indj;indk;indl];
                    %Geral_rest_simpl
                    Geral
                    
                    %teste1 = Feval(x0,A1,b1,c1,2, mystr,1,r1,ro1);
                    %teste1 = teste1(1:2,1:2);
                    %teste1 = teste1+ [0,0;0,lamb1/((1+x0(2)^2)^(3/2))];
                    %teste2 = Feval(x0,A2,b2,c2,2, mystr,2,r2,ro2);
                    %teste2 = teste2(1:2,1:2);
                    %teste2 = teste2+ [0,0;0,lamb2/((1+x0(5)^2)^(3/2))];
                    %teste=[eigs(teste1),eigs(teste2)];
                    %disp(teste);
                    
                    teste1 = Feval(x0,A1,b1,c1,2, mystr,1,r1,ro1);
                    teste2 = Feval(x0,A2,b2,c2,2, mystr,2,r2,ro2);
                    teste=[eigs(teste1),eigs(teste2)];
                    disp(teste);
                    
                    salvastr = strcat('.\Dados_apl1\newton',num2str(comeco+retr*indi),num2str(comeco+retr*indj),num2str(comeco+retr*indk),num2str(comeco+retr*indl));
                    save(salvastr,'x0','k','erro');
                    %pause
                    
                    
                end
            end
        end
    end
    
    
else
    
    
    x0=[-1;-1]; comeco=-2; retr = 1; %x0=[-1;-1;-1;-1];
    %Geral
    %salvastr = strcat('.\Dados_apl1\newton',num2str(-3),num2str(-3),num2str(-3),num2str(-3));
    %save(salvastr,'x0','k','erro');
    fazertotal = 1;
    %for indi=1:total
       % for indj=1:total
       while(fazertotal<100)
                    clear x0;
                    %disp([indi,indj]);
                    %x0=comeco*[1;1;0;1;1;0]+retr*[indi;indj;0;indk;indl;0];
                    %x0=comeco*[1;1]+retr*[indi;indj];
                    x0 = 10*[rand(1);rand(1)];
                    disp('x0: '); disp(x0);
                    %pause
                    %Geral_rest_simpl
                    Geral
                    
                    teste1 = Feval(x0,A1,b1,c1,2, mystr,1,r1,ro1);
                    teste2 = Feval(x0,A2,b2,c2,2, mystr,2,r2,ro2);
                    teste=[eigs(teste1),eigs(teste2)];
                    disp(teste);
                    if (eigs(teste1)<0)||(eigs(teste2)<0)
                       pause 
                    end
                    fazertotal=fazertotal+1;
       end %end while
   
       % end
    %end
    
    
end


%quadfile = matfile(mystr);
%A1 = quadfile.A1;