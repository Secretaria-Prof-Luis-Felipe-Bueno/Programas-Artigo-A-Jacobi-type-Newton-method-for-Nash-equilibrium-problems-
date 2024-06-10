function [y1,y2] = gradrest(x0,mystr,qual,bdarest1,bdarest2,rest1,rest2,H1,H2,lamb1,lamb2)
%retorna o par pra ficar mais facil no alg principal
if nargin<8
    H1 = 0; H2=0; lamb1=0; lamb2=0;
end

if strcmp(mystr,'apalfa1')
         r1 = 1; r2=1;
        ro1 = 1; ro2=1;
        [A1,A2,b1,b2,c1,c2,~,~,~,~,~] = montadados(mystr,3,3); 
end

n=length(rest1);
m=length(rest2);
numbless=1;
x0c = [2.6;0;-134;0];                                                                                                                                                 
%test_script

if qual==0 %gradlam
    y1 = -bdarest1 + rest1*x0(1:n);
    y2 = -bdarest2 + rest2*x0(n+1:n+m);
    if (strcmp(mystr,'tenta1'))||(strcmp(mystr,'exmisto'))%||(strcmp(mystr,'aplbuenovdd'))
        y1 = x0(1) + sqrt(1+x0(2)^2);
        y2 = x0(3) + sqrt(1+x0(4)^2);
    end
    if (strcmp(mystr,'aplbuenovdd'))
        y1 = x0(1) + sqrt(1+x0(2)^2) - 1.1727;
        y2 = x0(3) + sqrt(1+x0(4)^2) - 1.3994;
    end
    if (strcmp(mystr,'aplbuenovddiq'))
        %y1 = x0(1) + sqrt(1+x0(2)^2+ numbless)  + x0(3);
        %y2 = x0(4) + sqrt(1+x0(5)^2+ numbless)  + x0(6);
        y1 = (x0(1)-x0c(1))^2 + (x0(2)-x0c(2))^2 + x0(3)-numbless;
        y2 = (x0(3)-x0c(3))^2 + (x0(4)-x0c(4))^2 + x0(6)-numbless;
    end
    if strcmp(mystr,'apalfa1')     
            y1 = Feval(x0,A1,b1,c1,1,mystr,1,r1,ro1);
            y1=y1(3);
            y2 = Feval(x0,A2,b2,c2,1,mystr,2,r2,ro2);
            y2=y2(3);
    end
    
end
if qual==1 %rest
    if (strcmp(mystr,'tenta1'))||(strcmp(mystr,'exmisto'))||(strcmp(mystr,'aplbuenovdd'))
        y1 = [1, x0(2)/sqrt(1+x0(2)^2)];
        y2 = [1, x0(4)/sqrt(1+x0(4)^2)];
    else
        if strcmp(mystr,'aplbuenovddiq')
            %y1 = [1, x0(2)/sqrt(1+x0(2)^2), 1];
            %y2 = [1, x0(4)/sqrt(1+x0(4)^2), 1];
            y1 = [2*(x0(1)-x0c(1)),2*(x0(2)-x0c(2)),1];
            y2 = [2*(x0(3)-x0c(3)),2*(x0(4)-x0c(4)),1];
        else
            y1 = rest1;
            y2 = rest2;
        end
        if strcmp(mystr,'apalfa1')
            y1 = Feval(x0,A1,b1,c1,2,mystr,1,r1,ro1);
            %y1 = y1(3,1:2);
            y1 = y1(3,:);
            y2 = Feval(x0,A2,b2,c2,2,mystr,2,r2,ro2);
            %y2 = y2(3,1:2);
            y2 = y2(3,:);
        end
    end
end
if qual==2 %hess L
    if (strcmp(mystr,'tenta1'))||(strcmp(mystr,'exmisto'))||(strcmp(mystr,'aplbuenovdd'))
        y1 = H1 + [0,0;0,lamb1/((1+x0(2)^2)^(3/2))];
        y2 = H2 + [0,0;0,lamb2/((1+x0(4)^2)^(3/2))];
    else
        if strcmp(mystr,'aplbuenovddiq')
            %y1 = H1 + [0,0,0;0,lamb1/((1+x0(2)^2)^(3/2)),0;0,0,0];
            %y2 = H2 + [0,0,0;0,lamb2/((1+x0(5)^2)^(3/2)),0;0,0,0];
            y1 = H1 + [2*lamb1,0,2*lamb1;0,0,0;0,0,0];
            y2 = H1 + [2*lamb2,0,2*lamb2;0,0,0;0,0,0];
        else
            y1=H1;
            y2=H2;
        end
    end
end
if qual==3 %hess mista 
     [A1,A2,b1,b2,c1,c2,~,~,~,~,~] = montadados(mystr,3,3);
    if (strcmp(mystr,'tenta1'))||(strcmp(mystr,'exmisto'))
        lildim = length(b1);
        y1 = zeros(lildim,length(x0)-lildim);
        lildim = length(b2);
        y2 = zeros(lildim,length(x0)-lildim);
        
        %y1 = H1 + [0,0;0,lamb1/((1+x0(2)^2)^(3/2))];
        %y2 = H2 + [0,0;0,lamb2/((1+x0(4)^2)^(3/2))];
    else
        lildim = length(b1);
        y1 = zeros(lildim,length(x0)-lildim);
        lildim = length(b2);
        y2 = zeros(lildim,length(x0)-lildim);
    end
end

end

