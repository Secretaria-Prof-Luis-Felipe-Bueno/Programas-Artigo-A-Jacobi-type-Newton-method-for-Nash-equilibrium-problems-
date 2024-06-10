function [y] = jacaplbueno(A,b,qual,x0)


lildim=length(b);
        cols = length(A(1,:));
        if qual==1
            valx = x0(1:lildim);
            valy = x0(lildim+1:end);
        else %ai qual=2
            valx = x0(end-lildim+1:end);
            valy = x0(1:end-lildim);
        end
        
y = zeros(lildim,1);
           for j=1:cols
            den = norm(valx - A(:,j))^2 + norm(valx-valy)^2;
            y = y + (2*den*(valx-A(:,j)) - 2*norm(valx-A(:,j))^2*(2*valx-(valy+A(:,j))))/den^2;
           end


end

