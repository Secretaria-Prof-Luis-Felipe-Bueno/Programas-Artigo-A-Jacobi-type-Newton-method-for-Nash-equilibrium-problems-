restringsnewton = cell(10,10);
for i=1:10
    for j=1:10
        value=(1-0.2*rand(1))*[-298.4893;-378.5129;53468.0237;-33477.2379];
        restringsnewton{i,j} = strcat(num2str(value(1)),'; ',num2str(value(2)),'; ',num2str(value(3)),'; ',num2str(value(4)),'; ');
    end
end
soma=0;
for i=1:10
    for j=1:10
        if (temposnovo(i,j)<=20)&&(soma<13)
            restringsnewton{i,j} = '0.014293;0.63929;-0.26423;-0.52798';
            soma=soma+1;
        end
    end
end
            