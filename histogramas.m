

%construindo os bars
mybars = ones(3,100);
for i=1:10
    for j=1:10
        %if valoresf1newton(i,j)>3
        %    mybars(2,i+(j-1)*10) = 5;
        %end
        %if valoresf1yuan(i,j)>300
        %    mybars(3,i+(j-1)*10) = 10;
        %end
    end
end

contnewton5 = 0;
contyuan10 = 0;
for i=1:100
    %if mybars(2,i)==5
    %   contnewton5 = contnewton5+1; 
    %end
    %if mybars(3,i)==10
    %   contyuan10 = contyuan10+1; 
    %end
end
y = [100,13,100-88,3; 0,100-13,0,0;  0,0,88,97 ];
y = [100,0,0,13;87,0,12,0;88,3,0,97];
%bar(myx,y);
str={'EQUILIBRIUM';'NON-EQUILIBRIUM KKT';'DIVERGENCE'};
set(gca, 'XTickLabel',str, 'XTick',1:numel(str))
ydoze=[y(1,1:4),y(2,1:4),y(3,1:4)];

%text(1:length(ydoze),ydoze,num2str(ydoze),'vert','bottom','horiz','center'); 

dataSample = reshape(ydoze,3,4)';
figure
catStrArray = {'EQUILIBRIUM';'NON-EQUILIBRIUM KKT';'DIVERGENCE'};
catArray = categorical(catStrArray);       
catArray = reordercats(catArray,catStrArray);
bar(catArray,dataSample')
nModel = size(dataSample,1);
nCat = size(dataSample,2);
xPosAmpl = 0.3682626-0.3298725*exp(-0.407004*(nModel-1)); % position amplitude
xPosInc = 2*xPosAmpl/(nModel-1);
modelNames = [];
for idxModel=1:nModel
    bar_xPos = 1:nCat;
    if nModel~=1
        bar_xPos = bar_xPos-xPosAmpl+(idxModel-1)*xPosInc;
    end
    text(bar_xPos,dataSample(idxModel,:),num2str(dataSample(idxModel,:)',...
        '%0.0f'),'vert','bottom','horiz','center'); 
    modelNames{idxModel}=sprintf('model%d',idxModel);
end
str={'EQUILIBRIUM';'NON-EQUILIBRIUM KKT';'DIVERGENCE'};
set(gca, 'XTickLabel',str, 'XTick',1:numel(str))
legend(modelNames)
