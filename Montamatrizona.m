%%%%%%%%%%
%testepaper
%%%%%%%%%

%1 x 1

bigkstruct = {'nosso' 'dois' 'intern' 'jacobi' 'yuan';
0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 
0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 
0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 
0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 
0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0 
};
bigtstruct = bigkstruct; bigestruct = bigkstruct;

%for i=1:30
   filestr = strcat('.\Salvar\resultsproblemayuan1130'); 
   matrizdavez = matfile(filestr);
for i=1:30
   bigtstruct{i+1,1} = matrizdavez.tvec(i,1); 
   bigkstruct{i+1,1} = matrizdavez.kvec(i,1); 
   bigestruct{i+1,1} = matrizdavez.evec(i,1); 
   bigtstruct{i+1,2} = matrizdavez.tvecdois(i,1); 
   bigkstruct{i+1,2} = matrizdavez.kvecdois(i,1);
   bigestruct{i+1,2} = matrizdavez.evecdois(i,1);
   bigtstruct{i+1,3} = matrizdavez.tvecintern(i,1); 
   bigkstruct{i+1,3} = matrizdavez.kvecintern(i,1);
   bigestruct{i+1,3} = matrizdavez.evecintern(i,1);
   bigtstruct{i+1,4} = matrizdavez.tvecjacobi(i,1); 
   bigkstruct{i+1,4} = matrizdavez.kvecjacobi(i,1); 
   bigestruct{i+1,4} = matrizdavez.evecjacobi(i,1); 
   bigtstruct{i+1,5} = matrizdavez.tvecyuan(i,1); 
   bigkstruct{i+1,5} = matrizdavez.kvecyuan(i,1);
   bigestruct{i+1,5} = matrizdavez.evecyuan(i,1);
   if bigkstruct{i+1,5} > 100
       bigkstruct{i+1,5} = 1000;
   end
end

% Plot_Perfomance_Profiles(bigkstruct)
% Plot_Perfomance_Profiles(bigtstruct)
 Plot_Perfomance_Profiles(bigestruct)