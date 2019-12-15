%% bersih-bersih
clear
clc

%% membuka file
file_id=fopen('vertices_bangka.txt','r');

% koordinat titik-titik pada pantai Pulau Bangka
% satuan angka adalah meter
vertices_bangka=transp(fscanf(file_id,'%f %f',[2 Inf]));

% tutup file
fclose(file_id);

% nomor titik pada pantai Pulau Bangka
faces_bangka=horzcat(1:32,1);

%% perhitungan luas
prod1=0;
prod2=0;
iter=length(vertices_bangka);
for i=1:iter
    if i==iter
        j=1;
    else
        j=i+1;
    end
    
    % ubah jarak dari meter ke km
    prod1=prod1+(vertices_bangka(i,1)*vertices_bangka(j,2)/1e+6); 
    prod2=prod2+(vertices_bangka(i,2)*vertices_bangka(j,1)/1e+6);
end

% luas Pulau Bangka dalam kilometer persegi
area_bangka=0.5*(prod1-prod2);
