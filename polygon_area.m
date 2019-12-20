%% bersih-bersih
clear
clc

%% membuka file
file_id=fopen('vertices_bangka.txt','r');

% koordinat titik-titik pada pantai Pulau Bangka
% satuan angka adalah meter
vertices_bangka=transp(fscanf(file_id,'%f %f',[2 Inf]));
% rubah koordinat universal ke koordinat lokal
vertices_bangka(:,1)=vertices_bangka(:,1)-min(vertices_bangka(:,1));
vertices_bangka(:,2)=vertices_bangka(:,2)-min(vertices_bangka(:,2));

% tutup file
fclose(file_id);

% nomor titik pada pantai Pulau Bangka
faces_bangka=horzcat(1:32,1);

%% perhitungan luas Pulau Bangka dalam kilometer persegi
prod1=1e-6*sum(vertices_bangka(:,1).*circshift(vertices_bangka(:,2),1));
prod2=1e-6*sum(vertices_bangka(:,2).*circshift(vertices_bangka(:,1),1));
area_bangka=0.5*(prod1-prod2);

%% plotting
patch('Faces',faces_bangka,'Vertices',vertices_bangka./1e+3)
hold on
scatter(vertices_bangka(:,1)./1e+3,vertices_bangka(:,2)./1e+3)
txtplot=cellstr(num2str(transp(1:32)));
text((vertices_bangka(:,1)./1e+3)+1, (vertices_bangka(:,2)./1e+3)+1, txtplot);
title({'Pulau Bangka',['luasnya ',num2str(area_bangka,'%.0f'),' km^2.']})
hold off
