similarities = load('hnRNPA1_similarities.dat');

figure(1); clf; hold on
my_map = flipud(colormap);
colormap(my_map);
image(similarities,'CDataMapping','scaled');
curr_axis = gca;
set(curr_axis,'YDir','normal');
colorbar;
axis tight
xlabel('Ensemble number');
ylabel('Ensemble number');
axis equal
title('Comparison of hnRNPA1 ensembles');
set(gca,'FontSize',12);