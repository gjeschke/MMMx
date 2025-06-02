function h = plot_pae(pae,options)

if ~exist("options",'var')
    options.title = '';
end

if ~isfield(options,'quantity')
    options.quantity = 'PAE';
end

if ~isfield(options,'colors')
    options.colors = 256;
end

if isfield(options,'clipping') && options.clipping
    diff = max(max(round(pae)-pae));
    if diff == 0
        clipping = 31;
        options.colors = 32;
    else
        diff = max(max(round(4*pae)-4*pae));
        if diff == 0
            clipping = 31.75;
            options.colors = 128;
        else
            clipping = max(max(ceil(pae)));
            options.colors = clipping;
        end
    end
else
    options.clipping = false;
end

h = figure;

image(pae,'CDataMapping','scaled');

mymap = parula(options.colors);
curr_axis = gca;
curr_axis.YDir = 'normal';
if options.clipping
    mymap(options.colors,:) = [1,1,1]; % clipped values should be shown white
    curr_axis.CLim = [0,clipping];
end
colormap(mymap);
c = colorbar;
ylabel(c, sprintf('%s [Å]',options.quantity));
title(options.title);
axis tight
xlabel('Residue number');
ylabel('Residue number');
axis equal
set(gca,'FontSize',12);
