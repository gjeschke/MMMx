clear options
options.altlocs = false;
options.Bfactor = false;
options.element = false;
% id = cx_get_pdb('1a6m',options);

address = '(A)4.CG:B';
[xyz,exceptions1] = cx_get_location(id,address,'coor');
fprintf(1,'Coordinates of %s : %8.3f,%8.3f,%8.3f ?\n',address,xyz);
xyz = xyz + [1,-2,-4];

[attribute,exceptions] = cx_set_location(id,address,'coor',xyz);

fprintf(1,'Response: %s\n',attribute);

[attribute,exceptions2] = cx_get_location(id,address,'coor');
fprintf(1,'Coordinates of %s were changed to: %8.3f,%8.3f,%8.3f ?\n',address,xyz);
