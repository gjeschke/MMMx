function conformer_clicked(cb,~)
% conformer_clicked(cb,~)
%
% Function executed when a user clicks on a marker that represents one
% conformer in an ensemble

fprintf(1,'--- Conformer %i with population %.4f ---\n',cb.UserData.conformer,cb.UserData.population);
if isfield(cb.UserData,'Rg')
    fprintf(1,'Radius of gyration: %.1f Å\n',cb.UserData.Rg);
end
if isfield(cb.UserData,'asphericity')
    fprintf(1,'Asphericity: %.3f\n',cb.UserData.asphericity);
end
