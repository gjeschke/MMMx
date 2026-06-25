function GO_clicked(cb,~)
% gobject_clicked(cb,eventdata)
%
% Function executed when a user clicks on a sphere that represents one
% protein in the proteome

web(cb.UserData.url, '-browser');
fprintf(1,'--- GO %i: %s ---\n',cb.UserData.id,cb.UserData.term);
fprintf(1,'%i proteins assigned\n',cb.UserData.proteins);
fprintf(1,'fIDR   mean %.3f, std. dev. %.3f\n',cb.UserData.mean_p_IDR,cb.UserData.std_p_IDR);
fprintf(1,'RMSF   mean %.1f Å, std. dev. %.1f Å\n',cb.UserData.mean_RMSF,cb.UserData.std_RMSF);
fprintf(1,'SSP    mean %.3f, std. dev. %.3f\n',cb.UserData.mean_SSP,cb.UserData.std_SSP);
fprintf(1,'IFRs   mean %.1f, std. dev. %.1f\n',cb.UserData.mean_nd,cb.UserData.std_nd);
fprintf(1,'length mean %i, std. dev. %i\n',round(cb.UserData.mean_n),round(cb.UserData.std_n));

