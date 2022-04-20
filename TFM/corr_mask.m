function [m_corr]=corr_mask(regx,regy,mask)
%corrige le décalge sur une image prise AVANT dissolution des cellules

deltax = round(regx);
deltay = round(regy);

%correction du décalage
if (deltax>0) && (deltay>0)
        %n_corr=nonstressed(deltay+1:size(nonstressed,1),deltax+1:size(nonstressed,2));
       m_corr=mask(1:size(mask,1)-deltay,1:size(mask,2)-deltax);

%         s_corr=stressed(1:size(stressed,1)-deltay,1:size(stressed,2)-deltax);
%         p_corr=widefield(1:size(stressed,1)-deltay,1:size(stressed,2)-deltax);
%         t_corr=pattern(1:size(stressed,1)-deltay,1:size(stressed,2)-deltax);
    elseif (deltax>0) && (deltay<=0)
        %n_corr=nonstressed(1:size(nonstressed,1)-abs(deltay),deltax+1:size(nonstressed,2));
        m_corr=mask(abs(deltay)+1:size(mask,1),1:size(mask,2)-deltax);
        %s_corr=stressed(abs(deltay)+1:size(stressed,1),1:size(stressed,2)-deltax);
        %p_corr=widefield(abs(deltay)+1:size(stressed,1),1:size(stressed,2)-deltax);
        %t_corr=pattern(abs(deltay)+1:size(stressed,1),1:size(stressed,2)-deltax);
    elseif (deltax<=0) && (deltay>0)
        %n_corr=nonstressed(deltay+1:size(nonstressed,1),1:size(nonstressed,2)-abs(deltax));
        m_corr=mask(1:size(mask,1)-deltay,abs(deltax)+1:size(mask,2));

        %s_corr=stressed(1:size(stressed,1)-deltay,abs(deltax)+1:size(stressed,2));
        %p_corr=widefield(1:size(stressed,1)-deltay,abs(deltax)+1:size(stressed,2));
        %t_corr=pattern(1:size(stressed,1)-deltay,abs(deltax)+1:size(stressed,2));
    else
        %s_corr=stressed(abs(deltay)+1:size(stressed,1),abs(deltax)+1:size(stressed,2));
        %n_corr=nonstressed(1:size(nonstressed,1)-abs(deltay),1:size(nonstressed,2)-abs(deltax));
        m_corr=mask(abs(deltay)+1:size(mask,1),abs(deltax)+1:size(mask,2));

         %p_corr=widefield(abs(deltay)+1:size(stressed,1),abs(deltax)+1:size(stressed,2));
         %t_corr=pattern(abs(deltay)+1:size(stressed,1),abs(deltax)+1:size(stressed,2));
end