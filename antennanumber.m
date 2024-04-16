function Nt = antennanumber(Delta_d,r,phi)

tt = atand(Delta_d/2/r);  % tt必为正
tt(find(tt<0)) = 180+tt(find(tt<0)); % -90~90转换为0~180

Nt = floor(50.5/( tt*sind(phi) ) );

% if Nt > 15
%     Nt = 15;
% end



