function theta = arctand(x,y)

theta = atand(y./x);            % 真实散射点角度状态存储矩阵 
theta(find(theta<0)) = 180+theta(find(theta<0)); % -90~90转换为0~180