function theta = arctand(x,y)

theta = atand(y./x);            % ��ʵɢ���Ƕ�״̬�洢���� 
theta(find(theta<0)) = 180+theta(find(theta<0)); % -90~90ת��Ϊ0~180