function x = Cell_4p(nelx,nely,x_val)

% the parameter value can be from 0 to 1
[n,m]=find(x_val>1);
x_val(n,m)=1;

kv=x_val(1);
kh=x_val(2);
kd1=x_val(3);
kd2=x_val(4);

tv=kv*nelx/2;
th=kh*nelx/2;
td1=kd1*nelx/2; % is equal to the thicnkess *sqrt(2)/2
td2=kd2*nelx/2; % is equal to the thicnkess *sqrt(2)/2

xv=zeros(nely,nelx);
xh=zeros(nely,nelx);
xd1=zeros(nely,nelx);
xd2=zeros(nely,nelx);

for elx=1:nelx
    for ely=1:nely
        
        % Vertical beams
        if elx<=tv
            xv(ely,elx)=1;
        elseif elx<tv+1 && elx>tv
            xv(ely,elx)=tv-elx+1;
        end           
        
        % Horizontal beams
        if ely<=th
            xh(ely,elx)=1;
        elseif ely<th+1 && ely>th
            xh(ely,elx)=th-ely+1;
        end
        
        % Diagonal beam up-left corner to down-right corner        
        if td1<1
           if ely==elx
               xd1(ely,elx)=2*td1*(1-td1/2);
           elseif abs(ely-elx)==1 || nelx-abs(ely-elx)==1
               xd1(ely,elx)=1/2*(td1)^2;
           end
        else
            if ely==elx || abs(ely-elx)+1<=td1 || nelx+1-abs(ely-elx)<=td1
               xd1(ely,elx)=1;
            elseif abs(ely-elx)-1<=td1 && abs(ely-elx)>td1
               xd1(ely,elx)=1/2*(td1-abs(ely-elx)+1)^2;
            elseif nelx-abs(ely-elx)-1<=td1 && nelx-abs(ely-elx)>td1
               xd1(ely,elx)=1/2*(td1-nelx+abs(ely-elx)+1)^2;
            elseif abs(ely-elx)<=td1 && abs(ely-elx)+1>td1 
               xd1(ely,elx)=1/2+(td1-abs(ely-elx))*(1-td1+abs(ely-elx))+1/2*(td1-abs(ely-elx))^2;
            elseif nelx-abs(ely-elx)<=td1 && nelx-abs(ely-elx)+1>td1
               xd1(ely,elx)=1/2+(td1-nelx+abs(ely-elx))*(1-td1+nelx-abs(ely-elx))+1/2*(td1-nelx+abs(ely-elx))^2;
            end
        end
        if ely-elx==nely/2 || elx-ely==nelx/2
            xd1(ely,elx)=2*xd1(ely,elx);
            if xd1(ely,elx)>1
                xd1(ely,elx)=1;
            end
        end
        
        % Diagonal beam up-right corner to down-left corner to fliplr
        if td2<1
           if ely==elx
               xd2(ely,elx)=2*td2*(1-td2/2);
           elseif abs(ely-elx)==1 || nelx-abs(ely-elx)==1
               xd2(ely,elx)=1/2*(td2)^2;
           end
        else
            if ely==elx || abs(ely-elx)+1<=td2 || nelx+1-abs(ely-elx)<=td2
               xd2(ely,elx)=1;
            elseif abs(ely-elx)-1<=td2 && abs(ely-elx)>td2 
               xd2(ely,elx)=1/2*(td2-abs(ely-elx)+1)^2;
            elseif nelx-abs(ely-elx)-1<=td2 && nelx-abs(ely-elx)>td2
               xd2(ely,elx)=1/2*(td2-nelx+abs(ely-elx)+1)^2;
            elseif abs(ely-elx)<=td2 && abs(ely-elx)+1>td2 
               xd2(ely,elx)=1/2+(td2-abs(ely-elx))*(1-td2+abs(ely-elx))+1/2*(td2-abs(ely-elx))^2;
            elseif nelx-abs(ely-elx)<=td2 && nelx-abs(ely-elx)+1>td2
               xd2(ely,elx)=1/2+(td2-nelx+abs(ely-elx))*(1-td2+nelx-abs(ely-elx))+1/2*(td2-nelx+abs(ely-elx))^2;
            end
        end
        if ely-elx==nely/2 || elx-ely==nelx/2
            xd2(ely,elx)=2*xd2(ely,elx);
            if xd2(ely,elx)>1
                xd2(ely,elx)=1;
            end
        end
    end
end

if (-1)^nelx==-1
    xv(:,nelx/2+1/2)=2*xv(:,nelx/2+1/2);
end
if (-1)^nely==-1
    xh(nely/2+1/2,:)=2*xh(nely/2+1/2,:);
end

xv = max(xv,fliplr(xv));
xh = max(xh,flipud(xh));
xd2 = fliplr(xd2);

x=max(xd2,max(xd1,max(xv,xh)));

% figure
% colormap(gray); imagesc(1-x); caxis([0 1]); axis equal; axis off; drawnow;

end
