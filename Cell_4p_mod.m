function x = Cell_4p_mod(nelx,nely,x_val)

% the parameter value can be from 0 to 1
[n,m]=find(x_val>1);
x_val(n,m)=1;

kv=x_val(1);
kh=x_val(2);
kd1=x_val(3);
kd2=x_val(4);

tv=kv*nelx/2;
th=kh*nelx/2;
td1=kd1*sqrt(2)*nelx/4; % is equal to the thicnkess *sqrt(2)
td2=kd2*sqrt(2)*nelx/4; % is equal to the thicnkess *sqrt(2)

xv=zeros(nely,nelx);
xh=zeros(nely,nelx);
xd1=zeros(nely,nelx);
xd1c1=zeros(nely,nelx);
xd1c2=zeros(nely,nelx);
xd2=zeros(nely,nelx);
xd2c1=zeros(nely,nelx);
xd2c2=zeros(nely,nelx);

D=@(x,y,kx,ky) 1/sqrt(2)*abs((x-kx)-(y-ky));

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
        
        v=[ely-1,elx;
            ely-1,elx-1;
            ely,elx-1;
            ely,elx];
        
        D1=D(v(1,2),v(1,1),0,0);
        D2=D(v(2,2),v(2,1),0,0);
        D3=D(v(3,2),v(3,1),0,0);
        D4=D(v(4,2),v(4,1),0,0);
        Dvec=[D1,D2,D3,D4];
        n_zeros=sum(Dvec==0);        
        Dmin=min(Dvec);
        Dmax=max(Dvec);
        
        D1c1=D(v(1,2),v(1,1),0,nely);
        D2c1=D(v(2,2),v(2,1),0,nely);
        D3c1=D(v(3,2),v(3,1),0,nely);
        D4c1=D(v(4,2),v(4,1),0,nely);
        Dvecc1=[D1c1,D2c1,D3c1,D4c1];
        nc1_zeros=sum(Dvecc1==0);        
        Dminc1=min(Dvecc1);
        Dmaxc1=max(Dvecc1);
        
        D1c2=D(v(1,2),v(1,1),nelx,0);
        D2c2=D(v(2,2),v(2,1),nelx,0);
        D3c2=D(v(3,2),v(3,1),nelx,0);
        D4c2=D(v(4,2),v(4,1),nelx,0);
        Dvecc2=[D1c2,D2c2,D3c2,D4c2];
        nc2_zeros=sum(Dvecc2==0);        
        Dminc2=min(Dvecc2);
        Dmaxc2=max(Dvecc2);
        
        % Diagonal beam up-left corner to down-right corner     
        if td1<sqrt(2)/2
            if n_zeros==2
                xd1(ely,elx)=2*td1*(sqrt(2)-td1);
            elseif n_zeros==1
                xd1(ely,elx)=td1^2;
            end
        else
            if td1>=Dmax
                xd1(ely,elx)=1;
            elseif  td1>Dmin && td1<=Dmin+sqrt(2)/2
                xd1(ely,elx)=(td1-Dmin)^2;
            elseif td1>Dmin+sqrt(2)/2 && td1<Dmax
                xd1(ely,elx)=1-(Dmax-td1)^2;
            end
        end
        if td1<sqrt(2)/2 && nc1_zeros==1
            xd1c1(ely,elx)=td1^2;
        else
            if td1>=Dmaxc1
                xd1c1(ely,elx)=1;           
            elseif  td1>Dminc1 && td1<=Dminc1+sqrt(2)/2
                xd1c1(ely,elx)=(td1-Dminc1)^2;       
            elseif td1>Dminc1+sqrt(2)/2 && td1<Dmaxc1
                xd1c1(ely,elx)=1-(Dmaxc1-td1)^2;            
            end
        end
        if td1<sqrt(2)/2 && nc2_zeros==1           
            xd1c2(ely,elx)=td1^2;            
        else
            if td1>=Dmaxc2
                xd1c2(ely,elx)=1;            
            elseif td1>Dminc2 && td1<=Dminc2+sqrt(2)/2
                xd1c2(ely,elx)=(td1-Dminc2)^2;            
            elseif td1>Dminc2+sqrt(2)/2 && td1<Dmaxc2
                xd1c2(ely,elx)=1-(Dmaxc2-td1)^2;
            end
        end        
        xd1(ely,elx)=xd1(ely,elx)+xd1c1(ely,elx)+xd1c2(ely,elx);
        if  xd1(ely,elx)>1
             xd1(ely,elx)=1;
        end
        
        % Diagonal beam up-right corner to down-left corner to fliplr
        if td2<sqrt(2)/2
            if n_zeros==2
                xd2(ely,elx)=2*td2*(sqrt(2)-td2);
            elseif n_zeros==1
                xd2(ely,elx)=td2^2;
            end
        else
            if td2>=Dmax
                xd2(ely,elx)=1;
            elseif  td2>Dmin && td2<=Dmin+sqrt(2)/2
                xd2(ely,elx)=(td2-Dmin)^2;
            elseif td2>Dmin+sqrt(2)/2 && td2<Dmax
                xd2(ely,elx)=1-(Dmax-td2)^2;
            end
        end
        if td2<sqrt(2)/2 && nc1_zeros==1
            xd2c1(ely,elx)=td2^2;
        else
            if td2>=Dmaxc1
                xd2c1(ely,elx)=1;           
            elseif  td2>Dminc1 && td2<=Dminc1+sqrt(2)/2
                xd2c1(ely,elx)=(td2-Dminc1)^2;       
            elseif td2>Dminc1+sqrt(2)/2 && td2<Dmaxc1
                xd2c1(ely,elx)=1-(Dmaxc1-td2)^2;            
            end
        end
        if td2<sqrt(2)/2 && nc2_zeros==1           
            xd2c2(ely,elx)=td2^2;            
        else
            if td2>=Dmaxc2
                xd2c2(ely,elx)=1;            
            elseif td2>Dminc2 && td2<=Dminc2+sqrt(2)/2
                xd2c2(ely,elx)=(td2-Dminc2)^2;            
            elseif td2>Dminc2+sqrt(2)/2 && td2<Dmaxc2
                xd2c2(ely,elx)=1-(Dmaxc2-td2)^2;
            end
        end        
        xd2(ely,elx)=xd2(ely,elx)+xd2c1(ely,elx)+xd2c2(ely,elx);
        if  xd2(ely,elx)>1
             xd2(ely,elx)=1;
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

figure
colormap(gray); imagesc(1-x); caxis([0 1]); axis equal; axis off; drawnow;

end
