%% PERIODIC MATERIAL MICROSTRUCTURE DESIGN WITH 4 PARAMETERS
function [tens,obj,micro,par]=unitCell_4p(nelx,nely,density,penal,angle,cubicity,k_val)
%density : 0 for void, 1 for full material
%angle : 0 for 0 rad, 1 for pi/4 rads
%cubicity : 0 for only one privileged direction, 1 for cubic material
%the x_val has to be an horizontal vector, its components are related to
%the thicknesses of the beams that made the unit cell defined by Cell_4p
volfrac=density;
cubicity=sqrt(cubicity);
angle=angle*pi/4;
xval=k_val'+1e-6;
[~,szv]=size(k_val);
if szv==1
    error('k_val has to be 1x4')
end

%% MATERIAL PROPERTIES
E0=1;
Emin=1e-9;
nu=0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1, nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);

%% PERIODIC BOUNDARY CONDITIONS
e0 = eye(3);
ufixed = zeros(8,3);
U = zeros(2*(nely+1)*(nelx+1),3);
alldofs = (1:2*(nely+1)*(nelx+1));
n1 = [nodenrs(end,[1,end]),nodenrs(1,[end,1])];
d1 = reshape([(2*n1-1);2*n1],1,8);
n3 = [nodenrs(2:end-1,1)',nodenrs(end,2:end-1)];
d3 = reshape([(2*n3-1);2*n3],1,2*(nelx+nely-2));
n4 = [nodenrs(2:end-1,end)',nodenrs(1,2:end-1)];
d4 = reshape([(2*n4-1);2*n4],1,2*(nelx+nely-2));
d2 = setdiff(alldofs,[d1,d3,d4]);
   
for j = 1:3
    ufixed(3:4,j) = [e0(1,j),e0(3,j)/2;e0(3,j)/2,e0(2,j)]*[nelx;0];
    ufixed(7:8,j) = [e0(1,j),e0(3,j)/2;e0(3,j)/2,e0(2,j)]*[0;nely];
    ufixed(5:6,j) = ufixed(3:4,j)+ufixed(7:8,j);
end

wfixed = [repmat(ufixed(3:4,:),nely-1,1); repmat(ufixed(7:8,:),nelx-1,1)];

%% INITIALIZE ITERATION
qe= cell(3,3);
Q=zeros(3,3);
Qd=zeros(3,3);
change = 1;
loop = 0;
loop_after_min=0;
c_opt=0;

%% START ITERATION
while change >= 1e-5 && loop < 150 && loop_after_min < 50
    loop = loop+1;
    loop_after_min=loop_after_min+1;
    xPhys=Cell_4p(nelx,nely,xval);
    
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    Kr = [K(d2,d2), K(d2,d3)+K(d2,d4);K(d3,d2)+K(d4,d2), K(d3,d3)+K(d3,d4)+K(d4,d3)+K(d4,d4)];
    U(d1,:) = ufixed;
    U([d2,d3],:) = Kr\(-[K(d2,d1);K(d3,d1)+K(d4,d1)]*ufixed-[K(d2,d4); K(d3,d4)+K(d4,d4)]*wfixed);
    U(d4,:) = U(d3,:)+wfixed;

%   Plot def
%     for i=1:3
%         figure()
%         plot_def2D(U(:,i),nelx,nely)
%     end

    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    for i = 1:3
        for j = 1:3
            U1 = U(:,i); U2 = U(:,j);
            qe{i,j} = reshape(sum((U1(edofMat)*KE).*U2(edofMat),2),nely,nelx)/(nelx*nely);
            Q(i,j) = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*qe{i,j}));
        end
    end
    Q2=rotateTensorMatrix(Q,angle);
    c = -(1-0.5*cubicity)*Q2(1,1)-0.5*cubicity*Q2(2,2);
 
    % SENSITIVITIES
    h=1e-5;
    dc=zeros(size(xval));
    dv=zeros(size(xval));
    
    for ii=1:szv
        xdval=xval;
        xdval(ii)=xdval(ii)+h;
        xd=Cell_4p(nelx,nely,xdval);
        sK = reshape(KE(:)*(Emin+xd(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
        K = sparse(iK,jK,sK); K = (K+K')/2;
        Kr = [K(d2,d2), K(d2,d3)+K(d2,d4);K(d3,d2)+K(d4,d2), K(d3,d3)+K(d3,d4)+K(d4,d3)+K(d4,d4)];
        U(d1,:) = ufixed;
        U([d2,d3],:) = Kr\(-[K(d2,d1);K(d3,d1)+K(d4,d1)]*ufixed-[K(d2,d4); K(d3,d4)+K(d4,d4)]*wfixed);
        U(d4,:) = U(d3,:)+wfixed;    
        
        for i = 1:3
            for j = 1:3
                U1 = U(:,i); U2 = U(:,j);
                qe{i,j} = reshape(sum((U1(edofMat)*KE).*U2(edofMat),2),nely,nelx)/(nelx*nely);
                Qd(i,j) = sum(sum((Emin+xd.^penal*(E0-Emin)).*qe{i,j}));
            end
        end
        Q2d=rotateTensorMatrix(Qd,angle);
        c_d = -(1-0.5*cubicity)*Q2d(1,1)-0.5*cubicity*Q2d(2,2);
        dc(ii)=(c_d-c)/h;
        dv(ii)=(mean2(xd)-mean2(xPhys))/h;
    end
    
    %% OPTIMALITY CRITERIA UPDATE DESIGN VARIABLES
    l1 = 0; l2 = 1e9; move = density/4;
    while (l2-l1)/(l1+l2) > 1e-4
        lmid = 0.5*(l2+l1);
        xval_new = max(0,max(xval-move,min(1,min(xval+move,xval.*sqrt(-dc./dv/lmid)))));
        xPhys_OC=Cell_4p(nelx,nely,xval_new);
        if sum(xPhys_OC(:)) > volfrac*nelx*nely
            l1 = lmid;
        else
            l2 = lmid;
        end        
    end
    change = max(abs(xval_new-xval));
    xval_old=xval;
    xval=xval_new;
    
    %% Optimal Result
    if mean(xPhys(:))<=volfrac && c<c_opt && loop>1
        c_opt=c;
        obj=c;
        tens=Q(:);
        par=xval_old;
        micro=xPhys;
        loop_opt=loop; 
        loop_after_min=0;
    elseif exist('loop_opt','var')==0
        obj=c;
        tens=Q(:);
        par=xval_old;
        micro=xPhys;
    end  
    %% PRINT RESULTS
    fprintf(' It.:%4i  Obj.:%10.5f  Vol.:%6.3f  ch.:%8.5f\n',loop,c, mean(xPhys(:)),change);
    colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end
%% PRINT THE BEST RESULT
if exist('loop_opt','var')==1
    fprintf('\n The best result is  It.:%4i  Obj.:%10.5f  Vol.:%6.3f\n',loop_opt,obj, mean(micro(:)));
else
    fprintf('\n The best result is  It.:%4i  Obj.:%10.5f  Vol.:%6.3f\n',loop,obj, mean(micro(:)));
end
colormap(gray); imagesc(1-micro); caxis([0 1]); axis equal; axis off; drawnow;
