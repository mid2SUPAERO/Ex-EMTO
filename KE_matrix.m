% CALCULATION OF STIFFNESS MATRIX AND ITS FOUR PARTIAL DERIVATIVES
% The stiffness matrix and derivatives depends on the cell microstructure, defined by the
% macro-variables
function [KE, KE_dxdens, KE_dxcos, KE_dxsin, KE_dxcub] = KE_matrix(xdens, xcos, xsin, xcub)
    global B database sig;  
    %% get point's tensor and derivatives
    xdensd=xdens+0.01; %xdens+delta used to get the partial derivative approximation
    negdifdens=1;
    if xdensd>1 %if right partial derivative isn't accessible, get left partial derivative
        xdensd=xdensd-0.02;
        negdifdens=-1;
    end
    xcubd=xcub+0.01; %xcub+delta used to get the partial derivatives
    negdifcub=1;
    if xcubd>1 %if right partial derivative isn't accessible, get left partial derivative
        xcubd=xcubd-0.02;
        negdifcub=-1;
    end
    %derive xor from xcos and xsin. 
    cosalpha=2*xcos-1; sinalpha=2*xsin-1;
    xor=atan(sinalpha/cosalpha)/pi; %Here, xor is in [-1,1], representing an orientation angle in [-pi,pi]
    %put the orientation angle in [0,pi]
    if xor<0
        xor=xor+1;
    end 
    xord=xor+0.01; %xcor+delta used to get the partial derivatives
    negdifor=1;
    if xord>1 %if right partial derivative isn't accessible, get left partial derivative
        xord=xord-0.02;
        negdifor=-1;
    end
    %macrovariables for xi, xi', xi'' and xi'''
    xdensv=[xdens, xdensd, xdens, xdens];
    xorv=[xor, xor, xord, xor];
    xcubv=[xcub, xcub, xcub, xcubd];
    d11=[];
    d12=[];
    d13=[];
    d22=[];
    d23=[];
    d33=[];
    % GET ELASTICITY TENSORS FROM DATABASE METAMODEL FOR Xi, Xi', Xi'' and Xi'''
    for i = 1:4
        % find points in matrix within 3 kernel radii
        xdenslim1=max(1,round((xdensv(i)-3*sig)*(size(database,4)-1)+1));
        xdenslim2=min(size(database,4),round((xdensv(i)+3*sig)*(size(database,4)-1)+1));
        xorlim1=max(1,round((xorv(i)-3*sig)*(size(database,3)-1)+1));
        xorlim2=min(size(database,3),round((xorv(i)+3*sig)*(size(database,3)-1)+1));
        xcublim1=max(1,round((xcubv(i)-3*sig)*(size(database,2)-1)+1));
        xcublim2=min(size(database,2),round((xcubv(i)+3*sig)*(size(database,2)-1)+1));

        Msim=database(:,xcublim1:xcublim2,xorlim1:xorlim2,xdenslim1:xdenslim2);

        %get distance of each point
        distance=sqrt((Msim(1,:,:,:)-xdensv(i)).^2+(Msim(2,:,:,:)-xorv(i)).^2+(Msim(3,:,:,:)-xcubv(i)).^2);
        %Nadaraya-Watson kernel-weighted average
        gaussFactor = (1/((2*pi*sig^2)^(3/2)))*exp(-((distance).^2)./(2*sig^2));
        gaussFactors = repmat(gaussFactor,6,1,1,1);
        totalGauss=sum(sum(sum(gaussFactor)));
        pointGauss=sum(sum(sum(Msim(4:end,:,:,:).*gaussFactors,2), 3), 4)/totalGauss;


        d11 = [d11,pointGauss(1)];
        d12 = [d12,pointGauss(2)];
        d13 = [d13,pointGauss(3)];
        d22 = [d22,pointGauss(4)];
        d23 = [d23,pointGauss(5)];
        d33 = [d33,pointGauss(6)];
    end
    %Approximate partial derivatives
    d11ddens=negdifdens*(d11(2)-d11(1))/0.01;
    d12ddens=negdifdens*(d12(2)-d12(1))/0.01;
    d13ddens=negdifdens*(d13(2)-d13(1))/0.01;
    d22ddens=negdifdens*(d22(2)-d22(1))/0.01;
    d23ddens=negdifdens*(d23(2)-d23(1))/0.01;
    d33ddens=negdifdens*(d33(2)-d33(1))/0.01;
    
    d11dor=negdifor*(d11(3)-d11(1))/0.01;
    d12dor=negdifor*(d12(3)-d12(1))/0.01;
    d13dor=negdifor*(d13(3)-d13(1))/0.01;
    d22dor=negdifor*(d22(3)-d22(1))/0.01;
    d23dor=negdifor*(d23(3)-d23(1))/0.01;
    d33dor=negdifor*(d33(3)-d33(1))/0.01;
    
    d11dcub=negdifcub*(d11(4)-d11(1))/0.01;
    d12dcub=negdifcub*(d12(4)-d12(1))/0.01;
    d13dcub=negdifcub*(d13(4)-d13(1))/0.01;
    d22dcub=negdifcub*(d22(4)-d22(1))/0.01;
    d23dcub=negdifcub*(d23(4)-d23(1))/0.01;
    d33dcub=negdifcub*(d33(4)-d33(1))/0.01;

    D = [d11(1), d12(1), d13(1);
        d12(1), d22(1), d23(1);
        d13(1), d23(1), d33(1)];
    D_dxdens = [d11ddens, d12ddens, d13ddens;
        d12ddens, d22ddens, d23ddens;
        d13ddens, d23ddens, d33ddens];
    D_dxor = [d11dor, d12dor, d13dor;
        d12dor, d22dor, d23dor;
        d13dor, d23dor, d33dor];
    D_dxcub = [d11dcub, d12dcub, d13dcub;
        d12dcub, d22dcub, d23dcub;
        d13dcub, d23dcub, d33dcub];
    
    %DERIVE STIFFNESS MATRIX FROM ELASTICITY TENSOR
    %get gauss points weights and positions
    x = [-sqrt(3)/3, sqrt(3)/3]; % gauss point positions
    w = [1,1]; % gauss point weights
    ke = zeros(8);
    ke_dxdens = zeros(8);
    ke_dxor = zeros(8);
    ke_dxcub = zeros(8);
    [i, j] = meshgrid(1:2, 1:2);
    for m = 1:4
        ke = ke + w(i(m))*w(j(m))*B(x(i(m)), x(j(m)))' * D * B(x(i(m)), x(j(m)));
        ke_dxdens = ke_dxdens + w(i(m))*w(j(m))*B(x(i(m)), x(j(m)))' * D_dxdens * B(x(i(m)), x(j(m)));
        ke_dxor = ke_dxor + w(i(m))*w(j(m))*B(x(i(m)), x(j(m)))' * D_dxor * B(x(i(m)), x(j(m)));
        ke_dxcub = ke_dxcub + w(i(m))*w(j(m))*B(x(i(m)), x(j(m)))' * D_dxcub * B(x(i(m)), x(j(m)));
    end
    ke_dxor=ke_dxor/pi;
    ke_dxcosalpha=ke_dxor*(-sinalpha/(cosalpha^2+sinalpha^2));
    ke_dxcos=2*ke_dxcosalpha;
    ke_dxsinalpha=ke_dxor*(cosalpha/(cosalpha^2+sinalpha^2));
    ke_dxsin=2*ke_dxsinalpha;
    KE = ke; KE_dxdens = ke_dxdens; KE_dxcos = ke_dxcos; KE_dxsin = ke_dxsin; KE_dxcub = ke_dxcub;

end