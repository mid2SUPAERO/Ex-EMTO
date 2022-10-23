function plot_def2D(U,nelx,nely)

%     [x_in,y_in]=meshgrid(0:nely,0:nelx);
%     [x_in,y_in]=meshgrid(0.5:nely+0.5,0.5:nelx+0.5);
    [x_in,y_in]=meshgrid(0:nely,nelx:-1:0);

    x_fin(:)=x_in(:)+U(1:2:end);
%     y_fin(:)=y_in(:)-U(2:2:end);
    y_fin(:)=y_in(:)+U(2:2:end);
    
    plot(x_fin(:),y_fin(:),'ro','linewidth',1.2);
    plot(x_fin(1),y_fin(1),'r*',x_fin(nely+1),y_fin(nely+1),'r*',x_fin((nely+1)*(nelx+1)),y_fin((nely+1)*(nelx+1)),'r*',x_fin((nely+1)*(nelx+1)-nelx),y_fin((nely+1)*(nelx+1)-nely),'r*','linewidth',1.2);
%     set(gca, 'YDir','reverse');
    hold on
    plot(x_in(:),y_in(:),'bo','linewidth',1.2)

    axis equal; axis tight; axis on; xlabel('x'); ylabel('y');

end