function Draw_Barra(Element,Node,c)

nel = size(Element,1);
nnodel = size(Element,2);
xx = zeros(nnodel,1);
yy = zeros(nnodel,1);
figure
hold on
for k = 1:nel;
    for i = 1:nnodel;
        xx(i) = Node(Element(k,i),1);
        yy(i) = Node(Element(k,i),2);
        text(xx(i),yy(i),num2str(Element(k,i)),'VerticalAlignment','bottom','Color','r','FontSize',16);
    end
    plot(xx,yy,c)
    plot(xx,yy,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',4)
    x = mean(Node(Element(k,:),1));
    y = mean(Node(Element(k,:),2));
    text(x,y,num2str(k),'EdgeColor','k','Color','k','FontSize',16);
    
end

grid on
axis equal
