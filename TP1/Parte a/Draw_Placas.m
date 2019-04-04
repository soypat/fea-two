% Mesh plotter
function Draw(Element,Position,color)
hold off
nel=size(Element,1);
ndofel=2*size(Element,2);
daspect([1 1 1])
% axis ([-1,11,-1,11]);
hold on
for iel=1:nel
    xx=zeros(5,1); yy=zeros(5,1);
    for i=1:4
        xx(i)=Position(Element(iel,i),1); yy(i)=Position(Element(iel,i),2); zz(i)=Position(Element(iel,i),3);
    end
    xx(5)=xx(1); yy(5)=yy(1); zz(5)=zz(1);
    plot3(xx,yy,zz,color)
end
grid on