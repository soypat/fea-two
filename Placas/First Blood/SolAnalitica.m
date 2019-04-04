%
N=100;
xgv = 0:(pi)/(N-1):pi;
Z = zeros(length(xgv));
for i=1:length(xgv)
   for j = 1:length(xgv)
       Z(i,j) = sin(xgv(i)).*sin(xgv(j));
   end
end

surf(xgv,xgv,Z)