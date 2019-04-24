close all
plot(Neleck,timerk)
hold on
plot(Nelec4,timer4)
plot(Nelec8,timer8)
legend('Kirchoff','Q4',"Q8")
title('Tiempo de cálculo para varios elementos')
xlabel('Número de elementos')
ylabel('Tiempo [s]')

figure
semilogx(Neleck,errveck)
hold on
semilogx(Nelec4,errvec4)
semilogx(Nelec8,errvec8)
legend('Kirchoff','Q4',"Q8")
title('Error de varios elementos')
xlabel('Número de elementos')
ylabel('Error máximo absoluto [mm]')

figure
loglog(Neleck,errveck)
hold on
loglog(Nelec4,errvec4)
loglog(Nelec8,errvec8)
legend('Kirchoff','Q4',"Q8")
title('Error de varios elementos')
xlabel('Número de elementos')
ylabel('Error máximo absoluto [mm]')