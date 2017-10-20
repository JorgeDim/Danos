clear;
prefijos={...
...% '001' ...% (0)eta=0       , (0)E(u),    (1)    div(u)
...% '002' ...% (0)eta=0       , (0)E(u),    (2)abs(div(u))
...% '003' ...% (0)eta=0       , (0)E(u),    (3)max(div(u),0)
...% '011' ...% (0)eta=0       , (1)E(udot), (1)    div(u)
...% '012' ...% (0)eta=0       , (1)E(udot), (2)abs(div(u))
...% '013' ...% (0)eta=0       , (1)E(udot), (3)max(div(u),0)
...% '101' ...% (1)eta=articulo, (0)E(u),    (1)    div(u)
...% '102' ...% (1)eta=articulo, (0)E(u),    (2)abs(div(u))
...% '103' ...% (1)eta=articulo, (0)E(u),    (3)max(div(u),0)
...% '111' ...% (1)eta=articulo, (1)E(udot), (1)    div(u)
...% '112' ...% (1)eta=articulo, (1)E(udot), (2)abs(div(u))
...% '113' ...% (1)eta=articulo, (1)E(udot), (3)max(div(u),0)
...% '201' ...% (2)eta=MGrande , (0)E(u),    (1)    div(u)
     '202' ...% (2)eta=MGrande , (0)E(u),    (2)abs(div(u))
...% '203' ...% (2)eta=MGrande , (0)E(u),    (3)max(div(u),0)
     '211' ...% (2)eta=MGrande , (1)E(udot), (1)    div(u)
...% '212' ...% (2)eta=MGrande , (1)E(udot), (2)abs(div(u))
...% '213' ...% (2)eta=MGrande , (1)E(udot), (3)max(div(u),0)
    };
%%
np=length(prefijos);
for i=1:np
    d{i}=load([prefijos{i} '_extremealpha.dat']);
    dat=load([prefijos{i} '_datos.dat']);
    ac{i}=dat(1);MGrande{i}=dat(2);MGrande2{i}=dat(3);mu0{i}=dat(4);test=0;
end

figure(1);clf;
sp1=subplot(4,2,[1 2]);
hold on
Tmax=0;
ymax=0;
for i=1:np
    dd=d{i};
    Tmax=max(Tmax,d{i}(end,1));
    ymax=max([ymax,d{i}(end,2),d{i}(end,3)]);
    plot(dd(:,1), dd(:,1)*0+ac{i},'r','Linewidth',2)
    if i==1
        plot(dd(:,1)+i/1000, dd(:,2:3))
    else
        plot(dd(:,1)+i/1000, dd(:,2:3),'+')
    end
end
title('Min, max \alpha(t)')
xlabel('Time');
grid on
axis([0 Tmax 0 ymax*1.1]);

sp2=subplot(4,2,[3 4]);%figure(2);clf;
hold on

ymin=0;
ymax=0;
for i=1:np
    dd=d{i};
    ymax=max([ymax;d{i}(:,4);d{i}(:,5)]);
    ymin=min([ymin;d{i}(:,4);d{i}(:,5)]);
    if i==1
        plot(dd(:,1)+i/1000, dd(:,4:5))
    else
        plot(dd(:,1)+i/1000, dd(:,4:5),'+')
    end
end
title('Min, max z(t)')
xlabel('Time');
grid on
axis([0 Tmax ymin*1.1 ymax]);

sp3=subplot(4,1,3);hold on

ymin=0;
ymax=0;
for i=1:np
    ymax=max([ymax;d{i}(:,6)]);
    ymin=min([ymin;d{i}(:,6)]);
    if i==1
        plot(d{i}(:,1), d{i}(:,6),'.-')
    else
        plot(d{i}(:,1), d{i}(:,6),'+-')
    end
end
grid on
axis([0 Tmax ymin ymax]);
title('N° de Iteraciones');
%xlabel('Time');
sp4=subplot(4,1,4);hold on

ymin=0;
ymax=0;
for i=1:np
    ymax=max([ymax;d{i}(:,7)]);
    ymin=min([ymin;d{i}(:,7)]);
    if i==1
        plot(d{i}(:,1), d{i}(:,7),'.-')
    else
        plot(d{i}(:,1), d{i}(:,7),'+-')
    end
end
grid on
axis([0 Tmax ymin ymax]);
title('Residuo');


set(sp1,'Position',[0.0300 0.69 0.96 0.27])
set(sp2,'Position',[0.0300 0.35 0.96 0.27])
set(sp3,'Position',[0.0300 0.17 0.96 0.10])
set(sp4,'Position',[0.0300 0.03 0.96 0.10])


%FuncionIndicatriz
