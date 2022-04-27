% clear; clc; close all;

%% Visualize - Vg
%% Data load
Vd=1.00; dt=0.10;

FILENAME = sprintf('./data/Homework14_Transient/Homework14_Vd_%.2fV_dt_%.2fs.mat',Vd,dt);
load(FILENAME)

%% visulaize
Visual_solution_vector_Vg=zeros(size(Vertex,1),3);
for ii=1:size(Table_Jaco,1)
    Visual_solution_vector_Vg(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector_Vg_saved(ii,Newton_Vg);
end

% figure(1); % hole
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_Vg(:,3), 'EdgeColor','black','FaceColor','interp');
% title('hole')
% xlim([0 50*1e-9])
% ylim([0 14*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar
% 
% figure(2); % electron
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_Vg(:,2), 'EdgeColor','black','FaceColor','interp');
% title('elctron')
% xlim([0 50*1e-9])
% ylim([0 14*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar
% 
% figure(3); % potential
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_Vg(:,1), 'EdgeColor','black','FaceColor','interp');
% title('Initial potential')
% hold on
% patch('Faces',edge_interface,'Vertices',Vertex(:,[1 2]),'EdgeColor','black','FaceColor','none','LineWidth',5);
% hold on
% patch('Faces',Element_source_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
% hold on
% patch('Faces',Element_Dirichlet_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','g','FaceColor','none','LineWidth',5)
% hold off
% xlim([0 50*1e-9])
% ylim([0 14*1e-9])
% xlabel('X');
% ylabel('Y');
% colorbar
% 
% figure(4); % mesh 모양
% patch('Faces',Element,'Vertices',Vertex(:,[1 2]),'EdgeColor','black','FaceColor','none','LineWidth',1);
% title('Mesh')
% xlim([0 50*1e-9])
% ylim([0 14*1e-9])
% xlabel('X');
% ylabel('Y');
% hold on
% patch('Faces',edge_interface,'Vertices',Vertex(:,[1 2]),'EdgeColor','green','FaceColor','none','LineWidth',3);
% hold on
% patch('Faces',Element_si_source,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
% hold on
% patch('Faces',Element_si_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
% hold on
% patch('Faces',Element_si_channel,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','red','LineWidth',0.5)
% hold on
% patch('Faces',Element_ox1,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
% hold on
% patch('Faces',Element_ox2,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
% hold on
% patch('Faces',Element_ox2,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
% hold on
% patch('Faces',Element_source_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
% hold on
% patch('Faces',Element_Dirichlet_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','g','FaceColor','none','LineWidth',5)
% hold off
% 
% figure(5) % 점으로 3D
% plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_Vg(:,1),'*')
% xlim([0 50*1e-9])
% ylim([0 14*1e-9])

figure(8) % current
x=time;
y=I;
plot(x,y,'-o')
title('Drain Current vs time')
xlabel('Time [us]');
ylabel('Drain Current [V]');
xlim([0 5])

figure(7); % gate
x=time;
y=Vg_save;
plot(x,y, '-o')
title('Gate Voltage vs time')
xlabel('Time [us]');
ylabel('Gate Voltage [V]');
xlim([0 5])
ylim([0 1.2])


plot(x_010,I_010,'--r', x_005,I_005,'-.b', x_001,I_001,':m')
title('Drain Current vs time')
xlabel('Time [us]');
ylabel('Drain Current [V]');
legend('Time step: 0.10us','Time step: 0.05us', 'Time step: 0.01us')
xlim([0 5])

% figure(6);
% subplot(1,2,1)
% x=time;
% y=Vg_save;
% plot(x,y,'-o')
% title('Gate Voltage vs time')
% xlabel('Time [us]');
% ylabel('Gate Voltage [V]');
% xlim([0 5])
% 
% subplot(1,2,2)
% x=time;
% y=I;
% plot(x,y,'-o')
% title('Drain Current vs time')
% xlabel('Time [us]');
% ylabel('Drain Current [V]');
% xlim([0 5])