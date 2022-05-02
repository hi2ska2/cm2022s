% clear;
clc; close all;

%% Visualize - Vg
%% Data load

% 열기 기능을 이용해서 직접 load 하자.
% 볼 time만 설정하자.
Type='Vg';

%% visulaize %%
%% DD visual
if strcmp('DD',Type)==1
    Visual_solution_vector_DD=zeros(size(Vertex,1),3);
    for ii=1:size(Table_Jaco,1)
        Visual_solution_vector_DD(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector_DD_saved(ii,Newton_DD);
    end

    figure(1) % 점으로 3D
    plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_DD(:,1),'*')
    xlim([0 50*1e-9])

    figure(2); % mesh 모양
    patch('Faces',Element,'Vertices',Vertex(:,[1 2]),'EdgeColor','black','FaceColor','none','LineWidth',1);
    title('Mesh')
    xlim([0 50*1e-9])
    xlabel('X');
    ylabel('Y');
    hold on
    patch('Faces',edge_interface,'Vertices',Vertex(:,[1 2]),'EdgeColor','green','FaceColor','none','LineWidth',3);
    hold on
    patch('Faces',Element_si_source,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
    hold on
    patch('Faces',Element_si_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
    hold on
    patch('Faces',Element_si_channel,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','red','LineWidth',0.5)
    hold on
    patch('Faces',Element_ox1,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
    hold on
    patch('Faces',Element_ox2,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
    hold on
    patch('Faces',Element_ox2,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
    hold on
    patch('Faces',Element_source_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
    hold on
    patch('Faces',Element_Dirichlet_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','g','FaceColor','none','LineWidth',5)
    hold off

    figure(3); % potential
    patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_DD(:,1), 'EdgeColor','black','FaceColor','interp');
    title('Initial potential')
    hold on
    patch('Faces',edge_interface,'Vertices',Vertex(:,[1 2]),'EdgeColor','black','FaceColor','none','LineWidth',5);
    hold on
    patch('Faces',Element_source_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
    hold on
    patch('Faces',Element_Dirichlet_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','g','FaceColor','none','LineWidth',5)
    hold off
    xlim([0 50*1e-9])
    xlabel('X');
    ylabel('Y');
    colorbar

    figure(4); % electron
    patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_DD(:,2), 'EdgeColor','black','FaceColor','interp');
    title('elctron')
    xlim([0 50*1e-9])
    xlabel('X');
    ylabel('Y');
    colorbar

    figure(5); % hole
    patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_DD(:,3), 'EdgeColor','black','FaceColor','interp');
    title('hole')
    xlim([0 50*1e-9])
    xlabel('X');
    ylabel('Y');
    colorbar

    %% Vd visual
elseif strcmp('Vd',Type)==1
    Visual_solution_vector_Vd=zeros(size(Vertex,1),3);
    for ii=1:size(Table_Jaco,1)
        Visual_solution_vector_Vd(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector_Vd_saved(ii,Newton_Vd);
    end

    figure(1) % 점으로 3D
    plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_Vd(:,1),'*')
    xlim([0 50*1e-9])

    figure(2); % mesh 모양
    patch('Faces',Element,'Vertices',Vertex(:,[1 2]),'EdgeColor','black','FaceColor','none','LineWidth',1);
    title('Mesh')
    xlim([0 50*1e-9])
    xlabel('X');
    ylabel('Y');
    hold on
    patch('Faces',edge_interface,'Vertices',Vertex(:,[1 2]),'EdgeColor','green','FaceColor','none','LineWidth',3);
    hold on
    patch('Faces',Element_si_source,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
    hold on
    patch('Faces',Element_si_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
    hold on
    patch('Faces',Element_si_channel,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','red','LineWidth',0.5)
    hold on
    patch('Faces',Element_ox1,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
    hold on
    patch('Faces',Element_ox2,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
    hold on
    patch('Faces',Element_ox2,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
    hold on
    patch('Faces',Element_source_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
    hold on
    patch('Faces',Element_Dirichlet_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','g','FaceColor','none','LineWidth',5)
    hold off

    figure(3); % potential
    patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_Vd(:,1), 'EdgeColor','black','FaceColor','interp');
    title('Initial potential')
    hold on
    patch('Faces',edge_interface,'Vertices',Vertex(:,[1 2]),'EdgeColor','black','FaceColor','none','LineWidth',5);
    hold on
    patch('Faces',Element_source_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
    hold on
    patch('Faces',Element_Dirichlet_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','g','FaceColor','none','LineWidth',5)
    hold off
    xlim([0 50*1e-9])
    xlabel('X');
    ylabel('Y');
    colorbar

    figure(4); % electron
    patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_Vd(:,2), 'EdgeColor','black','FaceColor','interp');
    title('elctron')
    xlim([0 50*1e-9])
    xlabel('X');
    ylabel('Y');
    colorbar

    figure(5); % hole
    patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_Vd(:,3), 'EdgeColor','black','FaceColor','interp');
    title('hole')
    xlim([0 50*1e-9])
    xlabel('X');
    ylabel('Y');
    colorbar

    %% Vg visual
elseif strcmp('Vg',Type)==1



    prompt = 'Please enter the time? ';
    TIME=input(prompt);

%     for ii=1:length(time)
%         if abs(time(ii,1)-time_input)<=1e-15
%             TIME=ii;
%         end
%     end


    Visual_solution_vector_Vg=zeros(size(Vertex,1),3);
    for ii=1:size(Table_Jaco,1)
        Visual_solution_vector_Vg(Table_Jaco(ii,2),Table_Jaco(ii,3))=solution_vector_Vg_saved(ii,TIME);
    end

    figure(1); % hole
    patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_Vg(:,3), 'EdgeColor','black','FaceColor','interp');
    title('hole')
    xlabel('X');
    ylabel('Y');
    colorbar

    figure(2); % electron
    patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_Vg(:,2), 'EdgeColor','black','FaceColor','interp');
    title('elctron')
    xlabel('X');
    ylabel('Y');
    colorbar

    figure(3); % potential
    patch('Faces',Element,'Vertices',Vertex(:,[1 2]), 'FaceVertexCData',Visual_solution_vector_Vg(:,1), 'EdgeColor','black','FaceColor','interp');
    title('Initial potential')
    hold on
    patch('Faces',edge_interface,'Vertices',Vertex(:,[1 2]),'EdgeColor','black','FaceColor','none','LineWidth',5);
    hold on
    patch('Faces',Element_source_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
    hold on
    patch('Faces',Element_Dirichlet_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','g','FaceColor','none','LineWidth',5)
    hold off
    xlabel('X');
    ylabel('Y');
    colorbar

    figure(4); % mesh 모양
    patch('Faces',Element,'Vertices',Vertex(:,[1 2]),'EdgeColor','black','FaceColor','none','LineWidth',1);
    title('Mesh')
    xlabel('X');
    ylabel('Y');
    hold on
    patch('Faces',edge_interface,'Vertices',Vertex(:,[1 2]),'EdgeColor','green','FaceColor','none','LineWidth',3);
    hold on
    patch('Faces',Element_si_source,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
    hold on
    patch('Faces',Element_si_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','yellow','LineWidth',0.5)
    hold on
    patch('Faces',Element_si_channel,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','red','LineWidth',0.5)
    hold on
    patch('Faces',Element_ox1,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
    hold on
    patch('Faces',Element_ox2,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
    hold on
    patch('Faces',Element_ox2,'Vertices',Vertex(:,[1 2]), 'EdgeColor','black','FaceColor','blue','LineWidth',0.5)
    hold on
    patch('Faces',Element_source_drain,'Vertices',Vertex(:,[1 2]), 'EdgeColor','m','FaceColor','none','LineWidth',5)
    hold on
    patch('Faces',Element_Dirichlet_BC,'Vertices',Vertex(:,[1 2]), 'EdgeColor','g','FaceColor','none','LineWidth',5)
    hold off

    figure(5) % 점으로 3D
    plot3(Vertex(:,1),Vertex(:,2),Visual_solution_vector_Vg(:,1),'*')

    figure(8) % current
    plot(time, I )
    title('Drain Current vs time')
    xlabel('Time [s]');
    ylabel('Drain Current [A]');
    ylim([0 0.06])


    figure(7); % gate
    plot(time, Vg)
    title('Gate Voltage vs time')
    xlabel('Time [s]');
    ylabel('Gate Voltage [V]');

    figure(9) % current
    y1=I;
    yyaxis right
    plot(time, y1 )
    title('Drain Current vs time')
    xlabel('Time [s]');
    ylabel('Drain Current [A]');
    

    hold on
    % gate
    y2=Vg;
    yyaxis left
    plot(time, y2)
    title('Gate Voltage vs time')
    xlabel('Time [s]');
    ylabel('Gate Voltage [V]');
    hold off


end