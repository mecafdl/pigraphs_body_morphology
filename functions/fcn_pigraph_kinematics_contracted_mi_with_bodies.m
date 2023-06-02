function [G_cntr, G_mst, c_kin] = fcn_pigraph_kinematics_contracted_mi_with_bodies(PLOT_MST, ...
            W_cntr, qIndex, wIndex, myLayout_type, SHOW_CLUSTERS, constPar)

rng('default')
interpreter = constPar.plotting.interpreter;
fcolor      = rand(constPar.nob,3);
% myLayout_type = 'Force';
if PLOT_MST == 0
    warning('Section of code not ready')  
    return
elseif PLOT_MST == 1
%     close all
    N_joints = constPar.noj;
    N_bodies = constPar.nob;
 
    figure
    %tiledlayout(1,1, 'Padding', 'none', 'TileSpacing', 'tight')
    set(gcf,'color','w');

    % Clean inter-joint connections   
%     W_cntr(1:N_joints,1:N_joints) = zeros(N_joints);
    G_cntr         = graph(W_cntr);
    LWidths        = 5*G_cntr.Edges.Weight/max(G_cntr.Edges.Weight);
    p_G              = plot(G_cntr, 'LineStyle','-', 'LineWidth', LWidths, 'Layout',myLayout_type,'EdgeColor',[0.5 0.5 0.5]);
%     p              = plot(G_cntr, 'LineStyle',':', 'LineWidth',2,'EdgeColor',[0.7 0.7 0.7],'Layout',myLayout_type);
    p_G.EdgeAlpha    = 1;
    p_G.Interpreter  = interpreter;    
    p_G.MarkerSize   = 20;


    labels = repmat({''},(N_joints + N_bodies),1);%cell(123,1);
    for i=1:N_joints
        if N_joints>20
            labels{i} = num2str(qIndex(i));
        else
            labels{i} = strcat('$\dot{q}_{',num2str(qIndex(i)),'}$');
        end
    end
    rng('default')
    subindices = ['0ABCDEFGHIJKLMNOPQRSTUVWXYZ'];
    % subindices = subindices(randperm(numel(subindices),numel(subindices)));
    for i=1:N_bodies
        %labels{N_joints + i} = strcat('$\omega_{',num2str(wIndex(i)),'}$');
%         labels{N_joints + i} = strcat('$\omega_{',subindices(wIndex(i)),'}$');
        if N_bodies>20
            labels{N_joints + i} = subindices(i);
        else
            labels{N_joints + i} = strcat('$\omega_{',subindices(i),'}$');
        end
    end

    p_G.NodeColor    = [0.7 0.7 0.7];
    if constPar.useNodeLabels == 1
        p_G.NodeLabel    = labels;
    else
        p_G.NodeLabel    = {};
    end

    p_G.MarkerSize   = 2.5*max(LWidths);%7;%15;
%     p.NodeColor    = [0.3 0.3 0.3];
    p_G.NodeFontSize = 8;%15%30;   
    p_G.EdgeCData    = G_cntr.Edges.Weight;%nonzeros(triu(W_cntr)');
    highlight(p_G, 1:N_joints, 'Marker', 'o')
    highlight(p_G, (N_joints+1):size(W_cntr,2), 'Marker', '^')
    colormap(flipud(autumn))     
%     colormap(flipud(gray))     
%     colormap('hot')     
    axis equal
      
    % Plot max. spanning tree ---------------------------------------------
    hold on
    [G_mst,~]  = minspantree(graph(-W_cntr,'upper'));
    G_mst.Edges.Weight = abs(G_mst.Edges.Weight);
    W_mst = full(G_mst.adjacency('weighted')); 
    
%     G_mst = digraph(triu(W_mst));
%     p_mst = plot(G_mst, 'LineWidth', 4, 'EdgeColor',[0, 1, 0],'EdgeAlpha',0.3,'Layout',myLayout_type) ;
%     p_mst.ArrowSize = 20;
%     p_mst.ArrowPosition = 0.8;


    p_mst = plot(G_mst, 'LineWidth', 1.75*max(LWidths), 'EdgeColor',[0, 1, 0],'EdgeAlpha',0.5,'Layout',myLayout_type) ;


    p_mst.NodeLabel    = {};
    p_mst.XData = p_G.XData;
    p_mst.YData = p_G.YData;
    p_mst.ZData = p_G.ZData;
    p_mst.NodeColor    = [1 1 1];
    %Returns handles to the patch and line objects    
    chi=get(gca, 'Children');
    %Reverse the stacking order so that the patch overlays the line
    set(gca, 'Children',flipud(chi))
    
%     layout(p2,'layered','Sources',19)
    layout(p_mst,myLayout_type)
    p_G.XData = p_mst.XData;
    p_G.YData = p_mst.YData;
    p_G.ZData = p_mst.ZData;

    if isfield(constPar,'imu_global_pos') && isfield(constPar,'joint_global_pos')
%         p_G.XData = [constPar.joint_global_pos(1,:),constPar.imu_global_pos(1,:)];
%         p_G.ZData = 0*[constPar.joint_global_pos(2,:),constPar.imu_global_pos(2,:)];
%         p_G.YData = [constPar.joint_global_pos(3,:),constPar.imu_global_pos(3,:)];
% 
%         p_mst.XData = [constPar.joint_global_pos(1,:),constPar.imu_global_pos(1,:)];
%         p_mst.ZData = 0*[constPar.joint_global_pos(2,:),constPar.imu_global_pos(2,:)];
%         p_mst.YData = [constPar.joint_global_pos(3,:),constPar.imu_global_pos(3,:)];        

%         p_G.XData = [1.25*constPar.imu_global_pos(1,2:end),constPar.imu_global_pos(1,:)];
%         p_G.ZData = 0*[constPar.imu_global_pos(2,2:end),constPar.imu_global_pos(2,:)];
%         p_G.YData = [1.25*constPar.imu_global_pos(3,2:end),constPar.imu_global_pos(3,:)];
% 
%         p_mst.XData = [1.25*constPar.imu_global_pos(1,2:end),constPar.imu_global_pos(1,:)];
%         p_mst.ZData = 0*[constPar.imu_global_pos(2,2:end),constPar.imu_global_pos(2,:)];
%         p_mst.YData = [1.25*constPar.imu_global_pos(3,2:end),constPar.imu_global_pos(3,:)];        

%         dist      = 0.02;    
%         tmp = [2:11,12*ones(1,7),19:26];
%         p_G.XData = [dist*sign(constPar.imu_global_pos(1,tmp)) + constPar.imu_global_pos(1,2:end),constPar.imu_global_pos(1,:)];
%         p_G.ZData = zeros(1,constPar.noj + constPar.nob);
%         p_G.YData = [dist*sign(constPar.imu_global_pos(3,2:end)) + constPar.imu_global_pos(3,2:end),constPar.imu_global_pos(3,:)];
% 
%         p_mst.XData = [dist*sign(constPar.imu_global_pos(1,tmp)) + constPar.imu_global_pos(1,2:end),constPar.imu_global_pos(1,:)];
%         p_mst.ZData = zeros(1,constPar.noj + constPar.nob);
%         p_mst.YData = [dist*sign(constPar.imu_global_pos(3,2:end)) + constPar.imu_global_pos(3,2:end),constPar.imu_global_pos(3,:)];
        
%         p_G.XData = [dist*sign(constPar.imu_global_pos(1,2:end)) + constPar.imu_global_pos(1,2:end),constPar.imu_global_pos(1,:)];
%         p_G.ZData = zeros(1,constPar.noj + constPar.nob);
%         p_G.YData = [dist*sign(constPar.imu_global_pos(3,2:end)) + constPar.imu_global_pos(3,2:end),constPar.imu_global_pos(3,:)];
% 
%         p_mst.XData = [dist*sign(constPar.imu_global_pos(1,2:end)) + constPar.imu_global_pos(1,2:end),constPar.imu_global_pos(1,:)];
%         p_mst.ZData = zeros(1,constPar.noj + constPar.nob);
%         p_mst.YData = [dist*sign(constPar.imu_global_pos(3,2:end)) + constPar.imu_global_pos(3,2:end),constPar.imu_global_pos(3,:)];

    end
%     return
% =========================================================================

%     p_G.XData = [constPar.joint_global_pos(1,:),constPar.imu_global_pos(1,:)];
%     p_G.YData = [constPar.joint_global_pos(2,:),constPar.imu_global_pos(2,:)];
%     p_G.ZData = [constPar.joint_global_pos(3,:),constPar.imu_global_pos(3,:)];
%     
%     p_mst.XData = [constPar.joint_global_pos(1,:),constPar.imu_global_pos(1,:)];
%     p_mst.YData = [constPar.joint_global_pos(2,:),constPar.imu_global_pos(2,:)];
%     p_mst.ZData = [constPar.joint_global_pos(3,:),constPar.imu_global_pos(3,:)];
    
%     axis equal
% =========================================================================
    
    if(SHOW_CLUSTERS)
        % Cluster definition using spectral clustering ------------------------
        N_clusters = N_bodies;
        rng('default')
        cluster_order = randperm(N_clusters,N_clusters);
        Aux     = full(G_mst.adjacency);
        C_kin    = [Aux(1:constPar.noj,constPar.noj+1:end);
                   eye(constPar.nob, constPar.nob)].*cluster_order;
        c_kin    = [sum(C_kin(1:constPar.noj,:),2);
                   transpose(sum(C_kin(constPar.noj+1:end,:),1))];
        plt           = gscatter(p_G.XData',p_G.YData',c_kin,fcolor,'o',1);

        % Color nodes according to clusters
        for node = 1:numel(c_kin)
%             plt(j).MarkerFaceColor = plt(j).Color;
            if node <= constPar.noj
                highlight(p_G, node, 'NodeColor',fcolor(node + 1,:))
            else
                highlight(p_G, node, 'NodeColor',fcolor(node - constPar.noj,:))
            end
                
%                 highlight(p_G, node, 'NodeColor', plt(c_kin(node)).Color)
%             if j<25
%                 plt(j).Marker ='o';
%             else
%                 plt(j).Marker = '^';
%             end
        end
        delete(findobj('type','legend'))
        clear dummy
        dummy      = zeros(N_clusters,1);
        leg_labels = cell(N_clusters,1);
        colormat   = zeros(N_clusters,3);
        for cluster = 1:N_clusters
            colormat(cluster,:) = plt(cluster).Color;
            % highlight(p,find(c_kin==cluster),'NodeColor', colormat(cluster,:));
            dummy(cluster) = plot(NaN,NaN,'o','Color','k');
            set(dummy(cluster), 'MarkerFaceColor', colormat(cluster,:),'MarkerSize',10)    
            leg_labels{cluster} = ['$b_{',num2str(cluster),'}$'];
            [R(cluster),C(cluster,:),~] = ExactMinBoundCircle([p_G.XData(find(c_kin==cluster))',p_G.YData(find(c_kin==cluster))']);
            [~]     =  VisualizeBoundCircle( ...
                            [p_G.XData(find(c_kin==cluster))', p_G.YData(find(c_kin==cluster))'], ...
                            1.2*R(cluster) + 1*0.15, ...
                            C(cluster,:), ...
                            colormat(cluster,:), ...
                            0.01, ...
                            gca);
            
            
            tst = find(c_kin==cluster);
            tst = tst(tst>constPar.noj) - constPar.noj;
            
            % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            %tx(cluster) = text(C(cluster,1)-(R(cluster)+0.15),C(cluster,2)+(R(cluster)+0.15),['$b_{',num2str(cluster),'}$'],'Interpreter','latex','FontSize',30,'Color','blue');
%             tx(cluster) = text(C(cluster,1)-(R(cluster)+0.15),C(cluster,2)+(R(cluster)+0.15),['$b_{',num2str(cluster),'}$'],'Interpreter','latex','FontSize',30,'Color',[0.9 0.9 0.9 0.1]);
%             tx(cluster) = text(C(cluster,1)-(R(cluster)+0.15),C(cluster,2)+(R(cluster)+0.15),['$b_{',num2str(cluster_order(tst)),'}$'],'Interpreter',interpreter,'FontSize',30,'Color',[0.7 0.7 0.7 0.3]);
%             for group = 1:N_clusters
%                 if cluster == group
%                     continue
%                 else
%                     A_top(cluster,group) = norm(W_mst(c_kin==cluster,c_kin==group),'fro');
%     %                 A_top(cluster,group) = norm(W_mst(find(c_kin_ref==cluster),find(c_kin_ref==group)),'fro');
%                 end
%             end      
            % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
        end    
        delete(plt)
        chi=get(gca, 'Children');
        %Reverse the stacking order so that the patch overlays the line
    %     set(gca, 'Children',flipud(chi))
    %     set(gca, 'Children',[chi(77); chi(78); chi(1:76)])
%         set(gca, 'Children',[chi(29); chi(30); chi(1:28)])
        set(gca, 'Children',[chi(end-1); chi(end); chi(1:end-2)])
    else
        c_kin = [];
    end    
%     if(SHOW_CLUSTERS)
%         % Cluster definition using spectral clustering ------------------------
%         N_clusters = N_bodies;
%         rng('default')
%         cluster_order = randperm(N_clusters,N_clusters);
%         Aux     = full(G_mst.adjacency);
%         C_kin    = [Aux(1:constPar.noj,constPar.noj+1:end);
%                    eye(constPar.nob, constPar.nob)].*cluster_order;
%         c_kin    = [sum(C_kin(1:constPar.noj,:),2);
%                    transpose(sum(C_kin(constPar.noj+1:end,:),1))];
%         plt           = gscatter(p_G.XData',p_G.YData',c_kin,fcolor,'o',3);
%         delete(findobj('type','legend'))
%         clear dummy
%         dummy      = zeros(N_clusters,1);
%         leg_labels = cell(N_clusters,1);
%         colormat   = zeros(N_clusters,3);
%         for cluster = 1:N_clusters
%             colormat(cluster,:) = plt(cluster).Color;
%             % highlight(p,find(c_kin==cluster),'NodeColor', colormat(cluster,:));
%             dummy(cluster) = plot(NaN,NaN,'o','Color','k');
%             set(dummy(cluster), 'MarkerFaceColor', colormat(cluster,:),'MarkerSize',10)    
%             leg_labels{cluster} = ['$b_{',num2str(cluster),'}$'];
%             [R(cluster),C(cluster,:),~] = ExactMinBoundCircle([p_G.XData(find(c_kin==cluster))',p_G.YData(find(c_kin==cluster))']);
%             [~]     =  VisualizeBoundCircle( ...
%                             [p_G.XData(find(c_kin==cluster))', p_G.YData(find(c_kin==cluster))'], ...
%                             R(cluster) + 0.15, ...
%                             C(cluster,:), ...
%                             colormat(cluster,:), ...
%                             0.01, ...
%                             gca);
%             
%             
%             tst = find(c_kin==cluster);
%             tst = tst(tst>constPar.noj) - constPar.noj;
%             
%             % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%             %tx(cluster) = text(C(cluster,1)-(R(cluster)+0.15),C(cluster,2)+(R(cluster)+0.15),['$b_{',num2str(cluster),'}$'],'Interpreter','latex','FontSize',30,'Color','blue');
% %             tx(cluster) = text(C(cluster,1)-(R(cluster)+0.15),C(cluster,2)+(R(cluster)+0.15),['$b_{',num2str(cluster),'}$'],'Interpreter','latex','FontSize',30,'Color',[0.9 0.9 0.9 0.1]);
%             tx(cluster) = text(C(cluster,1)-(R(cluster)+0.15),C(cluster,2)+(R(cluster)+0.15),['$b_{',num2str(cluster_order(tst)),'}$'],'Interpreter',interpreter,'FontSize',30,'Color',[0.7 0.7 0.7 0.3]);
%             for group = 1:N_clusters
%                 if cluster == group
%                     continue
%                 else
%                     A_top(cluster,group) = norm(W_mst(c_kin==cluster,c_kin==group),'fro');
%     %                 A_top(cluster,group) = norm(W_mst(find(c_kin_ref==cluster),find(c_kin_ref==group)),'fro');
%                 end
%             end      
%             % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
%         end    
%         delete(plt)
%         chi=get(gca, 'Children');
%         %Reverse the stacking order so that the patch overlays the line
%     %     set(gca, 'Children',flipud(chi))
%     %     set(gca, 'Children',[chi(77); chi(78); chi(1:76)])
% %         set(gca, 'Children',[chi(29); chi(30); chi(1:28)])
%         set(gca, 'Children',[chi(end-1); chi(end); chi(1:end-2)])
%     else
%         c_kin = [];
%     end

% #########################################################################
%     blabels = repmat({''},(N_links),1);
%     for i=1:N_links
%         blabels{i} = strcat('$b_{',num2str(i-1),'}$');
%     end    
%     A_top(A_top~=0) = 1;
%     p_top = plot(graph(A_top),'layout','layered');
%     p_top.EdgeAlpha = 1;
%     p_top.LineWidth = 3;
%     p_top.EdgeColor = 'k';
%     p_top.MarkerSize   = 20;
%     p_top.NodeColor    = colormat(1:19,:);%'k';%[0.7 0.7 0.7];
%     p_top.NodeFontSize = 30;
%     p_top.NodeLabel = blabels;
%     p_top.Interpreter  = 'latex';
% [~,idx] = sort(p_top.XData);   
%      p_top.XData = [0,-2,-1,-2,-2,-4,-6,-2,-1,-2, 2, 1, 2, 2, 4, 6, 2, 1, 2];
%      p_top.YData = [0, 1, 3, 5, 0, 0, 0,-1,-3,-5, 1, 3, 5, 0, 0, 0,-1,-3,-5];
%      p_top.ZData = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
% #########################################################################

%     %Returns handles to the patch and line objects ========================
%     chi=get(gca, 'Children');
%     
%     %Reverse the stacking order so that the patch overlays the line
%     set(gca, 'Children',flipud(chi))
%     p_aux1 = scatter(p.XData(1:N_joints)' , p.YData(1:N_joints) ,300,'marker','o','MarkerEdgeColor', 'k','LineWidth', 1);    
%     p_aux2 = scatter(p.XData(N_joints+1:end)', p.YData(N_joints+1:end),300,'marker','^','MarkerEdgeColor', 'k','LineWidth', 1);
%     
%     for cluster = 1:N_clusters
%         tx(cluster) = text(C(cluster,1)-(R(cluster)+0.15), ...
%                            C(cluster,2)+(R(cluster)+0.15), ...
%                            ['$b_{',num2str(cluster),'}$'], ...
%                            'Interpreter','latex', ...
%                            'FontSize',30, ...
%                            'Color','blue', ...
%                            'FontWeight', 'bold', ...
%                            'EdgeColor', 'none', ...
%                            'BackgroundColor', 'none');
%      
%     end 

    % Set colorbar --------------------------------------------------------
    % cb = colorbar;
    % cb.Location = 'westoutside';
    % set(cb,'YTick',[])
    % ylabel(cb,'Edge weight','FontSize',20,'Rotation',90);
    
    % Set legends ---------------------------------------------------------
    %leg = legend(dummy, leg_labels, 'Interpreter','latex','NumColumns',2);    
    %leg.Location = 'northeast';  
    %leg.Box      = 'off';
    %colormap(flipud(gray))     
    
    axis off
end

dummy_plot = zeros(2, 1);
dummy_plot(1) = plot(NaN,NaN,'-g','LineWidth',3);
dummy_plot(2) = plot(NaN,NaN,'--k','LineWidth',3);
legend(dummy_plot, 'MST','Body');

% dummy_plot    = zeros(1, 1);
% dummy_plot(1) = plot(NaN,NaN,'-g','LineWidth',3);
% % dummy_plot(2) = plot(NaN,NaN,'--k','LineWidth',3);
% legend(dummy_plot, 'MST');

fig = gcf;
ax  = gca;
plt = [p_G, p_mst];
leg = findobj(gcf, 'Type', 'Legend');
% Arguments to the function below:
% (fig, ax, plt, leg, tx, text_width, k_scaling, k_width_height)
% delete(leg)
% delete(cb)
if ~exist('tx','var')
    tx = [];
end
% fcn_scrpt_prepare_graph(fig, ax, plt, leg, tx, 1*8.8, 4, 1)

fcn_scrpt_prepare_graph_science_std(fig, ax, plt, leg, tx, 18.3/2, 4, 1)

% layout(p,'force3','WeightEffect','inverse','UseGravity',false)


% h = zeros(3, 1);
% h(1) = plot(NaN,NaN,'or');
% h(2) = plot(NaN,NaN,'ob');
% h(3) = plot(NaN,NaN,'ok');



RUN_THIS = 0;
if SHOW_CLUSTERS && RUN_THIS == 1
figure('Color',[1 1 1])
    blabels = repmat({''},(N_bodies),1);
    for i=1:N_bodies
        blabels{i} = strcat('$b_{',num2str(i),'}$');
    end   
    A_aux = 0*A_top;
    A_aux(A_top~=0) = A_top(A_top~=0);
    G_aux = graph(A_aux);
    [baseDegree,baseIndex] = max(sum(A_aux,2));
    LWidths        = 3*G_aux.Edges.Weight/max(G_aux.Edges.Weight);
    p_top = plot(graph(A_aux),'layout','layered','LineWidth', LWidths);
    layout(p_top,'layered','Sources',baseIndex)
    p_top.EdgeAlpha = 1;
    p_top.LineWidth = 3;
    p_top.EdgeColor = [0.7 0.7 0.7];
    p_top.MarkerSize   = 20;
    p_top.NodeColor    = colormat;%([2;c_kin(1:18)],:);%colormat(1:19,:);%'k';%[0.7 0.7 0.7];
    p_top.NodeFontSize = 30;
    p_top.NodeLabel = blabels;
    p_top.Interpreter  = 'latex';
    % highlight(p_top, 1:N_links, 'Marker', 'd')
    layout(p_top,'force')%,'WeightEffect','inverse')
    [~,idx] = sort(p_top.XData);   
    axis off
    axis equal
end
