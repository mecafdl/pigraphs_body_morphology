function fcn_phantomx_display_learned_morphology(lambda_hat, p_R_c, pxmark4, constPar, mi_graphs, parents, children, fig)


%% ************************************************************************
%                      DISPLAY LEARNED KINEMATICS                         *
% *************************************************************************

% clc
% close all

% Display pose
q_phantomx_plot = pxmark4.homeConfiguration;%0*q_poppy(:, point);

% Parameters
gamma_hat         = lambda_hat(1:7,:);
rho_hat           = lambda_hat(8:13,:);    

% Initialization
sp_r_sc     = NaN(3,constPar.noj);
sp_R_sc_hat = NaN(3,3,constPar.noj);
w_T_sc      = NaN(4,4,constPar.nol);
w_T_sp      = NaN(4,4,constPar.nol);
w_T_sc_hat  = NaN(4,4,constPar.nol);
%w_R_sc_hat  = NaN(3,3,constPar.nol);
sp_T_sc_hat = NaN(4,4,constPar.nol);

% Get (actual) transformation from pelivs to world frame
% w_T_sc(:,:,1)      = getTransform(poppy_fb, q_poppy_plot, ...
%                          poppy_fb.Bodies{constPar.imuBodyIndices(constPar.reference_body_order(1))}.Name); 
% w_T_sc(:,:,1)      = getTransform(poppy, q_phantomx_plot(constPar.reference_joint_order), ...
%                          poppy.Bodies{constPar.imuBodyIndices(constPar.reference_body_order(1))}.Name); 
w_T_sc(:,:,1)      = getTransform(pxmark4, q_phantomx_plot, ...
                         pxmark4.Bodies{constPar.imuBodyIndices(constPar.reference_body_order(1))}.Name); 
%w_R_sc_hat(:,:,1)     = w_T_sc(1:3,1:3,1);
w_T_sc_hat(:,:,1)  = w_T_sc(:,:,1);
sp_R_sc_hat(:,:,1) = w_T_sc(1:3,1:3,1);
sp_T_sc_hat(:,:,1) = w_T_sc(:,:,1);

% trplot(rt2tr(s0_R_sj_hat(:,:,1), s0_T_sj(1:3,4,1)),...
%     'rgb','frame',suffixes{1},'arrow','width',0.1,'length', 0.05)   

% Show robot
clear alpha
% figure('Color','w')
% rb = show(poppy, q_poppy_plot(constPar.reference_joint_order), 'PreservePlot', false, 'Frames','off');
% alpha(rb,0.01)
% hold on    
% [~] = fcn_pigraph_kinematics_directed(triu(adjacency(mi_graphs.G_kin_mst,'weighted')), constPar, fig); hold on;
rb = show(pxmark4, q_phantomx_plot, 'PreservePlot', false, 'Frames','off');
alpha(rb,0.05)
hold on    

% Show world frame
trplot(eye(4),'rgb', 'frame',['W'],'length', 0.025,'thick',3)

% Sphere representing the pelvis IMU --------------------------------------
sensor_radius = 0.005;
frame_length  = 0.03;
fcn_plot_imu_with_frame(sensor_radius, frame_length, w_T_sc(:,:,1))

for j = 1:constPar.noj
    p = parents(j);
    c = children(j);
    
    % Display body name ---------------------------------------------------   
%     disp(poppy_fb.Bodies{constPar.imuBodyIndices(constPar.reference_body_order(c))}.Name)
%     joint_index = find(strcmp(jointNames, ...
%         poppy_fb.Bodies{imuBodyIndices(constPar.reference_body_order(c))}.Parent.Joint.Name));    

%     disp(poppy.Bodies{constPar.imuBodyIndices(constPar.reference_body_order(c))}.Name)


%     joint_index = find(strcmp(poppy_properties.jointNames, ...
%         poppy.Bodies{constPar.imuBodyIndices(constPar.reference_body_order(c))}.Parent.Joint.Name));    



    % Estimated transformation matrix -------------------------------------
    sp_R_sc_hat(:,:,c) = p_R_c(q_phantomx_plot((c-1)), gamma_hat(:,c-1));
    sp_r_sc(:,c-1)     = -rho_hat(1:3,c-1) + sp_R_sc_hat(:,:,c)*rho_hat(4:6,c-1);
    sp_T_sc_hat(:,:,c) = rt2tr(sp_R_sc_hat(:,:,c), sp_r_sc(:,c-1));
    w_T_sc_hat(:,:,c)  = w_T_sc_hat(:,:,p)*sp_T_sc_hat(:,:,c);

    % Plot IMU with CS ----------------------------------------------------
    fcn_plot_imu_with_frame(sensor_radius, frame_length, w_T_sc_hat(:,:,c))
   
    % Indetified joint rotation axes plotted at IMU location --------------
%     jointAxis = 0.05*gamma_hat(5:7,c-1);
%     rotAxis   = w_T_sc_hat(1:3,1:3,c)*jointAxis;
%     quiver3(w_T_sc_hat(1,4,c), w_T_sc_hat(2,4,c), w_T_sc_hat(3,4,c), rotAxis(1), rotAxis(2), rotAxis(3), 1, 'k:', 'LineWidth',2); 
        
    % Axis j in sensor j position vector ----------------------------------
    point_start     = w_T_sc_hat(1:3,4,c);
    point_end_jnt_p = w_T_sc_hat(:,:,c)*[(-rho_hat(4:6,c-1));1];
    v1              = [point_start(1),point_start(2),point_start(3)];
    v2              = [point_end_jnt_p(1),point_end_jnt_p(2),point_end_jnt_p(3)];
    arrow3(v1,v2,['m:' '3'],0.1,0.1) 

    % Axis j+1 in sensor j position vector --------------------------------
    point_start     = w_T_sc_hat(1:3,4,p);
    point_end_jnt_c = w_T_sc_hat(:,:,p)*[(-rho_hat(1:3,c-1));1];
    v1              = [point_start(1),point_start(2),point_start(3)];
    v2              = [point_end_jnt_c(1),point_end_jnt_c(2),point_end_jnt_c(3)];
    arrow3(v1,v2,['k:' '3'],0.1,0.1)           

    % Axis translated to identified joint center point --------------------
    jointAxis    = 0.05*gamma_hat(5:7,c-1);
    rotAxisStart = w_T_sc_hat(:,:,c)*[(-rho_hat(4:6,c-1));1];
    rotAxisEnd   = w_T_sc_hat(:,:,c)*[(-rho_hat(4:6,c-1) + jointAxis);1];
    v1=[rotAxisStart(1),rotAxisStart(2),rotAxisStart(3)];
    v2=[rotAxisEnd(1),rotAxisEnd(2),rotAxisEnd(3)];
    v=[v2;v1];
%     plot3(v(:,1),v(:,2),v(:,3),'y-','LineWidth',6)
    arrow3(v1,v2,['_b-' '3'],0.5,0.4)
%     text(1.01*rotAxisEnd(1),1.01*rotAxisEnd(2),1.01*rotAxisEnd(3), ...
%          strcat('$\zeta_{',constPar.suffixes(c),'}$'), ...
%          'interpreter','latex','FontSize',20)
end

axis equal
axis off

% /home/diaz/Dropbox/PhD_WORK/matlab_irt/phantomx/morphology_learning/figures
% view([-330 45])

view([0 0])
SAVE_FIG = 0;
if SAVE_FIG == 1
    export_fig(fullfile(fcn_get_matlab_files_dir,'poppy','morphology_learning','figures','phantomx_body_structure'),'-pdf')
    close(gcf);
end
