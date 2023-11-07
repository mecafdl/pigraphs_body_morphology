
clc
% constPar.displayJointConfiguration =  [0 0 0 0 0 deg2rad(120) 0]';
constPar.displayJointConfiguration =  q_panda(:,test_points(i));
h = figure('color','w');
% fcn_poppy_display_learned_morphology(mean(lambda_rams(:,:,end-1000:end),3), p_R_c, poppy, constPar, mi_graphs, parents, children, h)
fcn_franka_display_learned_morphology(lambda_hat_online_float, p_R_c, panda, constPar, mi_graphs, parents, children, MOVING_BASE, h)

%% ************************************************************************
%     Compute sensor positions errors
clc
close all
        
err_rot = zeros(constPar.noj,numel(test_points),3); 
err_loc = zeros(constPar.noj,numel(test_points),3); 


for selector = 1:3

    switch selector
        case 1 % offline+floating base
            gamma_hat = lambda_hat_offline_float(1:7,:);
            rho_hat   = lambda_hat_offline_float(8:13,:);
    
        case 2 % online+floating base
            gamma_hat = lambda_hat_online_float(1:7,:);
            rho_hat   = lambda_hat_online_float(8:13,:);
    
        case 3 % offline+fixed base
            gamma_hat = lambda_hat_offline_fix(1:7,:);
            rho_hat   = lambda_hat_offline_fix(8:13,:);
    end
    
    j_T_sj      = NaN(4,4,constPar.noj);
    s0_T_sj     = NaN(4,4,constPar.noj);
    si_T_sj     = NaN(4,4,constPar.noj);
    si_R_sj     = NaN(3,3,constPar.noj);
    w_T_j       = NaN(4,4,constPar.noj);
    i_T_j       = NaN(4,4,constPar.noj);
    
    si_T_sj_hat = NaN(4,4,constPar.noj);
    %s0_T_sj_hat = NaN(4,4,constPar.noj);
    si_R_sj_hat = NaN(3,3,constPar.noj);
    %si_r_sj     = NaN(3,constPar.noj);
   
    cprintf('*Yellow', ['>> Computing errors for ', num2str(numel(test_points)), ' different configurations\n'])
    for i = 1:numel(test_points)
        disp(['Point ' num2str(i),'/',num2str(numel(test_points))])
        for j=1:constPar.noj
            %disp(panda.Bodies{j}.Name)
            w_T_j(:,:,j)   = getTransform(panda,q_panda(:,test_points(i)),panda.Bodies{j}.Name);
            j_T_sj(:,:,j)  = rt2tr(j_R_sj(:,:,j),j_p_s(1:3,j));
            s0_T_sj(:,:,j) = w_T_j(:,:,j)*j_T_sj(:,:,j);
            if j == 1
                i_T_j(:,:,j)   = w_T_j(:,:,j);
                si_T_sj(:,:,j) = s0_T_sj(:,:,j);
            else
                i_T_j(:,:,j)   = w_T_j(:,:,j-1)\w_T_j(:,:,j);
                si_T_sj(:,:,j) = s0_T_sj(:,:,j-1)\s0_T_sj(:,:,j);
            end
            si_R_sj(:,:,j)     = si_T_sj(1:3,1:3,j);
            
            % Sensor-to-sensor rotation error
            si_R_sj_hat(:,:,j) = p_R_c(q_panda(j,test_points(i)),gamma_hat(:,j)); 
            Q = si_R_sj(:,:,j)*transpose(si_R_sj_hat(:,:,j));
            err_rot(j,i,selector) = acos((trace(Q)-1)/2);        
    
            % Estimated transformation matrix
            si_r_sj            = (-rho_hat(1:3,j)) + si_R_sj_hat(:,:,j)*(rho_hat(4:6,j));
            si_T_sj_hat(:,:,j) = rt2tr(si_R_sj_hat(:,:,j), si_r_sj);
    
            % Sensor-to-sensor location error
            err_loc(j,i,selector) = norm(si_T_sj(1:3,4,j)  - si_T_sj_hat(1:3,4,j));
        end
    end
end
cprintf('*Green', '>> Errors computed!\n')

%% Bar plot of rotation errors ============================================
pause(1)
% close all
clc
interpreter = 'latex';
figure
mean_err_rot = squeeze(mean(err_rot,2));

% error = [ mean_err_rot,mean_err_loc];

h = bar(mean_err_rot);
set(gca,'YScale','log')
% ylabel('$\tilde{\boldmath{r}}_j$','interpreter','latex')
ylabel('$\bar{\delta}$','interpreter',interpreter)
xaxisproperties= get(gca, 'XAxis');
labels = cell(constPar.noj,1);

% for i=1:constPar.noj
%     labels{i} = ['$^{s_{',num2str(i-1),'}}\hat{R}_{s_{',num2str(i),'}}$'];
% end
subindices = ['0ABCDEFGHIJKLMNOPQRSTUVWXYZ'];
for i=1:constPar.noj
    p = i;
    c = i+1;    
    labels{i} = ['$s_{',subindices(p),'}, s_{',subindices(c),'}$'];
end

set(gca, 'XTick', 1:18)
xticklabels(labels)
xaxisproperties.TickLabelInterpreter = interpreter; % latex for x-axis
leg = legend('offline (floating)', 'online (floating)', 'offline (fixed)');
set(leg,'interpreter','latex')

fig = gcf;
ax  = gca;
plt = h;
fcn_scrpt_prepare_graph_science_std(fig, ax, plt, leg, [], 18.3/2, 3, 0.5)
% grid minor
xtickangle(90)
% tightfig(fig)
% ax = gca;
% ax.XAxis.FontSize = 20;
set(leg,'location','northeast')
leg.Orientation = 'vertical';
leg.Interpreter = interpreter;
% ylim([min(mean_err_rot,[],'all'),1E-1])
SAVE_FIG = 0;
if SAVE_FIG == 1
%     export_fig(fullfile(fileparts(matlab.desktop.editor.getActiveFilename),'figures','panda_morphology_errors_real'),'-pdf')
%     close(gcf);
    fig = gcf;           % generate a figure
    tightfig(fig)
    fig.Units = 'centimeters';        % set figure units to cm
    % f.Position = [1203 646 478 174];
    fig.PaperUnits = 'centimeters';   % set pdf printing paper units to cm
    fig.PaperSize = fig.Position(3:4);  % assign to the pdf printing paper the size of the figure
    print(fig,'-dpdf','-r600','-painters',fullfile(fileparts(matlab.desktop.editor.getActiveFilename),'figures','panda_morphology_rot_errors_simulated'));
    disp('Printing done!')
    % close(fig)
end

%% Bar plot of location errors ============================================
pause(1)
close all
clc
interpreter = 'latex';
figure
mean_err_loc = squeeze(mean(err_loc,2));

h = bar(mean_err_loc);
set(gca,'YScale','log')
% ylabel('$\tilde{\boldmath{r}}_j$','interpreter','latex')
ylabel('$\bar{\tilde{r}}$','interpreter',interpreter)
xaxisproperties= get(gca, 'XAxis');
labels = cell(constPar.noj,1);

% for i=1:constPar.noj
%     labels{i} = ['$^{s_{',num2str(i-1),'}}\hat{R}_{s_{',num2str(i),'}}$'];
% end
subindices = ['0ABCDEFGHIJKLMNOPQRSTUVWXYZ'];
for i=1:constPar.noj
    p = i;
    c = i+1;    
    labels{i} = ['$s_{',subindices(p),'}, s_{',subindices(c),'}$'];
end

set(gca, 'XTick', 1:18)
xticklabels(labels)
xaxisproperties.TickLabelInterpreter = interpreter; % latex for x-axis
leg = legend('offline (floating)', 'online (floating)', 'offline (fixed)');
set(leg,'interpreter','latex')

fig = gcf;
ax  = gca;
plt = h;
fcn_scrpt_prepare_graph_science_std(fig, ax, plt, leg, [], 18.3/2, 3, 0.5)
% grid minor
xtickangle(90)
% tightfig(fig)
% ax = gca;
% ax.XAxis.FontSize = 20;
set(leg,'location','northeast')
leg.Orientation = 'vertical';
leg.Interpreter = interpreter;
% ylim([min(mean_err_rot,[],'all'),1E-1])
SAVE_FIG = 0;
if SAVE_FIG == 1
%     export_fig(fullfile(fileparts(matlab.desktop.editor.getActiveFilename),'figures','panda_morphology_errors_real'),'-pdf')
%     close(gcf);
    fig = gcf;           % generate a figure
    tightfig(fig)
    fig.Units = 'centimeters';        % set figure units to cm
    % f.Position = [1203 646 478 174];
    fig.PaperUnits = 'centimeters';   % set pdf printing paper units to cm
    fig.PaperSize = fig.Position(3:4);  % assign to the pdf printing paper the size of the figure
    print(fig,'-dpdf','-r600','-painters',fullfile(fileparts(matlab.desktop.editor.getActiveFilename),'figures','panda_morphology_loc_errors_simulated'));
    disp('Printing done!')
    % close(fig)
end
%% ************************************************************************
%                    PROTO KINEMATICS (COMPARISON)                        *
% *************************************************************************

pause(1)

figure('Color','w')
       
r_hat       = rho_hat;

j_T_sj      = NaN(4,4,constPar.noj);
s0_T_sj     = NaN(4,4,constPar.noj);
si_T_sj     = NaN(4,4,constPar.noj);
si_R_sj     = NaN(3,3,constPar.noj);
w_T_j       = NaN(4,4,constPar.noj);
i_T_j       = NaN(4,4,constPar.noj);

si_T_sj_hat = NaN(4,4,constPar.noj);
s0_T_sj_hat = NaN(4,4,constPar.noj);
si_R_sj_hat = NaN(3,3,constPar.noj);
si_r_sj     = NaN(3,constPar.noj);


q_test = [0 0 0 0 0 deg2rad(120) 0];

% Sphere representing the sensors
[x,y,z] = sphere;
% Scale to desire radius.
radius = 0.01;
x = x * radius;
y = y * radius;
z = z * radius; 

% err_rot = zeros(constPar.noj,numel(test_points)); 
for i = 1
    for j=1:constPar.noj
        w_T_j(:,:,j)   = getTransform(panda, q_test, panda.Bodies{j}.Name);            
        j_T_sj(:,:,j)  = rt2tr(j_R_sj(:,:,j),j_p_s(1:3,j));
        s0_T_sj(:,:,j) = w_T_j(:,:,j)*j_T_sj(:,:,j);
        if j == 1
            i_T_j(:,:,j)   = w_T_j(:,:,j);
            si_T_sj(:,:,j) = s0_T_sj(:,:,j);
        else
            i_T_j(:,:,j)   = w_T_j(:,:,j-1)\w_T_j(:,:,j);
            si_T_sj(:,:,j) = s0_T_sj(:,:,j-1)\s0_T_sj(:,:,j);
        end
        si_R_sj(:,:,j)     = si_T_sj(1:3,1:3,j);
        si_R_sj_hat(:,:,j) = p_R_c( q_test(j), gamma_hat(:,j)); 

        %======================================================================        
        cprintf('*yellow', '>> Rotation error in...\n')
        disp(['JOINT: ' num2str(j)])
        Q = si_R_sj(:,:,j)*transpose(si_R_sj_hat(:,:,j));
        disp(acos((trace(Q)-1)/2))
        %====================================================================== 
    end
end

i = 1;
for j = 1:constPar.noj
    if j < constPar.noj
        si_r_sj(:,j+1)   = (-r_hat(1:3,j+1)) + si_R_sj_hat(:,:,j+1)*(r_hat(4:6,j+1));
    end

    % Estimated transformation from sensor frame to sensor frame ----------
    if j <= 4 && MOVING_BASE == 0
        si_T_sj_hat(:,:,j) = si_T_sj(:,:,j);
    elseif j == 1 && MOVING_BASE == 1
        si_T_sj_hat(:,:,j) = rt2tr(si_R_sj_hat(:,:,j), ...
                                   (-r_hat(1:3,1)) + si_R_sj_hat(:,:,1)*(r_hat(4:6,1)));
    else
        si_T_sj_hat(:,:,j) = rt2tr(si_R_sj_hat(:,:,j), ...
                                   si_r_sj(:,j));
    end

    % Estimated transformation from sensor frame to world frame -----------
    T = eye(4);
    for i=1:j
        T = T*si_T_sj_hat(:,:,i);
    end
    s0_T_sj_hat(:,:,j) = T;     
end

% -------------------------------------------------------------------------
i = 1;

rb = show(panda, q_test, 'PreservePlot', false, 'Frames','off');
alpha(rb,0.1)
hold on
%trplot(eye(4),'rgb','framelabel',['s' num2str(0)],'arrow','width',0.1,'length', 0.05)
trplot(eye(4),'rgb', 'notext','length', 0.05,'thick',3)
for j=1:constPar.noj    
    if j>0
        % Plot indentified IMU coordinate system ------------------------------
        %trplot(s0_T_sj_hat(:,:,j),'rgb','framelabel',['s' num2str(j)],'arrow','width',0.1,'length', 0.05)
        trplot(s0_T_sj_hat(:,:,j),'rgb', 'notext','length', 0.05,'thick',3)
        % Sphere representing the IMU -----------------------------------------
        sph = surf(x+s0_T_sj_hat(1,4,j),y+s0_T_sj_hat(2,4,j),z+s0_T_sj_hat(3,4,j));    
        set(sph,'FaceColor',[1 0 0], ...
          'FaceAlpha',1.0,'FaceLighting','gouraud','EdgeColor','none')        

        % Position vector of axis j in sensor sj frame ------------------------
        point_start     = s0_T_sj_hat(1:3,4,j);
        point_end_jnt_p = s0_T_sj_hat(:,:,j)*[(-r_hat(4:6,j));1];
        v1              = [point_start(1),point_start(2),point_start(3)];
        v2              = [point_end_jnt_p(1),point_end_jnt_p(2),point_end_jnt_p(3)];
        arrow3(v1,v2,['m',':','3.0'],0.2,0.3) 

        % Axis translated to identified joint center point
        jointAxis    = 0.1*gamma_hat(5:7,j);
        rotAxisStart = s0_T_sj_hat(:,:,j)*[(-r_hat(4:6,j));1];
        rotAxisEnd   = s0_T_sj_hat(:,:,j)*[(-r_hat(4:6,j) + jointAxis);1];
        v1=[rotAxisStart(1),rotAxisStart(2),rotAxisStart(3)];
        v2=[rotAxisEnd(1),rotAxisEnd(2),rotAxisEnd(3)];
        v=[v2;v1];
        plot3(v(:,1),v(:,2),v(:,3),'y-','LineWidth',6)
        arrow3(v1,v2,['b-' '4'],0.4,0.5)
        %text(rotAxisEnd(1),rotAxisEnd(2),rotAxisEnd(3),['j' num2str(j)])

        if j < nlinks
            % Axis j+1 in sensor j position vector --------------------------------
            point_start     = s0_T_sj_hat(1:3,4,j);
            point_end_jnt_c = s0_T_sj_hat(:,:,j)*[(-r_hat(1:3,j+1));1];
            v1              = [point_start(1),point_start(2),point_start(3)];
            v2              = [point_end_jnt_c(1),point_end_jnt_c(2),point_end_jnt_c(3)];
            arrow3(v1,v2,['k',':','3.0'],0.2,0.3)    
        end  
    end
end
xlim([-0.3 0.3])
ylim([-0.3 0.3])
zlim([0 1.5])
view([30 30])
axis equal
axis off

view([30 30])
if SAVE_FIG == 1
    export_fig('./figures/panda_morphology_experiment_single_front','-pdf')
    close(gcf);
end


RUN_THIS = 0;

if RUN_THIS == 1
% *************************************************************************
%                    PROTO KINEMATICS (COMPARISON)                        *
% *************************************************************************                                
     
close all
clc
factor = 1;
figure('Color','w')
for f = 1%1:3
    if f<=2
%         switch f
%             case 1
%                 MOVING_BASE = 1;
%                 zeta_hat = lambda_hat_offline_float(1:7,:);
%                 r_hat    = lambda_hat_offline_float(8:13,:);                
%             case 2
%                 MOVING_BASE = 1;
%                 zeta_hat = lambda_rams(1:7,:,end);
%                 r_hat    = lambda_rams(8:13,:,end);
% %             case 3
% %                 MOVING_BASE = 0;
% %                 zeta_hat = lambda_hat_offline_fix(1:7,:);
% %                 r_hat    = lambda_hat_offline_fix(8:13,:);
%         end
        
        gamma_hat = lambda_hat_offline_fix(1:7,:);
        r_hat    = lambda_hat_offline_fix(8:13,:);        

        [x,y,z] = sphere;
        % Scale to desire radius.
        radius = 0.01;
        x = x * radius;
        y = y * radius;
        z = z * radius; 

        
        
        j_T_sj      = NaN(4,4,constPar.noj);
        w_T_j       = NaN(4,4,constPar.noj);
        i_T_j       = NaN(4,4,constPar.noj);
        for j=1:constPar.noj
            w_T_j(:,:,j)   = getTransform(panda,q_ref(:,i)',panda.Bodies{j}.Name);            
            j_T_sj(:,:,j)  = rt2tr(j_R_sj(:,:,j),j_p_s(1:3,j));
            s0_T_sj(:,:,j) = w_T_j(:,:,j)*j_T_sj(:,:,j);
            if j == 1
                i_T_j(:,:,j)       = w_T_j(:,:,j);
                si_T_sj(:,:,j)     = s0_T_sj(:,:,j);
                
            else
                i_T_j(:,:,j)       = w_T_j(:,:,j-1)\w_T_j(:,:,j);
                si_T_sj(:,:,j)     = s0_T_sj(:,:,j-1)\s0_T_sj(:,:,j);
            end
            si_R_sj(:,:,j)     = si_T_sj(1:3,1:3,j);
        end
%T_i_j_aux    = getTransform(panda,q_ref(:,i)',panda.Bodies{j}.Name,panda.Bodies{j-1}.Name);
% Transformation from frame j to frame i = j-1 ----------------
i_T_j              = w_T_j(:,:,j-1)\w_T_j(:,:,j);
i_R_j              = i_T_j(1:3,1:3);
si_R_j             = transpose(j_R_sj(:,:,j-1))*i_R_j;
si_R_sj(:,:,j)     = si_R_j*j_R_sj(:,:,j);
si_R_sj_hat(:,:,j) = p_R_c(q_ref(j,i),gamma_hat(:,j));
s0_T_sj(:,:,j)     = w_T_j(:,:,j)*j_T_sj(:,:,j);                
        
        
        
        si_T_sj_hat = NaN(4,4,constPar.noj);
        s0_T_sj_hat = NaN(4,4,constPar.noj);
        si_r_sj     = NaN(3,constPar.noj);

        
        
        % Transformation from sensor frame to joint frame -----------------
        j_T_sj(:,:,j) = rt2tr(j_R_sj(:,:,j),j_p_s(1:3,j));
        if j == 1
            w_T_j(:,:,j)       = getTransform(panda,q_ref(:,i)',panda.Bodies{j}.Name);
            s0_T_sj(:,:,j)     = w_T_j(:,:,j)*j_T_sj;
            si_R_sj(:,:,j)     = s0_T_sj(1:3,1:3,j);
            si_R_sj_hat(:,:,j) = p_R_c(q_ref(j,i),gamma_hat(:,j));            
        else
            w_T_j(:,:,j)       = getTransform(panda,q_ref(:,i)',panda.Bodies{j}.Name);
            %T_i_j_aux    = getTransform(panda,q_ref(:,i)',panda.Bodies{j}.Name,panda.Bodies{j-1}.Name);
            % Transformation from frame j to frame i = j-1 ----------------
            i_T_j              = w_T_j(:,:,j-1)\w_T_j(:,:,j);
            i_R_j              = i_T_j(1:3,1:3);
            si_R_j             = transpose(j_R_sj(:,:,j-1))*i_R_j;
            si_R_sj(:,:,j)     = si_R_j*j_R_sj(:,:,j);
            si_R_sj_hat(:,:,j) = p_R_c(q_ref(j,i),gamma_hat(:,j));
            s0_T_sj(:,:,j)     = w_T_j(:,:,j)*j_T_sj(:,:,j);        
        end
        
        
        
        % Get all sensor-to-sensor rotation matrices before plotting
        % position vectors
        for j=1:constPar.noj
            % Transformation from frame j to world frame ------------------
            w_T_j(:,:,j) = getTransform(panda,q_ref(:,i)',panda.Bodies{j}.Name);
            % Transformation from frame j to frame i = j-1 ----------------
            i_T_j(:,:,j) = w_T_j(:,:,j-1)\w_T_j(:,:,j);            
            
            if MOVING_BASE == 0 && j == 1

                w_T_j              = getTransform(panda,q_ref(:,1)',panda.Bodies{j}.Name);
                si_R_sj_hat(:,:,j) = w_T_j(1:3,1:3)*j_R_sj(:,:,j);

                %si_R_sj(:,:,j)     = s0_T_sj(1:3,1:3,j);
                si_R_sj_hat(:,:,j)     = si_R_sj(1:3,1:3,j);

            else
                si_R_sj_hat(:,:,j) = p_R_c(q_ref(j,1),gamma_hat(:,j));
            end
        end     
        
        
        subplot(1,3,f)
        rb = show(panda,[q_ref(:,1)'], 'PreservePlot', false, 'Frames','off');
        alpha(rb,0.1)
        hold on
        trplot(eye(4),'rgb','framelabel',['s' num2str(0)],'arrow','width',0.1,'length', 0.05)
        for j = 1:constPar.noj

            if j<nlinks
                si_r_sj(:,j+1)   = (-r_hat(1:3,j+1)) + si_R_sj_hat(:,:,j+1)*(r_hat(4:6,j+1));
            end    

            if j == 1
                if MOVING_BASE == 0
                    w_T_sj             = w_T_j*rt2tr(j_R_sj(:,:,j),j_p_s(:,j)); 
                    si_T_sj_hat(:,:,j) = w_T_sj;         
        %             w_T_j2 = getTransform(panda,q_ref(:,1)',panda.Bodies{2}.Name);
        %             si_T_sj_hat(:,:,j) = rt2tr(w_T_j2(1:3,1:3), w_T_s1);
                else
                    si_T_sj_hat(:,:,j) = rt2tr(si_R_sj_hat(:,:,j), (-r_hat(1:3,1)) + si_R_sj_hat(:,:,1)*(r_hat(4:6,1)));
                end
            else
                si_T_sj_hat(:,:,j) = rt2tr(si_R_sj_hat(:,:,j), si_r_sj(:,j));
            end

            % Transformation from sensor frame to world frame ---------------------
            T = eye(4);
            for i=1:j
                T = T*si_T_sj_hat(:,:,i);
            end
            s0_T_sj_hat(:,:,j) = T;

            % Plot indentified IMU coordinate system ------------------------------
            trplot(T,'rgb','framelabel',['s' num2str(j)],'arrow','width',0.1,'length', 0.05)
        %     trplot(getTransform(panda,q_ref(:,1)',panda.Bodies{j}.Name),'frame',[num2str(j)],'width',0.1,'length', 0.05)    

            % Sphere representing the IMU -----------------------------------------
            sph = surf(x+s0_T_sj_hat(1,4,j),y+s0_T_sj_hat(2,4,j),z+s0_T_sj_hat(3,4,j));    
            set(sph,'FaceColor',[1 0 0], ...
              'FaceAlpha',1.0,'FaceLighting','gouraud','EdgeColor','none')        
            
            % Position vector of axis j in sensor sj frame ------------------------
            point_start     = s0_T_sj_hat(1:3,4,j);
            point_end_jnt_p = s0_T_sj_hat(:,:,j)*[(-r_hat(4:6,j));1];
            v1              = [point_start(1),point_start(2),point_start(3)];
            v2              = [point_end_jnt_p(1),point_end_jnt_p(2),point_end_jnt_p(3)];
            arrow3(v1,v2,['m','--','2.0'],0.2,0.3) 

            % Axis translated to identified joint center point
            jointAxis    = 0.1*gamma_hat(5:7,j);
            rotAxisStart = s0_T_sj_hat(:,:,j)*[(-r_hat(4:6,j));1];
            rotAxisEnd   = s0_T_sj_hat(:,:,j)*[(-r_hat(4:6,j) + jointAxis);1];
            v1=[rotAxisStart(1),rotAxisStart(2),rotAxisStart(3)];
            v2=[rotAxisEnd(1),rotAxisEnd(2),rotAxisEnd(3)];
            v=[v2;v1];
            plot3(v(:,1),v(:,2),v(:,3),'y-','LineWidth',6)
            arrow3(v1,v2,['b','-','3'],0.2,0.3)  
            text(rotAxisEnd(1),rotAxisEnd(2),rotAxisEnd(3),['j' num2str(j)])

            if j < nlinks
                % Axis j+1 in sensor j position vector --------------------------------
                point_start     = s0_T_sj_hat(1:3,4,j);
                point_end_jnt_c = s0_T_sj_hat(:,:,j)*[(-r_hat(1:3,j+1));1];
                v1              = [point_start(1),point_start(2),point_start(3)];
                v2              = [point_end_jnt_c(1),point_end_jnt_c(2),point_end_jnt_c(3)];
                arrow3(v1,v2,['r','--','2.0'],0.2,0.3)    
            end
        end
        xlim([-0.3 0.3])
        ylim([-0.3 0.3])
        zlim([0 1.5])
        view([30 30])
        axis equal
        axis off
        zoom(factor)
    else
        subplot(1,3,f)
gamma_hat = lambda_hat_offline_fix(1:7,:);
r_hat    = lambda_hat_offline_fix(8:13,:);

clear alpha

% q_saveaux = q_franka(:,1);
% q_ref(:,1) = deg2rad(randi(360,7,1))';% [0 0 0 0 0 deg2rad(60) 0];%zeros(7,1);
q_ref(:,1) = [0 0 0 0 0 deg2rad(120) 0];%zeros(7,1);
[x,y,z] = sphere;
% Scale to desire radius.
radius = 0.006;
x = x * radius;
y = y * radius;
z = z * radius;   

si_R_sj     = NaN(3,3,constPar.noj);
si_R_sj_hat = NaN(3,3,constPar.noj);
for i = 1
    rb = show(panda,[q_ref(:,i)'], 'PreservePlot', false, 'Frames','off');
    alpha(rb,0.1)
    hold on    
    trplot(eye(4),'rgb','frame',['s' num2str(0)],'arrow','width',0.1,'length', 0.05)        
    if i ==1000
        test = 1;
    end
    if mod(i,1000) == 0
        disp(['Iteration: ' num2str(i) '/' num2str(samples)])
    end
    w_T_j             = zeros(4,4,nlinks);
    w_T_j(4,4,:)      = 1;
    w_T_j(1:3,1:3,1)  = eye(3);

    s0_T_sj            = zeros(4,4,nlinks);
    s0_T_sj(4,4,:)     = 1;
    s0_T_sj(1:3,1:3,1) = eye(3);    
    for j=1:nlinks
        % Transformation from sensor frame to joint frame -----------------
        j_T_sj = rt2tr(j_R_sj(:,:,j),j_p_s(1:3,j));
        if j == 1
            w_T_j(:,:,j)       = getTransform(panda,q_ref(:,i)',panda.Bodies{j}.Name);
            s0_T_sj(:,:,j)     = w_T_j(:,:,j)*j_T_sj;
            si_R_sj(:,:,j)     = s0_T_sj(1:3,1:3,j);
            si_R_sj_hat(:,:,j) = p_R_c(q_ref(j,i),gamma_hat(:,j));            
        else
            w_T_j(:,:,j)       = getTransform(panda,q_ref(:,i)',panda.Bodies{j}.Name);
            %T_i_j_aux    = getTransform(panda,q_ref(:,i)',panda.Bodies{j}.Name,panda.Bodies{j-1}.Name);
            % Transformation from frame j to frame i = j-1 ----------------
            i_T_j              = w_T_j(:,:,j-1)\w_T_j(:,:,j);
            i_R_j              = i_T_j(1:3,1:3);
            si_R_j             = transpose(j_R_sj(:,:,j-1))*i_R_j;
            si_R_sj(:,:,j)     = si_R_j*j_R_sj(:,:,j);
            si_R_sj_hat(:,:,j) = p_R_c(q_ref(j,i),gamma_hat(:,j));
%             T_0_j(:,:,j) = getTransform(panda,q_ref(:,i)',panda.Bodies{j}.Name);
            s0_T_sj(:,:,j) = w_T_j(:,:,j)*j_T_sj;
        end
% %===========================================        
% cprintf('*yellow', '>> Rotation error in...\n')
% disp(['JOINT: ' num2str(j)])
% disp(si_R_sj(:,:,j)  - si_R_sj_hat(:,:,j))
% %===========================================        
%===========================================        
cprintf('*yellow', '>> Rotation error in...\n')
disp(['JOINT: ' num2str(j)])
disp(si_R_sj(:,:,j)*transpose(si_R_sj_hat(:,:,j)))
disp(norm(si_R_sj(:,:,j)*transpose(si_R_sj_hat(:,:,j)) - (si_R_sj(:,:,j)*transpose(si_R_sj_hat(:,:,j))).*eye(3),'fro'))
%===========================================               
        % CS of the IMUs ---------------------------
        %trplot(T_0_j(:,:,j),'rgb','frame',num2str(j),'arrow','width',0.1,'length', 0.05)
        trplot(s0_T_sj(:,:,j),'rgb','framelabel',['s' num2str(j)],'arrow','width',0.1,'length', 0.05)        
        
        % Sphere representing the IMU -------------------------------------
        sph = surf(x+s0_T_sj(1,4,j),y+s0_T_sj(2,4,j),z+s0_T_sj(3,4,j));    
        set(sph,'FaceColor',[1 0 0], ...
          'FaceAlpha',1.0,'FaceLighting','gouraud','EdgeColor','none')  
      
        
        if(exist('r_hat','var'))            
            % Axis j in sensor j position vector --------------------------
            point_start     = s0_T_sj(1:3,4,j);
            point_end_jnt_p = s0_T_sj(:,:,j)*[(-r_hat(4:6,j));1];
            v1              = [point_start(1),point_start(2),point_start(3)];
            v2              = [point_end_jnt_p(1),point_end_jnt_p(2),point_end_jnt_p(3)];
            arrow3(v1,v2,['m--' '2'],0.2,0.3)
            
            % Axis translated to identified joint center point
            rotAxisStart = s0_T_sj(:,:,j)*[(-r_hat(4:6,j));1];
            rotAxisEnd   = s0_T_sj(:,:,j)*[(-r_hat(4:6,j) + jointAxis);1];
            v1=[rotAxisStart(1),rotAxisStart(2),rotAxisStart(3)];
            v2=[rotAxisEnd(1),rotAxisEnd(2),rotAxisEnd(3)];
            v=[v2;v1];
            plot3(v(:,1),v(:,2),v(:,3),'y-','LineWidth',6)
            arrow3(v1,v2,['b-.' '3'],0.2,0.3)
            text(rotAxisEnd(1),rotAxisEnd(2),rotAxisEnd(3),['j' num2str(j)])
            % quiver3(rotAxisStart(1),rotAxisStart(2),rotAxisStart(3),rotAxis(1),rotAxis(2),rotAxis(3),1,'k-.','LineWidth',2);

            if j<7
                % Axis j+1 in sensor j position vector --------------------
                point_start     = s0_T_sj(1:3,4,j);
                point_end_jnt_c = s0_T_sj(:,:,j)*[(-r_hat(1:3,j+1));1];
                v1              = [point_start(1),point_start(2),point_start(3)];
                v2              = [point_end_jnt_c(1),point_end_jnt_c(2),point_end_jnt_c(3)];
                arrow3(v1,v2,['r--' '2'],0.2,0.3)
            end
        end
    end
end
        xlim([-0.3 0.3])
        ylim([-0.3 0.3])
        zlim([0 1.5])
        view([30 30])
        axis equal
        axis off
        zoom(factor)
    end
end

if SAVE_FIG == 1
    export_fig('./figures/panda_morphology_comparison','-pdf')
    close(gcf);
end

end