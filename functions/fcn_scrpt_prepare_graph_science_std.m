function fcn_scrpt_prepare_graph_science_std(fig, ax, plt, leg, tx, text_width, k_scaling, k_width_height)
    %% IEEE Standard Figure Configuration - Version 1.0

    % run this code before the plot command

    %%
    % According to the standard of IEEE Transactions and Journals: 

    % Times New Roman is the suggested font in labels. 

    % For a singlepart figure, labels should be in 8 to 10 points,
    % whereas for a multipart figure, labels should be in 8 points.

    % Width: column width: 8.8 cm; page width: 18.1 cm.

    %% width & hight of the figure
%     k_scaling = 5;          % scaling factor of the figure

    % (You need to plot a figure which has a width of (8.8 * k_scaling)
    % in MATLAB, so that when you paste it into your paper, the width will be
    % scalled down to 8.8 cm  which can guarantee a preferred clearness.
    % k_width_height = 2;      % width:hight ratio of the figure
    % width = 8.8 * k_scaling;
    % height = width / k_width_height;
    
%     text_width     = 8.4;%13.97;
%     k_width_height = 0.6;      % width:hight ratio of the figure    
%     k_width_height = 1;      % width:hight ratio of the figure    
    width          = text_width * k_scaling;
    height         = width * k_width_height;

    %% figure margins
    top    = 0.5;  % normalized top margin
    bottom = 3;	   % normalized bottom margin
    left   = 3.5;  % normalized left margin
    right  = 1;    % normalized right margin

    %% set default figure configurations

    % Figure config
    if exist('fig','var')
        set(fig,'Units','centimeters');
        set(fig,'Position',[0 0 width height]); %[1.7727 0 35.1896 23.3098]);
%         set(fig,'Position',[1.7727 0 35.1896 23.3098]); %);
        set(fig,'PaperPositionMode','auto');
        set(fig,'color','w');
    end

    % Plot config
    if exist('plt','var')
        for i = 1:numel(plt)
            if strcmp(get(plt(i),'Type'),'graphplot')
                set(plt(i),'NodeFontSize',10*k_scaling)
            else
                set(plt(i),'LineWidth',0.5*k_scaling);
            end
        end
    end

    %Axes config
    if exist('ax','var') && ~isempty(ax)
        %set(ax,'LineWidth',0.25*k_scaling);
        set(ax,'LineWidth',0.75*k_scaling);
        set(ax,'GridLineStyle',':');
        set(ax,'XGrid','on');
        set(ax,'YGrid','on');
        set(ax,'ZGrid','on');
        set(ax,'FontName','Helvetica');
        set(ax,'FontSize',10*k_scaling);
        set(ax,'Units','normalized');
        set(ax,'LineStyleOrder','-|--');
        set(ax,'TickDir','in');  
        set(ax,'Units','normalized');
    %     set(ax,'Position',[left/width bottom/hight (width-left-right)/width  (hight-bottom-top)/hight]);    
    %     set(ax,'Position',[.12 .2 .75 .7]);   
    %     set(ax,'Position',[.12 .2 .75 .6]);
    end

    % Text cofig
    if exist('tx','var')
        set(tx,'FontName', 'Helvetica' );
        set(tx,'FontSize',10*k_scaling);
    end

    % Legend config
    if exist('leg','var')
        set(leg,'Interpreter','latex')
        set(leg,'FontName',  'Helvetica' );
        set(leg,'FontSize',10*k_scaling);
        set(leg,'Location','southeast');
        set(leg,'Box','off');
        set(leg,'Orientation','vertical');
    end
end
