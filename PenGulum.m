function PenGulum
%%% Graphic user interface for the pendulum test model

%% Init
% Create handle hGui and Figure
hGui.Fig = figure('Name','PenGulum, a Gui for pendulum simulations',...
    'NumberTitle','off','Resize','off','Position',[100 50 900 650]);
data = guihandles(hGui.Fig);

% Create axis
hGui.Axis = axes('Parent',hGui.Fig,'Position',[.1 .45 .58 .5]);

% Set axis properties
xlabel(hGui.Axis,'time [s]');
ylabel(hGui.Axis,'\theta [^o]')
set(hGui.Axis,'xlim',[0 12],'ylim',[-140 20],'XTick',0:2:12,...
    'YTick',-120:20:0);

%% Controls
% Create parameters controls (text, edit text and slider)
controlParameters = {'Mass (kg)','Length (m)','Inertia (kgm^2)',...
    'Baseline tone (Nm)','Damping (Nm.s)','Duration (s)','Feedback gain','Derivative feedback gain',...
    'Starting Angle'};
for i = 1:size(controlParameters,2)
    hGui.controlText(i) = uicontrol(hGui.Fig,'Style','text',...
        'String',controlParameters(i),'Position',[630 620-(70*(i-1)) 130 20]);
    
    hGui.controlEdit(i) = uicontrol(hGui.Fig,'Style','edit',...
        'Position',[630 600-(70*(i-1)) 130 20],'Tag',num2str(i),...
        'Callback',@changeValue2);
    
    hGui.controlSlide(i) = uicontrol(hGui.Fig,'Style','slider',...
        'Position',[630 580-(70*(i-1)) 130 20],...
        'value',0, 'min',0, 'max',1,'SliderStep',[0.01 0.1],...
        'Tag',num2str(i),'Callback',@changeValue);
end

hGui.controlText(10) = uicontrol(hGui.Fig,'Style','text',...
    'String','Scale','Position',[765 620-(70) 130 20]);

hGui.controlEdit(10) = uicontrol(hGui.Fig,'Style','edit',...
    'Position',[765 600-(70) 130 20],'Tag','10',...
    'Callback',@changeValue2);

hGui.controlSlide(10) = uicontrol(hGui.Fig,'Style','slider',...
    'Position',[765 580-(70) 130 20],...
    'value',0, 'min',0, 'max',1,'SliderStep',[0.01 0.1],...
    'Tag','10','Callback',@changeValue);

hGui.controlText(11) = uicontrol(hGui.Fig,'Style','text',...
    'String','Scale','Position',[765 620-(140) 130 20]);

hGui.controlEdit(11) = uicontrol(hGui.Fig,'Style','edit',...
    'Position',[765 600-(140) 130 20],'Tag','11',...
    'Callback',@changeValue2);

hGui.controlSlide(11) = uicontrol(hGui.Fig,'Style','slider',...
    'Position',[765 580-(140) 130 20],...
    'value',0, 'min',10, 'max',1000,'SliderStep',[1 1],...
    'Tag','11','Callback',@changeValue);

% Model selection
hGui.modelFlags = uibuttongroup('Position',[0.04 0.151 .49 .21],...
    'SelectionChangedFcn',@modelSelection);
hGui.flag1 = uicontrol(hGui.modelFlags,'Style','radiobutton',...
    'Position',[10 70 300 30],'Tag','1',...
    'String','State derivative of pendulum model with SRS');
hGui.flag2 = uicontrol(hGui.modelFlags,'Style','radiobutton',...
    'Position',[10 40 380 30],'Tag','2',...
    'String','State derivative of pendulum model with SRS and force-related reflexes');
hGui.flag3 = uicontrol(hGui.modelFlags,'Style','radiobutton',...
    'Position',[10 10 380 30],'Tag','3',...
    'String','State derivative of pendulum model with SRS and velocity-related reflexes');

% Autofit algorithm
hGui.algFlags = uibuttongroup('Position',[0.04 0.028 .49 .12],...
    'SelectionChangedFcn',@algSelection);
hGui.flag4 = uicontrol(hGui.algFlags,'Style','radiobutton',...
    'Position',[10 40 300 30],'Tag','0',...
    'String','Fmincon');
hGui.flag5 = uicontrol(hGui.algFlags,'Style','radiobutton',...
    'Position',[10 10 380 30],'Tag','3',...
    'String','BADS');

%% Buttons
% Load button
hGui.Load = uicontrol('style','pushbutton','string','Load recorded data',...
    'Position',[485 70 130 40],'Callback',@loaddata);

% Save button
hGui.Save = uicontrol('style','pushbutton','string','Save figure',...
    'Position',[485 20 130 40]);

% Auto-Optimization button
hGui.autoFit = uicontrol('style','pushbutton','string','Autofit',...
    'Position',[485 120 130 40],'Enable','off','Callback',@autofit);

% Update button
hGui.Update = uicontrol('style','pushbutton','string','Update plot',...
    'Position',[485 170 130 40],'Callback',@simulation);

%%
%%
% Plot handles
hold on
data.simulationPlot = plot(0,0);
data.recordedPlot = plot(0,0);

% Init parameters
data.model.flag = 1;
data.model.alg = 0;
data.model.m = 4.9575;
data.model.lc = 0.2556;
data.model.I = 0.4455;
data.model.d = 0.15;
data.model.knee_r_range = [-1.9199,0];
data.model.klim = 3;
data.model.kM = 10;
data.model.kF = 0.00001;
data.model.kdF = 0.00001;
data.model.delta_theta_crit = 0.0262;
data.model.dur = 12;
data.model.Tb = 0.95;
data.model.q0 = 0;
data.model.scale = 1;
data.model.Fs = 120;

data.model.vectorID = {'m','lc','I','Tb','d','dur','kF','kdF','q0','scale','Fs'};

data.model.vector = [data.model.m data.model.lc data.model.I data.model.Tb...
    data.model.d data.model.dur data.model.kF data.model.kdF data.model.q0...
    data.model.scale data.model.Fs];
data.model.lbounds = [0.1 0.1 0.1 0 0 1 0.00001 0.00001 -90 0.1 10];
data.model.ubounds = [10 2 2 10 3 20 3 1 10 2 1000];

data.simulation.q = [];
data.simulation.time = [];
data.recorded.q = [];
data.recorded.time = [];

%Update values on gui
for i = 1:size(hGui.controlEdit,2)
    set(hGui.controlEdit(i),'String',num2str(data.model.vector(i)));
    set(hGui.controlSlide(i),'Value',data.model.vector(i),...
        'min',data.model.lbounds(i), 'max',data.model.ubounds(i),...
        'SliderStep',[0.01 0.1]);
end

% update data
guidata(hGui.Fig,data);
simulation(hGui.Fig,[]);

%% Callbacks
    %% Update values
    function changeValue(hObject,event)
        handle = guidata(hObject);
        set(findobj('Style','edit','Tag',get(hObject,'Tag')),...
            'String',num2str(get(hObject,'Value')));
        handle.model.vector(str2double(get(hObject,'Tag'))) = get(hObject,'Value');
        
        if str2double(get(hObject,'Tag'))==10            
            handle.model.vector(1) = (get(hObject,'Value')^2)*data.model.m;
            set(findobj('Style','edit','Tag','1'),'String',...
                handle.model.vector(1));
            
            handle.model.vector(2) = (get(hObject,'Value'))*data.model.lc;
            set(findobj('Style','edit','Tag','2'),'String',...
                handle.model.vector(2));
            
            handle.model.vector(3) = (get(hObject,'Value')^4)*data.model.I;     
            set(findobj('Style','edit','Tag','3'),'String',...
                handle.model.vector(3));                  
        end
        
        guidata(hObject, handle)
        simulation(hObject,event);
    end

    function changeValue2(hObject,event)
        handle = guidata(hObject);
        set(findobj('Style','slider','Tag',get(hObject,'Tag')),...
            'Value',str2double(get(hObject,'String')));
        handle.model.vector(str2double(get(hObject,'Tag'))) = str2double(get(hObject,'String'));
        
        if str2double(get(hObject,'Tag'))==10            
            handle.model.vector(1) = (str2double(get(hObject,'String'))^2)*data.model.m;
            set(findobj('Style','edit','Tag','1'),'String',...
                handle.model.vector(1));
            
            handle.model.vector(2) = (str2double(get(hObject,'String')))*data.model.lc;
            set(findobj('Style','edit','Tag','2'),'String',...
                handle.model.vector(2));
            
            handle.model.vector(3) = (str2double(get(hObject,'String'))^4)*data.model.I;     
            set(findobj('Style','edit','Tag','3'),'String',...
                handle.model.vector(3));                  
        end
        
        guidata(hObject, handle)
        simulation(hObject,event);
    end
    
    %% Select Model
    function modelSelection(hObject,event)
        handle = guidata(hObject);
        handle.model.flag = str2double(event.NewValue.Tag);
        guidata(hObject, handle)
    end

    %% Select Algorithm
    function algSelection(hObject,event)
        handle = guidata(hObject);
        handle.model.alg = str2double(event.NewValue.Tag);
        guidata(hObject, handle)
    end

    %% Load recorded data
    function loaddata(hObject,~)
        [basename, folder] = uigetfile('*.mat');
        if basename == 0
        else
            fullfilename = fullfile(folder,basename);
            load(fullfilename,'q','time')
            handle = guidata(hObject);
            set(handle.recordedPlot,'xdata',time,'ydata',q+5,...
                'Color',[0.5 0.5 0.5],'Linewidth',2)
            handle.recorded.q = q+5;
            handle.recorded.time = time;            
            xlabel(gca,'time [s]');
            ylabel(gca,'\theta [^o]')
            set(gca,'xlim',[0 12],'ylim',[-140 20],'XTick',0:2:12,...
                'YTick',-120:20:0);
            set(findobj('Style','pushbutton','string','Autofit'),...
                'Enable','on')
            set(findobj('Style','edit','Tag','11'),'String',...
                num2str(1/(time(2)-time(1))));
            handle.model.vector(11) = 1/(time(2)-time(1));
            guidata(hObject, handle)
        end
    end
    
    %% Autofit Data
    function autofit(hObject,~)
        handle = guidata(hObject);
        set(gcf, 'pointer', 'watch')
        drawnow
        
        lower_bound = 0; % Lower bound of gains and time delay
        upper_bound = 3;
        gain = [handle.model.vector(strcmp(cellstr(handle.model.vectorID),'kF'))
            handle.model.vector(strcmp(cellstr(handle.model.vectorID),'kdF'))];
        glb = [0.00001 0.00001];
        gub = [2 0.5];
        Tb = handle.model.vector(strcmp(cellstr(handle.model.vectorID),'Tb'));
                
        q = handle.recorded.q;
        time = handle.recorded.time;
        options = [];
        alg = handle.model.alg;
        if alg == 3
            if exist('bads','file')~=2
                warndlg('BADS not installed, please download it from https://github.com/lacerbi/bads')
                return
            end
        end
        flags = handle.model.flag+alg;
        optimize_options = optimset('display','off','algorithm',...
    'interior-point','TolX',1e-9,'TolFun',1e-8,'MaxFunEvals',10000); 

        switch flags
            %fmincon
            case 1
                [opt_gains,~] = fmincon(@(datav)costFunction(datav,q,time,flags,hObject),...
                    Tb,[],[],[],[],lower_bound, upper_bound,[],optimize_options);
            case 2
                [opt_gains,~] = fmincon(@(datav)costFunction(datav,q,time,flags,hObject),...
                    [Tb gain'],[],[],[],[],[lower_bound glb],[upper_bound gub],[],optimize_options);
                handle.model.vector(strcmp(cellstr(handle.model.vectorID),'kF')) = opt_gains(2);
                handle.model.vector(strcmp(cellstr(handle.model.vectorID),'kdF')) = opt_gains(3);
            case 3
                [opt_gains,~] = fmincon(@(datav)costFunction(datav,q,time,flags,hObject),...
                    [Tb gain'],[],[],[],[],[lower_bound glb], [upper_bound gub],[],optimize_options);
                handle.model.vector(strcmp(cellstr(handle.model.vectorID),'kF')) = opt_gains(2);
                handle.model.vector(strcmp(cellstr(handle.model.vectorID),'kdF')) = opt_gains(3);
            %bads
            case 4
                [opt_gains,~] = bads(@(datav)costFunction(datav,q,time,flags,hObject),...
                                    Tb,lower_bound, upper_bound,[],[],[],options);
            case 5
                [opt_gains,~] = bads(@(datav)costFunction(datav,q,time,flags,hObject),...
                    [Tb gain'],[lower_bound glb],[upper_bound gub],[],[],[],options);
                handle.model.vector(strcmp(cellstr(handle.model.vectorID),'kF')) = opt_gains(2);
                handle.model.vector(strcmp(cellstr(handle.model.vectorID),'kdF')) = opt_gains(3);
            case 6
                [opt_gains,~] = bads(@(datav)costFunction(datav,q,time,flags,hObject),...
                    [Tb' gain],[lower_bound glb], [upper_bound gub],[],[],[],options);
                handle.model.vector(strcmp(cellstr(handle.model.vectorID),'kF')) = opt_gains(2);
                handle.model.vector(strcmp(cellstr(handle.model.vectorID),'kdF')) = opt_gains(3);
        end
        
        handle.model.vector(strcmp(cellstr(handle.model.vectorID),'Tb'))= opt_gains(1);        
        handle.model.vector(strcmp(cellstr(handle.model.vectorID),'dur')) = time(end);
        handle.model.vector(strcmp(cellstr(handle.model.vectorID),'q0')) = q(1);  
        
        for ii = 1:size(handle.model.vector,2)
            set(findobj('Style','edit','Tag',num2str(ii)),...
                'String',num2str(handle.model.vector(ii)));
            set(findobj('Style','slider','Tag',num2str(ii)),...
                'Value',handle.model.vector(ii));
        end
        set(gcf, 'pointer', 'arrow')
        guidata(hObject, handle)
        simulation(hObject);
    end

    %% Run model simulation
    function simulation(hObject,~)        
        handle = guidata(hObject);
        set(gcf, 'pointer', 'watch')
        drawnow

        load inputdata.mat inputdata;
        inputdata.mass = handle.model.vector(strcmp(cellstr(handle.model.vectorID),'m'));
        inputdata.lc = handle.model.vector(strcmp(cellstr(handle.model.vectorID),'lc'));
        inputdata.I = handle.model.vector(strcmp(cellstr(handle.model.vectorID),'I'));
        inputdata.d = handle.model.vector(strcmp(cellstr(handle.model.vectorID),'d'));
        inputdata.kF = handle.model.vector(strcmp(cellstr(handle.model.vectorID),'kF'));
        inputdata.fdF = handle.model.vector(strcmp(cellstr(handle.model.vectorID),'kdF'));
        inputdata.kl = handle.model.vector(strcmp(cellstr(handle.model.vectorID),'kdF'));
        inputdata.kv = handle.model.vector(strcmp(cellstr(handle.model.vectorID),'kF'));
        dur = handle.model.vector(strcmp(cellstr(handle.model.vectorID),'dur'));
        inputdata.Tb = handle.model.vector(strcmp(cellstr(handle.model.vectorID),'Tb'));
        q0 = handle.model.vector(strcmp(cellstr(handle.model.vectorID),'q0'));
        x0 = [q0/180*pi 0 0]';
        tspan = [0 dur]';
        inputdata.theta0 = x0(1);
        options = [];
        lags = 0.05;
        Fs = handle.model.vector(strcmp(cellstr(handle.model.vectorID),'Fs'));
        
        global hist params;
        
        switch handle.model.flag
            case {1,4}
                hist = 1;
                [tM,sim] = ode15s(@pendulumStateDerivative_SRS, tspan, x0, options, inputdata);
                set(handle.simulationPlot,'xdata',tM,'ydata',sim(:,1)*180/pi,'Color',[0 0 0],'Linewidth',2);
            case {2,5}
                x0 = [q0/180*pi 0 0 0];
                params = inputdata;
                hist = 1;
                solls = ddensd(@pendulumStateDerivative_SRS_Ffb, lags, lags, x0, tspan);
                sim = deval(solls,tspan(1):1/Fs:tspan(2))';
                %sim = interp1(solls.x,solls.y(1,:)',tspan(1):0.01:tspan(2))';
                set(handle.simulationPlot,'xdata',tspan(1):1/Fs:tspan(2),'ydata',...
                    sim(:,1)*180/pi,'Color',[0 0 0],'Linewidth',2);
            case {3,6}
                x0 = [q0/180*pi 0 0 0];
                params = inputdata;
                hist = 1;
                solls_v = ddensd(@pendulumStateDerivative_SRS_vfb, lags, lags, x0, tspan);
                sim = deval(solls_v,tspan(1):1/Fs:tspan(2))';
                %sim = interp1(solls_v.x,solls_v.y(1,:)',tspan(1):1/Fs:tspan(2))';
                set(handle.simulationPlot,'xdata',tspan(1):1/Fs:tspan(2),'ydata',...
                    sim(:,1)*180/pi,'Color',[0 0 0],'Linewidth',2);
        end
        xlabel(gca,'time [s]');
        ylabel(gca,'\theta [^o]')
        set(gca,'xlim',[0 12],'ylim',[-140 20],'XTick',0:2:12,...
            'YTick',-120:20:0);
        handle.simulation.q = sim(:,1)*180/pi;
        handle.simulation.time = [];
        set(gcf, 'pointer', 'arrow')
        guidata(hObject, handle)
    end

    %% Cost Function for Autofit
    function J = costFunction(datav,recorded,time,flags,hObject)
        handle = guidata(hObject);
        
        load inputdata.mat inputdata;
        inputdata.mass = handle.model.vector(strcmp(cellstr(handle.model.vectorID),'m'));
        inputdata.lc = handle.model.vector(strcmp(cellstr(handle.model.vectorID),'lc'));
        inputdata.I = handle.model.vector(strcmp(cellstr(handle.model.vectorID),'I'));
        inputdata.d = handle.model.vector(strcmp(cellstr(handle.model.vectorID),'d'));
        
        x0 = [recorded(1)/180*pi 0 0 0]';
        tspan = time;
        inputdata.theta0 = x0(1);
        options = [];
        lags = 0.05;
        
        global hist params;
        
        inputdata.Tb = datav(1);
        
        hist = 1;
        params = inputdata;
        
        switch flags
            case {1,4}
                [~,sim] = ode15s(@pendulumStateDerivative_SRS, tspan, x0(1:3), options, inputdata);
                
            case {2,5}
                inputdata.kF = datav(2);
                inputdata.fdF = datav(3);
                solls = ddensd(@pendulumStateDerivative_SRS_Ffb, lags, lags, x0, tspan);
                sim = deval(solls,tspan)';
                %sim = interp1(solls.x,solls.y(1,:)',tspan)';
            case {3,6}
                inputdata.kl = datav(3);
                inputdata.kv = datav(2);
                solls_v = ddensd(@pendulumStateDerivative_SRS_vfb, lags, lags, x0, tspan);
                sim = deval(solls_v,tspan)';
                %sim = interp1(solls_v.x,solls_v.y(1,:)',tspan)';                
        end
        
        error = recorded - (sim(:,1)'*180/pi);
        J = sum(error.^2)/length(error); %Calculate sum of squared errors
    end
end











