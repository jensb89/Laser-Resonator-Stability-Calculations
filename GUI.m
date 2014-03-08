function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 07-Mar-2014 14:55:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

set(handles.edit1,'String',num2str(3.2E-3)); %dK
set(handles.edit2,'String',num2str(0.1)); %R
set(handles.edit3, 'String', num2str(1.84)); %L_ges

warning('off','MATLAB:plot:IgnoreImaginaryXYPart');

% Update handles structure
guidata(hObject, handles);

pushbutton1_Callback(hObject, eventdata, handles)

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;





% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dK = str2double(get(handles.edit1,'String'));
R = str2double(get(handles.edit2,'String'));
L_ges = str2double(get(handles.edit3,'String'));

res = Resonator(dK,1.76,atan(1.76),R,L_ges,800E-9,'s');

theta_ges=res.getTheta();
fprintf('Theta_ges = %1.2fdeg\n',theta_ges/pi*180)

etaL = get(handles.slider1,'Value');
etaTheta = get(handles.slider4,'Value');

res.setArmLengthDist(etaL)%63); %Arm1 = val*(L1+L2); Arm2 = (1-val)*(L1+L2)
res.setThetaDist(etaTheta)%0.63); %Theta1 = val*theta_ges, theta2 = (1-val)*theta_ges
res.setRelCrystalPos(0.5); %Crystal in the middle of the cavity

L1 = res.L1;
L2 = res.L2;
Theta1 = res.theta1*180/pi;
Theta2 = res.theta2*180/pi;

set(handles.text4,'String',['L1 = ',num2str(L1),' L2 = ', num2str(L2)]);
set(handles.text5,'String',['Theta1 = ',num2str(Theta1),' Theta2 = ', num2str(Theta2)]);

delta = linspace(0,10E-3); 

for i=1:length(delta)
    res.setDelta(delta(i)); %deviation from Cavity distance 2*f+dK
    
    [~,~,ws_arm1(i),ws_arm2(i),Ss_arm1(i),Ss_arm2(i)] = res.calcRoundTrip();
    res.setPolarisation('p');
    [~,~,wp_arm1(i),wp_arm2(i),Sp_arm1(i),Sp_arm2(i)] = res.calcRoundTrip();
    res.setPolarisation('s');
    
end

beam_profile = get(handles.checkbox1,'Value');
beam_profile_pump = get(handles.checkbox2,'Value');

%x = (res.R+delta)/1E-3;
x=delta/1E-3;

axes(handles.axes1);
if beam_profile
    plot(x,ws_arm1/1E-3,x,wp_arm1/1E-3,x,ws_arm2/1E-3,x,wp_arm2/1E-3)
    legend('S','P','S','P');
else
    plot(x,ws_arm1/1E-3,x,wp_arm1/1E-3);
    legend('S','P');
    title(['Arm1 = ',num2str(res.L1),'mm']);
end
xlabel('Delta / mm');
ylabel('beam diameter / mm');


axes(handles.axes2);
if ~beam_profile
    plot(x,ws_arm2/1E-3,x,wp_arm2/1E-3);
    xlabel('Delta / mm');
    ylabel('beam diameter / mm');
    legend('S','P');
    title(['Arm2 = ',num2str(res.L2),'mm']);
else
    val = str2double(get(handles.edit4,'String'));
    res.setDelta(val);
    res.calcRoundTrip();
    res.setPolarisation('s')
    [zs,ws] = res.beamInsideCrystal();
    res.setPolarisation('p')
    res.calcRoundTrip();
    [zp,wp] = res.beamInsideCrystal();
    plot(zs,ws,'b',zp,wp,'g'); 
    xlabel('Distance inside Crystal / m'); ylabel('beam diameter / m');
    crystal_end = res.dK/(cos(asin(sin(res.thetaK)/res.nK)));
    hold on 
    plot(crystal_end*ones(1,length(zs)),linspace(0,max(wp),length(zs)),'r');
    legend('S','P');
    hold off
    
    if beam_profile_pump
        hold on
        L3 = get(handles.slider6,'Val');
        res.setPolarisation('s')
        [zpumps,wpumps] = res.PumpBeamInsideCrystal(532E-9,2.1E-3,-1,10E-2,9E-3,11E-2,L3);
        res.setPolarisation('p')
        [zpumpp,wpumpp] = res.PumpBeamInsideCrystal(532E-9,2.1E-3,-1,10E-2,9E-3,11E-2,L3);
        plot(zpumps,wpumps,'b-.',zpumpp,wpumpp,'g-.');
        hold off
    end
        
end



% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = get(hObject,'Value');
set(handles.text7,'String',num2str(val));
pushbutton1_Callback(hObject, eventdata, handles);

% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

val = get(hObject,'Value')
set(handles.text8,'String',num2str(val));
pushbutton1_Callback(hObject, eventdata, handles);

% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value')
    set(handles.slider6,'Visible','on');
    set(handles.text11,'Visible','on');
    set(handles.text12,'Visible','on');
else
    set(handles.slider6,'Visible','off');
    set(handles.text11,'Visible','of');
    set(handles.text12,'Visible','off');
end

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value')
    set(handles.text6,'Visible','on');
    set(handles.edit4,'Visible','on');
    set(handles.checkbox2,'Visible','on');
else
    set(handles.text6,'Visible','off');
    set(handles.edit4,'Visible','off');
    set(handles.checkbox2,'Value',0);
    checkbox2_Callback(hObject, eventdata, handles);
    set(handles.checkbox2,'Visible','off');
end


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = get(hObject,'Value');
set(handles.text12,'String',num2str(val));


%____________________________________________________________________
%
%


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end








% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
