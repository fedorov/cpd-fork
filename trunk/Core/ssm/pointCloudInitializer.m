function varargout = pointCloudInitializer(varargin)
% POINTCLOUDINITIALIZER MATLAB code for pointCloudInitializer.fig
%      POINTCLOUDINITIALIZER, by itself, creates a new POINTCLOUDINITIALIZER or raises the existing
%      singleton*.
%
%      H = POINTCLOUDINITIALIZER returns the handle to a new POINTCLOUDINITIALIZER or the handle to
%      the existing singleton*.
%
%      POINTCLOUDINITIALIZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POINTCLOUDINITIALIZER.M with the given input arguments.
%
%      POINTCLOUDINITIALIZER('Property','Value',...) creates a new POINTCLOUDINITIALIZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pointCloudInitializer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pointCloudInitializer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pointCloudInitializer

% Last Modified by GUIDE v2.5 07-Nov-2013 18:39:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pointCloudInitializer_OpeningFcn, ...
                   'gui_OutputFcn',  @pointCloudInitializer_OutputFcn, ...
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


% --- Executes just before pointCloudInitializer is made visible.
function pointCloudInitializer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pointCloudInitializer (see VARARGIN)

% Choose default command line output for pointCloudInitializer
set(gcf, 'Toolbar', 'figure');
handles.output = hObject;
handles.X = varargin{1};
handles.Xf = varargin{2};
handles.Y = varargin{3};
handles.Yf = varargin{4};
handles.tx = 0.0;
handles.ty = 0.0;
handles.tz = 0.0;
handles.rx = 0.0;
handles.ry = 0.0;
handles.rz = 0.0;
handles.s = 1.0;
handles.Rx = eye(3);
handles.Ry = eye(3);
handles.Rz = eye(3);
handles.R =  eye(3);
handles.t = zeros(1,3);

handles.txscale = 100.0;
handles.tyscale = 100.0;
handles.tzscale = 100.0;

handles.rxscale = pi;
handles.ryscale = pi;
handles.rzscale = pi;

set(handles.slider1,'Value',0.5);
set(handles.slider3,'Value',0.5);
set(handles.slider4,'Value',0.5);
set(handles.slider5,'Value',0.5);
set(handles.slider6,'Value',0.5);
set(handles.slider7,'Value',0.5);
set(handles.slider8,'Value',1.0);

UpdateRigidTransform(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pointCloudInitializer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pointCloudInitializer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% Translation along X-axis in mms
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.tx = handles.txscale*(get(hObject,'Value') - 0.5);
UpdateRigidTransform(hObject, eventdata, handles)


function UpdateRigidTransform(hObject, eventdata, handles)
handles.Rx(2,2) = cos(handles.rx); 
handles.Rx(2,3) = -sin(handles.rx); 
handles.Rx(3,2) = sin(handles.rx); 
handles.Rx(3,3)= cos(handles.rx);
handles.Ry(1,1) = cos(handles.ry);
handles.Ry(1,3) = sin(handles.ry);
handles.Ry(3,1) = -sin(handles.ry);
handles.Ry(3,3)= cos(handles.ry);
handles.Rz(1,1) = cos(handles.rz);
handles.Rz(1,2) = -sin(handles.rz);
handles.Rz(2,1) = sin(handles.rz);
handles.Rz(2,2)= cos(handles.rz);
handles.R = handles.Rz*handles.Ry*handles.Rx;
handles.t = [handles.tx handles.ty handles.tz];
Yt = bsxfun(@plus, handles.s*handles.R*(handles.Y'), handles.t')';
cla;
patch('Vertices',handles.X,'Faces',handles.Xf,'FaceColor','r','FaceAlpha',.5);
hold on;
patch('Vertices',Yt,'Faces',handles.Yf,'FaceColor','b','FaceAlpha',.5);
hold off;
disp('Rotation: '), [handles.rx handles.ry handles.rz]
disp('Translation: '), [handles.tx handles.ty handles.tz]
guidata(hObject, handles);


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
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.ty = handles.tyscale*(get(hObject,'Value') - 0.5);
UpdateRigidTransform(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.tz = handles.tzscale*(get(hObject,'Value') - 0.5);
UpdateRigidTransform(hObject, eventdata, handles);


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
handles.rx = handles.rxscale*(get(hObject,'Value') - 0.5);
UpdateRigidTransform(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.ry = handles.ryscale*(get(hObject,'Value') - 0.5);
UpdateRigidTransform(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider7_Callback(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.rz = handles.rzscale*(get(hObject,'Value') - 0.5);
UpdateRigidTransform(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function slider7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider8_Callback(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.s = get(hObject,'Value');
UpdateRigidTransform(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function slider8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
