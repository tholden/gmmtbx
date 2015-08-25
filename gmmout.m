function varargout = gmmout(varargin)
% GMMOUT M-file for gmmout.fig
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gmmout_OpeningFcn, ...
                   'gui_OutputFcn',  @gmmout_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
     gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

function gmmout_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);
set(handles.values, 'String', varargin{3});
set(handles.conint, 'String', num2str(varargin{4}));
set(handles.jt, 'String', sprintf('%6.4g',varargin{5}));
set(handles.pjt, 'String', sprintf('%6.3g',varargin{6}));         

function varargout = gmmout_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;
function values_CreateFcn(hObject, eventdata, handles)
function values_Callback(hObject, eventdata, handles)
function conint_CreateFcn(hObject, eventdata, handles)
function conint_Callback(hObject, eventdata, handles)