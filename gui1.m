function varargout = gui1(varargin)
% GUI1 M-file for gui1.fig
%      GUI1, by itself, creates a new GUI1 or raises the existing
%      singleton*.
%
%      H = GUI1 returns the handle to a new GUI1 or the handle to
%      the existing singleton*.
%
%      GUI1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI1.M with the given input arguments.
%
%

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui1_OpeningFcn, ...
                   'gui_OutputFcn',  @gui1_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before gui1 is made visible.
function gui1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui1 (see VARARGIN)

% Choose default command line output for gui1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

a=ones(256,256);
axes(handles.one);
imshow(a);
axes(handles.two);
imshow(a);
axes(handles.three);
imshow(a);


% --- Outputs from this function are returned to the command line.
function varargout = gui1_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile(('*.tif;*.jpg;*.bmp'), 'Pick an Image');
 
    if isequal(filename,0) | isequal(pathname,0)
        
       warndlg('User pressed cancel')
    else
        
        
        a=imread(filename);
        a=imresize(a,[256,256]);
        axes(handles.one);
        imshow(a);

        handles.filename=filename;
        handles.a=a;
        % Update handles structure
guidata(hObject, handles);
        
    end
 
    
                                                             
% --- Executes on button press in encoding.
function wavelets_Callback(hObject, eventdata, handles)
% hObject    handle to wavelets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% contents= (get(handles.popupmenu1,'Value'));
contents= get(handles.filter,'Value')
filename=handles.filename;
X=imread(filename);
X=imresize(X,[256,256]);
[r c p]=size(X);
if p==3
    
X=X(:,:,1);
end

if contents==1
filter = 'haar';
X=double(X);
sX = size(X);
% record start time


% perform single-level decomposition of X. 
[cA1,cH1,cV1,cD1] = dwt2(X, filter);

% record wavelet transform time and output
% dwttime = toc;
% fprintf('\nDWT time:    %6.3f seconds\n', dwttime);


% put it into the tree structure
dec2d = [... 
        cA1,     cH1;     ... 
        cV1,     cD1      ... 
        ];
 
 axes(handles.two);
 imshow(dec2d,[]);
 title('Transformed Image');
handles.dec2d=dec2d;
 guidata(hObject, handles);

else 
 
     x=double(X);
[r,c]=size(x);

% ROW TRANSFORM

for i=1:r
for j=1:c/2
k=j-1;
if j==1
lo(i,j)=(2*x(i,2*j-1)+6*x(i,2*j)+2*x(i,2*j+1)-x(i,2*j+2)+2)/4;
hi(i,j)=(2*x(i,2*k+1)-x(i,2*k+2))/2;
elseif j==c/2
lo(i,j)=(-x(i,2*j-2)+2*x(i,2*j-1)+6*x(i,2*j)+2)/4;
hi(i,j)=(-x(i,2*k)+2*x(i,2*k+1)-x(i,2*k+2))/2;
else
lo(i,j)=(-x(i,2*j-2)+2*x(i,2*j-1)+6*x(i,2*j)+2*x(i,2*j+1)-x(i,2*j+2)+2)/4;
hi(i,j)=(-x(i,2*k)+2*x(i,2*k+1)-x(i,2*k+2))/2;
end
j=j+1;
end
i=i+1;
end
lo=floor(lo);
hi=floor(hi);

% COLUMN TRANSFORM
for j=1:c/2
for i=1:r/2
k=i-1;
if i==1
a(i,j)=(2*lo(2*i-1,j)+6*lo(2*i,j)+2*lo(2*i+1,j)-lo(2*i+2,j)+2)/4;
h(i,j)=(2*lo(2*k+1,j)-lo(2*k+2,j))/2;
v(i,j)=(2*hi(2*i-1,j)+6*hi(2*i,j)+2*hi(2*i+1,j)-hi(2*i+2,j)+2)/4;
d(i,j)=(2*hi(2*k+1,j)-hi(2*k+2,j))/2; 
elseif i==c/2
a(i,j)=(-lo(2*i-2,j)+2*lo(2*i-1,j)+6*lo(2*i,j)+2)/4;
h(i,j)=(-lo(2*k,j)+2*lo(2*k+1,j)-lo(2*k+2,j))/2;
v(i,j)=(-hi(2*i-2,j)+2*hi(2*i-1,j)+6*hi(2*i,j)+2)/4;
d(i,j)=(-hi(2*k,j)+2*hi(2*k+1,j)-hi(2*k+2,j))/2;
else
a(i,j)=(-lo(2*i-2,j)+2*lo(2*i-1,j)+6*lo(2*i,j)+2*lo(2*i+1,j)-lo(2*i+2,j)+2)/4;
h(i,j)=(-lo(2*k,j)+2*lo(2*k+1,j)-lo(2*k+2,j))/2;
v(i,j)=(-hi(2*i-2,j)+2*hi(2*i-1,j)+6*hi(2*i,j)+2*hi(2*i+1,j)-hi(2*i+2,j)+2)/4;
d(i,j)=(-hi(2*k,j)+2*hi(2*k+1,j)-hi(2*k+2,j))/2;
end
i=i+1;
end
j=j+1;
end

a=floor(a);
h=floor(h);
v=floor(v);
d=floor(d);
dec2d=[a,h;v,d];
axes(handles.two);
imshow(dec2d,[]);
handles.dec2d=dec2d;
title('Transformed Image');
end
 guidata(hObject, handles);

% --- Executes on button press in inverse_wave.
function encoding_Callback(hObject, eventdata, handles)
% hObject    handle to encoding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dec =   handles.dec2d;
 guidata(hObject, handles);
global true;
global false;
true  =  1;
false =  0;
% % parameters
bitbudget = 1000000;
% X = imread(filename);

  
% sX = size(X);


% round all coefficients
% bits1=2;
%  dec2d = round(dec2d);   
%  dec2d=quant1(dec2d,bits1);

% reset start time
tic;

% perform SPIHT compression where encoded contains output and bits contains
% the amount of bits used.
[encoded bits] = cSPIHT(dec, 1, bitbudget);
enc_time = toc;

fprintf('Encoding time: %6.3f seconds\n', enc_time);

save enc_time enc_time;

save encoded encoded
save bits bits


save sX

warndlg('Encoding process completed');


% --- Executes on button press in decoding.
function decoding_Callback(hObject, eventdata, handles)
% hObject    handle to decoding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~exist ('encoded.mat')

    msgbox ('start encoding first');
else



load encoded 
load bits
load  sX


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % perform inverse
 tic;
[decoded level] = dSPIHT(encoded, bits);
dec_time = toc;
nobits =bits;
handles.decoded=decoded;
guidata(hObject, handles);
save nobits nobits;
save dec_time dec_time;
 % record cSPIHT time and output

fprintf('Decoding time: %6.3f seconds\n',dec_time);

warndlg('Decoding process completed');
end
% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes on button press in validate.
function validate_Callback(hObject, eventdata, handles)
% hObject    handle to validate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load PSNR;
load CR;
load MSE;
load enc_time;
load dec_time;
load InputFILESIZE;
load nobits;

PSNR=num2str(PSNR);
CR=num2str(CR);
set(handles.psnr,'String',PSNR);
set(handles.cr,'String',CR);
set(handles.enc_time,'String',enc_time);
set(handles.dec_time,'String',dec_time);
set(handles.orig_size,'String',InputFILESIZE);
set(handles.comp_size,'String',nobits);

fprintf('\nMSE: %7.2f \nPSNR: %9.7f dB', MSE, PSNR);
fprintf(' \nCR: %9.7f dB',  CR);
% --- Executes during object creation, after setting all properties.
function psnr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to psnr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function psnr_Callback(hObject, eventdata, handles)
% hObject    handle to psnr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of psnr as text
%        str2double(get(hObject,'String')) returns contents of psnr as a double


% --- Executes during object creation, after setting all properties.
function cr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function cr_Callback(hObject, eventdata, handles)
% hObject    handle to cr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cr as text
%        str2double(get(hObject,'String')) returns contents of cr as a double


% --- Executes during object creation, after setting all properties.
function enc_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enc_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function enc_time_Callback(hObject, eventdata, handles)
% hObject    handle to enc_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enc_time as text
%        str2double(get(hObject,'String')) returns contents of enc_time as a double


% --- Executes during object creation, after setting all properties.
function dec_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dec_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function dec_time_Callback(hObject, eventdata, handles)
% hObject    handle to dec_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dec_time as text
%        str2double(get(hObject,'String')) returns contents of dec_time as a double


% --- Executes on button press in wavelets.

function inverse_wave_Callback(hObject, eventdata, handles)
% hObject    handle to inverse_wave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents= get(handles.filter,'Value')
% load sX;
decoded=handles.decoded;
sX=size(decoded);
load nobits;
% filename=handles.filename;
% B=imread(filename);
B=handles.a;
B=imresize(B,[256,256]);
[r c p]=size(B);

if p==3
   
X=B(:,:,1);
image1=B(:,:,2);
image2=B(:,:,3);
else
    X=B;
end
X=double(X);
cA1 = decoded(1:(sX(1)/2), 1:(sX(1)/2));
cH1 = decoded(1:(sX(1)/2), (sX(1)/2 + 1):sX(1));
cV1 = decoded((sX(1)/2 + 1):sX(1), 1:(sX(1)/2));
cD1 = decoded((sX(1)/2 + 1):sX(1), (sX(1)/2 + 1):sX(1));

% cD1 =  zeros(1,4096);
% cD1 = reshape(cH1,[64 64]);

% reset start time
% tic;
if contents==1
% reconstruct image from wavelet coefficients
dec = idwt2(cA1,cH1,cV1,cD1,'haar');
dec=uint8(dec);
% record IDWT time and output
% idwttime = toc;
% fprintf('IDWT time:   %6.3f seconds\n', idwttime);

% output total times
% fprintf('\nTotal encoding time:   %6.3f seconds\n', dwttime + cspihttime);
% fprintf('\nTotal decoding time: %6.3f seconds\n', idwttime + dspihttime);


% calculate Mean Square Error and PSNR
% dec=double(dec);
% title('reconstructed image');
% figure,
% subplot(1,2,1);imshow(X,[]);
% title('original image');
% subplot(1,2,2);imshow(dec,[]);
% title('reconstructed image');
if p==3
Out(:,:,1)=dec;
Out(:,:,2)=image1;
Out(:,:,3)=image2;
else
    Out=dec;
end
axes(handles.three);
imshow(Out,[]);
n=size(B);
M=n(1);
N=n(2);
InputFILESIZE = M * N *16*3;
CR = InputFILESIZE / nobits ;
% Out=double(Out);
dec=double(dec);
MSE = sum(sum((X-dec).^2))/(M*N);
% MSE = sum(sum((X-dec).^2))/(size(X,1))/(size(X,2));
PSNR = 10*log10(255*255/MSE);
save CR CR;
save PSNR PSNR;
save MSE MSE;
save InputFILESIZE InputFILESIZE ;
save nobits nobits;
else
    

    dec = IDBMW(cA1,cH1,cV1,cD1);
dec=uint8(dec);
if p==3
Out(:,:,1)=dec;
Out(:,:,2)=image1;
Out(:,:,3)=image2;
else
    Out=dec;
end
axes(handles.three);
imshow(Out,[]);
n=size(B);
M=n(1);
N=n(2);
InputFILESIZE = M * N *16*3;
CR = InputFILESIZE / nobits ;
% Out=double(Out);
dec=double(dec);
X=double(X);
MSE = sum(sum((X-dec).^2))/(M*N);
% MSE = sum(sum((X-dec).^2))/(size(X,1))/(size(X,2));
PSNR = 10*log10(255*255/MSE);
save CR CR;
save PSNR PSNR;
save MSE MSE;
save InputFILESIZE InputFILESIZE ;
save nobits nobits;

end



% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete *.mat
PSNR = 0;
CR=0;
enc_time=0;
dec_time=0;
a=ones(256,256);
axes(handles.one);
imshow(a);
axes(handles.two);
imshow(a);
axes(handles.three);
imshow(a);

PSNR=num2str(PSNR);
CR=num2str(CR);
set(handles.psnr,'String',PSNR);
set(handles.cr,'String',CR);
set(handles.enc_time,'String',enc_time);
set(handles.dec_time,'String',dec_time);
warndlg('Files erased succesfully');



% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close all
clear all
close gui1
% --- Executes during object creation, after setting all properties.
function filter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in filter.
function filter_Callback(hObject, eventdata, handles)
% hObject    handle to filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns filter contents as cell array
%        contents{get(hObject,'Value')} returns selected item from filter


% --- Executes during object creation, after setting all properties.
function orig_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to orig_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function orig_size_Callback(hObject, eventdata, handles)
% hObject    handle to orig_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of orig_size as text
%        str2double(get(hObject,'String')) returns contents of orig_size as a double


% --- Executes during object creation, after setting all properties.
function comp_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to comp_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function comp_size_Callback(hObject, eventdata, handles)
% hObject    handle to comp_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of comp_size as text
%        str2double(get(hObject,'String')) returns contents of comp_size as a double


