%% driver file
close all;
if 1
S0_path = 'C:\Users\Peihong\Desktop\Data\FaceWarehouse_Data_0\Tester_1\Blendshape\shape_0.obj';
S1_path = 'C:\Users\Peihong\Desktop\Data\FaceWarehouse_Data_0\Tester_1\Blendshape\shape_22.obj';
T0_path = 'C:\Users\Peihong\Desktop\Data\FaceWarehouse_Data_0\Tester_12\Blendshape\shape_0.obj';
T1_path = 'C:\Users\Peihong\Desktop\Data\FaceWarehouse_Data_0\Tester_12\Blendshape\shape_22.obj';
%T0_path = 'C:\Users\Peihong\Desktop\Data\FaceWarehouse_Data_0\Tester_106\Blendshape\shape_0.obj';
%T1_path = 'C:\Users\Peihong\Desktop\Data\FaceWarehouse_Data_0\Tester_106\Blendshape\shape_22.obj';
else    
%S0_path = 'horse-poses/horse-01.obj';
%S1_path = 'horse-poses/horse-03.obj';
S0_path = 'cube.obj';
S1_path = 'cube1.obj';
T0_path = S0_path;
T1_path = S1_path;
end

S0 = read_wobj(S0_path); S0.faces = S0.objects.data.vertices;
S1 = read_wobj(S1_path); S1.faces = S1.objects.data.vertices;
T0 = read_wobj(T0_path); T0.faces = T0.objects.data.vertices;
T1 = read_wobj(T1_path); T1.faces = T1.objects.data.vertices;

doTransfer;