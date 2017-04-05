clc;
clear;

global markers % marker data array NFx(NM*3)
global markers_world
global numFrames % number of frames
global numMarkers % number of markers
global NF %total number of frames

NF = 2000;

%Add Paths
currentPath = path();
addpath((currentPath));

%% Part 2a

%Run mdrplot to convert data
data = mrdplot_convert('d00121');

%Parameters
duration = 10;
dt = 0.01;
numMarkers = 8; % number of markers
numFrames = duration/dt;

%Get Lander Position
lander_position = data(1:2000,16:18);
rotz_90 = [0 -1 0; 1 0 0; 0 0 1];

%Get Marker Data
for i=1:numFrames
    markers_world(i,:) = [data(i,1) data(i,57:80)]; 
end

%assignment marker position
for i=1:numFrames
    markers(i,:) = [data(i,1) data(i,81:104)];
end

%Use optimization to find initial COM and Velocity
options = optimset('MaxFunEvals',1000000);
for i = 1:6
    p0(i) = 0;
end
[answer,fval,exitflag]=fminunc(@criterion,p0,options);
initialCOM_world = answer(1:3);
initialCOM = initialCOM_world - lander_position(1,:);
initialVel = answer(4:6);

%COM(t)
com(1,:) = initialCOM;
com_full(1,:) = initialCOM;
com_world(1,:) = initialCOM_world;
com_world_full(1,:) = initialCOM_world;
for i=2:numFrames
   com(i,:) = initialCOM + initialVel*dt*(i-1);   
   com_world(i,:) = initialCOM_world + initialVel*dt*(i-1);
end
for i=2:NF
    com_full(i,:) = initialCOM + initialVel*dt*(i-1);
    com_world_full(i,:) = initialCOM_world + initialVel*dt*(i-1); %(lander_position(i,:) * rotz_90)
end

% %Check COM
% figure();
% hold on;
% plot(markers(:,1),com_world(:,1));
% plot(markers(:,1),com_world(:,2));
% plot(markers(:,1),com_world(:,3));
% hold off;
% title('COM of Alien Artifact'); 
% legend('x','y','z');

%Find orientation of alien artifact over time
q0 = [1 0 0 0]; %initial quaternion orientation
for i=1:numFrames
    %get orthonormal basis from markers
    x_axis = markers_world(i,14:16)-markers_world(i,2:4);
    x_axis = x_axis./norm(x_axis);
    y_axis = markers_world(i,8:10)-markers_world(i,2:4);
    y_axis = y_axis./norm(y_axis);
    z_axis = markers_world(i,5:7)-markers_world(i,2:4);
    z_axis = z_axis./norm(z_axis);
    
    %make rotation matrix
    rotm_world(:,:,i) = [x_axis' y_axis' z_axis'];
    
    %convert rotation matrix to assigment coordinates
    rotm(:,:,i) = rotm_world(:,:,i)* rotz_90;
    
    %convert to quaternion
    q(i,:) = rotm2quat(rotm(:,:,i));
    if i>1
        if norm(q(i,:)-q(i-1,:))>1
            q(i,:) = -q(i,:);
        end
    end
    q_world(i,:) = rotm2quat(rotm_world(:,:,i));
    if i>1
        if norm(q_world(i,:)-q_world(i-1,:))>1
            q_world(i,:) = -q_world(i,:);
        end
    end
end

%Write part 2a results to dat file
part2a = [com q];
save('problem_2_0.dat','part2a','-ascii');

% %Check quaternions
% figure();
% hold on;
% plot(markers(:,1),q_world(:,1));
% plot(markers(:,1),q_world(:,2));
% plot(markers(:,1),q_world(:,3));
% plot(markers(:,1),q_world(:,4));
% hold off;
% title('Orientation of Alien Artifact (Quaternions)');
% legend('q scalar','q x','q y','q z');

% %Check basis vectors
% figure();
% hold on;
% for i=1:numMarkers
%     scatter3(markers(1,(i-1)*3+1+1), markers(1,(i-1)*3+2+1), markers(1,(i-1)*3+3+1));
% end
% plot3([markers(1,2) markers(1,2)+x_axis(1)], [markers(1,3) markers(1,3)+x_axis(2)], [markers(1,4) markers(1,4)+x_axis(3)]);
% plot3([markers(1,2) markers(1,2)+y_axis(1)], [markers(1,3) markers(1,3)+y_axis(2)], [markers(1,4) markers(1,4)+y_axis(3)]);
% plot3([markers(1,2) markers(1,2)+z_axis(1)], [markers(1,3) markers(1,3)+z_axis(2)], [markers(1,4) markers(1,4)+z_axis(3)]);
% hold off;


%% Part 2B
clear wx wy wz; clc;

%Angular velocites from rotation matrices
for i=1:numFrames-1
    %finite difference
    drotm(:,:,1) = (rotm(:,:,i+1)-rotm(:,:,i))/dt;
    drotm_world(:,:,1) = (rotm_world(:,:,i+1)-rotm_world(:,:,i))/dt;
    
    %skew symmetric angular velocity tensor
    W(:,:,1) = transpose((rotm(:,:,i))) * drotm(:,:,1);
    W_world(:,:,1) = transpose((rotm_world(:,:,i))) * drotm_world(:,:,1);
    
    %angular velocities in assignment coordinates
    wx(i) = W(3,2,1);
    wy(i) = W(1,3,1);
    wz(i) = W(2,1,1);
    
    %angular velocities in world coordinates
    wx_world(i) = W_world(3,2,1);
    wy_world(i) = W_world(1,3,1);
    wz_world(i) = W_world(2,1,1);
end
wx = [wx wx(i)];
wy = [wy wy(i)];
wz = [wz wz(i)];
wx_world = [wx_world wx_world(i)];
wy_world = [wy_world wy_world(i)];
wz_world = [wz_world wz_world(i)];

%Write part 2b results to dat file
part2b = [wx' wy' wz'];
save('problem_2_1.dat','part2b','-ascii');

% %Check angular velocities
% figure();
% hold on;
% plot(markers(1:1000,1),wx_world(:));
% plot(markers(1:1000,1),wy_world(:));
% plot(markers(1:1000,1),wz_world(:));
% hold off;
% title('Angular Velocities of Alien Artifact');
% legend('omega x','omega y','omega z');

%% Part 2c

%Angular Accelerations
%finite differences
for i=1:numFrames-2
   ax(i) =  (wx(i+1)-wx(i))/dt;
   ay(i) =  (wy(i+1)-wy(i))/dt;
   az(i) =  (wz(i+1)-wz(i))/dt;
   
   ax_world(i) =  (wx_world(i+1)-wx_world(i))/dt;
   ay_world(i) =  (wy_world(i+1)-wy_world(i))/dt;
   az_world(i) =  (wz_world(i+1)-wz_world(i))/dt;
end
ax = smooth([ax ax(i) ax(i)],10);
ay = smooth([ay ay(i) ay(i)],10);
az = smooth([az az(i) az(i)],10);
ax_world = smooth([ax_world ax_world(i) ax_world(i)],10);
ay_world = smooth([ay_world ay_world(i) ay_world(i)],10);
az_world = smooth([az_world az_world(i) az_world(i)],10);

%Write part 2c results to text file
part2c = [ax ay az];
save('problem_2_2.dat','part2c','-ascii');

% %Check angular accelerations
% figure();
% hold on;
% plot(markers(1:1000,1),ax_world(:));
% plot(markers(1:1000,1),ay_world(:));
% plot(markers(1:1000,1),az_world(:));
% hold off;
% title('Angular Accelerations of Alien Artifact');
% legend('alpha x','alpha y','alpha z');

%% Part 2d

% Factorize t=I*w_dot +w X Iw
% syms I11 I12 I13 I22 I23 I33;
% syms a_x a_y a_z;
% syms w_x w_y w_z;
% w_dot = [a_x; a_y; a_z];
% w = [w_x; w_y; w_z];
% I = [I11 I12 I13; I12 I22 I23; I13 I23 I33];
% equation = simplify(I*w_dot + cross(w,I*w));

%Solve for Inertia Tensor
A_tot = [];
for i=1:50:1000
    w_dot = part2c(i,:)';
    w = part2b(i,:)';
    A = [w_dot(1) w_dot(2)-w(1)*w(3) w_dot(3)+w(1)*w(3) -w(2)*w(3) w(2)^2-w(3)^2 w(2)*w(3);...
         w(1)*w(3) w_dot(1)+w(2)*w(3) w(3)^2-w(1)^2 w_dot(2) w_dot(3)-w(1)*w(2) -w(1)*w(3);...
         -w(1)*w(2) w(1)^2 w_dot(1)-w(2)*w(3) w(1)*w(2) w_dot(2)+w(1)*w(3) w_dot(3)];
    A_tot = [A_tot;A];
end
[U,S,V] = svd(A_tot);

%Build Inertia Matrix
I11 = V(1,6);
I12 = V(2,6);
I13 = V(3,6);
I21 = I12;
I22 = V(4,6);
I23 = V(5,6);
I31 = I13;
I32 = I23;
I33 = V(6,6);

I = [I11 I12 I13; I21 I22 I23; I31 I32 I33];

%Principal Axis
I_principle = I;
I_principle(2,1) = 0;
I_principle(1,3) = 0;
I_principle(3,2) = 0;
I_principle(1,2) = I_principle(2,1);
I_principle(3,1) = I_principle(1,3);
I_principle(2,3) = I_principle(3,2);
I_principle(1,1) = I(2,2);
I_principle(2,2) = I(1,1);

rot_I_to_I_prince = I_principle/I;

%% Part 2e
clear w; clear S;
%Simulation
%initialize stuff
w_x_world(1) = wx_world(end);
w_y_world(1) = wy_world(end);
w_z_world(1) = wz_world(end);
w_world(1,:) = [w_x_world(1) w_y_world(1) w_z_world(1)];
alpha_x_world(1) = ax_world(end);
alpha_y_world(1) = ay_world(end);
alpha_z_world(1) = az_world(end);
R_world(:,:,1) = rotm_world(:,:,end);

w_x(1) = wx(end);
w_y(1) = wy(end);
w_z(1) = wz(end);
w(1,:) = [w_x(1) w_y(1) w_z(1)];
alpha_x(1) = ax(end);
alpha_y(1) = ay(end);
alpha_z(1) = az(end);
R(:,:,1) = rotm(:,:,end);

%future COM
for i=2:numFrames
   com_future(i,:) = com(end,:) + initialVel*dt*(i-1);
end

%simulate forward for 1000 time steps (10 seconds)
for i=1:1000
    %generate skew symmetric matrix S
    S(:,:,i) = [0 -w_z(i) w_y(i);...
                w_z(i) 0 -w_x(i);...
               -w_y(i) w_x(i) 0];
           
    S_world(:,:,i) = [0 -w_z_world(i) w_y_world(i);...
                     w_z_world(i) 0 -w_x_world(i);...
                     -w_y_world(i) w_x_world(i) 0];

    %find derivative of rotation matrix
    Rdot(:,:,i) = R(:,:,i)*S(:,:,i);
    Rdot_world(:,:,i) = R_world(:,:,i)*S_world(:,:,i);

    %euler integrate rotations forward
    R(:,:,i+1) = R(:,:,i)+Rdot(:,:,i)*dt;
    w_x(i+1) = w_x(i)+alpha_x(i)*dt;
    w_y(i+1) = w_y(i)+alpha_y(i)*dt;
    w_z(i+1) = w_z(i)+alpha_z(i)*dt;
    w(i+1,:) = [w_x(i+1) w_y(i+1) w_z(i+1)]; 
    
    R_world(:,:,i+1) = R_world(:,:,i)+Rdot_world(:,:,i)*dt;
    w_x_world(i+1) = w_x_world(i)+alpha_x_world(i)*dt;
    w_y_world(i+1) = w_y_world(i)+alpha_y_world(i)*dt;
    w_z_world(i+1) = w_z_world(i)+alpha_z_world(i)*dt;
    w_world(i+1,:) = [w_x_world(i+1) w_y_world(i+1) w_z_world(i+1)]; 
    
    %solve torque equations for new alpha
    alpha(i,:) = [alpha_x(i) alpha_y(i) alpha_z(i)];
    alpha(i+1,:) = -(cross(w(i+1,:)', I_principle*w(i+1,:)'))'/I_principle;
    alpha_x(i+1) = alpha(i+1,1);
    alpha_y(i+1) = alpha(i+1,2);
    alpha_z(i+1) = alpha(i+1,3);
    
    alpha_world(i,:) = [alpha_x_world(i) alpha_y_world(i) alpha_z_world(i)];
    alpha_world(i+1,:) = -(cross(w_world(i+1,:)', I_principle*w_world(i+1,:)'))'/I_principle;
    alpha_x_world(i+1) = alpha_world(i+1,1);
    alpha_y_world(i+1) = alpha_world(i+1,2);
    alpha_z_world(i+1) = alpha_world(i+1,3);
end
t = linspace(0, 20, 2000);

for i=1:1000
   quat(i,:) = rotm2quat(R(:,:,i));
end
for i=1:1000
   quat_world(i,:) = rotm2quat(R_world(:,:,i));
end
for i=1:1000
    if i>1
        if norm(quat(i,:)-quat(i-1,:))>1
            quat(i,:) = -quat(i,:);
        end
    end
end
for i=1:1000
    if i>1
        if norm(quat_world(i,:)-quat_world(i-1,:))>1
            quat_world(i,:) = -quat_world(i,:);
        end
    end
end

q_full =  [q; quat];
wx_full = [wx'; w_x(1:1000)'];
wy_full = [wy'; w_y(1:1000)'];
wz_full = [wz'; w_z(1:1000)'];
ax_full = [ax; alpha_x(1:1000)'];
ay_full = [ay; alpha_y(1:1000)'];
az_full = [az; alpha_z(1:1000)'];

q_full_world =  [q_world; quat_world];
wx_full_world = [wx_world'; w_x_world(1:1000)'];
wy_full_world = [wy_world'; w_y_world(1:1000)'];
wz_full_world = [wz_world'; w_z_world(1:1000)'];
ax_full_world = [ax_world; alpha_x_world(1:1000)'];
ay_full_world = [ay_world; alpha_y_world(1:1000)'];
az_full_world = [az_world; alpha_z_world(1:1000)'];

for i=1:length(q_full_world)
    if i>1
        if norm(q_full_world(i,:)-q_full_world(i-1,:))>1
            q_full_world(i,:) = -q_full_world(i,:);
        end
    end
end
for i=1:length(q_full)
    if i>1
        if norm(q_full(i,:)-q_full(i-1,:))>1
            q_full(i,:) = -q_full(i,:);
        end
    end
end


%Check COM
figure();
hold on;
plot(t,com_world_full(:,1));
plot(t,com_world_full(:,2));
plot(t,com_world_full(:,3));
hold off;
title('COM');
legend('x','y','z');

%Check quaternions
figure();
hold on;
plot(t,q_full_world(:,1));
plot(t,q_full_world(:,2));
plot(t,q_full_world(:,3));
plot(t,q_full_world(:,4));
hold off;
title('Orientation');
legend('q scalar','q x','q y','q z');

%Check angular velocities
figure();
hold on;
plot(t,wx_full_world(:));
plot(t,wy_full_world(:));
plot(t,wz_full_world(:));
hold off;
title('Angular Velocities');
legend('omega x','omega y','omega z');


%Check angular accelerations
figure();
hold on;
plot(t,ax_full_world(:));
plot(t,ay_full_world(:));
plot(t,az_full_world(:));
hold off;
title('Angular Accelerations');
legend('alpha x','alpha y','alpha z');


%Write part 2e results to text file
part2e = [com_future quat];
save('problem_2_3.dat','part2e','-ascii');


%% Part 3

%Get Rotation Matrices
R(:,:,end)=[];
R_world(:,:,end)=[];
for i=1:NF
%     R_full(:,:,i) = quat2rotm(q_full_world(i,:));
    if i<=1000
        R_full(:,:,i) = rotm_world(:,:,i);
    else
        R_full(:,:,i) = R_world(:,:,i-1000);
    end
end

%Lander Trajectory
com_trajectory = com_world_full;
for i = 1:NF
    if i < 500
        offset = [-0.5 6.2 0];
    else
        offset = [-0.5 3.2 0];
    end
    com_trajectory(i,:) = com_trajectory(i,:) + offset*transpose(R_full(:,:,i));
 end
    
%Save data for simulation control
fid = fopen('q_scalar.txt','wt');
fprintf(fid,'q_scalar\n');
for i =1:NF
    fprintf(fid, '%f\n',q_full(i,1));
end
fclose(fid);

fid = fopen('q_x.txt','wt');
fprintf(fid,'q_x\n');
for i =1:NF
    fprintf(fid, '%f\n',q_full(i,2));
end
fclose(fid);

fid = fopen('q_y.txt','wt');
fprintf(fid,'q_y\n');
for i =1:NF
    fprintf(fid, '%f\n',q_full(i,3));
end
fclose(fid);

fid = fopen('q_z.txt','wt');
fprintf(fid,'q_z\n');
for i =1:NF
    fprintf(fid, '%f\n',q_full(i,4));
end
fclose(fid);


fid = fopen('COM_X.txt','wt');
fprintf(fid,'COM_X\n');
for i =1:NF
    fprintf(fid, '%f\n',com_trajectory(i,1));
end
fclose(fid);

fid = fopen('COM_Y.txt','wt');
fprintf(fid,'COM_Y\n');
for i =1:NF
    fprintf(fid, '%f\n',com_trajectory(i,2));
end
fclose(fid);

fid = fopen('COM_Z.txt','wt');
fprintf(fid,'COM_Z\n');
for i =1:NF
    fprintf(fid, '%f\n',com_trajectory(i,3));
end
fclose(fid);

fid = fopen('W_X.txt','wt');
fprintf(fid,'W_X\n');
for i =1:NF
    fprintf(fid, '%f\n',wx_full(i));
end
fclose(fid);

fid = fopen('W_Y.txt','wt');
fprintf(fid,'W_Y\n');
for i =1: NF
    fprintf(fid, '%f\n',wy_full(i));
end
fclose(fid);

fid = fopen('W_Z.txt','wt');
fprintf(fid,'W_Z\n');
for i =1: NF
    fprintf(fid, '%f\n',wz_full(i));
end
fclose(fid);
