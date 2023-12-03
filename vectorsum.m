function [ ll,dsi, osi, vs_angle,idx ,nidx, o1dx, o2dx, os_angle, global_osi] = vectorsum( slist )
num_dir = length(slist);
theta0=0:(360/num_dir):(360-(360/num_dir)); theta=theta0*2*pi/360;
vs=zeros(length(theta),2);
resp = slist;
slist=slist/sum(slist);
for i = 1:8
    vs(i,1)=cos(theta(i)).*slist(i);
    vs(i,2)=sin(theta(i)).*slist(i);
end
vs=sum(vs);  
vs_angle = atan2(vs(2),vs(1));
vs_angle = vs_angle*(180/pi);
if vs_angle < 0;
    vs_angle = vs_angle + 360;
end

ll=sqrt(sum(vs.^2)); 
if vs(2)>0
    degrees=acos(vs(1)/ll);
elseif vs(2)<=0;
    degrees=pi+acos(-vs(1)/ll);
end
degrees=degrees*360/(2*pi);
% asign theta to 
[~ ,idx]=min(abs(degrees-[theta0 360]));
if idx==9
    idx=1;
end 
degrees=theta0(idx);
nidx=mod(idx+4,8);

if nidx==0
    nidx=8;
end 

o1dx = idx + 2;
if o1dx == 9
    o1dx = 1;
elseif o1dx == 10
    o1dx = 2;
end

if o1dx <= 4
    o2dx = o1dx +4;
elseif o1dx > 4
    o2dx = o1dx - 4;
end
    
dsi=(resp(idx)-resp(nidx))/(resp(idx)+resp(nidx));


% osi = ((slist(idx)+slist(nidx))-(slist(o1dx)+slist(o2dx)))/(slist(idx)+slist(nidx)+ slist(o1dx)-slist(o2dx));
% 
% o_vs1=zeros(length(theta),2);
% o_vs2=zeros(length(theta),2);
% slist1 = slist;
% slist2 = slist;
% slist1(6:8,:) = 0;
% slist2(2:4,:) = 0;
% for i = 1:length(theta)
%     o_vs1(i,1)=cos(theta(i)).*slist1(i);
%     o_vs1(i,2)=sin(theta(i)).*slist1(i);
%     o_vs2(i,1)=cos(theta(i)).*slist2(i);
%     o_vs2(i,2)=sin(theta(i)).*slist2(i);
% end
% 
% o_vs1=sum(o_vs1); 
% o_vs2=sum(o_vs2);
% os_angle1 = atan2(o_vs1(2),o_vs1(1));
% os_angle1 = os_angle1*(180/pi);
% 
% os_angle2 = atan2(o_vs2(2),o_vs2(1));
% os_angle2 = os_angle2*(180/pi) + 180;
% 
% os_angle = (os_angle1 + os_angle2)/2;

for q = 1:num_dir
    complex_numerator(1,q) = resp(q)*exp(2i*theta(q));
end

complex_phase = sum(complex_numerator)/sum(resp);

global_osi = abs(complex_phase);

os_angle = 0.5*(180/pi*angle(complex_phase));

if os_angle < 0
    os_angle = os_angle + 180;
end

diff = abs(theta0-os_angle);
op1_index = find(diff == min(diff));
o_pref1 = theta0(op1_index);
o_pref2 = o_pref1 + 180;
if o_pref2 >= 360
    o_pref2 = o_pref2 - 360;
end
op2_index = find(theta0 == o_pref2);
o_ortho1 = o_pref1 + 90;
oo1_index = find(theta0 == o_ortho1);
o_ortho2 = o_ortho1 + 180;
if o_ortho2 >= 360
    o_ortho2 = o_ortho2 - 360;
end
oo2_index = find(theta0 == o_ortho2);

r_pref = (resp(op1_index) + resp(op2_index))/2;
r_ortho = (resp(oo1_index) + resp(oo2_index))/2;
osi=(r_pref - r_ortho)/(r_pref + r_ortho);

end

