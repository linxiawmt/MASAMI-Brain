% ===== 书 1.1 例 1.1 · Table 1.2 Plot_sinc.m =====
%Plot_sinc.m  Plotting of sinc(x) function
% 绘制 sinc(x)=sin(pi*x)/(pi*x)，零点在 x=±1,±2,±3,...
% 注：x=0 处是 0/0，MATLAB 会得到 NaN，plot 只是在峰值处漏画一个点；
%     如需平滑峰值，可把上一行改为 sinc=sin(pi*x+eps)./(pi*x+eps) 或对 x==0 补 1。
x=-5.5:0.01:5.5;
sinc=sin(pi*x)./(pi*x);
plot(x,sinc)
axis([-5.5 5.5 -0.3 1.1])
grid on
xlabel('x')
ylabel('sinc (x)')
