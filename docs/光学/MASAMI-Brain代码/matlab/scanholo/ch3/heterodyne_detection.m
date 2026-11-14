% ===== 书 3.2 节 · 光学外差 + 电子多路检测（式 3.2-1 ~ 3.2-4） =====
% 时域仿真：信息光 psi_p = A exp(j*phi) exp(j*w0 t) 与参考光 B exp[j(w0+Omega)t]
% 在光电探测器上拍频，电流 (3.2-1)：
%     i(t) = D [ A^2 + B^2 + 2AB cos(Omega t - phi) ]
% 电子学链路（图 3.6）：BPF@Omega 提取外差电流 -> 乘 cos(Omega t)+LPF 得同相分量 i_c (3.2-3a)
% -> 乘 sin(Omega t)+LPF 得正交分量 i_s (3.2-3b) -> 复数加法 i_c + j*i_s = A exp(j*phi) (3.2-4)
%
% 差异提示：
% 1) 理想带通/低通用 FFT 选带与时间平均实现；实际系统用模拟滤波器。
% 2) 本例无傅里叶约定问题（只做频域选带）；涉及书式变换的脚本见同目录其他文件。

clear
A = 0.7; phi = 1.3; B = 1.0;     % 信息光幅度/相位（真值）、参考光幅度（实数）
Omega = 2*pi*1.0e6;              % 外差（拍频）频率 1 MHz
fs = 64e6;                       % 采样率
T = 8e-6;                        % 总时长 8 us（8 个拍频周期）
t = 0:1/fs:T-1/fs;

% ---- 光电探测电流（式 3.2-1，取探测器面积 D=1） ----
i_t = A^2 + B^2 + 2*A*B*cos(Omega*t - phi);     % DC(基带) + AC(外差)

% ---- BPF @ Omega：提取外差电流 i_Omega ∝ 2AB cos(Omega t - phi) ----
I_f = fft(i_t);
f_axis = (0:length(t)-1)*fs/length(t);
band = f_axis > 0.5e6 & f_axis < 1.5e6;         % 理想带通（保留 Omega，滤掉直流基带）
iOmega = real(ifft(I_f.*band));

% ---- lock-in 解调（式 3.2-2 / 3.2-3）：乘 cos/sin 后低通（时间平均） ----
seg = round(0.5*length(t)):length(t);           % 丢弃滤波器建立段
ic = 2*mean(iOmega(seg).*cos(Omega*t(seg)));    % i_c = A cos(phi)   (3.2-3a)
is_ = 2*mean(iOmega(seg).*sin(Omega*t(seg)));   % i_s = A sin(phi)   (3.2-3b)
A_rec = hypot(ic, is_)/2;                       % 幅度（含常数因子 2，归一化后）
phi_rec = atan2(is_, ic);                       % 相位
fprintf('恢复幅度 A = %.6f（真值 %.2f），误差 %.1e\n', A_rec, A, abs(A_rec-A));
fprintf('恢复相位 phi = %.6f rad（真值 %.2f），误差 %.1e\n', phi_rec, phi, abs(phi_rec-phi));
fprintf('式(3.2-4)：i_c + j*i_s ∝ A exp(j*phi) —— 幅度与相位完整恢复\n');

% ---- 可视化（对应笔记图 sim-heterodyne.png） ----
figure(1)
subplot(2,2,1); plot(t*1e6, i_t, 'LineWidth', 0.8); grid on
title('i(t) = D[A^2+B^2+2AB cos(Omega t - phi)]  (Eq. 3.2-1)'); xlabel('t (us)')
subplot(2,2,2); plot(f_axis/1e6, abs(I_f), 'LineWidth', 0.8); xlim([0 4]); grid on
title('spectrum: DC(baseband) + heterodyne line at \Omega'); xlabel('f (MHz)')
subplot(2,2,3); plot(t*1e6, iOmega, 'LineWidth', 0.8); grid on
title('after BPF @ \Omega:  i_\Omega = 2AB cos(\Omega t - \phi)'); xlabel('t (us)')
subplot(2,2,4); plot(t*1e6, iOmega.*cos(Omega*t), 'LineWidth', 0.6); hold on
plot(t*1e6, iOmega.*sin(Omega*t), 'LineWidth', 0.6)
yline(ic/2, 'r--'); yline(-is_/2, 'b:'); grid on
title('lock-in detection (Eq. 3.2-2 / 3.2-3)'); xlabel('t (us)')
legend('i_\Omega \times cos\Omegat', 'i_\Omega \times sin\Omegat', ...
       sprintf('LPF -> i_c/2 = %.3f', ic/2), sprintf('LPF -> -i_s/2 = %.3f', -is_/2), ...
       'Location', 'best')
