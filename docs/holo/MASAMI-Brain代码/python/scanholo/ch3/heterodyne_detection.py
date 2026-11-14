# ===== 书 3.2 节 · 光学外差 + 电子多路检测（式 3.2-1 ~ 3.2-4） =====
# 时域仿真：信息光 psi_p = A exp(j*phi) exp(j*w0 t) 与参考光 B exp[j(w0+Omega)t]
# 在光电探测器上拍频，电流 (3.2-1)：
#     i(t) = D [ A^2 + B^2 + 2AB cos(Omega t - phi) ]
# 电子学链路（图 3.6）：BPF@Omega 提取外差电流 -> 上通道乘 cos(Omega t) + LPF 得同相分量 i_c (3.2-3a)
# -> 下通道乘 sin(Omega t) + LPF 得正交分量 i_s (3.2-3b) -> 复数加法 i_c + j*i_s = A exp(j*phi) (3.2-4)
#
# 验证结论（本笔记实跑）：A_rec = 0.700000（误差 2e-16）、phi_rec = 1.300000 rad（误差 2e-15）。
# 即外差电流完整保留信息光的幅度与相位；lock-in 解调可无损恢复。
#
# MATLAB / Python 差异提示：
# 1) MATLAB 的 FFT 前向核是 e^{-j...}，numpy 同为 e^{-j...}，此处只做频域选带，无约定问题。
# 2) 理想带通/低通直接用 FFT 选带与时间平均实现；实际系统用模拟滤波器，会引入纹波与延迟。

import numpy as np
import matplotlib.pyplot as plt

# ----- 1. 参数 -----
A, phi, B = 0.7, 1.3, 1.0        # 信息光幅度/相位（真值）、参考光幅度（实数）
Omega = 2*np.pi*1.0e6            # 外差（拍频）频率 1 MHz
fs = 64e6                        # 采样率
T = 8e-6                         # 总时长 8 us（8 个拍频周期）
t = np.arange(0, T, 1/fs)

# ----- 2. 光电探测电流（式 3.2-1，取探测器面积 D=1） -----
i_t = A**2 + B**2 + 2*A*B*np.cos(Omega*t - phi)      # DC(基带) + AC(外差)

# ----- 3. BPF @ Omega：提取外差电流 i_Omega ∝ 2AB cos(Omega t - phi) -----
I_f = np.fft.rfft(i_t)
f_axis = np.fft.rfftfreq(len(t), 1/fs)
band = (f_axis > 0.5e6) & (f_axis < 1.5e6)           # 理想带通（保留 Omega，滤掉直流基带）
iOmega = np.fft.irfft(np.where(band, I_f, 0), len(t))

# ----- 4. lock-in 解调（式 3.2-2 / 3.2-3）：乘 cos/sin 后低通（时间平均） -----
seg = slice(int(0.5*len(t)), None)                   # 丢弃滤波器建立段
ic = 2*np.mean(iOmega[seg]*np.cos(Omega*t[seg]))     # i_c = A cos(phi)   (3.2-3a)
is_ = 2*np.mean(iOmega[seg]*np.sin(Omega*t[seg]))    # i_s = A sin(phi)   (3.2-3b)
A_rec = np.hypot(ic, is_)/2                          # 幅度（含常数因子 2，归一化后）
phi_rec = np.arctan2(is_, ic)                        # 相位
print(f'恢复幅度 A = {A_rec:.6f}（真值 {A}），误差 {abs(A_rec-A):.1e}')
print(f'恢复相位 phi = {phi_rec:.6f} rad（真值 {phi}），误差 {abs(phi_rec-phi):.1e}')
print('式(3.2-4)：i_c + j*i_s ∝ A exp(j*phi) —— 幅度与相位完整恢复')

# ----- 5. 可视化（对应笔记图 sim-heterodyne.png） -----
fig, ax = plt.subplots(2, 2, figsize=(12, 7))
ax[0,0].plot(t*1e6, i_t, lw=0.8)
ax[0,0].set_title('i(t) = D[A$^2$+B$^2$+2AB cos($\\Omega$t-$\\phi$)]   (Eq. 3.2-1)')
ax[0,0].set_xlabel('t (us)'); ax[0,0].set_ylabel('i(t)')
ax[0,1].plot(f_axis/1e6, np.abs(I_f), lw=0.8); ax[0,1].set_xlim(0, 4)
ax[0,1].set_title('spectrum: DC(baseband) + heterodyne line at $\\Omega$')
ax[0,1].set_xlabel('f (MHz)')
ax[1,0].plot(t*1e6, iOmega, lw=0.8)
ax[1,0].set_title('after BPF @ $\\Omega$:  i$_\\Omega$ = 2AB cos($\\Omega$t-$\\phi$)')
ax[1,0].set_xlabel('t (us)')
ax[1,1].plot(t*1e6, iOmega*np.cos(Omega*t), lw=0.6, label='i$_\\Omega\\times$cos$\\Omega$t')
ax[1,1].plot(t*1e6, iOmega*np.sin(Omega*t), lw=0.6, alpha=0.6, label='i$_\\Omega\\times$syn$\\Omega$t')
ax[1,1].axhline(ic/2, color='r', ls='--', label=f'LPF -> i_c/2 = {ic/2:.3f}')
ax[1,1].axhline(-is_/2, color='b', ls=':', label=f'LPF -> -i_s/2 = {-is_/2:.3f}')
ax[1,1].set_title('lock-in detection (Eq. 3.2-2 / 3.2-3)')
ax[1,1].set_xlabel('t (us)'); ax[1,1].legend(fontsize=8)
plt.tight_layout(); plt.show()
