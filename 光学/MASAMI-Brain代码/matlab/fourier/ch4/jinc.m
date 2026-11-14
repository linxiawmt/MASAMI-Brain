function[out]=jinc(x);
% ===== 书 4.5 节：jinc 函数 M 文件（逐行忠实 + 中文注释） =====
% 书 4.5 节在 New M-file（命名为 "jinc"）里建立的函数（原书第 1-14 行）：
%   1  function[out]=jinc(x);
%   5  % J1(2*pi*x)/x
%   9  mask=(x~=0);
%   11 out=pi*ones(size(x));
%   13 out(mask)=besselj(1,2*pi*x(mask))./(x(mask));
%   14 end
% 要点（书 4.5 节原话与讲解）：
%   - 本函数计算 J1(2*pi*x)/x；用掩码（mask）办法避开 x=0 的除零：
%     第 9 行 mask 取 x 非零元素（~= 即 ≠），第 11 行 out 先全部初始化为
%     pi——它是 x=0 处的函数值（J1(u)≈u/2 -> J1(2*pi*x)/x -> pi），
%     第 13 行用逻辑索引 out(mask)、x(mask) 只对非零元素求值。
%   - 掩码写法看似绕远，但允许输入 x 是向量或数组（全书一贯的向量化风格）。
%   - besselj(1,...) 是第一类 1 阶贝塞尔函数；Python 对应 scipy.special.j1。
%   - 书特意提醒："not all jinc functions are the same"——本书的 jinc(x)=
%     J1(2*pi*x)/x，jinc(0)=pi，第一个零点在 x=0.61（J1 零点 1.22*pi）。
%     与公共模块 python/fourier/common/base.py 的 jinc(rho, R=1.0) 完全一致。
% MATLAB/Python 差异：logical 索引 out(mask) 对应 numpy 的布尔掩码 x[mask]；
%   MATLAB 的 pi 是常量，numpy 用 np.pi。
%
% ------------- 以下为书原代码（第 1-14 行） -------------
%
% jinc function
%
% J1(2*pi*x)/x
% divide by zero fix
%
% locate non-zero elements of x
mask=(x~=0);
% initialize output with pi (value for x=0)
out=pi*ones(size(x));
% compute output values for all other x
out(mask)=besselj(1,2*pi*x(mask))./(x(mask));
end
