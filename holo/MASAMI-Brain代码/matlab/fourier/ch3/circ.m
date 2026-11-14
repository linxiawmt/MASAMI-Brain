function out = circ(r)
% ===== 书习题 3.3：圆函数 M 文件（仿书 3.1 节 rect 的编码风格） =====
% 习题 3.3(a) 要求自建 circle function M 文件。书没给标准答案，这里按
% rect.m 的同一风格编写：输入径向坐标数组 r，输出 0/1。
% 调用方式：circ(sqrt(X.^2+Y.^2)/R) -> 半径 R 的圆盘。
% 半径约定警告（承接第 1 章勘误 1）：表 1.2 的 circ(r) 是"半径 1/2"圆盘，
% 表 1.3 的变换对却按"半径 = 分母"（Goodman 约定）计算。本函数用调用式
% circ(r/R) 把半径明确为 R，频谱就是可靠的锚点 R*J1(2*pi*R*rho)/rho
% （jinc 形），峰值 = 圆盘面积 pi*R^2（中心纵坐标定理）。
% 边界 r==1 取 1（与书 rect 的"单点跳变"风格一致）；Python 参考实现
% python/fourier/common/base.py 的 circle() 用严格小于——两者只差恰好
% 落在圆周上的少数样本，对频谱的影响可忽略（ch3_exercises 有数值验证）。
out = (r <= 1);
end
