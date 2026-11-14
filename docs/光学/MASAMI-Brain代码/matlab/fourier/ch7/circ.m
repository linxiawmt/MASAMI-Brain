function[out]=circ(r)
% ===== 附录 C.1:circ 圆盘函数(逐行忠实 + 中文注释) =====
% 书附录 C.1 列出的函数(原书第 1-8 行)。
% 约定:宗量 r 已经"除好半径"——circ(r)=1 当 |r|<=1。
% 因此 P(x,y)=circ(sqrt(x^2+y^2)/wXP)(式 7.11)画的是半径 wXP 的圆瞳
% (Goodman 半径约定);注意表 1.2 的 circ 印刷定义对应"宗量 = r/(2*半径)",
% 两套约定差 2 倍半径(ch1 勘误)。Python 侧等价写法是公共模块的
% base.circle(X, Y, R)(直接传半径),见 python/fourier/common/base.py。
% ------------- 以下为书原代码(附录 C.1 第 1-8 行) -------------
%
% circle function
%
% evaluates circ(r)
% note: returns odd number of samples for diameter
%
out=abs(r)<=1;
end
