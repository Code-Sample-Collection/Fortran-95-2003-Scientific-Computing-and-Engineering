# Fortran 95/2003 Scientific Computing and Engineering

> 《Fortran 95/2003科学计算与工程》

NOTE: 目前中文编码为 GBK！

## [配书光盘说明](readme.txt)
```
《Fortran 95/2003科学计算与工程》配书光盘说明

光盘内容：

	1.各章代码
	2.readme.txt

光盘运行环境：

    本光盘中范例程序的运行环境在书中附录A说明。


注意：本书配套光盘中的文件，仅用于学习和练习时使用，未经许可不能用于任何商业行为。

```


## 目录
```
封面 1
扉页 2
内容简介 3
版权页 3
前言 4
引子 6
目录 10

第1章　矩阵分解与线性方程组的直接方法 14
1.1 三角方程组 14
1.2 高斯消去法 20
1.3 选主元消去法 25
1.4 Crout分解 30
1.5 Doolittle 分解 33
1.6 LU 分解法计算线性方程组 37
1.7 追赶法计算三对角方程 41
1.8 对称正定阵的乔里斯基（Cholesky）分解 46
1.9 用Cholesky分解计算对称正定方程 49
1.10 行列式的计算 53
1.11 矩阵方程的计算 56
1.12 逆矩阵的计算 63
1.13 线性方程组解的迭代改进 68
    本章小结 74

第2章　解线性方程组的迭代方法 75
2.1 Jacobi迭代法 75
2.2 Gauss-Seidel 迭代法 79
2.3 逐次超松弛迭代法 83
2.4 Richardson 同步迭代法 88
2.5 广义Richardson 迭代法 92
2.6 Jacobi 超松弛迭代法 95
2.7 最速下降法 99
2.8 共轭梯度法 105
    本章小结 112

第3章　最小二乘与数据拟合 113
3.1 Cholesky分解法计算最小二乘 113
3.2 Householder 镜像变换之QR分解 119
3.3 修正的Gram-Schimdt 正交化方法的QR分解 126
3.4 QR分解法计算最小二乘问题 130
3.5 最小二乘曲线拟合 136
    本章小结 142

第4章　矩阵特征值及特征向量 143
4.1 幂法计算主特征值及其特征向量 143
4.2 幂法2范数单位化方法 147
4.3 Rayleigh 加速方法 152
4.4 修正的Rayleigh 加速方法 157
4.5 QR分解方法求全部特征值 162
    本章小结 166

第5章　非线性方程求根 167
5.1 Bolzano 二分法 168
5.2 Picard 迭代法 173
5.3 Aitken加速与Steffensen 迭代方法 178
5.4 Newton-Raphson 迭代法 184
5.5 重根时的迭代改进 189
5.6 割线法 195
5.7 多重迭代法 199
5.8 4 阶收敛多重迭代法 204
5.9 开普勒方程的计算 209
    本章小结 214

第6章　非线性方程组的数值方法 215
6.1 牛顿迭代法 215
6.2 简化牛顿法 221
6.3 拟牛顿之Broyden 方法 228
6.4 Broyden 第二公式计算非线性方程组 237
6.5 DFP方法 247
6.6 BFS方法 256
6.7 拓展收敛域之数值延拓法 266
6.8 拓展收敛域之参数微分法 277
    本章小结 287

第7章　插值法 288
7.1 拉格朗日插值 288
7.2 牛顿插值法 292
7.3 Hermite 插值 296
7.4 三次样条插值之固支条件 300
7.5 三次样条插值之自然边界条件 308
7.6 三次样条之周期边界条件 315
7.7 反插值 324
7.8 第一类标准B样条 328
7.9 第二类标准B样条 336
7.10 第三类标准B样条 343
    本章小结 351

第8章　数值微分 352
8.1 简单的中点公式 352
8.2 三点公式法 355
8.3 五点公式法 358
8.4 Richardson 外推方法 361
8.5 数值微分应用范例—雷达跟踪微分求速 364
    本章小结 368

第9章　数值积分 369
9.1 复合梯形求积法 369
9.2 复合Simpson积分 373
9.3 自动变步长Simpson方法 377
9.4 复合高阶Newton-Cotes 方法 382
9.5 Romberg 积分方法 386
9.6 Gauss-Legendre 积分 390
9.7 Gauss-Laguerre 方法计算反常积分 395
9.8 Gauss-Hermite 方法计算反常积分 399
9.9 复合高斯积分法 403
9.10 变步长高斯积分方法 407
9.11 重积分的数值方法 412
    本章小结 416

第10章　常见的特殊函数计算 417
10.1 Gamma 函数 417
10.2 不完全Gamma 函数及其互补函数 420
10.3 Beta函数及卡方分布函数 426
10.4 误差函数、余误差函数 及标准正态分布表的制作 431
10.5 第一类整数阶贝塞尔函数 440
10.6 第二类整数阶贝塞尔函数 448
    本章小结 457

第11章　常微分方程（组）的数值方法 458
11.1 经典龙格-库塔方法 458
11.2 Gill方法 465
11.3 Rung-Kutta方法计算微分方程组 468
11.4 Adams-Bashforth 3 步三阶方法 472
11.5 Adams-Bashforth 4 步四阶方法 478
11.6 三阶Adams 预测校正方法（PECE） 483
11.7 四阶Adams 预测校正方法（PECE） 489
    本章小结 494

第12章　应用范例 495
12.1 航天器轨道外推 495
12.2 卫星三位置矢量的Gibbs定初轨方法 502
12.3 空间导航基本原理 506
12.4 计算机辅助设计中的Bézier 样条曲线 515
12.6 人体生理周期预测 518
    本章小结 524

正文结束 524
附录A 集成开发环境介绍 525
附录B 程序调试方法 533
附录C 代码编辑器Ultra Edit 539
参考文献 541
```
