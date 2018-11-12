基于180820性能提升版本，采用最新的解析解+Bisection来计算居民决策路径，允许负资产（能够刻画外债带来的过度消费）

ExtraScripts提供测试脚本test_Analytical.jl，证明了采用解析解可以保证同一套数据下最优路径上不同起点决策结果完全一致

另外：
1. 更新了NumAlgo模块，新增Bisection()函数，高精度二分法求根
2. 新增Policy_Analytical, Policy_Analytical_Retired模块，PolicyFunctions（DP，值函数求解）has been deprecated
3. 新增New_DepartableUtilityFunction系列文档(PART 1,2)给出解析解的推导、算法和代码中的符号约定



by Tianhao Zhao
2018-9-10