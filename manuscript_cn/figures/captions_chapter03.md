# 第3章图注与来源

## 图3.1：共同更新规则与独立物理读数

文件：`chapter03_update_and_readout_cn.{pdf,svg,png}`。

完整状态 $Z=(q_\ell,p_\ell,R,P_R)$ 包含格点场与共享内部时钟。图中展示纯包络：共同势能中的同一逆对数平方响应同时进入场力与时钟反馈。可分离 Hamilton 量的完整数值步采用 $K(h/2)\to D(h)\to K(h/2)$：每一次 kick 在该子步同一配置上同时求出所有耦合力，第一次使用漂移前的配置，第二次使用漂移后的配置。物理读数通过另行指定的读出与校准规则建立。下部区分计算依赖与待定义的物理因果关系；共享的均匀 $R$ 是当前全局约化，不是已经建立的局域驱动场。$h$ 是可加密验证的数值步长，不是已证实的宇宙基本时间节拍。

本图为新增的方法示意，不包含新实验、拟合或观测数据。依据 [现有耦合系统推导的第9节](../../experiments/log_clock_coupling_v01/reports/derivation_and_limits.md) 绘制，未声称数值积分精确保能、光锥已被证明或数据依赖等于物理瞬时作用。

复现命令：`python manuscript_cn/figures/make_chapter03_figures.py`。
提供矢量 PDF、SVG 与 240 dpi PNG；沿用既有中文字体和配色。
