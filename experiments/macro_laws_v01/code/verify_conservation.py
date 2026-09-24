#!/usr/bin/env python3
"""Finite, self-contained checks of Chapter 4 energy/momentum identities.

Run from any working directory. The protocol fixes inputs and thresholds before
the first execution. This code does not integrate trajectories or fit data.
"""
from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SEED = 20260925
BALANCE_TOL = 5e-12
NORMALIZATION_TOL = 5e-13
WAVE_TOL = 5e-13
COMPLEX_STEP = 1e-25


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(residual, terms, threshold):
    absolute = float(np.max(np.abs(np.asarray(residual))))
    scale = max(1.0, *(float(np.max(np.abs(np.asarray(t)))) for t in terms))
    return {"absolute_residual": absolute, "scale": scale,
            "scaled_residual": absolute / scale, "threshold": threshold,
            "passed": bool(np.isfinite(absolute) and absolute / scale <= threshold)}


def symbolic_checks():
    checks = []

    def record(name, expression):
        residual = sp.simplify(sp.expand(expression))
        checks.append({"name": name, "residual": str(residual),
                       "passed": bool(residual == 0)})

    a, v, m2, dm2, g = sp.symbols("a v m2 dm2 g", nonzero=True)
    ql, q, qr, ul, u, ur = sp.symbols("ql q qr ul u ur")
    e = (u**2 / 2 + m2*q**2 / 2 + g*q**4 / 4
         + v**2*((qr-q)**2 + (q-ql)**2)/(4*a**2))
    acc = v**2*(qr-2*q+ql)/a**2 - m2*q - g*q**3
    edot = (sp.diff(e,u)*acc + sp.diff(e,ql)*ul + sp.diff(e,q)*u
            + sp.diff(e,qr)*ur + sp.diff(e,m2)*dm2)
    jl = -v**2*(ul+u)*(q-ql)/(2*a)
    jr = -v**2*(u+ur)*(qr-q)/(2*a)
    record("finite_a_local_energy", edot + (jr-jl)/a - dm2*q**2/2)

    # A separate five-site periodic Hamiltonian, with a three-site block.
    nq = 5
    qs = sp.symbols(f"q0:{nq}")
    us = sp.symbols(f"u0:{nq}")
    site_e = [us[i]**2/2 + m2*qs[i]**2/2 + g*qs[i]**4/4
              + v**2*((qs[(i+1)%nq]-qs[i])**2
                       +(qs[i]-qs[(i-1)%nq])**2)/(4*a**2)
              for i in range(nq)]
    hfield = sum(a*us[i]**2/2 + a*m2*qs[i]**2/2 + a*g*qs[i]**4/4
                 + v**2*(qs[(i+1)%nq]-qs[i])**2/(2*a)
                 for i in range(nq))
    record("finite_a_density_normalization", a*sum(site_e)-hfield)
    accelerations = [-sp.diff(hfield,qs[i])/a for i in range(nq)]
    eb = a*sum(site_e[1:4])
    deb = (sum(sp.diff(eb,qs[i])*us[i]
               + sp.diff(eb,us[i])*accelerations[i] for i in range(nq))
           + sp.diff(eb,m2)*dm2)
    flux = [-v**2*(us[(i+1)%nq]+us[i])*(qs[(i+1)%nq]-qs[i])/(2*a)
            for i in range(nq)]
    record("finite_a_block_flux", deb-flux[0]+flux[3]
           - a*dm2*sum(qs[i]**2 for i in range(1,4))/2)

    x, t = sp.symbols("x t")
    f = sp.Function("phi")(x,t)
    mass = sp.Function("M2")(t)
    potential = mass*f**2/2 + g*f**4/4
    ft, fx = sp.diff(f,t), sp.diff(f,x)
    acc_cont = v**2*sp.diff(f,x,2)-mass*f-g*f**3
    energy = ft**2/2 + v**2*fx**2/2 + potential
    current = -v**2*ft*fx
    record("continuum_energy_with_source",
           (sp.diff(energy,t)+sp.diff(current,x)-sp.diff(mass,t)*f**2/2)
           .subs(sp.diff(f,t,2),acc_cont))
    momentum = -ft*fx
    stress = ft**2/2 + v**2*fx**2/2 - potential
    record("continuum_momentum_uniform_clock",
           (sp.diff(momentum,t)+sp.diff(stress,x))
           .subs(sp.diff(f,t,2),acc_cont))
    inertia, clock_velocity, q2, mass_prime = sp.symbols("I U Q2 Mprime")
    clock_acc = -mass_prime*q2/(2*inertia)
    record("integrated_field_plus_clock_energy",
           mass_prime*clock_velocity*q2/2
           + inertia*clock_velocity*clock_acc)

    # Candidate extension only: G(R) replaces the constant v^2.
    r = sp.symbols("R")
    gradnorm, pr = sp.symbols("Qgrad PR")
    gr, mr = sp.Function("G")(r), sp.Function("M2R")(r)
    extended_h = pr**2/(2*inertia)+gr*gradnorm/2+mr*q2/2
    record("candidate_G_clock_force",
           -sp.diff(extended_h,r)+sp.diff(gr,r)*gradnorm/2
           +sp.diff(mr,r)*q2/2)
    gt = sp.Function("Gtime")(t)
    extended_acc = gt*sp.diff(f,x,2)-mass*f-g*f**3
    extended_energy = ft**2/2+gt*fx**2/2+potential
    extended_flux = -gt*ft*fx
    record("candidate_G_continuum_energy",
           (sp.diff(extended_energy,t)+sp.diff(extended_flux,x)
            -sp.diff(gt,t)*fx**2/2-sp.diff(mass,t)*f**2/2)
           .subs(sp.diff(f,t,2),extended_acc))
    extended_stress = ft**2/2+gt*fx**2/2-potential
    record("candidate_G_continuum_momentum",
           (sp.diff(momentum,t)+sp.diff(extended_stress,x))
           .subs(sp.diff(f,t,2),extended_acc))
    return checks


def mass_squared(r):
    return 1.0 + 2.0*(np.log(r/np.exp(-2.0))**-2 - 0.25)


def mass_prime(r):
    return -4.0/(r*np.log(r/np.exp(-2.0))**3)


def hamiltonian(y, n, a, inertia=7.5, v=0.9, g=0.2):
    q, p, r, pr = y[:n], y[n:2*n], y[-2], y[-1]
    difference = np.roll(q,-1)-q
    return (pr*pr/(2*inertia) + np.sum(p*p)/(2*a)
            + v*v*np.sum(difference*difference)/(2*a)
            + a*np.sum(mass_squared(r)*q*q/2+g*q**4/4))


def complex_gradient(y, n, a):
    gradient = np.empty(y.size)
    for j in range(y.size):
        perturbed = y.astype(complex)
        perturbed[j] += 1j*COMPLEX_STEP
        gradient[j] = np.imag(hamiltonian(perturbed,n,a))/COMPLEX_STEP
    return gradient


def assembled_gradient(y, n, a, inertia=7.5, v=0.9, g=0.2):
    """Separate explicit Hamilton gradient assembled bond by bond."""
    q, p, r, pr = y[:n], y[n:2*n], y[-2], y[-1]
    grad = np.zeros_like(y)
    grad[:n] = a*(mass_squared(r)*q+g*q**3)
    for i in range(n):
        j = (i+1)%n
        bond_derivative = v*v*(q[j]-q[i])/a
        grad[i] -= bond_derivative
        grad[j] += bond_derivative
    grad[n:2*n] = p/a
    grad[-2] = a*mass_prime(r)*np.sum(q*q)/2
    grad[-1] = pr/inertia
    return grad


def lattice_checks():
    rng = np.random.Generator(np.random.PCG64(SEED))
    cases = []
    inertia, v, g, r, udot_clock = 7.5, 0.9, 0.2, 1.4, 0.13
    for n in (7,16,31):
        for a in (0.03125,0.2,0.7):
            for realization in range(4):
                q, u = 0.4*rng.normal(size=n), 0.3*rng.normal(size=n)
                p, pr = a*u, inertia*udot_clock
                y = np.r_[q,p,r,pr]
                grad = complex_gradient(y,n,a)
                qdot, pdot = grad[n:2*n], -grad[:n]
                rdot, prdot = grad[-1], -grad[-2]
                plus, minus = np.roll(q,-1)-q, q-np.roll(q,1)
                density = (u*u/2+mass_squared(r)*q*q/2+g*q**4/4
                           +v*v*(plus*plus+minus*minus)/(4*a*a))
                energy_dot = (u*pdot/a+(mass_squared(r)*q+g*q**3)*qdot
                              +mass_prime(r)*rdot*q*q/2
                              +v*v*(plus*(np.roll(qdot,-1)-qdot)
                                     +minus*(qdot-np.roll(qdot,1)))/(2*a*a))
                flux = -v*v*(np.roll(qdot,-1)+qdot)*plus/(2*a)
                divergence = (flux-np.roll(flux,1))/a
                source = mass_prime(r)*rdot*q*q/2
                field_energy = float(hamiltonian(y,n,a)-pr*pr/(2*inertia))
                field_dot, clock_dot = a*np.sum(energy_dot), pr*prdot/inertia
                entry = {
                    "N": n, "a": a, "realization": realization,
                    "state_sha256": hashlib.sha256(y.astype('<f8').tobytes()).hexdigest(),
                    "rhs_gradient": check(grad-assembled_gradient(y,n,a),
                                          [grad,assembled_gradient(y,n,a)],BALANCE_TOL),
                    "density_normalization": check(a*np.sum(density)-field_energy,
                                                    [a*np.sum(density),field_energy],
                                                    NORMALIZATION_TOL),
                    "local_energy": check(energy_dot+divergence-source,
                                          [energy_dot,divergence,source],BALANCE_TOL),
                    "total_energy": check(field_dot+clock_dot,
                                          [field_dot,clock_dot],BALANCE_TOL),
                    "blocks": [],
                }
                for start,length in ((0,1),(n-2,n//2),(1,n)):
                    indices = (start+np.arange(length))%n
                    change = a*np.sum(energy_dot[indices])
                    left, right = flux[(start-1)%n], flux[(start+length-1)%n]
                    production = a*np.sum(source[indices])
                    result = check(change-left+right-production,
                                   [change,left,right,production],BALANCE_TOL)
                    entry["blocks"].append({"start":start,"length":length,**result})
                cases.append(entry)
    return cases


def plane_wave_moments(amplitude, k, v, phase_count, directions):
    phase = 2*np.pi*(np.arange(phase_count)+0.5)/phase_count
    sine = np.sin(phase)
    wavevectors = k*np.asarray(directions,dtype=float)
    # Product average over direction and independent uniformly sampled phase.
    velocity = np.broadcast_to(amplitude*v*k*sine,(len(directions),phase_count))
    gradient = -amplitude*wavevectors[:,None,:]*sine[None,:,None]
    gradient_squared = np.sum(gradient*gradient,axis=-1)
    energy = velocity*velocity/2+v*v*gradient_squared/2
    stress = (v*v*np.einsum('dti,dtj->dtij',gradient,gradient)
              +np.eye(3)[None,None,:,:]
              *(velocity*velocity/2-v*v*gradient_squared/2)[:,:,None,None])
    momentum = -velocity[:,:,None]*gradient
    return float(np.mean(energy)), np.mean(stress,axis=(0,1)), np.mean(momentum,axis=(0,1))


def wave_checks():
    directions = [[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]]
    cases = []
    for phase_count in (16,32,64):
        for amplitude in (0.4,1.2):
            for k in (0.3,2.1):
                for v in (0.7,1.7):
                    rho,stress,momentum = plane_wave_moments(amplitude,k,v,phase_count,directions)
                    expected_energy = amplitude*amplitude*v*v*k*k/2
                    cases.append({
                        "phase_count":phase_count,"amplitude":amplitude,"k":k,"v":v,
                        "mean_energy":rho,"mean_stress":stress.tolist(),
                        "mean_momentum":momentum.tolist(),
                        "energy_average":check(rho-expected_energy,[rho,expected_energy],WAVE_TOL),
                        "isotropic_stress":check(stress-np.eye(3)*rho/3,[stress,rho],WAVE_TOL),
                        "zero_mean_momentum":check(momentum,[rho/v],WAVE_TOL),
                    })
    rho,stress,momentum = plane_wave_moments(1.2,2.1,1.7,64,[[1,0,0]])
    anisotropic = {
        "directions":[[1,0,0]],"phase_count":64,"amplitude":1.2,"k":2.1,"v":1.7,
        "mean_energy":rho,"mean_stress":stress.tolist(),"mean_momentum":momentum.tolist(),
        "directional_stress":check(stress-np.diag([rho,0,0]),[stress,rho],WAVE_TOL),
        "max_deviation_from_rho_over_3_relative_to_rho":float(np.max(np.abs(stress-np.eye(3)*rho/3))/rho),
    }
    return cases,anisotropic


def collect_results(node, path=""):
    records = []
    if isinstance(node,dict):
        if "passed" in node:
            records.append({"path":path,**node})
        else:
            for key,value in node.items():
                records.extend(collect_results(value,f"{path}/{key}"))
    elif isinstance(node,list):
        for i,value in enumerate(node):
            records.extend(collect_results(value,f"{path}/{i}"))
    return records


def write_report(data):
    def worst(cases, key):
        return max(case[key]["scaled_residual"] for case in cases)
    lattice, wave = data["lattice_cases"], data["wave_cases"]
    block_max = max(block["scaled_residual"] for case in lattice for block in case["blocks"])
    failure_lines = ("无。全部预设门限均通过。" if not data["failures"] else
                     "\n".join(f"- `{failure['path']}`：{failure}" for failure in data["failures"]))
    symbolic_rows = "\n".join(f"| `{row['name']}` | `{row['residual']}` | {'通过' if row['passed'] else '失败'} |"
                              for row in data["symbolic_checks"])
    text = f"""# 第4章宏观收支规律的独立有限核验

执行时间：{data['metadata']['run_utc']}。本次执行共{data['summary']['total_checks']}项预定检查，{data['summary']['passed_checks']}项通过，{data['summary']['failed_checks']}项失败。检查不是新观测拟合或宇宙模拟，也没有计算新的时间演化轨迹。

## 输入和复算

在仓库根目录运行 `python experiments/macro_laws_v01/code/verify_conservation.py`。脚本只使用其内置输入及NumPy、SymPy。完整分例记录见 `../results/conservation_checks.json`，执行前固定的约定见 `../protocol.md`。

- 协议SHA256：`{data['metadata']['protocol_sha256']}`。
- 源码SHA256：`{data['metadata']['script_sha256']}`。
- Python {data['metadata']['python']}；NumPy {data['metadata']['numpy']}；SymPy {data['metadata']['sympy']}。
- 随机数：PCG64，固定种子{SEED}。每个随机状态另保存SHA256。

## 符号结果

| 恒等式 | 化简残差 | 状态 |
| --- | --- | --- |
{symbolic_rows}

有限格距能量密度采用对称分配键能，实际格点能量为密度乘以a。5点周期链中独立选取3点区域，内部键通量相消；连续场的时间质量项产生能量源，而均匀时钟不产生空间动量源。总场能与自由时钟能的变化精确抵消。

推广G(R)的三项符号检查仅说明：若在共同势能中另引入该梯度系数，则必须同时加入其时钟反作用。它不声称已运行的CMB代理实现了这个共同Hamilton。

## 非均匀格点数值结果

N取7、16、31，格距取0.03125、0.2、0.7，每组4个随机状态，共36例；每例检查3个区域，共108个区域，含跨周期边界区域。独立RHS通过完整Hamilton的复步长梯度获得，再与逐键装配的解析梯度比较。

| 数值检查 | 最大尺度化残差 | 预定门限 |
| --- | --- | --- |
| 两种Hamilton梯度 | {worst(lattice,'rhs_gradient'):.6e} | {BALANCE_TOL:.1e} |
| 密度乘a的能量归一化 | {worst(lattice,'density_normalization'):.6e} | {NORMALIZATION_TOL:.1e} |
| 每个格点的局部能量平衡 | {worst(lattice,'local_energy'):.6e} | {BALANCE_TOL:.1e} |
| 全场加时钟能量平衡 | {worst(lattice,'total_energy'):.6e} | {BALANCE_TOL:.1e} |
| 区域边界通量收支 | {block_max:.6e} | {BALANCE_TOL:.1e} |

尺度为相应等式各项绝对值的最大值与1取大值；JSON同时记录绝对残差和尺度。这里检验的是非均匀状态上ODE右端的代数收支，不代表KDK按有限步长推进时原Hamilton逐步不变。

## 人工各向同性波动集合

24例使用三组相位网格和两组幅度、波数及速度。六轴等权方向给定各向同性二阶矩；每个方向与独立均匀相位做乘积平均。

| 检查 | 最大尺度化残差 | 预定门限 |
| --- | --- | --- |
| 平均能量A²v²k²/2 | {worst(wave,'energy_average'):.6e} | {WAVE_TOL:.1e} |
| 平均应力ρδij/3 | {worst(wave,'isotropic_stress'):.6e} | {WAVE_TOL:.1e} |
| 平均动量为零 | {worst(wave,'zero_mean_momentum'):.6e} | {WAVE_TOL:.1e} |

仅沿正x方向的对照得到Πxx=ρ、Πyy=Πzz=0，其与各向同性ρ/3张量的最大差异除以ρ为{data['anisotropic_control']['max_deviation_from_rho_over_3_relative_to_rho']:.12g}。该对照通过的是定向波自身的应力关系，不能据它声称各向同性。三维ρ/3关系依赖额外的三维推广、线性无质量近似及人为指定的各向同性统计集合；本轮没有证明初始场会自行各向同性化或热化。

## 全部失败记录

{failure_lines}

## 结论边界

这些检查验证公式、积分权重及有限数值实现彼此一致，不建立真实粒子、光子、暗物质或宇宙膨胀模型。当前共享R仍是全局自由度：场局部能量方程含时钟源项，只有积分后与时钟能相消；它不自动给出相对论性的局部总应力能量张量。连续场动量恒等式也不能直接替代一般非线性有限格点的连续平移动量定理。对候选G(R)只做代数检查，未改变现有CMB计算或完成物理单位标定。
"""
    (ROOT/"reports/conservation_checks_cn.md").write_text(text,encoding="utf-8")


def main():
    data = {
        "metadata":{
            "run_utc":datetime.now(timezone.utc).isoformat(),
            "protocol_sha256":digest(ROOT/"protocol.md"),
            "script_sha256":digest(Path(__file__).resolve()),
            "python":platform.python_version(),"numpy":np.__version__,"sympy":sp.__version__,
            "random_generator":"PCG64","seed":SEED,"complex_step":COMPLEX_STEP,
            "scope":"finite algebraic and deterministic numerical checks; no trajectory integration or data fit",
        },
        "symbolic_checks":symbolic_checks(),
        "lattice_cases":lattice_checks(),
    }
    data["wave_cases"],data["anisotropic_control"] = wave_checks()
    records = collect_results(data)
    data["failures"] = [record for record in records if not record["passed"]]
    data["summary"] = {"total_checks":len(records),
                       "passed_checks":sum(record["passed"] for record in records),
                       "failed_checks":len(data["failures"]),
                       "symbolic_checks":len(data["symbolic_checks"]),
                       "lattice_states":len(data["lattice_cases"]),
                       "lattice_blocks":sum(len(case["blocks"]) for case in data["lattice_cases"]),
                       "isotropic_wave_cases":len(data["wave_cases"]),
                       "anisotropic_controls":1}
    (ROOT/"results/conservation_checks.json").write_text(
        json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    write_report(data)
    print(json.dumps(data["summary"],ensure_ascii=False))
    raise SystemExit(1 if data["failures"] else 0)


if __name__ == "__main__":
    main()
