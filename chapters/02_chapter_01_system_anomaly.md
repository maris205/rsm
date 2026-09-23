<!-- Faithful format conversion; scientific claims are unreviewed source text. Source: ../ori_paper/merged_document_v2_with_latex.docx; Pandoc blocks [144, 185). -->

# Chapter 1: System Anomaly (The Bug) --- Two New Dark Clouds Over the Physics Sky

**Abstract:**

High-precision observational data from the James Webb Space Telescope (JWST) has finally confirmed a $5\sigma$ discrepancy in the determination of the Hubble Constant between the early universe (CMB) and the late universe (SH0ES). This "Hubble Tension," alongside the Vacuum Energy Catastrophe---often termed the "worst theoretical prediction in history"---constitutes two new dark clouds hovering over the Standard Cosmological Model ($\Lambda$CDM). This chapter proposes, from the perspective of Simulation Engineering, that these anomalies are not measurement errors or the effects of unknown fields, but rather inevitable intrinsic drifts of a discrete computational system on an asymptotic time scale ($t \rightarrow \infty$). This is not merely a crisis of physics; it is a "Fault Alarm" for the universe as a computational system.

## 1.1 Introduction: Webb's "Death Sentence"

In 1900, Lord Kelvin delivered a famous lecture at the Royal Institution, declaring that the edifice of physics was complete, and that only "two dark clouds" floated in the sky: Ether drift and the Ultraviolet Catastrophe. These two clouds ultimately sparked the revolutions of Relativity and Quantum Mechanics.

One hundred and twenty-six years later, history is humorously repeating itself.

Before the James Webb Space Telescope (JWST) launched in 2021, a sense of optimism pervaded the physics community. The mainstream view held that the "Hubble Tension" plaguing cosmology for years---the inconsistency between two different methods of measuring the universe\'s expansion rate---was simply because our rulers weren\'t precise enough. Everyone hoped that Webb, with its unparalleled infrared resolution, would correct observational errors regarding Cepheid variables, thereby restoring the Standard Model to perfection.

However, data from 2024 to 2026 handed the Standard Model a "death sentence."

Webb did not smooth out the difference; it confirmed it. It clearly distinguished Cepheid variables from the crowded galactic background, eliminating photometric contamination. The results showed: Adam Riess's team was correct. The universe today is indeed expanding faster than the Planck satellite predicted.

This $5\sigma$ confidence discrepancy (implying a probability of error of only one in a million) formally declared that physics has entered the phase of **"System Anomaly."**

## 1.2 The First Cloud: Hubble Tension --- Server Clock Drift

Let us open this "trouble ticket" that is driving physicists mad:

**Early Universe Measurement (Server Startup Phase):** Based on the Planck satellite\'s observation of the Cosmic Microwave Background (CMB). Extrapolating to today using the $\Lambda$CDM model, the Hubble Constant $H_{0}$ should be $67.4 \pm 0.5$ km/s/Mpc.

**Late Universe Measurement (Server Running to Date):** Based on Webb's actual measurements of the current universe, the result is $73.0 \pm 1.0$ km/s/Mpc.

The gap is as high as 9%. To balance the books, theoretical physicists attempt to introduce new "fields" or "particles." But in the eyes of a simulation engineer, no new particles are needed. This is a classic phenomenon of **Sampling Clock Drift**.

Imagine a massive discrete event simulation:

During the system\'s "Cold Start," the load is extremely low, and the Global Clock is very precise.

After running for 13.8 billion years, the number of entities in the system explodes (entropy increase), and the computational load reaches its peak. Although the logical Tick remains $+ 1$, due to the accumulation of processing latency in computational nodes, the sampling step at the physical layer ($\Delta t$) has subtly elongated.

The Hubble Tension is not space expanding faster, but the time benchmark defining "speed" becoming slower. It is as if a video player was quietly adjusted from 1.0x speed to 1.1x speed; to the observer, everything on the screen appears to have accelerated.

## 1.3 The Second Cloud: Vacuum Energy Catastrophe --- Floating Point Truncation Error

If the Hubble Tension merely gives physicists a headache, the "Vacuum Energy Catastrophe" causes them to break down.

**Theoretical Prediction:** According to Quantum Field Theory, the vacuum is filled with quantum fluctuations. The zero-point energy density after superimposing all modes is astoundingly large.

**Actual Observation:** We require an extremely tiny Cosmological Constant ($\Lambda$) to drive expansion.

**The Error:** The theoretical value is $10^{120}$ times larger than the observed value.

Physicists are helpless against this, resorting to the "Anthropic Principle." But within the Riemann Simulation Architecture, this is inevitable **Precision Truncation**.

In computer graphics, this is called LOD (Level of Detail) filtering. Quantum Field Theory calculates the energy at the underlying Lattice Level, which is full of high-frequency noise; meanwhile, gravity runs at the macroscopic Render Level. To prevent numerical overflow, the simulation system performs a forced **Zero-point Normalization** when outputting data across layers.

That $10^{120}$ energy didn\'t disappear; it was simply discarded by the system as "invalid pointer overhead." The remaining $10^{- 120}$ residue (i.e., Dark Energy) is merely the **Rounding Error** left behind by the system\'s floating-point operations.

## 1.4 Engineer\'s Diagnosis: Performance Degradation from Excessive Uptime

When we examine these two dark clouds, we do not see the collapse of physical laws, but rather a **Kernel Panic** log issued by a distributed system that has been running for 13.8 billion years.

**\[System Log Snapshot\]**

Plaintext

\[LOG\_ERROR\] Global\_Clock jitter detected: 67.4 vs 73.0.\
\[LOG\_WARN\] Precision loss at Planck level: 10\^-120 residual noise found.\
\[STATUS\] Kernel: Riemann\_Kernel\_v1.0; Uptime: 4.35e17 seconds.

**Diagnostic Conclusion:**

**Hubble Constant Inconsistency:** Caused by a drop in sampling frequency due to increased total system load.

**Vacuum Energy Catastrophe:** Caused by numerical precision truncation when mapping from the microscopic lattice to the macroscopic render.

**Root Cause:** The **Uptime** (continuous running time) of the Universe System is too long, and the system is approaching its initially defined precision limits.

The drift of these parameters is not chaotic or disordered; it follows strict number-theoretic distribution laws. This indicates that behind the anomalies, a lower-level code logic is controlling this decay.

## 1.5 Conclusion

Physicists are attempting to apply the 100th patch (introducing new scalar fields) to this server. What we need to do, however, is open the console and check the underlying source code.

These two dark clouds are not there to obscure our vision; they are breaches leading to the truth. In the next chapter, we will pass through these Bugs to explore the very first line of code that generated this universe.

**Next --- Chapter 2: The Source Code.**
