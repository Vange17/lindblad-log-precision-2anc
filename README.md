# lindblad-log-precision-2anc
    The Matlab program is compatible with ​​MATLAB R2022b​​ and is designed to simulate a 5-qubit transverse-field Ising model with ​​periodic boundary conditions​​, where the first qubit is subject to ​​Amplitude Damping (AD) noise​​. The dataset corresponding to this program is: miti_Ising_OQS4_sp5n500gam1p5rho1O1_NMnoshots. This program is adapted from the paper "​​Lindbladian simulation with logarithmic precision scaling via two ancillas​​" by Wenjun Yu, Xiaogang Li, Qi Zhao, and Xiao Yuan, published in Physical Review Letters.

## Supplementary variance and standard errors

The [variance supplement](variance/README.md) provides independent trajectory
data, reproducible variance estimates and standard errors of the Monte Carlo
mean for the five-qubit benchmark. The original script and dataset are preserved.

## Thesis reference version

This fork records the numerical version used in the revised thesis. It retains
the upstream MATLAB source and archived means from revision
`db4e49f9ab20aaf907803dc788f71206efbdfd6a`, and adds the supplementary trajectory
data, variance/SEM analysis, regression checks, and
[reproducible thesis figures](figures/README.md).

Figure 4.4 labels the rotation-related curves as estimated T counts: their
archived values include a factor of 66 per Rz rotation. The figure documentation
specifies the synthesis convention and the pre-synthesis CNOT-count scope;
the numerical values and plot layout are unchanged.

Figure 4.5 keeps the original 5,000,000-trajectory mean estimates and adds SEM
estimates obtained from 10,000 supplementary trajectories. The new sampler
corrects a per-step phase initialization; its SEM estimates apply to the
historical means only if the estimator distributions match. This is not
established by the archived data. See the variance supplement for details.

For reproducible citation, use the full commit identifier recorded in the
thesis bibliography, rather than the moving `main` branch. The upstream
project remains [xiaogangli1169/lindblad-log-precision-2anc](https://github.com/xiaogangli1169/lindblad-log-precision-2anc);
this fork preserves its attribution and history.
