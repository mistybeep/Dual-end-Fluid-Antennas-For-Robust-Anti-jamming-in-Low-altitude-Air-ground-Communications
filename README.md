# Dual-end Fluid Antennas For Robust Anti-jamming in Low-altitude Air-ground Communications

This repository contains the official MATLAB simulation code for the paper:  
**"Dual-end Fluid Antennas For Robust Anti-jamming in Low-altitude Air-ground Communications"**, to be published in **IEEE Journal of Selected Topics in Signal Processing (JSTSP)**, 2026.

---

## 📖 Project Overview

In low-altitude Air-to-Ground (A2G) networks, communications are highly vulnerable to malicious jamming. This research proposes a **Dual-end Fluid Antenna (FA)** framework where antenna elements can freely move within a confined space at both the transmitter and the receiver. By optimizing the positions of these "fluid" elements, the system can autonomously navigate to spatial locations that maximize signal strength while nullifying interference.

### System Model
<p align="center">
  <img src="https://youke.xn--y7xa690gmna.cn/s1/2026/02/10/698ad09c76048.webp" width="50%" alt="System Model">
</p>

---

## ✨ Main Innovation: `find_pos.m`

The core contribution of this repository is the `find_pos.m` function. This script implements the primary optimization logic described in the paper.

---

## 🚀 Getting Started

### Prerequisites
* MATLAB R2023b or higher.

### Usage
1.  Clone this repository to your local machine.
2.  Open MATLAB and set the project folder as the current directory.
3.  Run the main simulation script:
    ```matlab
    run('Func/JamPow.m')
    ```

---

## 📜 Citation

If you use this code or our findings in your research, please cite our JSTSP 2026 paper:

```bibtex
@ARTICLE{DualEndFluidAntenna2026JSTSP,
  author={Yifan Guo, Junshan Luo, Fanggang Wang, Haiyang Ding, Shilian Wang, and Zhenhai Xu},
  journal={IEEE Journal of Selected Topics in Signal Processing}, 
  title={Dual-end Fluid Antennas For Robust Anti-jamming in Low-altitude Air-ground Communications}, 
  year={2026},
  volume={},
  number={},
  pages={1-16},
  doi={10.1109/JSTSP.2026.3666647}
}
