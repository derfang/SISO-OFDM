# Project 1: Adaptive OFDM System with Equalization Analysis

## 1. Project Overview
This project implements a complete end-to-end **Orthogonal Frequency Division Multiplexing (OFDM)** communication link using MATLAB. The simulation models a realistic wireless environment characterized by **Rayleigh Fading** and multipath interference.

The core objective is to analyze and compare the performance of two fundamental channel equalization algorithms:
1.  **Zero-Forcing (ZF)**
2.  **Minimum Mean Square Error (MMSE)**

Additionally, the system features an **Adaptive Equalizer** that dynamically switches strategies based on the estimated Signal-to-Noise Ratio (SNR), optimizing the trade-off between computational complexity and Bit Error Rate (BER) performance.

## 2. Key Features
* **End-to-End Simulation:** Implements full Tx/Rx chains including QAM mapping, IFFT/FFT processing, and Cyclic Prefix (CP) insertion/removal.
* **High-Performance Computing:** Utilizes MATLAB's `parfor` (Parallel Computing Toolbox) to accelerate Monte Carlo simulations, processing over **5,000 OFDM packets** for statistical accuracy.
* **Adaptive Logic:** Automates the decision-making process, switching algorithms at a defined threshold (15 dB) to balance performance and complexity.
* **Visual Analytics:** Generates professional-grade visualization of Phase Rotation and Constellation recovery.

## 3. System Architecture

### A. Transmitter
* **Modulation:** 16-QAM (Quadrature Amplitude Modulation).
* **OFDM Parameters:** 64 Subcarriers, 16-sample Cyclic Prefix.
* **Operation:** Bits $\rightarrow$ Symbol Map $\rightarrow$ IFFT $\rightarrow$ CP Insertion.

### B. Channel Model
* **Physics:** 4-tap multipath channel with Rayleigh distributed fading coefficients.
* **Noise:** Complex Additive White Gaussian Noise (AWGN).
* **Equation:**
    $$y = h * x + n$$

### C. Receiver & Equalization
The receiver recovers the transmitted signal $\hat{X}$ using Channel State Information (CSI) $H$ and the received signal $Y$.

**1. Zero Forcing (ZF)**
Inverts the channel frequency response. While computationally simple, it suffers from "Noise Amplification" when the channel gain $|H|$ is small (deep fade).
$$\hat{X}_{ZF} = \frac{Y[k]}{H[k]}$$

**2. Minimum Mean Square Error (MMSE)**
Minimizes the total error by accounting for the noise variance $\sigma^2_{noise}$. This prevents division by zero during deep fades.

$$\hat{X}_{MMSE} = \frac{H[k]^* \cdot Y[k]}{|H[k]|^2 + \sigma^2_{noise}}$$

**3. Adaptive Logic**
The system monitors SNR and switches modes:
* **Low SNR (< 15dB):** Uses MMSE to prioritize signal integrity.
* **High SNR (> 15dB):** Uses ZF to prioritize lower computational complexity.

## 4. Simulation Results

### Figure 1: BER Waterfall Curve
![BER Waterfall Curve](<BER Performance.jpg>)
![BER Waterfall Curve zoomed](<BER Performance zoom.jpg>)

* **Observation:** The MMSE algorithm (Blue) consistently outperforms Zero Forcing (Red), typically providing a **3-4 dB coding gain**.
* **Analysis:** The Zero Forcing curve flattens out at high SNR, indicating an "Error Floor" caused by unrecoverable spectral nulls.
* **Adaptive Performance:** The Black dashed line tracks the optimal performance path, verifying the switching logic.

### Figure 2: Constellation Analysis
![Constellation](<Constellation Analysis.jpg>)

* **Raw Signal (Left):** Visualizes the "Donut" effect caused by phase rotation in the complex channel.
* **Zero Forcing (Center):** Shows the impact of noise amplification—points are scattered ("exploded") far from their ideal grid locations.
* **MMSE (Right):** Demonstrates superior noise suppression, resulting in tight, clean symbol clusters.

## 5. How to Run
1.  Ensure you have MATLAB installed with the **Parallel Computing Toolbox** (optional but recommended for speed).
2.  Open `Project1_Final.m`.
3.  Run the script.
4.  The console will display progress for each SNR step:
    ```text
    --- Starting Final Simulation ---
    1. Calculating BER Waterfall (Parallel)...
    Completed 0 dB
    Completed 2 dB
    ...
    ```
5.  The script will output two figures: **BER Performance** and **Constellation Analysis**.

## 6. Future Work
* **MIMO Implementation:** Upgrade the system to 2x2 MIMO using Alamouti Space-Time Block Coding (STBC) for diversity gain.
* **Channel Estimation:** Replace "Perfect CSI" with Least Squares (LS) or MMSE pilot-based estimation.