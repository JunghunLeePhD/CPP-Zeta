# **CPP-Zeta**

**CPP-Zeta** is a C++ project designed to numerically compute and visualize the zeros of the **Riemann Zeta Function** on the critical line.

This repository provides an efficient C++ implementation of advanced algorithms—including the **Euler-Maclaurin**summation and the **Riemann-Siegel** formula—to evaluate the Zeta function and the **Hardy $Z$-function**. It generates image frames visualizing the behavior of the function as t increases, which can be compiled into animations.


## **🎥 Demos**

![Hardy $Z$-function and $\zeta$-function](/assets/CPP-Zeta.gif)

### Hardy $Z$-function with the Euler-Maclaurin method
![Hardy $Z$-function by EM](/assets/hardyEM.gif)

### Hardy $Z$-function with the Riemann-Siegel method
![Hardy $Z$-function by RS](/assets/hardyRS.gif)

### Riemann $\zeta$-function with the Euler-Maclaurin method
![Riemann $\zeta$-function by RS](/assets/zetaEM.gif)

## **📖 Mathematical Background**


The Riemann Zeta function is a function of a complex variable $s$, 
defined in the half-plane containing complex numbers with real part greater than 1 by the Dirichlet series:

$$
    \zeta (s) 
    =
    \sum_{n = 1}^{\infty} n^{-s}
$$

To find non-trivial zeros on the critical line $Re(s)=1/2$, 
we analyze the **Hardy $Z$-function**, $Z(t)$. 
The zeros of $Z(t)$ correspond exactly to the zeros of the Riemann Zeta function on the critical line.

The definition used in this project is:

$$
    Z (t) = e^{i \theta(t)} \zeta \left(\frac{1}{2} + it \right)
$$

Where the Riemann-Siegel theta function $\theta(t)$ is defined as:

$$
    \theta(t) = \mathrm{arg}\left(\Gamma \left(\frac{1}{4} + \frac{i t}{2} \right) \right) - \frac{t}{2} \ln \pi
$$

### **Computational Methods Implemented**


To efficiently compute $\zeta(s)$ and $Z(t)$ for large $t > 0$, 
this library implements:

1. **Euler-Maclaurin Summation:** Used for high-precision evaluation at lower ranges of $t$.

2. **Riemann-Siegel Formula:** An asymptotic expansion allowing for much faster evaluation of $Z(t)$ at very high heights $t > 0$ on the critical line.


```
CPP-Zeta/
├── src/
│   └── main.cpp           # Main application logic
├── lib/
│   ├── zeta_euler.cpp     # Euler-Maclaurin implementation
│   ├── zeta_siegel.cpp    # Riemann-Siegel implementation
│   ├── hardy_z.cpp        # Hardy Z-function logic
│   ├── plotter.cpp        # Canvas/Graphing utilities
│   └── utils.cpp          # Complex number helpers & math utils
├── .devcontainer/         # VS Code DevContainer configuration
│   ├── devcontainer.json
│   └── Dockerfile
├── output/                # Generated frames and videos (created at runtime)
├── Makefile               # Build instructions
└── README.md              # Project documentation
```

## **🚀 Getting Started**

You can set up this project locally or use the provided VS Code DevContainer for an instant, pre-configured environment.

### **Option A: VS Code DevContainer (Recommended)**

This project is configured with a Development Container. If you use Visual Studio Code and Docker, you can open this repository in a container that comes pre-installed with **GCC**, **Make**, and **FFmpeg**.

1. Open the project in VS Code.

2. Click **"Reopen in Container"** when the popup appears (or run the command from the palette).

3. The environment will automatically install all necessary tools defined in `.devcontainer/Dockerfile`.

### **Option B: Local Installation**

If you prefer to run directly on your host machine, ensure you have the following prerequisites installed:

- **C++ Compiler:** GCC or Clang supporting C++17.

- **Make:** For building the project.

- **FFmpeg:** Required to convert the generated image frames into video files.


### **Installation**

1. **Clone the repository:**

    ```bash
    git clone https://github.com/JunghunLeePhD/CPP-Zeta.git
    cd CPP-Zeta
    ```

## **💻 Usage**

The usage consists of two steps: running the mathematical simulation to generate plot frames, and then compiling those frames into video files.

### **1. Run the Simulation**

Execute the following command to compile the project and run the simulation. This will generate numerical data and `.ppm` image frames inside the `output/` directory.

```bash
make run
```

*Note: Calculating high values of $t > 0$ using Riemann-Siegel may take time depending on the range specified in `*main.cpp*`.*

### **2. Generate Videos**

Once the simulation completes, use `ffmpeg` to create `.mp4` videos from the generated frames.

**Hardy Z-Function with the Euler-Maclaurin method:**

```bash
ffmpeg -framerate 300 -i output/frames_hardyEM/frame_%04d.ppm -c:v libx264 -pix_fmt yuv420p output/hardyEM.mp4
```

**Hardy Z-Function with the Euler-Maclaurin method:**

```bash
ffmpeg -framerate 300 -i output/frames_hardyRS/frame_%04d.ppm -c:v libx264 -pix_fmt yuv420p output/hardyRS.mp4
```

**Riemann Zeta-Function with the Euler-Maclaurin method:**

```bash
ffmpeg -framerate 300 -i output/frames_zeta/frame_%04d.ppm -c:v libx264 -pix_fmt yuv420p output/zeta.mp4
```

## **🧹 Cleanup**

To remove the generated frames and executable to save space:

```bash
make clean
```

## **🤝 Contributing**

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the Project

2. Create your Feature Branch (`git checkout -b feature/NewAlgorithm`)

3. Commit your Changes (`git commit -m 'Add Odlyzko-Schonhagen method'`)

4. Push to the Branch (`git push origin feature/NewAlgorithm`)

5. Open a Pull Request


## **📄 License**

Distributed under the MIT License. See `LICENSE` for more information.

Created by [*Junghun Lee, PhD*](https://www.google.com/search?q=https://github.com/JunghunLeePhD)
