# Simulating-the-Double-Slit-Experiment-with-Python

## Overview
This project simulates the **double-slit experiment** in the far-field (Fraunhofer diffraction) using Python.  
The aperture is modeled as two rectangular slits, and the far-field intensity distribution is computed using a **Fast Fourier Transform (FFT)**.  
The program includes an **interactive interface** with sliders to adjust physical parameters such as wavelength, slit width, separation, and screen distance.

This project demonstrates:
- Understanding of the **physics of diffraction and interference**.
- Use of **numerical methods** (FFT) for solving physics problems.
- Skills in **Python programming and data visualization**.

---

## Background Theory

### Single Slit Diffraction
Light passing through a single slit of width $$\(a\)$$ produces a diffraction pattern given by:

$$
I(\theta) \propto \left(\frac{\sin \beta}{\beta}\right)^2,
\quad
\beta = \frac{\pi a}{\lambda} \sin \theta
$$

where:
- $$\(\lambda\)$$ = wavelength of light,
- $$\(\theta\)$$ = angle on the screen.

---

### Double Slit Interference
Two slits of width $$\(a\)$$, separated by distance $$\(d\)$$, give an intensity pattern:

$$
I(\theta) \propto
\cos^2\!\left(\frac{\pi d}{\lambda} \sin \theta\right)
\times
\left(\frac{\sin \beta}{\beta}\right)^2
$$

- The **cosine term** produces fine interference fringes.
- The **sinc term** (single slit diffraction) acts as an envelope.

---

### Fraunhofer Approximation
In the far-field, the intensity is proportional to the squared modulus of the Fourier transform of the aperture function:

$$
I(x,y) \propto \left| \mathcal{F}\{A(x',y')\} \right|^2
$$

- $$\(A(x',y')\)$$: aperture transmission function (two slits).  
- Fourier transform implemented numerically with **FFT**.

---

## Methodology

### 1. Aperture Construction
- Define a 2D grid representing the aperture plane.
- Set pixels to `1` where light passes (two slits), `0` elsewhere.

### 2. FFT for Diffraction
- Compute 2D FFT of the aperture.
- Shift zero-frequency component to center (`fftshift`).
- Square modulus to obtain intensity distribution.

### 3. Mapping to Physical Units
- Spatial frequencies $$\(f_x, f_y\)$$ from FFT map to screen coordinates via:

$$
x = \lambda D f_x, \quad y = \lambda D f_y
$$

where $$\(D\)$$ = screen distance.

### 4. Visualization
- `imshow`: 2D intensity pattern.
- Lineout at $$\(y = 0\)$$: interference fringes.
- Sliders allow real-time adjustment of physical parameters.

---

## Results

### Simulation prediction for high wavelength compared to low wavelength

<img width="1755" height="821" alt="Screenshot 2025-09-25 023125" src="https://github.com/user-attachments/assets/73824c8b-d0a8-4bdd-85ac-d1402b38b9bb" />
<img width="1645" height="810" alt="Screenshot 2025-09-25 023206" src="https://github.com/user-attachments/assets/03094d09-11ac-428a-bcca-57ef68c3b20f" />

### Simulation prediction for a smaller slit separation distance compared to a larger separation distance
<img width="1671" height="805" alt="Screenshot 2025-09-25 023229" src="https://github.com/user-attachments/assets/0e4197e2-2d3a-43fb-bf53-9273c8103aed" />
<img width="1586" height="814" alt="Screenshot 2025-09-25 023245" src="https://github.com/user-attachments/assets/dc0a6cb5-c12a-48a7-bd7f-68dcd83a0e7e" />

### Simulation prediction for a smaller slit width compared to a larger slit width
<img width="1642" height="819" alt="Screenshot 2025-09-25 023300" src="https://github.com/user-attachments/assets/f7035ac5-9ca1-46ec-a350-4b0be41a9661" />
<img width="1716" height="815" alt="Screenshot 2025-09-25 023340" src="https://github.com/user-attachments/assets/54d2ef77-b774-4b17-963a-a0eed95fc492" />

### Simulation prediction for a smaller slit height compared to a larger slit height
<img width="1596" height="814" alt="Screenshot 2025-09-25 023402" src="https://github.com/user-attachments/assets/f86efffc-e7a3-4b73-8dc9-828269bb65a2" />
<img width="1679" height="806" alt="Screenshot 2025-09-25 023416" src="https://github.com/user-attachments/assets/98ebb6bc-54f9-403b-87e3-a845313fd95f" />


## Discussion

- **Observations:**
  - Longer wavelength → fringes spread out.  
  - Increasing slit separation → fringes become closer.  
  - Increasing slit width → diffraction envelope narrows.  
  - Increasing slit height → diffraction squeezes in the y plane. 

- **Limitations:**
  - Simulation assumes Fraunhofer (far-field) approximation.  
  - Aperture edges are ideal (no roughness).  
  - Resolution limited by grid size $$\(N\)$$.  

- **Extensions:**
  - Fresnel regime (near-field diffraction).  
  - Arbitrary apertures (gratings, circular holes).  
  - 3D volumetric visualizations.  


---

## Code

All code is implemented in **Python** using:
- `numpy` for numerical calculations.
- `matplotlib` for visualization and interactivity.

Key functions:
```python
make_aperture(N, L, slit_width, slit_height, separation)
fraunhofer_pattern(aperture)
map_screen_coords(N, dx, wavelength, D)
