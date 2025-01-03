[![Build Status](https://github.com/sapphimars/physkit/actions/workflows/test.yml/badge.svg)](https://github.com/sapphimars/physkit/actions/workflows/test.yml)
[![PyPI Version](https://img.shields.io/pypi/v/physkit)](https://pypi.org/project/physkit/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

<div align="center">
  <figure>
    <img src="https://github.com/sapphimars/physkit/blob/main/assets/physkit-logo-light.png?raw=true" alt="Light mode logo" width=240 />
  </figure>
</div>

###### Icon created by [astroastra](https://github.com/astroastrastudio)

# Physkit - A Python toolkit for constants, unit conversions, and equations

**Physkit** is a Python library for performing scientific computations, unit conversions, and working with physical constants and equations.  
It provides tools for astrophysical and physical calculations, supporting multiple unit systems and CLI functionalities.

---

## **License**

This project is licensed under the terms of the **GNU General Public License v3.0**.  
See the [LICENSE](https://github.com/sapphimars/physkit/blob/main/LICENSE) file for details.  

Copyright (C) 2024 [sapphimars](https://github.com/sapphimars)  

---

### **Third-Party Licenses**

This project makes use of the following third-party libraries:

- **[Matplotlib](https://matplotlib.org/)** - Licensed under the Matplotlib License Agreement. See [licenses/LICENSE-MATPLOTLIB.txt](https://github.com/sapphimars/physkit/blob/main/licenses/LICENSE-MATPLOTLIB.txt) or [directly on their repository](https://github.com/matplotlib/matplotlib/tree/main/LICENSE).
- **[Pint](https://github.com/hgrecco/pint)** - Licensed under the BSD 3-Clause License. See [licenses/LICENSE-PINT.txt](https://github.com/sapphimars/physkit/blob/main/licenses/LICENSE-PINT.txt) or [directly on their repository](https://github.com/hgrecco/pint/blob/master/LICENSE).

Note: Third-party libraries included in this project are licensed under their respective terms. See the `licenses/` directory for full details.

---

## **Installation**

### **Requirements**

- Python **3.12** or higher  
- Dependencies:
  - **Matplotlib** >= 3.5  
  - **Pint** >= 0.20  

---

### **Install via pip**
```bash
pip install physkit
```

---

### **For Development**
```bash
pdm install
```

To install test dependencies:
```bash
pdm install -G test
```

---

## **Usage Examples**

### **Access Constants**
```python
import physkit as pk
from physkit.constants import constants as csts

print(csts.G)  # Gravitational constant in SI
pk.set_default("cgs")
print(csts.G)  # Gravitational constant in CGS
```

---

### **Unit Conversion**
```python
from physkit.conversions import convert_unit

result = convert_unit(1, 'm', 'cm')
print(result)  # Outputs: 100.0
```

---

### **Equations**
```python
from physkit.equations import equations

# Default units are SI. Specify other units explicitly, like so:
mass = 1.0  # Solar mass
radius = equations.gravitational_radius(mass, 'M_sun', 'km') # input mass in solar masses
print(radius)  # Gravitational radius in km
```

---
### **Plot Styling**
This is just a quick way to make good looking plots simply.
```python
import physkit as pk
import matplotlib.pyplot as plt

x_data = [...]  # Example data for x and y
y_data = [...]

fig, ax = plt.subplots()
pk.plot_styler(x_data, y_data, ax=ax, title="test", 
               ylabel="y label", xlabel="x label", loglog=True)
plt.show()

```
---
### **Command Line Interface (CLI)**
```bash
physkit constant G --system cgs
```
```bash
physkit convert 1 m cm
```

---

## **Contributing**

Contributions are welcome!  

1. Fork the repository.  
2. Create a new branch (`git checkout -b feature-name`).  
3. Commit changes (`git commit -m "Add new feature"`).  
4. Push to your branch (`git push origin feature-name`).  
5. Open a Pull Request.

---

## **Issues**

If you encounter any issues or have suggestions, feel free to open an issue on [GitHub](https://github.com/sapphimars/physkit/issues).

---

## **Contact**

For inquiries, contact me via GitHub: [sapphimars](https://github.com/sapphimars)

