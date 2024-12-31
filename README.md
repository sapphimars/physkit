# Physkit - A Python toolkit for constants, unit conversions, and equations

Physkit is a Python library for performing scientific computations, unit conversions, and working with physical constants and equations. It provides tools for astrophysical and physical calculations, supporting multiple unit systems and CLI functionalities.

---

## **License**

This project is licensed under the terms of the **GNU General Public License v3.0**.  
See the [LICENSE](./LICENSE) file for details.  

Copyright (C) 2024 [sapphimars](https://github.com/sapphimars)  

---

### **Third-Party Licenses**

This project makes use of the following third-party libraries:

- **[Matplotlib](https://matplotlib.org/)** - Licensed under the Matplotlib License Agreement. See [licenses/LICENSE-MATPLOTLIB.txt](./licenses/LICENSE-MATPLOTLIB.txt).
- **[Pint](https://github.com/hgrecco/pint)** - Licensed under the BSD 3-Clause License. See [licenses/LICENSE-PINT.txt](./licenses/LICENSE-PINT.txt).

Note: Third-party libraries included in this project are licensed under their respective terms. See the `licenses/` directory for full details.

---

## **Installation**

### **Requirements**

- Python 3.8 or higher
- Dependencies:
  - **Matplotlib** >= 3.5
  - **Pint** >= 0.20

Install dependencies via pip
```bash
pip install -r requirements.txt
```


---

### **Install via pip**
```bash
pip install physkit
```


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

### **Command Line Interface (CLI)**
In terminal:
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
