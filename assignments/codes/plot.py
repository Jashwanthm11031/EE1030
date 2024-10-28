import ctypes
import numpy as np
import matplotlib.pyplot as plt
import os

# Load the shared library
# Change 'libdata_generator.so' to 'data_generator.dll' on Windows
lib = ctypes.CDLL(os.path.abspath("libdata_generator.so"))

# Define the function argument and return types
lib.generate_sqrt_x_minus_1.argtypes = [
    ctypes.POINTER(ctypes.c_double),  # x_values
    ctypes.POINTER(ctypes.c_double),  # y_values
    ctypes.c_int                      # count
]

# Generate x values in Python
x_values = np.linspace(1, 10, 100)
y_values = np.zeros_like(x_values)

# Call the C function to fill y_values
lib.generate_sqrt_x_minus_1(
    x_values.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    y_values.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    len(x_values)
)

# Plot the data
plt.plot(x_values, y_values, label=r'$\sqrt{x - 1}$', color='blue')
plt.axvline(x=1, color='red', linestyle='--', label='x = 1')
plt.axvline(x=5, color='green', linestyle='--', label='x = 5')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Plot of $\sqrt{x - 1}$ generated by C and plotted in Python')
plt.legend()
plt.grid(True)
plt.show()

