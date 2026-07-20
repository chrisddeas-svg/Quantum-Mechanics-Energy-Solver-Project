{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "87b0efae-9943-455b-b2de-764c7b4cf784",
   "metadata": {},
   "outputs": [],
   "source": [
    "import numpy as np\n",
    "import matplotlib.pyplot as plt\n",
    "import scipy.linalg as la\n",
    "import re "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "48438910-55a4-4350-ae4c-98036dd0375b",
   "metadata": {},
   "outputs": [],
   "source": [
    "print(\"This solver operates in Atomic units, hence for any prior calculations required, assume mass and the reduced planck constant both equal 1.\") "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "68a0312d-e4e7-4391-bab1-5ea1d7e67710",
   "metadata": {},
   "outputs": [],
   "source": [
    "x_max = float(input(\"maximum x value of graph, x_max = \"))\n",
    "x_min = float(input(\"minimum x value of graph, x_min =\"))\n",
    "Nx = int(input(\"Number of data points, (the more data points the more precise the result):\"))\n",
    "\n",
    "x = np.linspace (x_min, x_max, Nx)\n",
    "dx = x[1] - x[0]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "1e608cfd-01fa-4a0f-9429-44a5052d0f49",
   "metadata": {},
   "outputs": [],
   "source": [
    "constant_V = input(\"V(x)= ?, (two examples: 2 * sin(x) + delta, x ** 2 + sin (2 * x))\")\n",
    "\n",
    "alpha = 0.0\n",
    "is_delta = False\n",
    "if \"delta\" in constant_V.lower():\n",
    "    alpha = float(input(\"delta function detected, alpha (strength) value?:\"))\n",
    "    is_delta = True \n",
    "    constant_V = constant_V.lower().replace(\"delta\", \"0\") \n",
    "\n",
    "Math_Functions = {\"x\", \"m\", \"sin\", \"cos\", \"tan\", \"exp\", \"pi\", \"sqrt\", \"log\", \"sinh\", \"cosh\", \"abs\", \"tanh\", \"arcsin\", \"arccos\", \"arctan\"} \n",
    "eval_env = {\"x\": x, \"m\": 1.0, \"pi\": np.pi, \"e\": np.e} \n",
    "for f in Math_Functions: \n",
    "    if f not in {\"x\",\"m\", \"pi\"}: \n",
    "        eval_env[f] = getattr(np, f, None) \n",
    "\n",
    "words = set(re.findall(r'\\b[a-zA-Z_][a-zA-Z0-9_]*\\b', constant_V)) \n",
    "undetected_variables = [w for w in words if w not in Math_Functions and w != \"e\"] \n",
    "for var in undetected_variables: \n",
    "    user_value = input(f\"detected variable '{var}'. Please enter its value: \") \n",
    "    eval_env[var] = float(user_value) \n",
    "\n",
    "V_numbers= eval(constant_V, {\"__builtins__\": None}, eval_env) \n",
    "if isinstance(V_numbers, (int, float)):\n",
    "    V_numbers= np.full(Nx, float(V_numbers)) \n",
    "\n",
    "if is_delta:\n",
    "    center_idx=Nx//2\n",
    "    V_numbers[center_idx] += (alpha/dx)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "3d304911-0c47-4315-abf2-d2b350bad8ca",
   "metadata": {},
   "outputs": [],
   "source": [
    "main_diag = np.full(Nx, -2.0) \n",
    "off_diag = np.full(Nx - 1, 1.0) \n",
    "\n",
    "Second_Derivative_Matrix = np.diag(main_diag) + np.diag(off_diag, 1) + np.diag(off_diag, -1)\n",
    "\n",
    "Kinetic_Multiplier = -1 / (2 * (dx**2))\n",
    "\n",
    "T = Kinetic_Multiplier * Second_Derivative_Matrix"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "2946c525-22d3-4f9e-8b79-0083bc16398e",
   "metadata": {},
   "outputs": [],
   "source": [
    "if np.isscalar(V_numbers):\n",
    "    V_numbers = np.full_like(x, V_numbers)\n",
    "\n",
    "V_Matrix = np.diag(V_numbers)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "8e946349-4d07-4edd-ab7f-852e8e757a16",
   "metadata": {},
   "outputs": [],
   "source": [
    "Hamiltonian_Matrix = T + V_Matrix "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "15979824-ceff-48a4-961d-e6f472d54817",
   "metadata": {},
   "outputs": [],
   "source": [
    "print(\"Choose boundary condition of particle:\")\n",
    "print(\"1 = Infinite Wall (Dirichlet)\")\n",
    "print(\"2 = Flat Slope (Neumann)\")\n",
    "Choice = input(\"Enter 1 or 2: \") \n",
    "\n",
    "if Choice == \"1\":\n",
    "    Hamiltonian_Matrix[0, :] = 0\n",
    "    Hamiltonian_Matrix[0, 0] = 1e9\n",
    "    Hamiltonian_Matrix[-1, :] = 0\n",
    "    Hamiltonian_Matrix[-1, -1] = 1e9\n",
    "\n",
    "elif Choice == \"2\":\n",
    "    Hamiltonian_Matrix[0, 0] = Hamiltonian_Matrix[0, 0] + Kinetic_Multiplier\n",
    "    Hamiltonian_Matrix[-1, -1] = Hamiltonian_Matrix[-1, -1] + Kinetic_Multiplier\n",
    "\n",
    "else: \n",
    "    print(\"Invalid Choice\") "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "573b5727-1fe4-40b8-933b-62f16605a1b0",
   "metadata": {},
   "outputs": [],
   "source": [
    "energies, wavefunctions = la.eigh(Hamiltonian_Matrix) "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "67ec3fc0-761c-43c1-bd76-351628c76079",
   "metadata": {},
   "outputs": [],
   "source": [
    "E_0 = energies[0] \n",
    "psi_0 = wavefunctions[:,0] \n",
    "\n",
    "E_1 = energies[1] \n",
    "psi_1 = wavefunctions[:,1] \n",
    "\n",
    "E_2 = energies[2] \n",
    "psi_2 = wavefunctions[:,2] \n",
    "\n",
    "E_3 = energies[3] \n",
    "psi_3 = wavefunctions[:,3] \n",
    "\n",
    "E_4 = energies[4] \n",
    "psi_4 = wavefunctions[:,4] "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "56621b7a-7a71-4ce8-8dde-9e5a66586b2d",
   "metadata": {},
   "outputs": [],
   "source": [
    "plt.figure(figsize=(16, 9))\n",
    "\n",
    "plt.plot(x, psi_0, label=f\"Ground State (E0 = {E_0:.6f})\", color=\"blue\")\n",
    "\n",
    "\n",
    "plt.title(r\"Ground State Wavefunction ($\\psi$)\") \n",
    "plt.xlabel(\"Position (x)\")\n",
    "plt.ylabel(r\"$\\psi(x)$\") \n",
    "plt.grid(True, linestyle=\"--\", alpha=0.5)\n",
    "plt.legend()\n",
    "\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "c8a365fa-2d19-4414-a4d5-62410fd24a44",
   "metadata": {},
   "outputs": [],
   "source": [
    "plt.figure(figsize=(16, 9))\n",
    "plt.plot(x, psi_1, label=f\"1st Excited State (E1 = {E_1:.6f})\", color=\"red\")\n",
    "plt.title(r\"1st Excited State Wavefunction ($\\psi$)\") \n",
    "plt.xlabel(\"Position (x)\")\n",
    "plt.ylabel(r\"$\\psi(x)$\") \n",
    "plt.grid(True, linestyle=\"--\", alpha=0.5)\n",
    "plt.legend()\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "efe32ea5-62b3-4538-bc9f-7f14eb372956",
   "metadata": {},
   "outputs": [],
   "source": [
    "plt.figure(figsize=(16, 9))\n",
    "plt.plot(x, psi_2, label=f\"2nd Excited State (E2 = {E_2:.6f})\", color=\"green\")\n",
    "plt.title(r\"2nd Excited State Wavefunction  ($\\psi$)\") \n",
    "plt.xlabel(\"Position (x)\")\n",
    "plt.ylabel(r\"$\\psi(x)$\") \n",
    "plt.grid(True, linestyle=\"--\", alpha=0.5)\n",
    "plt.legend()\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "296a2c15-15fb-4a2e-9d46-d3cd7ecb214f",
   "metadata": {},
   "outputs": [],
   "source": [
    "plt.figure(figsize=(16, 9))\n",
    "plt.plot(x, psi_3, label=f\"3rd Excited State (E3 = {E_3:.6f})\", color=\"purple\")\n",
    "plt.title(r\"3rd Excited State Wavefunction  ($\\psi$)\") \n",
    "plt.xlabel(\"Position (x)\")\n",
    "plt.ylabel(r\"$\\psi(x)$\") \n",
    "plt.grid(True, linestyle=\"--\", alpha=0.5)\n",
    "plt.legend()\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "43a4441b-8a54-4ac3-bff2-67cafe38193d",
   "metadata": {},
   "outputs": [],
   "source": [
    "plt.figure(figsize=(16, 9))\n",
    "plt.plot(x, psi_4, label=f\"4th Excited State (E4 = {E_4:.6f})\", color=\"orange\")\n",
    "plt.title(r\"4th Excited State Wavefunction  ($\\psi$)\") \n",
    "plt.xlabel(\"Position (x)\")\n",
    "plt.ylabel(r\"$\\psi(x)$\") \n",
    "plt.grid(True, linestyle=\"--\", alpha=0.5)\n",
    "plt.legend()\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "b1a11384-11ac-4158-9bee-f205b851ee51",
   "metadata": {},
   "outputs": [],
   "source": [
    "plt.figure(figsize=(16, 9))\n",
    "\n",
    "plt.plot(x, psi_0, label=f\"Ground State (E0 = {E_0:.6f})\", color=\"blue\")\n",
    "plt.plot(x, psi_1, label=f\"1st Excited State (E1 = {E_1:.6f})\", color=\"red\")\n",
    "plt.plot(x, psi_2, label=f\"2nd Excited State (E2 = {E_2:.6f})\", color=\"green\")\n",
    "plt.plot(x, psi_3, label=f\"3rd Excited State (E3 = {E_3:.6f})\", color=\"purple\")\n",
    "plt.plot(x, psi_4, label=f\"4th Excited State (E4 = {E_4:.6f})\", color=\"orange\")\n",
    "\n",
    "plt.title(r\"Numerical Wavefunctions ($\\psi$)\") \n",
    "plt.xlabel(\"Position (x)\")\n",
    "plt.ylabel(r\"$\\psi(x)$\") \n",
    "plt.grid(True, linestyle=\"--\", alpha=0.5)\n",
    "plt.legend()\n",
    "\n",
    "plt.show()"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.14.6"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
