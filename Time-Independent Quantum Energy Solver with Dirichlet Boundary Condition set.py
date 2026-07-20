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
   "id": "68a0312d-e4e7-4391-bab1-5ea1d7e67710",
   "metadata": {},
   "outputs": [],
   "source": [
    "x_max = float(input(\"x_max:\"))\n",
    "x_min = float(input(\"x_min:\"))\n",
    "Nx = int(input(\"Number of data points:\"))\n",
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
    "constant_V = input(\"V(x) = \")\n",
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
    "print (\"\\nSuccess! V(x) numbers calculated completely.\") "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "bfe7006e-fe7e-43ed-ab9f-aefbc3466a6b",
   "metadata": {},
   "outputs": [],
   "source": [
    "Kinetic_Multiplier = -1 / (2 * (dx**2))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "3e8ee00c-7b7f-424a-a3b0-42f0d514fcd6",
   "metadata": {},
   "outputs": [],
   "source": [
    "main_diag = np.full(Nx, -2.0)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "3d304911-0c47-4315-abf2-d2b350bad8ca",
   "metadata": {},
   "outputs": [],
   "source": [
    "off_diag = np.full(Nx - 1, 1.0) "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "6c0b6807-9890-4a9a-ba22-f6691b46b65e",
   "metadata": {},
   "outputs": [],
   "source": [
    "Second_Derivative_Matrix = np.diag(main_diag) + np.diag(off_diag, 1) + np.diag(off_diag, -1)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "ba13518b-982c-4c3c-825e-2467e5329cb8",
   "metadata": {},
   "outputs": [],
   "source": [
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
   "id": "bc712b02-ee48-45ef-963f-3ba3efd0beee",
   "metadata": {},
   "outputs": [],
   "source": [
    "Hamiltonian_Matrix[0, :] = 0\n",
    "Hamiltonian_Matrix[0, 0] = 1e9\n",
    "\n",
    "Hamiltonian_Matrix[-1, :] = 0\n",
    "Hamiltonian_Matrix[-1, -1] = 1e9"
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
    "E_0 = energies[0] "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "8e59fdde-d8e3-4fa3-8d33-c2133cc73f48",
   "metadata": {},
   "outputs": [],
   "source": [
    "psi_0 = wavefunctions[:,0] "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "1c32dafb-6e93-4c26-b519-4bf0a7bf6033",
   "metadata": {},
   "outputs": [],
   "source": [
    "E_1 = energies[1] "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "d23ebe45-047e-49f2-bd4b-d432e0ed1237",
   "metadata": {},
   "outputs": [],
   "source": [
    "psi_1 = wavefunctions[:,1] "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "9b2626e6-641c-4e63-a7c4-126df8b1ab73",
   "metadata": {},
   "outputs": [],
   "source": [
    "E_2 = energies[2] "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "749d11f6-eefa-4d3b-9964-a3926ccebd4b",
   "metadata": {},
   "outputs": [],
   "source": [
    "psi_2 = wavefunctions[:,2] "
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
    "plt.plot(x, psi_0, label=f\"Ground State (E0 = {E_0:.4f})\", color=\"blue\")\n",
    "plt.plot(x, psi_1, label=f\"1st Excited State (E1 = {E_1:.4f})\", color=\"red\")\n",
    "plt.plot(x, psi_2, label=f\"2nd Excited State (E2 = {E_2:.4f})\", color=\"green\")\n",
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
