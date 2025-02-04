# Quantum Tunnelling Simulation

This repository provides a Python-based simulation of quantum tunnelling, an effect in quantum mechanics where particles can penetrate potential barriers that they classically shouldn't be able to cross. 

The simulation solves the time-independent Schrödinger equation for a particle in a potential, described by:

$$
-\frac{\hbar^2}{2m} \frac{d^2}{dx^2} \psi(x) + V(x)\psi(x) = E\psi(x)
$$

Where $$\ \psi(x) \$$ is the wavefunction, $$\ V(x)\$$ is the potential, $$\ E \$$ is the energy, and $$\ \hbar \$$ is the reduced Planck constant. The barrier is modeled as a square potential, where the particle encounters the barrier if $$\ E < V_0 \$$.

### Key Features:
- Numerical solution of the Schrödinger equation using finite difference methods.
- Visualization of the wavefunction’s evolution across a potential barrier.
- Interaction with the model through a Streamlit web app, enabling real-time modifications to barrier parameters (height, width) and particle energy.

### Streamlit App:
[Quantum Tunnel to the Streamlit app](https://conors-tunnelling-demo.streamlit.app/)
