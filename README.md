# Relativistic Boson Star Numerical Solutions

This repository contains a project dedicated to characterizing numerical solutions for relativistic boson stars, both interacting and non-interacting cases (Massive Boson Stars and Mini-Boson Stars).

---

## Overview

Compact objects like white dwarfs and neutron stars are part of the entities that populate the universe and play a role in astrophysical and cosmological observations. These objects belong to the family of compact fermionic objects that conform to the universe.

Parallel to these objects are **boson stars**, hypothetical compact objects made of spin $\(s = 0\)$ bosons, and **Proca stars**, compact objects made of spin $\(s = 1\)$ bosons. These hypothetical objects could account for dark matter: massive self-gravitating scalar (or vector) fields can form stable astrophysical (or subatomic) objects supported by self-gravity and Heisenberg’s uncertainty principle.

Although the cold dark matter model (CDM) within the standard cosmological model $\(\Lambda\)$ CDM explains large-scale structure so well, it faces issues on galactic scales:  

- CDM simulations predict cuspy density profiles at galactic centers, whereas rotation curves suggest smooth cores.  
- The predicted number of satellite galaxies around galactic halos is higher than observed.  

If dark matter is composed of boson particles in a Bose-Einstein condensate (as in boson stars), these problems might be mitigated. Boson and Proca star phenomenology could address the cusp and missing satellite problems.

---

## Scientific Context

A **compact object** is a collective arrangement of quantum particles forming a stable macroscopic entity. John Wheeler (1955) proposed stable solutions to the electromagnetic field under general relativity called **geons**, which are unstable. Conversely, **boson** and **Proca stars** are stable solutions to the Einstein-Klein-Gordon and Einstein-Proca systems, respectively.

These objects form a macroscopic condensate of integer spin particles characterized by a single wave function $\(\Psi(t, \vec{x})\)$. Boson stars may have masses comparable to neutron stars or even larger, making them compelling subjects of study and candidates for dark matter.

---

## Relativistic and Non-Relativistic Boson Stars

We study:

1. **Einstein-Klein-Gordon system** (relativistic boson stars)  

These systems can include a self-interaction term, parameterized by a dimensionless coupling constant $\(\lambda\)$. When $\(\lambda\)$ is negligible, the systems reduce to:

- **Mini-boson stars** (relativistic)  
- **Schrödinger-Poisson system** (non-relativistic)

If the mass is small, relativistic effects can be neglected, and the Newtonian approximation is valid. For masses approaching the Kaup limit $(\(M_{\text{Kaup}} = 0.633 M_{\text{Pl}}^2 / m\))$, relativistic effects are crucial. Beyond this mass, no equilibrium configurations exist, similar to fermionic stars.  

Self-interacting boson stars have a maximum mass $\(M_{\text{max}} = 0.06 \lambda M_{\text{Pl}}^3 / m^2\)$, comparable to fermion stars, making them particularly interesting.

---

## Astrophysical Relevance

Due to their mass, size, and stability, boson stars could be viable dark matter candidates. Observations via gravitational waves and lensing could shed light on their existence. Non-relativistic boson stars can be modeled as **Bose-Einstein condensates** with wave function $\(\psi(t, \vec{x})\)$, governed by:

- **Gross-Pitaevskii equation** (wave function evolution)  
- **Poisson equation** (Newtonian gravitational potential)

Equilibrium arises from the balance of self-interaction pressure, self-gravity, and quantum pressure.

--- 

# Purpose of the Repository

The main purpose of this repository is to present and analyze **numerical solutions to Relativistic Boson Stars**. We provide a clear methodology in the file `metodos.py`, and present results through Jupyter notebooks. In particular, the notebooks illustrate:

- Equilibrium configurations for different values of the self-interaction parameter $\(\lambda\)$.  
- Excited configurations $n= 1,2,...$ for mini-Boson Stars ($\lambda$).  
- Gravitational potential behavior for a mini-Boson Star.  
- Mass $\(M\)$ and $\(N \cdot m_{0}\)$ as a function of the initial radial profile $\sigma_{0}$.  
- Total mass $M_T$ as a function of the central radial profile $\sigma_{0}$.  
- Total mass  $M_T$ as a function of the radial profile $\sigma_{0}$ for different values of $\(\lambda\)$.  

This repository is intended as both a **research tool** and a **reference** for students and researchers interested in the numerical study of boson stars.

---
---

## Background and References

This project is based on and inspired by several key works in the field of boson star research, including:  

- *TOPICAL REVIEW: General Relativistic Boson Stars* – Franz E. Schunck, Eckehard W. Mielke  
- *Boson Stars and Oscillatons: A Review* – Luca Visinelli  
- *The Structure and Formation of Boson Stars* – Andrew R. Liddle, Mark S. Madsen, 1992  

Along with these and other related works, this repository builds on the established theoretical and numerical framework of boson stars.
The project is also part of the ongoing research and writing of the author’s **Pdh Thesis**, where the numerical characterization of relativistic boson stars plays a central role. In addition, the **numerical methodology** applied here is mainly based on:  

- *Linear stability of nonrelativistic self-interacting boson stars* – Emmanuel Chávez Nambo, Alberto Diez-Tejedor, Armando A. Roque, and Olivier Sarbach  
