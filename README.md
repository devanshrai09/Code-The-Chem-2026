## Virtual Heat Exchanger Simulator

A physics-based **virtual heat exchanger simulator** developed in Julia for studying heat transfer, temperature evolution, flow configurations, and basic equipment design.

The simulator models a cylindrical heat exchanger using a **segment-by-segment numerical approach**, allowing the temperature of the hot fluid, cold fluid, and heat-exchanger wall to be tracked along the equipment.

The project aims to go beyond a simple inlet–outlet heat-transfer calculation by providing insight into the **internal thermal behaviour of the equipment**.

---

##  Problem Statement

Heat exchangers are commonly analyzed using overall heat-transfer calculations that provide quantities such as heat-transfer rate and outlet temperatures.

However, these calculations can hide what happens **inside the equipment**.

In a real heat exchanger:

- Fluid temperatures change continuously along the exchanger.
- The local temperature driving force changes with position.
- Heat passes through multiple thermal-resistance layers.
- Parallel-flow and counter-current configurations behave differently.
- Equipment dimensions may need to be determined from a desired performance target.

This project develops a virtual simulator that models these effects computationally.

---

#  Features

## 1. Cylindrical Heat-Transfer Model

The simulator represents a cylindrical heat exchanger in which heat passes through three thermal-resistance mechanisms:

```text
Hot Fluid
    │
    │  Convection
    ▼
Inner Wall
    │
    │  Cylindrical Conduction
    ▼
Outer Wall
    │
    │  Convection
    ▼
Cold Fluid
