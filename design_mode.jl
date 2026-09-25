using Plots

include("physics.jl")


# ==========================================
# HEAT EXCHANGER DESIGN MODE
# ==========================================


# ------------------------------------------
# DESIGN INPUTS
# ------------------------------------------

hi = 1000.0
ho = 500.0

ri = 0.005
ro = 0.006

kw = 15.0


# ------------------------------------------
# FLUID CONDITIONS
# ------------------------------------------

Ti = 90.0       # Hot fluid inlet °C
To = 20.0       # Cold fluid inlet °C


# ------------------------------------------
# FLOW CONDITIONS
# ------------------------------------------

m_dot_inner = 0.01
Cp_inner = 4180.0

m_dot_outer = 0.015
Cp_outer = 4180.0


# ------------------------------------------
# TARGET DESIGN CONDITION
# ------------------------------------------

T_target = 60.0

println("==========================================")
println("        HEAT EXCHANGER DESIGN MODE")
println("==========================================")

println()
println("Hot fluid inlet temperature = ", Ti, " °C")
println("Cold fluid inlet temperature = ", To, " °C")
println("Target hot outlet temperature = ", T_target, " °C")


# ------------------------------------------
# HEAT CAPACITY RATES
# ------------------------------------------

Ci = heat_capacity_rate(
    m_dot_inner,
    Cp_inner
)

Co = heat_capacity_rate(
    m_dot_outer,
    Cp_outer
)


# ------------------------------------------
# NUMBER OF SEGMENTS
# ------------------------------------------

N = 100


# ------------------------------------------
# FIND REQUIRED LENGTH
# ------------------------------------------

L_required,
Ti_design,
To_design,
Twi_design,
Two_design,
Q_design,
iterations = find_required_length(
    Ti,
    To,
    T_target,
    hi,
    ho,
    ri,
    ro,
    kw,
    Ci,
    Co,
    N;

    L_low = 1,   # add low
    L_high = 20.0  #add high
)


# ==========================================
# DESIGN RESULTS
# ==========================================

println()
println("==========================================")
println("           DESIGN RESULTS")
println("==========================================")

println()

println(
    "Required heat exchanger length = ",
    L_required,
    " m"
)

println(
    "Target hot outlet temperature = ",
    T_target,
    " °C"
)

println(
    "Calculated hot outlet temperature = ",
    Ti_design[end],
    " °C"
)

println(
    "Cold fluid outlet temperature = ",
    To_design[end],
    " °C"
)

println(
    "Total heat transferred = ",
    sum(Q_design),
    " W"
)

println(
    "iterations required = ",
    iterations
)


# ==========================================
# DESIGN VALIDATION
# ==========================================

Q_hot =
    Ci *
    (
        Ti -
        Ti_design[end]
    )

Q_cold =
    Co *
    (
        To_design[end] -
        To
    )

Q_segments =
    sum(Q_design)

println()
println("==========================================")
println("         ENERGY BALANCE CHECK")
println("==========================================")

println()

println("Heat lost by hot fluid = ", Q_hot, " W")

println("Heat gained by cold fluid = ", Q_cold, " W")

println("Segment heat transfer = ", Q_segments, " W")

println()

println(
    "Hot-side error = ",
    abs(Q_hot - Q_segments),
    " W"
)

println(
    "Cold-side error = ",
    abs(Q_cold - Q_segments),
    " W"
)