using Plots
include("physics.jl")

hi = 1000.0
ho = 500.0

ri = 0.005      # 5 mm inner radius
ro = 0.006      # 6 mm outer radius

L = 20.0         # 20 m heat exchanger
kw = 15.0        # W/(m·K)

Ti = 90.0        # °C, hot water
To = 20.0        # °C, cold water

m_dot_inner = 0.01       # kg/s
Cp_inner = 4180.0       # J/(kg·K)

m_dot_outer = 0.015       # kg/s
Cp_outer = 4180.0       # J/(kg·K)

N = 200                 # 200 segments
Ri, Rwall, Ro = cylindrical_resistances(
    hi,
    ho,
    ri,
    ro,
    L,
    kw
)

println("Inner convection resistance = ", Ri)
println("Wall conduction resistance = ", Rwall)
println("Outer convection resistance = ", Ro)

Rtotal = total_thermal_resistance(Ri, Rwall, Ro)

println("Total thermal resistance = ", Rtotal, " K/W")

Ui, Uo = overall_heat_transfer_coefficient(
    Rtotal,
    ri,
    ro,
    L
)

println("Overall U based on inner area = ", Ui, " W/(m²·K)")
println("Overall U based on outer area = ", Uo, " W/(m²·K)")

Ai = 2 * π * ri * L
Ao = 2 * π * ro * L

println("Ui * Ai = ", Ui * Ai) #verification that UiAi=UoAo
println("Uo * Ao = ", Uo * Ao)  




Q, Twi, Two = heat_transfer_and_wall_temperatures(
    Ti,
    To,
    Ri,
    Rwall,
    Ro
)

println()
println("Heat transfer rate = ", Q, " W")
println("Inner wall temperature = ", Twi, " °C")
println("Outer wall temperature = ", Two, " °C")





Ci = heat_capacity_rate(
    m_dot_inner,
    Cp_inner
)

Co = heat_capacity_rate(
    m_dot_outer,
    Cp_outer
)

println()
println("Inner fluid heat capacity rate = ", Ci, " W/K")
println("Outer fluid heat capacity rate = ", Co, " W/K")


Ti_out, To_out = outlet_temperatures(
    Ti,
    To,
    Q,
    Ci,
    Co
)  #whenever this function is called then for a segment ti_out is different for ti_in cuz usne heat lose ya gain ki hai due to radial heat transfer

println()
println("Inner fluid outlet temperature = ", Ti_out, " °C")
println("Outer fluid outlet temperature = ", To_out, " °C")

Q_inner, Q_outer, error_inner, error_outer =
    energy_balance_check(
        Ti,
        Ti_out,
        To,
        To_out,
        Ci,
        Co,
        Q
    )

println()
println("----- Energy Balance Check -----")

println("Calculated Q       = ", Q, " W")
println("Heat lost by inner = ", Q_inner, " W")
println("Heat gained by outer = ", Q_outer, " W")

println()
println("Inner-side error = ", error_inner, " W")
println("Outer-side error = ", error_outer, " W")




ΔL = segment_length(L, N)

println()
println("Number of segments = ", N)
println("Length of each segment = ", ΔL, " m")




Ti_next, To_next, Twi, Two, Qsegment =
    simulate_segment(
        Ti,
        To,
        hi,
        ho,
        ri,
        ro,
        ΔL,
        kw,
        Ci,
        Co
    )

println()
println("----- One Segment Results -----")

println("Segment heat transfer = ", Qsegment, " W")

println("Inner fluid before = ", Ti, " °C")
println("Inner fluid after  = ", Ti_next, " °C")

println()

println("Outer fluid before = ", To, " °C")
println("Outer fluid after  = ", To_next, " °C")

println()

println("Inner wall temperature = ", Twi, " °C")
println("Outer wall temperature = ", Two, " °C")


Ti_profile,
To_profile,
Twi_profile,
Two_profile,
Q_profile = simulate_parallel_flow(
    Ti,
    To,
    hi,
    ho,
    ri,
    ro,
    L,
    kw,
    Ci,
    Co,
    N
)  
println()
println("----- Full Heat Exchanger Results -----")

println("Inner fluid inlet  = ", Ti_profile[1], " °C")
println("Inner fluid outlet = ", Ti_profile[end], " °C")

println()

println("Outer fluid inlet  = ", To_profile[1], " °C")
println("Outer fluid outlet = ", To_profile[end], " °C")

println()

println("Total heat transferred = ", sum(Q_profile), " W") #note that this value is smaller than the one we calculated 
 #assuming L as whole length because we thought about only the inlet and out temp is the drivng force but in reality ◬T decreases as we go across x direction
 #The segmented result is lower because the temperature difference is not actually constant throughout the exchanger.

 x = segment_positions(L, N)

println()
println("First position = ", x[1], " m")
println("Last position = ", x[end], " m")
println("Number of points = ", length(x))


Q_inner,
Q_outer,
Q_segments,
error_inner,
error_outer = validate_simulation(
    Ti_profile,
    To_profile,
    Ci,
    Co,
    Q_profile
)

println()
println("----- Full Simulation Validation -----")

println("Heat lost by inner fluid = ", Q_inner, " W")
println("Heat gained by outer fluid = ", Q_outer, " W")
println("Sum of segment heat transfer = ", Q_segments, " W")

println()
println("Inner error = ", error_inner, " W")
println("Outer error = ", error_outer, " W")


segment = 50

Twi = Twi_profile[segment]
Qsegment = Q_profile[segment]

r, T_radial = radial_temperature_profile(
    Twi,
    Qsegment,
    ri,
    ro,
    kw,
    ΔL,
    100
)

# ------------------------------------------
# COUNTER-CURRENT FLOW SIMULATION
# ------------------------------------------

Ti_counter,
To_counter,
Twi_counter,
Two_counter,
Q_counter = simulate_countercurrent_flow(
    Ti,
    To,
    hi,
    ho,
    ri,
    ro,
    L,
    kw,
    Ci,
    Co,
    N
)

println()
println("----- Counter-Current Results -----")

println("Inner fluid inlet  = ", Ti_counter[1], " °C")
println("Inner fluid outlet = ", Ti_counter[end], " °C")

println()

# Remember: cold fluid flows in the opposite direction

println("Outer fluid inlet  = ", To_counter[end], " °C")
println("Outer fluid outlet = ", To_counter[1], " °C")

println()

println(
    "Total heat transferred = ",
    sum(Q_counter),
    " W"
)
# ------------------------------------------
# COUNTER-CURRENT ENERGY BALANCE VALIDATION
# ------------------------------------------

Q_hot_counter =
    Ci * (
        Ti_counter[1] -
        Ti_counter[end]
    )

Q_cold_counter =
    Co * (
        To_counter[1] -
        To_counter[end]
    )

Q_segments_counter =
    sum(Q_counter)

error_hot_counter =
    abs(Q_hot_counter - Q_segments_counter)

error_cold_counter =
    abs(Q_cold_counter - Q_segments_counter)

println()
println("----- Counter-Current Energy Balance -----")

println("Heat lost by hot fluid = ",
    Q_hot_counter, " W")

println("Heat gained by cold fluid = ",
    Q_cold_counter, " W")

println("Sum of segment heat transfer = ",
    Q_segments_counter, " W")

println()

println("Hot-side error = ",
    error_hot_counter, " W")

println("Cold-side error = ",
    error_cold_counter, " W")

    # ------------------------------------------
# COUNTER-CURRENT TEMPERATURE DIFFERENCE CHECK
# ------------------------------------------

ΔT_counter =
    Ti_counter .- To_counter

println()
println("----- Counter-Current Temperature Difference Check -----")

println(
    "Minimum temperature difference = ",
    minimum(ΔT_counter),
    " °C"
)

println(
    "Maximum temperature difference = ",
    maximum(ΔT_counter),
    " °C"
)

# ------------------------------------------
# COUNTER-CURRENT AXIAL TEMPERATURE PROFILE
# ------------------------------------------

x_counter = segment_positions(L, N)

# Wall temperatures are calculated per segment,
# so assign each wall temperature to the midpoint
# of its corresponding segment.

x_wall_counter = [
    (x_counter[i] + x_counter[i + 1]) / 2
    for i in 1:N
]

p_counter = plot(
    x_counter,
    Ti_counter,

    label = "Hot Fluid",
    xlabel = "Axial Position (m)",
    ylabel = "Temperature (°C)",
    title = "Counter-Current Heat Exchanger Temperature Profile",

    linewidth = 2
)

plot!(
    p_counter,
    x_counter,
    To_counter,

    label = "Cold Fluid",
    linewidth = 2
)

plot!(
    p_counter,
    x_wall_counter,
    Twi_counter,

    label = "Inner Wall",
    linewidth = 2,
    linestyle = :dash
)

plot!(
    p_counter,
    x_wall_counter,
    Two_counter,

    label = "Outer Wall",
    linewidth = 2,
    linestyle = :dash
)

display(p_counter)

# ------------------------------------------
# PARALLEL VS COUNTER-CURRENT COMPARISON
# ------------------------------------------

# Axial positions
x_parallel = segment_positions(L, N)
x_counter = segment_positions(L, N)

p_comparison = plot(
    x_parallel,
    Ti_profile,

    label = "Hot Fluid - Parallel",
    xlabel = "Axial Position (m)",
    ylabel = "Temperature (°C)",
    title = "Parallel vs Counter-Current Heat Exchanger",

    linewidth = 2
)

# Parallel cold fluid
plot!(
    p_comparison,
    x_parallel,
    To_profile,

    label = "Cold Fluid - Parallel",
    linewidth = 2
)

# Counter-current hot fluid
plot!(
    p_comparison,
    x_counter,
    Ti_counter,

    label = "Hot Fluid - Counter",
    linewidth = 2,
    linestyle = :dash
)

# Counter-current cold fluid
plot!(
    p_comparison,
    x_counter,
    To_counter,

    label = "Cold Fluid - Counter",
    linewidth = 2,
    linestyle = :dash
)

display(p_comparison)

println()
println("----- FLOW CONFIGURATION COMPARISON -----")

Q_parallel = sum(Q_profile)
Q_countercurrent = sum(Q_counter)

println("Parallel-flow heat transfer = ",
    Q_parallel, " W")

println("Counter-current heat transfer = ",
    Q_countercurrent, " W")

improvement =
    (
        Q_countercurrent -
        Q_parallel
    ) / Q_parallel * 100

println()

println(
    "Counter-current improvement = ",
    improvement,
    " %"
)