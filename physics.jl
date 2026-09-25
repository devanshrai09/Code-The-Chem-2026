function cylindrical_resistances(hi, ho, ri, ro, L, kw)

    Ai = 2 * π * ri * L
    Ao = 2 * π * ro * L

    Ri = 1 / (hi * Ai)

    Rwall = log(ro / ri) / (2 * π * kw * L)

    Ro = 1 / (ho * Ao)

    return  abs(Ri), abs(Rwall), abs(Ro)
end

function total_thermal_resistance(Ri, Rwall, Ro)

    Rtotal = Ri + Rwall + Ro

    return Rtotal
end

function overall_heat_transfer_coefficient(Rtotal, ri, ro, L)

    # Inner and outer heat transfer areas
    Ai = 2 * π * ri * L
    Ao = 2 * π * ro * L

    # Overall heat transfer coefficients
    Ui = 1 / (Rtotal * Ai)
    Uo = 1 / (Rtotal * Ao)

    return Ui, Uo
end

function heat_transfer_and_wall_temperatures(
    Ti,
    To,
    Ri,
    Rwall,
    Ro
)

    # Total thermal resistance
    Rtotal = Ri + Rwall + Ro

    # Heat transfer rate
    Q = (Ti - To) / Rtotal

    # Inner wall temperature
    Twi = Ti - Q * Ri

    # Outer wall temperature
    Two = Twi - Q * Rwall

    return Q, Twi, Two
end

function heat_capacity_rate(m_dot, Cp) #to calculate C for that mass in L length tube

    C = m_dot * Cp

    return C
end

function outlet_temperatures(
    Ti_in,
    To_in,
    Q,
    Ci,
    Co
)

    # Inner fluid loses heat
    Ti_out = Ti_in - Q / Ci

    # Outer fluid gains heat
    To_out = To_in + Q / Co

    return Ti_out, To_out
end

function energy_balance_check(
    Ti_in,
    Ti_out,
    To_in,
    To_out,
    Ci,
    Co,
    Q
)

    # Heat lost by inner fluid
    Q_inner = Ci * (Ti_in - Ti_out)

    # Heat gained by outer fluid
    Q_outer = Co * (To_out - To_in)

    # Errors compared with calculated heat transfer
    error_inner = abs(Q - Q_inner)
    error_outer = abs(Q - Q_outer)

    return Q_inner, Q_outer, error_inner, error_outer
end


function segment_length(L, N) # dividing our heat exchanger in small segments

    ΔL = L / N

    return ΔL
end



function simulate_segment(
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

    # 1. Thermal resistances for this small segment
    Ri, Rwall, Ro = cylindrical_resistances(
        hi,
        ho,
        ri,
        ro,
        ΔL,
        kw
    )  #for each segment now you calculate resistance and therefore within that ◬L length

    # 2. Heat transfer and wall temperatures
    Q, Twi, Two = heat_transfer_and_wall_temperatures(
        Ti,
        To,
        Ri,
        Rwall,
        Ro
    ) #how much heat tranfer is happening from that segment and what are the wall temperatures

    # 3. Update fluid temperatures

    # Inner fluid loses heat
    Ti_next = Ti - Q / Ci       #for the next segment the Ti_out of previous is Ti_in samjhe?? YES SIR!!!!


    # Outer fluid gains heat
    To_next = To + Q / Co  # same for the temp_out!!
    

    return Ti_next, To_next, Twi, Two, Q
end


function simulate_parallel_flow(
    Ti_in,
    To_in,
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

    # Length of each segment
    ΔL = L / N

    # Arrays to store temperatures
    Ti_profile = zeros(N + 1)
    To_profile = zeros(N + 1)

    Twi_profile = zeros(N)
    Two_profile = zeros(N)

    Q_profile = zeros(N)

    # Initial temperatures
    Ti_profile[1] = Ti_in
    To_profile[1] = To_in

    # Simulate each segment
    for i in 1:N

        Ti_next, To_next, Twi, Two, Qsegment =
            simulate_segment(
                Ti_profile[i],
                To_profile[i],
                hi,
                ho,
                ri,
                ro,
                ΔL,
                kw,
                Ci,
                Co
            )

        # Store next fluid temperatures
        Ti_profile[i + 1] = Ti_next
        To_profile[i + 1] = To_next

        # Store wall temperatures
        Twi_profile[i] = Twi
        Two_profile[i] = Two

        # Store heat transferred in this segment
        Q_profile[i] = Qsegment
    end

    return Ti_profile,
           To_profile,
           Twi_profile,
           Two_profile,
           Q_profile
end

function segment_positions(L, N)

    x = collect(range(0, L, length=N + 1))

    return x
end

function validate_simulation(
    Ti_profile,
    To_profile,
    Ci,
    Co,
    Q_profile
)

    # Total heat lost by inner fluid
    Q_inner = Ci * (
        Ti_profile[1] - Ti_profile[end]
    )

    # Total heat gained by outer fluid
    Q_outer = Co * (
        To_profile[end] - To_profile[1]
    )

    # Sum of heat transferred by all segments
    Q_segments = sum(Q_profile)

    # Numerical differences
    error_inner = abs(Q_inner - Q_segments)
    error_outer = abs(Q_outer - Q_segments)

    return Q_inner, Q_outer, Q_segments,
           error_inner, error_outer
end



function radial_temperature_profile(
    Twi,
    Qsegment,
    ri,
    ro,
    kw,
    ΔL,
    num_points #i.e 100
)

    # Radial positions through the tube wall
    r = collect(range(ri, ro, length=num_points))

    # Temperature at each radial position
    T_radial = zeros(num_points)

    for i in 1:num_points

        T_radial[i] =
            Twi -
            Qsegment /
            (2 * π * kw * ΔL) *
            log(r[i] / ri)

    end

    return r, T_radial
end

function segment_midpoints(L, N)

    ΔL = L / N

    x_mid = [
        (i - 0.5) * ΔL
        for i in 1:N
    ]

    return x_mid
end



# COUNTER CURRENT FLOW
function simulate_countercurrent_flow(
    Ti_in,
    To_in,
    hi,
    ho,
    ri,
    ro,
    L,
    kw,
    Ci,
    Co,
    N;
    tolerance = 1e-6,
    max_iterations = 100
)

    # Length of each heat exchanger segment
    ΔL = L / N

    # ------------------------------------------
    # Guess the cold-fluid outlet temperature
    # ------------------------------------------

    lower_guess = To_in
    upper_guess = Ti_in

    Tc_out_guess = (lower_guess + upper_guess) / 2

    # ------------------------------------------
    # Bisection iteration
    # ------------------------------------------

    for iteration in 1:max_iterations

        # Temperature profiles
        Ti_profile = zeros(N + 1)
        To_profile = zeros(N + 1)

        # Wall temperature profiles
        Twi_profile = zeros(N)
        Two_profile = zeros(N)

        # Heat transfer in each segment
        Q_profile = zeros(N)

        # ------------------------------------------
        # Boundary conditions
        # ------------------------------------------

        # Hot fluid enters at x = 0
        Ti_profile[1] = Ti_in

        # Cold fluid leaves at x = 0
        # This is initially a guessed value
        To_profile[1] = Tc_out_guess

        # ------------------------------------------
        # Solve segment by segment
        # ------------------------------------------

        for i in 1:N

            # Resistance of one segment
            Ri, Rwall, Ro = cylindrical_resistances(
                hi,
                ho,
                ri,
                ro,
                ΔL,
                kw
            )

            # Heat transfer and wall temperatures
            Qsegment, Twi, Two =
                heat_transfer_and_wall_temperatures(
                    Ti_profile[i],
                    To_profile[i],
                    Ri,
                    Rwall,
                    Ro
                )

            # Hot fluid flows from x = 0 to x = L
            Ti_profile[i + 1] =
                Ti_profile[i] - Qsegment / Ci

            # Cold fluid physically flows from x = L to x = 0.
            # Therefore, when moving numerically from x = 0 to L,
            # its temperature decreases.
            To_profile[i + 1] =
                To_profile[i] - Qsegment / Co

            # Store results
            Twi_profile[i] = Twi
            Two_profile[i] = Two
            Q_profile[i] = Qsegment
        end

        # ------------------------------------------
        # Check cold-fluid inlet boundary condition
        # ------------------------------------------

        calculated_cold_inlet = To_profile[end]

        error = calculated_cold_inlet - To_in

        # Convergence check
        if abs(error) < tolerance

            println(
                "Counter-current solver converged in ",
                iteration,
                " iterations"
            )

            return (
                Ti_profile,
                To_profile,
                Twi_profile,
                Two_profile,
                Q_profile
            )
        end

        # ------------------------------------------
        # Bisection update
        # ------------------------------------------

        if error > 0
            upper_guess = Tc_out_guess
        else
            lower_guess = Tc_out_guess
        end

        Tc_out_guess =
            (lower_guess + upper_guess) / 2
    end

    error("Counter-current solver did not converge")

end



#-----------------------------------------------------
#                  DESIGN MODE
#-----------------------------------------------------

function find_required_length(
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
    L_low = 0.01,
    L_high = 5.0,
    tolerance = 1e-4,
    max_iterations = 100
)

    # ------------------------------------------
    # FUNCTION TO CALCULATE OUTLET TEMPERATURE
    # ERROR FOR A GIVEN LENGTH
    # ------------------------------------------

    function design_error(L)

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

        f =
            Ti_profile[end] -
            T_target

        return f
    end


    # ------------------------------------------
    # PHYSICAL FEASIBILITY CHECK
    # ------------------------------------------

    T_equilibrium =
        (
            Ci * Ti +
            Co * To
        ) /
        (
            Ci + Co
        )

    if T_target <= T_equilibrium

        error(
            "Target temperature is physically impossible for parallel flow. " *
            "The limiting outlet temperature is " *
            string(T_equilibrium) *
            " °C."
        )

    end


    # ------------------------------------------
    # INITIAL FUNCTION VALUES
    # ------------------------------------------

    f_low = design_error(L_low)
    f_high = design_error(L_high)


    println()
    println("----- DESIGN MODE BRACKET -----")

    println("L_low = ", L_low, " m")
    println("f(L_low) = ", f_low)

    println()

    println("L_high = ", L_high, " m")
    println("f(L_high) = ", f_high)


    # ------------------------------------------
    # CHECK THAT ROOT IS BRACKETED
    # ------------------------------------------

    if f_low * f_high > 0

        error(
            "Target solution is not bracketed. " *
            "Choose a smaller L_low or larger L_high."
        )

    end


    # ------------------------------------------
    # BISECTION METHOD
    # ------------------------------------------

    for iteration in 1:max_iterations

        L_mid =
            (L_low + L_high) / 2

        f_mid =
            design_error(L_mid)


        # --------------------------------------
        # CONVERGENCE CHECK
        # --------------------------------------

        if abs(f_mid) < tolerance

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
                L_mid,   #  see L mid is given as input here
                kw,
                Ci,
                Co,
                N
            )

            return (
                L_mid,
                Ti_profile,
                To_profile,
                Twi_profile,
                Two_profile,
                Q_profile,
                iteration
            )

        end


        # --------------------------------------
        # UPDATE BRACKET
        # --------------------------------------

        if f_low * f_mid < 0

            L_high = L_mid
            f_high = f_mid

        else

            L_low = L_mid
            f_low = f_mid

        end

    end


    error("Design mode did not converge.")

end