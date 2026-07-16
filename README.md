# The Death Correction

An interactive Shiny app that shows why an open population needs the "death correction" to hold relational duration on target, a companion to the SISMID "Network Modeling for Epidemics" (NME) course, Module 7.

In a closed population an edge ends only one way: the two partners break up, at a hazard of `1 / duration` per step. Open the population and an edge gains a second way to end, a partner departing (dying), which adds a hazard of about `2 * d.rate`. These are competing risks. Layer departures on a dissolution model that never knew about them and realized durations come out short. The death correction lowers the breakup hazard just enough that the two hazards together sum back to `1 / duration`, restoring the duration, which is why it makes edges more persistent.

Drag the departure-rate slider and watch the uncorrected realized duration fall away from target while the corrected one holds. An analytic core gives the exact numbers instantly; a simulation panel shows the same result with run-to-run variability (the mean age of active partnerships over time, as in the Module 7 tutorial).

Available at: https://epimodel.github.io/death-correction
