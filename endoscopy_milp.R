#  D3M Assignment-1 | Task 1 — GI Endoscopy MILP Model
#  Queen's University Belfast — Queen's Business School
#  Student: Harsh Vikram Singh | 
#--------
##  PURPOSE OF THIS SCRIPT
#--------
#  This script solves a Mixed-Integer Linear Programme (MILP)
#  to optimise the weekly configuration and scheduling of
#  GI (Gastrointestinal) endoscopy rooms across a 26-week
#  planning horizon. Each room can be set to:
#    - Diagnostic mode  (colonoscopy, gastroscopy etc.)
#    - Therapeutic mode (longer intervention procedures)
#    - Unavailable      (switched off / not activated)
#
#  The model minimises total cost = allocation cost + setup cost,
#  subject to meeting patient demand, clinician capacity limits,
#  and room configuration feasibility constraints.
#
#  EXACT METHOD USED
#  ------
#  The MILP is solved exactly using the HiGHS solver via the
#  ROI (R Optimisation Infrastructure) interface. HiGHS uses
#  Branch-and-Bound to find and certify the optimal integer
#  solution.

#  SECTION 0 — PACKAGE INSTALLATION----
# Install all required packages.Run this only once if the packages are not already installed.
# install.packages(c(
#  "ompr", "ompr.roi", "ROI", "ROI.plugin.highs",
#  "ggplot2", "dplyr", "tidyr", "patchwork", "scales"
# ))

# SECTION 1 — LIBRARIES-----
# Load all libraries needed for modelling, solving the optimisation, and producing plots.
library(ompr)             # Used to build the MILP optimisation model.
library(ompr.roi)         # Connects the ompr model to optimisation solvers.
library(ROI)              # General optimisation infrastructure in R.
library(ROI.plugin.highs) # HiGHS solver used to compute the optimal solution.
library(dplyr)            # For data manipulation and summarising results.
library(tidyr)            # For reshaping data between wide and long formats.
library(ggplot2)          # For creating visualisations of results.
library(patchwork)        # Allows combining multiple plots into one figure.
library(scales)           # Helps format axis labels such as comma-separated numbers.

# SECTION 2 — OUTPUT DIRECTORY + SOLVER SETTINGS----
# Setting up a folder to save results and pointing to HiGHS as our solver —
# it's free, fast, and handles MILP problems well.
OUT_DIR <- "outputs"
dir.create(OUT_DIR, showWarnings = FALSE)
# Define solver settings used when solving the MILP model.
SOLVER       <- "highs"  # HiGHS is an exact MILP solver using branch-and-bound.
VERBOSE_SOLVER <- TRUE   # Prints solver progress and optimality gap to console.
# Define file paths for saving results and solver log.
result_file <- file.path(OUT_DIR, "result_stage2.rds")
log_file    <- file.path(OUT_DIR, "highs_log_stage2.txt")

# SECTION 3 — INPUT DATA----
# All the raw numbers from the problem brief — rooms, weeks, demand, and capacity.
W <- 26;      weeks <- 1:W   # 26-week planning horizon
R <- 10;      rooms <- 1:R   # 10 endoscopy rooms
# Weekly patient demand vectors (length = W = 26)
P_diag  <- c(230,235,225,230,240,245,280,295,300,285,
             240,235,230,225,230,235,240,245,250,245,
             240,235,230,230,235,230)
P_thera <- c(145,145,150,150,155,160,168,176,182,188,
             190,192,188,195,187,182,178,172,168,165,
             162,158,155,152,150,148)
# Weekly clinician capacity (maximum total procedures allowed)
C <- c(520,515,525,520,530,525,490,480,485,500,
       510,515,520,525,530,525,520,480,470,475,
       500,510,515,520,525,530)
stopifnot(length(P_diag)==W, length(P_thera)==W, length(C)==W)
# Room parameters — grouped by size (3 small, 4 medium, 3 large)
B_diag  <- c(rep(60,3),  rep(120,4), rep(180,3))
B_thera <- c(rep(30,3),  rep(60,4),  rep(120,3))
S       <- c(rep(15000,3), rep(20000,4), rep(25000,3))
A       <- c(rep(180,3),   rep(200,4),   rep(220,3))

# SECTION 4 — BUILD THE MILP MODEL----
# This is the core of the assignment — building the MILP using ompr's pipe syntax.
# I found the |> chain approach cleaner than nesting, but took a while to get used to.
model_stage2 <- MIPModel() |>
  
  # DECISION VARIABLES
  # xD and xT are binary switches for room mode; hD and hT track actual hours scheduled.
  # y is the setup trigger — this is what drives the one-time activation costs.
  # xD[r,w]: binary — 1 if room r in DIAGNOSTIC mode in week w
  add_variable(xD[r,w], r=rooms, w=weeks, type="binary") |>
  # xT[r,w]: binary — 1 if room r in THERAPEUTIC mode in week w
  add_variable(xT[r,w], r=rooms, w=weeks, type="binary") |>
  # hD[r,w]: continuous >= 0 — diagnostic slots scheduled in r, w
  add_variable(hD[r,w], r=rooms, w=weeks, lb=0) |>
  # hT[r,w]: continuous >= 0 — therapeutic slots scheduled in r, w
  add_variable(hT[r,w], r=rooms, w=weeks, lb=0) |>
  # y[r,w]: binary — 1 if room r incurs SETUP COST in week w
  add_variable(y[r,w],  r=rooms, w=weeks, type="binary") |>
  
  # CONSTRAINT 1: MODE EXCLUSIVITY
  # Room cannot be in both modes simultaneously.
  add_constraint(xD[r,w] + xT[r,w] <= 1,
                 r=rooms, w=weeks, .show_progress_bar=FALSE) |>
  
  # CONSTRAINT 2: DIAGNOSTIC DEMAND (hard)
  # Total diagnostic slots across all rooms >= demand each week.
  add_constraint(sum_over(hD[r,w], r=rooms) >= P_diag[w],
                 w=weeks, .show_progress_bar=FALSE) |>
  
  # CONSTRAINT 3: THERAPEUTIC DEMAND (hard)
  # Total therapeutic slots across all rooms >= demand each week.
  add_constraint(sum_over(hT[r,w], r=rooms) >= P_thera[w],
                 w=weeks, .show_progress_bar=FALSE) |>
  
  # CONSTRAINT 4: THERAPEUTIC MODE RESTRICTION
  # Therapeutic slots only if room is in therapeutic mode.
  add_constraint(hT[r,w] <= B_thera[r]*xT[r,w],
                 r=rooms, w=weeks, .show_progress_bar=FALSE) |>
  
  # CONSTRAINT 5: ROOM CAPACITY
  # Total slots in room r ≤ configured capacity.
  add_constraint(hD[r,w] + hT[r,w] <= B_diag[r]*xD[r,w] + B_thera[r]*xT[r,w],
                 r=rooms, w=weeks, .show_progress_bar=FALSE) |>
  
  # CONSTRAINT 6: CLINICIAN CAPACITY
  # Total slots across ALL rooms ≤ available clinician hours.
  add_constraint(sum_over(hD[r,w]+hT[r,w], r=rooms) <= C[w],
                 w=weeks, .show_progress_bar=FALSE) |>
  
  # CONSTRAINT 7a: SETUP TRIGGER — WEEK 1
  # All rooms assumed unavailable at week 0.
  add_constraint(y[r,1] >= xD[r,1]+xT[r,1],
                 r=rooms, .show_progress_bar=FALSE) |>
  
  # CONSTRAINT 7b: SETUP TRIGGER — WEEKS 2–26
  # Setup cost triggered only when transitioning from off - on.
  add_constraint(y[r,w] >= (xD[r,w]+xT[r,w]) - (xD[r,w-1]+xT[r,w-1]),
                 r=rooms, w=2:W, .show_progress_bar=FALSE) |>
  
  # OBJECTIVE: MINIMISE TOTAL COST
  # Total cost = Allocation cost + Setup cost
  # Allocation: A[r] * capacity_activated (running cost per week)
  # Setup:      S[r] * y[r,w]            (one-time activation cost)
  set_objective(
    sum_over(A[r]*(B_diag[r]*xD[r,w] + B_thera[r]*xT[r,w]), r=rooms, w=weeks) +
      sum_over(S[r]*y[r,w], r=rooms, w=weeks),
    "min"
  )

# SECTION 5 — SOLVE THE MILP + DISPLAY KEY RESULTS----
# Solve the optimisation model using the HiGHS solver.
if (file.exists(result_file)) {
  result_stage2 <- readRDS(result_file)   # # Load previously saved optimal solution
} else {
  captured <- capture.output({
    result_stage2 <- solve_model(model_stage2,
                                 with_ROI(solver=SOLVER, verbose=VERBOSE_SOLVER))
  })
  writeLines(captured, con=log_file)      # Save solver output log for documentation of optimisation performance
  saveRDS(result_stage2, result_file)     # Save the solved model result so it can be reused later
}

#  Print key solver information so we can confirm the optimisation succeeded.
cat("\n======\n")
cat("SOLVER RESULTS\n")
cat("=======\n")

# Displays whether the solver successfully found an optimal solution
cat("Solver status:   ", result_stage2 |> solver_status(), "\n")

# Extract the minimum total cost obtained from the optimisation
obj_value <- result_stage2 |> objective_value()
# Print the optimal cost with readable formatting
cat("Optimal total cost: ", format(obj_value, big.mark=","), "\n")
# Stop execution if the solver did not reach an optimal solution
stopifnot(result_stage2 |> solver_status() == "success")

# SECTION 6 — EXTRACT SOLUTION TABLES----
# Extract the values of decision variables from the solved optimisation model.
xD_sol <- get_solution(result_stage2, xD[r,w]) |> rename(xD=value)
xT_sol <- get_solution(result_stage2, xT[r,w]) |> rename(xT=value)
hD_sol <- get_solution(result_stage2, hD[r,w]) |> rename(hD=value)
hT_sol <- get_solution(result_stage2, hT[r,w]) |> rename(hT=value)
y_sol  <- get_solution(result_stage2, y[r,w])  |> rename(setup=value)
# Merge all solution variables into one master table for easier analysis
sol <- xD_sol |>
  left_join(xT_sol, by=c("r","w")) |>
  left_join(hD_sol, by=c("r","w")) |>
  left_join(hT_sol, by=c("r","w")) |>
# Convert binary decisions into a readable room mode label
  mutate(mode = case_when(
    xD > 0.5 ~ "Diagnostic",
    xT > 0.5 ~ "Therapeutic",
    TRUE      ~ "Unavailable"
  ))

# Count how many weeks each room operates in each configuration.
cat("\n======\n")
cat("ROOM MODE FREQUENCY (how many weeks each room is in each mode)\n")
cat("=======\n")
mode_counts <- sol |> count(r, mode)
# Print summary table showing diagnostic/therapeutic usage of each room
print(mode_counts)  

# SECTION 7 — FEASIBILITY VALIDATION----
# Rolling up all the weekly totals and checking every constraint is satisfied —
# this is my proof that the solution is actually valid, not just "solver says so".
weekly_check <- sol |>
  group_by(w) |>
  summarise(
    diag_done  = sum(hD), # Total diagnostic procedures delivered
    thera_done = sum(hT), # Total therapeutic procedures delivered
    total_done = sum(hD+hT), # Total procedures performed
    .groups = "drop"
  ) |>
# Compare scheduled procedures with required demand and capacity
  mutate(
    diag_req    = P_diag[w],
    thera_req   = P_thera[w],
    clin_cap    = C[w],
    diag_slack  = diag_done  - diag_req,
    thera_slack = thera_done - thera_req,
    clin_slack  = clin_cap   - total_done
  )

# Using a tiny tolerance of 1e-6 to handle floating point rounding from the solver.
stopifnot(all(weekly_check$diag_slack  >= -1e-6))
stopifnot(all(weekly_check$thera_slack >= -1e-6))
stopifnot(all(weekly_check$clin_slack  >= -1e-6))

# Print full weekly validation table so we can inspect demand satisfaction.
cat("\n-----\n")
cat("WEEKLY FEASIBILITY CHECK (all slacks should be ≥ 0)\n")
cat("-------\n")
print(weekly_check)

# PRINT: SLACK SUMMARY STATISTICS
cat("\n--- Slack summary (min should be ≥ 0 for feasibility) ---\n")
print(summary(weekly_check[, c("diag_slack","thera_slack","clin_slack")]))

# PRINT: ALL/ANY DEMAND MET?
cat("\nAll diagnostic demand met:    ", all(weekly_check$diag_slack  >= -1e-6), "\n")
cat("All therapeutic demand met:   ", all(weekly_check$thera_slack >= -1e-6), "\n")
cat("All clinician limits respected:", all(weekly_check$clin_slack  >= -1e-6), "\n")

# Listing every setup event — useful for spotting if any room is being turned on/off
# too frequently, which would inflate costs unnecessarily.
setup_table <- y_sol |>
  filter(setup > 0.5) |>
  arrange(w, r) |>
  mutate(setup_cost = S[r]) |>
  select(w, r, setup_cost)

# Print all setup events to understand when rooms were switched on
cat("\n------\n")
cat("SETUP EVENTS (weeks where a room was activated from idle)\n")
cat("--------\n")
print(setup_table)

# Setup count per room
setup_summary <- y_sol |>
  filter(setup > 0.5) |>
  count(r, name="num_setups")

cat("\n--- Number of setup events per room ---\n")
print(setup_summary)

# SECTION 8 — COST BREAKDOWN----
# Computing allocation and setup costs separately each week so I can
# verify they add up to the solver's objective value — good reconciliation check.
weekly_cost <- sol |>
  group_by(w) |>
  summarise(
    alloc_cost = sum(A[r]*(B_diag[r]*xD + B_thera[r]*xT)),
    .groups = "drop"
  )
# Compute setup costs incurred when rooms are activated.
weekly_setup_cost <- y_sol |>
  group_by(w) |>
  summarise(setup_cost = sum(S[r]*setup), .groups="drop")
# Combine allocation cost and setup cost to compute total weekly spending.
cost_ts <- weekly_cost |>
  left_join(weekly_setup_cost, by="w") |>
  mutate(
    setup_cost   = ifelse(is.na(setup_cost), 0, setup_cost),
    total_weekly = alloc_cost + setup_cost,
    cumulative   = cumsum(total_weekly)
  )

# Verify that the calculated cost matches the optimisation objective value.
diff_obj <- sum(cost_ts$total_weekly) - obj_value

# Print a summary of the financial cost of the optimisation plan.
cat("\n======\n")
cat("COST SUMMARY\n")
cat("========\n")

total_setup_cost   <- sum(setup_table$setup_cost)
alloc_cost_implied <- obj_value - total_setup_cost

cat("Total cost (objective):   ", format(obj_value,          big.mark=","), "\n")
cat("Total allocation cost:    ", format(alloc_cost_implied,  big.mark=","), "\n")
cat("Total setup cost:         ", format(total_setup_cost,    big.mark=","), "\n")
cat("Cost reconciliation diff: ",  format(diff_obj,            scientific=FALSE), "\n")
# diff_obj should be 0 or extremely small (floating point)

# PRINT: WEEKLY COST TABLE
cat("\n--- Weekly cost breakdown (first 10 weeks shown) ---\n")
print(head(cost_ts, 10))   # Print first 10 rows as a preview

cat("\n--- Full 26-week cumulative cost ---\n")
print(cost_ts[, c("w","total_weekly","cumulative")])

# ── PRINT: ROOMS 8-9
# Checking rooms 8 and 9 specifically — these large rooms are expensive and
# I wanted to see if the solver was using them at all.
cat("\n--- Configuration of rooms 8 and 9 across all weeks ---\n")
print(
  sol |>
    filter(r %in% c(8, 9)) |>
    select(r, w, mode) |>
    arrange(r, w)
)

# SECTION 9 — PLOTS----
# PLOT 1: ROOM CONFIGURATION HEATMAP
# BUSINESS QUESTION: "What mode is each room in, across every week?"
p1 <- ggplot(sol, aes(x = w, y = factor(r), fill = mode)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_manual(values = c(
    "Diagnostic"  = "#4c78a8",   # Blue  = diagnostic mode
    "Therapeutic" = "#54a24b",   # Green = therapeutic mode
    "Unavailable" = "#e0e0e0"    # Grey  = room switched off
  )) +
  scale_x_continuous(breaks = seq(1, 26, by = 2)) +
  labs(
    title    = "Plot 1: Room Configuration Plan by Week",
    subtitle = "Blue = Diagnostic | Green = Therapeutic | Grey = Unavailable (idle capital)",
    x        = "Week",
    y        = "Room",
    fill     = "Mode"
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10, color = "grey40"),
    legend.position = "top"
  )

p1   # Explicit call required to display in console/viewer
cat("Plot 1 — Room Configuration Heatmap: generated\n")

# PLOT 2: ROOM ACTIVATION GANTT CHART
# BUSINESS QUESTION: "When exactly does each room switch on/off,
#                     and how long does each activation block last?"

# I wanted a Gantt-style view to see continuous blocks clearly — had to use
# cumsum() with lag() to detect mode changes, which Claude helped me figure out.
gantt_data <- sol |>
  arrange(r, w) |>
  group_by(r) |>
  mutate(
    # block_id increments every time the mode changes within a room
    block_id = cumsum(mode != lag(mode, default = first(mode))) + 1
  ) |>
  group_by(r, block_id, mode) |>
  summarise(
    week_start = min(w),
    week_end   = max(w) + 0.9,  # +0.9 fills the visual cell gap neatly
    .groups    = "drop"
  )

p2 <- ggplot(gantt_data, aes(
  xmin = week_start,
  xmax = week_end,
  ymin = r - 0.4,
  ymax = r + 0.4,
  fill = mode
)) +
  geom_rect(color = "white", linewidth = 0.4) +
  scale_fill_manual(values = c(
    "Diagnostic"  = "#4c78a8",
    "Therapeutic" = "#54a24b",
    "Unavailable" = "#e0e0e0"
  )) +
  scale_x_continuous(breaks = seq(1, 26, by = 2), limits = c(1, 27)) +
  scale_y_continuous(breaks = 1:R, labels = paste("Room", 1:R)) +
  labs(
    title    = "Plot 2: Room Activation Gantt Chart",
    subtitle = "Each bar = continuous block in that mode | Grey-to-colour transition = setup cost incurred",
    x        = "Week",
    y        = "",
    fill     = "Mode"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, color = "grey40"),
    legend.position  = "top",
    panel.grid.minor = element_blank()
  )

p2
cat("Plot 2 — Room Activation Gantt Chart: generated\n")

# PLOT 3: VALIDATION — SCHEDULED vs DEMAND + CLINICIAN
# BUSINESS QUESTION: "Does the plan meet every patient demand
#                     without exceeding what staff can deliver?"

# # Two-panel plot: Panel A checks demand is met, Panel B checks the clinician
# ceiling is never breached — together these are my visual proof of feasibility.

# Panel A: reshaping to long format so I can plot all four series in one ggplot call.
df_panel1 <- weekly_check |>
  select(w, diag_done, thera_done, diag_req, thera_req) |>
  pivot_longer(
    cols      = c(diag_done, thera_done, diag_req, thera_req),
    names_to  = "series",
    values_to = "hours"
  ) |>
  mutate(
    category = case_when(             # Category from column name prefix
      grepl("^diag",  series) ~ "Diagnostic",
      grepl("^thera", series) ~ "Therapeutic",
      TRUE ~ "Other"
    ),
    measure = case_when(              # Measure from column name suffix
      grepl("_done$", series) ~ "Scheduled",
      grepl("_req$",  series) ~ "Demand",
      TRUE ~ "Other"
    )
  )

p3a <- ggplot(df_panel1,
              aes(x = w, y = hours, color = category, linetype = measure)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.4) +
  scale_linetype_manual(values = c(
    "Demand"    = "dashed",    # Dashed = requirement (what patients need)
    "Scheduled" = "solid"      # Solid  = delivery   (what MILP scheduled)
  )) +
  scale_color_manual(values = c(
    "Diagnostic"  = "#4c78a8",
    "Therapeutic" = "#54a24b"
  )) +
  scale_x_continuous(breaks = seq(1, 26, by = 2)) +
  labs(
    title    = "Panel A: Scheduled vs Demand by Procedure Type",
    subtitle = "Solid line ≥ Dashed line every week = all patient demand met = plan is feasible",
    x        = "Week",
    y        = "Hours",
    color    = "Procedure type",
    linetype = "Series"
  ) +
  theme_minimal() +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

# Panel B: the shaded ribbon between scheduled and capacity is a clean way to show
# how much workforce slack exists — more ribbon means more unused clinician time.
df_panel2 <- weekly_check |>
  mutate(
    slack_lower = total_done,   # Bottom of ribbon = hours actually scheduled
    slack_upper = clin_cap      # Top of ribbon    = maximum available hours
  )

p3b <- ggplot(df_panel2, aes(x = w)) +
  geom_ribbon(aes(ymin = slack_lower, ymax = slack_upper),
              fill = "grey70", alpha = 0.35) +  # Shaded slack area
  geom_line(aes(y = clin_cap,   color = "Clinician capacity limit"),
            linewidth = 1.2, linetype = "dashed") +
  geom_line(aes(y = total_done, color = "Total scheduled"),
            linewidth = 1.1) +
  scale_color_manual(values = c(
    "Clinician capacity limit" = "#d63031",   # Red = hard limit
    "Total scheduled"          = "#0984e3"    # Blue = actual usage
  )) +
  scale_x_continuous(breaks = seq(1, 26, by = 2)) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    title    = "Panel B: Total Scheduled vs Clinician Availability",
    subtitle = "Grey ribbon = unused workforce slack | Red line = hard ceiling that must never be breached",
    x        = "Week",
    y        = "Total hours",
    color    = ""
  ) +
  theme_minimal() +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

# Stack Panel A above Panel B using patchwork operator /
p3 <- p3a / p3b +
  plot_annotation(
    title = "Plot 3: Feasibility Validation — Demand Met & Clinician Limits Respected",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

p3
cat("Plot 3 — Feasibility Validation: generated\n")

# PLOT 4: WEEKLY AVAILABLE BAY-HOURS BY CONFIGURATION
# BUSINESS QUESTION: "How much capacity did we activate to meet
#                     demand, and is the diagnostic/therapeutic
#                     split appropriate?"

# Stacked bar showing how much total capacity was switched on each week
# and how it's split between diagnostic and therapeutic modes.
capacity_week <- sol |>
  group_by(w) |>
  summarise(
    cap_diag  = sum(B_diag[r]  * xD),   # Total diagnostic bay-hours
    cap_thera = sum(B_thera[r] * xT),   # Total therapeutic bay-hours
    .groups   = "drop"
  ) |>
  pivot_longer(
    cols      = c(cap_diag, cap_thera),
    names_to  = "cap_type",
    values_to = "bay_hours"
  ) |>
  mutate(cap_type = recode(cap_type,
                           cap_diag  = "Diagnostic-configured bay-hours",
                           cap_thera = "Therapeutic-configured bay-hours"
  ))

p4 <- ggplot(capacity_week, aes(x = w, y = bay_hours, fill = cap_type)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c(
    "Diagnostic-configured bay-hours"  = "#4c78a8",
    "Therapeutic-configured bay-hours" = "#54a24b"
  )) +
  scale_x_continuous(breaks = seq(1, 26, by = 2)) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    title    = "Plot 4: Weekly Available Bay-Hours by Configuration",
    subtitle = "Total capacity activated by MILP each week — compare bar heights with Plot 3 demand lines",
    x        = "Week",
    y        = "Bay-hours configured",
    fill     = ""
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10, color = "grey40"),
    legend.position = "top"
  )

p4
cat("Plot 4 — Weekly Available Bay-Hours: generated\n")

# PLOT 5: PLANNED HOURS PER ROOM PER WEEK
# BUSINESS QUESTION: "Which rooms are doing real work week by week,
#                     and which are sitting half-empty?"

# Faceted bar chart — one panel per room — so you can quickly see which rooms are
# carrying the load and which are barely contributing across the 26 weeks.
hours_long <- sol |>
  select(r, w, hD, hT) |>
  pivot_longer(
    cols      = c(hD, hT),
    names_to  = "type",
    values_to = "hours"
  ) |>
  mutate(type = recode(type,
                       hD = "Diagnostic hours",
                       hT = "Therapeutic hours"
  ))

p5 <- ggplot(hours_long, aes(x = w, y = hours, fill = type)) +
  geom_col(width = 0.85) +
  facet_wrap(
    ~ r, ncol = 2,                   # One panel per room, 2-column layout
    labeller = labeller(r = function(x) paste("Room", x))  # "Room 1" labels
  ) +
  scale_fill_manual(values = c(
    "Diagnostic hours"  = "#4c78a8",
    "Therapeutic hours" = "#54a24b"
  )) +
  scale_x_continuous(breaks = seq(1, 26, by = 5)) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    title    = "Plot 5: Planned Procedure Hours per Room per Week",
    subtitle = "Small bars in active rooms = low utilisation | Tall bars = bottleneck risk",
    x        = "Week",
    y        = "Hours scheduled",
    fill     = "Procedure type"
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10, color = "grey40"),
    legend.position = "top",
    strip.text      = element_text(face = "bold")   # Bold room panel labels
  )

p5
cat("Plot 5 — Planned Hours per Room per Week: generated\n")

# PLOT 6 — ROOM UTILISATION RATE HEATMAP
# BUSINESS QUESTION: "Which rooms are being underutilised?"

# Computing utilisation as (actual / max possible) per room per week — grey tiles
# show idle rooms, and dark blue shows rooms running near full capacity.
util_data <- sol |>
  mutate(
    # Maximum possible hours depends on the active mode
    max_hours = case_when(
      xD > 0.5 ~ B_diag[r],   # Diagnostic mode capacity
      xT > 0.5 ~ B_thera[r],  # Therapeutic mode capacity
      TRUE      ~ NA_real_     # Unavailable — no capacity
    ),
    # Actual hours scheduled as a percentage of maximum
    util_rate = case_when(
      is.na(max_hours) ~ NA_real_,          # Room unavailable
      max_hours == 0   ~ NA_real_,          # Avoid divide-by-zero
      TRUE             ~ (hD + hT) / max_hours * 100  # % utilisation
    )
  )

p6 <- ggplot(util_data, aes(x = w, y = factor(r), fill = util_rate)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low      = "white",
    mid      = "#74b9ff",
    high     = "#0036a3",
    midpoint = 50,
    na.value = "grey90",      # Grey = room unavailable that week
    limits   = c(0, 100),
    name     = "Utilisation\nRate (%)"
  ) +
  scale_x_continuous(breaks = seq(1, 26, by = 2)) +
  labs(
    title    = "Plot 6: Room Utilisation Rate Heatmap",
    subtitle = "Grey = room unavailable | Dark blue = fully utilised | White = idle",
    x        = "Week",
    y        = "Room"
  ) +
  theme_minimal() +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    legend.position = "right"
  )

p6
cat("Plot 6 — Room Utilisation Rate: generated\n")

# PLOT 7 — ROOM CONTRIBUTION ANALYSIS (STACKED BAR)
# BUSINESS QUESTION: "Which rooms are doing the heavy lifting?"

# Ordering rooms by total contribution makes it immediately obvious which rooms
# are doing the heavy lifting and which are almost never used.
room_contribution <- sol |>
  group_by(r) |>
  summarise(
    total_diag  = sum(hD),
    total_thera = sum(hT),
    .groups     = "drop"
  ) |>
  pivot_longer(
    cols      = c(total_diag, total_thera),
    names_to  = "type",
    values_to = "total_hours"
  ) |>
  mutate(type = recode(type,
                       total_diag  = "Diagnostic",
                       total_thera = "Therapeutic"
  ))

# Compute total per room for ordering
room_totals <- room_contribution |>
  group_by(r) |>
  summarise(grand_total = sum(total_hours)) |>
  arrange(desc(grand_total))

p7 <- ggplot(room_contribution,
             aes(x = reorder(factor(r), -total_hours, sum),
                 y = total_hours, fill = type)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Diagnostic" = "#4c78a8", "Therapeutic" = "#54a24b")) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    title    = "Plot 7: Total 26-Week Procedure Contribution by Room",
    subtitle = "Rooms with near-zero bars are under-utilised assets",
    x        = "Room (ordered by total contribution)",
    y        = "Total procedures scheduled (26 weeks)",
    fill     = "Type"
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10, color = "grey40"),
    legend.position = "top"
  )

p7
cat("Plot 7 — Room Contribution Analysis: generated\n")

# PLOT 8 — DIAGNOSTIC vs THERAPEUTIC DEMAND BALANCE
# BUSINESS QUESTION: "Is our capacity mix aligned with demand mix?"

# Comparing demand mix with capacity mix as percentages — a gap between the lines
# means we allocated proportionally more than what patients actually needed.
balance_data <- weekly_check |>
  left_join(
    sol |>
      group_by(w) |>
      summarise(
        cap_diag  = sum(B_diag[r]  * xD),
        cap_thera = sum(B_thera[r] * xT),
        .groups   = "drop"
      ),
    by = "w"
  ) |>
  mutate(
    demand_diag_pct  = diag_req  / (diag_req + thera_req) * 100,
    demand_thera_pct = thera_req / (diag_req + thera_req) * 100,
    cap_total        = cap_diag + cap_thera,
    cap_diag_pct     = ifelse(cap_total > 0, cap_diag  / cap_total * 100, NA),
    cap_thera_pct    = ifelse(cap_total > 0, cap_thera / cap_total * 100, NA)
  )

p8 <- ggplot(balance_data, aes(x = w)) +
  geom_line(aes(y = demand_diag_pct, color = "Diagnostic demand %"),   linewidth = 1.1) +
  geom_line(aes(y = cap_diag_pct,    color = "Diagnostic capacity %"),  linewidth = 1.1, linetype = "dashed") +
  geom_ribbon(aes(ymin = demand_diag_pct, ymax = cap_diag_pct),
              fill = "#4c78a8", alpha = 0.12) +
  scale_color_manual(values = c(
    "Diagnostic demand %"   = "#4c78a8",
    "Diagnostic capacity %" = "#74b9ff"
  )) +
  scale_x_continuous(breaks = seq(1, 26, by = 2)) +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Plot 8: Diagnostic Share — Demand vs Capacity Mix",
    subtitle = "Solid = actual demand mix | Dashed = capacity mix allocated by MILP | Shaded = mismatch",
    x        = "Week",
    y        = "Diagnostic share (% of total)",
    color    = ""
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10, color = "grey40"),
    legend.position = "top"
  )

p8
cat("Plot 8 — Demand vs Capacity Balance: generated\n")

# PLOT 9: COST BREAKDOWN — WEEKLY + CUMULATIVE
# BUSINESS QUESTION: "What is the total financial cost of this
#                     plan, and where are the cost spikes coming from?"

# Two-panel cost chart: Panel A shows the week-by-week split, Panel B shows the
# running total climbing toward the final optimal value — good for budget reporting.
cost_long <- cost_ts |>
  select(w, alloc_cost, setup_cost, cumulative) |>
  pivot_longer(
    cols      = c(alloc_cost, setup_cost),
    names_to  = "cost_type",
    values_to = "weekly_cost"
  ) |>
  mutate(cost_type = recode(cost_type,
                            alloc_cost = "Allocation cost",   # Recurring weekly running cost
                            setup_cost = "Setup cost"         # One-time room activation cost
  ))

# Panel A: stacked bar — allocation + setup per week
p9a <- ggplot(cost_long, aes(x = w, y = weekly_cost, fill = cost_type)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c(
    "Allocation cost" = "#4c78a8",   # Blue   = steady running cost
    "Setup cost"      = "#f58518"    # Orange = one-time setup spike
  )) +
  scale_y_continuous(labels = scales::label_comma(prefix = "£")) +
  scale_x_continuous(breaks = seq(1, 26, by = 2)) +
  labs(
    title    = "Panel A: Weekly Cost Breakdown",
    subtitle = "Orange = room activation weeks (setup cost) | Blue = allocation (running cost)",
    x        = "Week",
    y        = "Weekly cost (£)",
    fill     = ""
  ) +
  theme_minimal() +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))

# Annotating the final total directly on the chart saves the reader from
# having to look it up in the console output separately.
p9b <- ggplot(cost_ts, aes(x = w, y = cumulative)) +
  geom_line(color = "#2d3436", linewidth = 1.1) +
  geom_point(color = "#2d3436", size = 1.6) +
  # Annotate the final total cost on the chart
  geom_text(
    data  = cost_ts |> filter(w == max(w)),
    aes(label = paste0("Total: £", format(round(cumulative), big.mark = ","))),
    vjust = -1, hjust = 1, size = 3.5, fontface = "bold", color = "#d63031"
  ) +
  scale_y_continuous(labels = scales::label_comma(prefix = "£")) +
  scale_x_continuous(breaks = seq(1, 26, by = 2)) +
  labs(
    title    = "Panel B: Cumulative Cost — Budget Commitment Curve",
    subtitle = "Final annotated value = MILP-certified minimum total cost across all 26 weeks",
    x        = "Week",
    y        = "Cumulative cost (£)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

# Stack Panel A above Panel B
p9 <- p9a / p9b +
  plot_annotation(
    title = "Plot 9: Cost Breakdown — Weekly Components + Cumulative Commitment",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

p9
cat("Plot 9 — Cost Breakdown (Weekly + Cumulative): generated\n")

# PLOT 10 — COST EFFICIENCY: COST PER PROCEDURE OVER TIME
# BUSINESS QUESTION: "Is the model getting more or less efficient?"

# Flagging setup weeks in red makes spikes self-explanatory — they're not
# inefficiency, just the one-time cost of activating a room that week.
cost_per_proc <- cost_ts |>
  left_join(weekly_check |> select(w, total_done), by = "w") |>
  mutate(
    cost_per_proc = total_weekly / total_done,
    has_setup     = setup_cost > 0  # Flag weeks with setup events
  )

p10 <- ggplot(cost_per_proc, aes(x = w, y = cost_per_proc)) +
  # Highlight setup weeks with vertical bands
  geom_col(
    data = cost_per_proc |> filter(has_setup),
    aes(x = w, y = Inf, fill = "Setup event week"),
    width = 0.8, alpha = 0.15, inherit.aes = FALSE
  ) +
  geom_line(color = "#0984e3", linewidth = 1.1) +
  geom_point(aes(color = has_setup), size = 2.5) +
  scale_color_manual(
    values = c("FALSE" = "#0984e3", "TRUE" = "#d63031"),
    labels = c("FALSE" = "Normal week", "TRUE" = "Setup event week"),
    name   = ""
  ) +
  scale_fill_manual(values = c("Setup event week" = "#d63031"), name = "") +
  scale_x_continuous(breaks = seq(1, 26, by = 2)) +
  scale_y_continuous(labels = scales::label_comma(prefix = "£")) +
  labs(
    title    = "Plot 10: Cost per Procedure over Time",
    subtitle = "Red points = weeks with room setup events (cost spikes) | Blue = normal operation",
    x        = "Week",
    y        = "Cost per procedure (£)"
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10, color = "grey40"),
    legend.position = "top"
  )

p10
cat("Plot 10 — Cost per Procedure: generated\n")

# PLOT 11 — CAPACITY UTILISATION RIBBON (TIGHTNESS ANALYSIS)
# BUSINESS QUESTION: "How close are we to hitting our limits?"

# Three lines + two ribbons in one chart is a lot — Claude helped me layer
# the geom_ribbon and geom_line calls correctly so legends don't clash.
ribbon_data <- weekly_check |>
  mutate(total_demand = diag_req + thera_req)

p11 <- ggplot(ribbon_data, aes(x = w)) +
  # Ribbon between scheduled and clinician cap = workforce slack
  geom_ribbon(
    aes(ymin = total_done, ymax = clin_cap, fill = "Workforce slack"),
    alpha = 0.25
  ) +
  # Ribbon between demand and scheduled = scheduling buffer
  geom_ribbon(
    aes(ymin = total_demand, ymax = total_done, fill = "Scheduling buffer"),
    alpha = 0.25
  ) +
  geom_line(aes(y = clin_cap,     color = "Clinician capacity limit"), linewidth = 1.2, linetype = "dashed") +
  geom_line(aes(y = total_done,   color = "Total scheduled"),          linewidth = 1.1) +
  geom_line(aes(y = total_demand, color = "Total demand"),             linewidth = 1.1, linetype = "dotted") +
  scale_fill_manual(values  = c("Workforce slack" = "#fdcb6e", "Scheduling buffer" = "#55efc4")) +
  scale_color_manual(values = c(
    "Clinician capacity limit" = "#d63031",
    "Total scheduled"          = "#0984e3",
    "Total demand"             = "#6c5ce7"
  )) +
  scale_x_continuous(breaks = seq(1, 26, by = 2)) +
  scale_y_continuous(labels = scales::label_comma()) +
  labs(
    title    = "Plot 11: Capacity Tightness — Demand vs Scheduled vs Clinician Limit",
    subtitle = "Yellow ribbon = workforce slack | Green ribbon = scheduling buffer",
    x        = "Week",
    y        = "Total hours",
    color    = "",
    fill     = ""
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10, color = "grey40"),
    legend.position = "bottom"
  )

p11
cat("Plot 11 — Capacity Tightness Ribbon: generated\n")

# PLOT 12 — WEEKLY DEMAND PRESSURE INDEX
# BUSINESS QUESTION: "Which weeks are highest risk for failure?"

# Expressing demand as a % of clinician capacity gives a clean risk signal —
# anything above 85% needs a contingency plan in the real operational setting.
pressure_data <- weekly_check |>
  mutate(
    total_demand   = diag_req + thera_req,
    pressure_index = total_demand / clin_cap * 100,
    risk_level     = case_when(
      pressure_index > 85 ~ "High Risk (>85%)",
      pressure_index > 70 ~ "Moderate (70-85%)",
      TRUE                ~ "Comfortable (<70%)"
    ),
    risk_level = factor(risk_level,
                        levels = c("Comfortable (<70%)", "Moderate (70-85%)", "High Risk (>85%)"))
  )

p12 <- ggplot(pressure_data, aes(x = w, y = pressure_index, fill = risk_level)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = 85, color = "#d63031", linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = 70, color = "#e17055", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = 26.5, y = 86.5, label = "High risk", color = "#d63031", size = 3, hjust = 1) +
  annotate("text", x = 26.5, y = 71.5, label = "Moderate",  color = "#e17055", size = 3, hjust = 1) +
  scale_fill_manual(values = c(
    "Comfortable (<70%)"   = "#55efc4",
    "Moderate (70-85%)"    = "#fdcb6e",
    "High Risk (>85%)"     = "#ff7675"
  )) +
  scale_x_continuous(breaks = seq(1, 26, by = 2)) +
  scale_y_continuous(limits = c(0, 105), labels = function(x) paste0(x, "%")) +
  labs(
    title    = "Plot 12: Weekly Demand Pressure Index",
    subtitle = "Pressure = total demand / clinician capacity | Red weeks need contingency planning",
    x        = "Week",
    y        = "Demand Pressure (%)",
    fill     = "Risk Level"
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10, color = "grey40"),
    legend.position = "top"
  )

p12
cat("Plot 12 — Demand Pressure Index: generated\n")

# COMBINED EXECUTIVE DASHBOARD

# I originally had these three plots as separate outputs, but combining them into
# a single dashboard layout felt more appropriate for a management presentation.
# The layout — (p12 / p10) | p7 — was suggested and refined with Claude's help,
# as getting patchwork to balance panel sizes correctly took a few attempts.
executive_dashboard <- (p12 / p10) | p7
executive_dashboard <- executive_dashboard +
  plot_annotation(
    title   = "GI Endoscopy Unit — Executive Planning Dashboard",
    subtitle = "Demand pressure | Cost efficiency | Room contribution",
    theme   = theme(
      plot.title    = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 11, color = "grey40")
    )
  )

executive_dashboard
cat("Executive Dashboard: generated\n")

# SECTION 10 — SAVE ALL OUTPUTS----
# Save all generated plots and tables to the output directory
# so they can be included in the assignment report.
ggsave(
  file.path(OUT_DIR, "Plot1_RoomConfigHeatmap.png"),
  plot   = p1,
  width  = 11, height = 5, dpi = 300
)
ggsave(
  file.path(OUT_DIR, "Plot2_RoomActivationGantt.png"),
  plot   = p2,
  width  = 12, height = 5, dpi = 300
)
ggsave(
  file.path(OUT_DIR, "Plot3_FeasibilityValidation.png"),
  plot   = p3,
  width  = 11, height = 9, dpi = 300   # Taller — two stacked panels
)
ggsave(
  file.path(OUT_DIR, "Plot4_WeeklyBayHours.png"),
  plot   = p4,
  width  = 11, height = 5, dpi = 300
)
ggsave(
  file.path(OUT_DIR, "Plot5_PlannedHoursPerRoom.png"),
  plot   = p5,
  width  = 11, height = 10, dpi = 300  # Taller — 10 facet panels
)
ggsave(
  file.path(OUT_DIR, "Plot6_UtilisationRateHeatmap.png"),
  plot   = p6,
  width  = 11, height = 5, dpi = 300
)
ggsave(
  file.path(OUT_DIR, "Plot7_RoomContribution.png"),
  plot   = p7,
  width  = 10, height = 5, dpi = 300
)
ggsave(
  file.path(OUT_DIR, "Plot8_DemandCapacityBalance.png"),
  plot   = p8,
  width  = 11, height = 5, dpi = 300
)
ggsave(
  file.path(OUT_DIR, "Plot9_CostBreakdown.png"),
  plot   = p9,
  width  = 11, height = 9, dpi = 300   # Taller — two stacked panels
)
ggsave(
  file.path(OUT_DIR, "Plot10_CostPerProcedure.png"),
  plot   = p10,
  width  = 11, height = 5, dpi = 300
)
ggsave(
  file.path(OUT_DIR, "Plot11_CapacityTightness.png"),
  plot   = p11,
  width  = 11, height = 5, dpi = 300
)
ggsave(
  file.path(OUT_DIR, "Plot12_DemandPressureIndex.png"),
  plot   = p12,
  width  = 11, height = 5, dpi = 300
)
ggsave(
  file.path(OUT_DIR, "Executive_Dashboard.png"),
  plot   = executive_dashboard,
  width  = 16, height = 10, dpi = 300  # Wider — multi-panel layout
)

# Save feasibility validation table showing all weekly demand constraints satisfied.
write.csv(weekly_check,
          file.path(OUT_DIR, "weekly_check.csv"),
          row.names = FALSE)

# Save weekly and cumulative cost time series.
write.csv(cost_ts,
          file.path(OUT_DIR, "cost_ts.csv"),
          row.names = FALSE)

# Save all room setup events and associated costs.
write.csv(setup_table,
          file.path(OUT_DIR, "setup_table.csv"),
          row.names = FALSE)

# Final message confirming that all outputs and solver logs were generated successfully.
message("All 12 plots + dashboard saved to: ", normalizePath(OUT_DIR))
message("Solver log (MILP gap evidence):     ", normalizePath(log_file))

# END OF SCRIPT

