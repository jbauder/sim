# --- SIMULATION FUNCTIONS -----------------------------------------------------
# General functions for simulation construction and analysis of simulation.
# ------------------------------------------------------------------------------

# CompileAllAgentsStepData Function --------------------------------------------

### Creates a compiled dataframe of all the agent$step_data objects within 
###  sim$agents$all
### Usage: CompileAllAgentsStepData(sim)
### Arguments: sim = a 'sim' list object
### Returns: a dataframe of all the step_data objects 
### Notes: 
### Blake Massey
### 2015.03.29

CompileAllAgentsStepData <- function(sim = sim) {
  all_step_data <- data.frame()
  all <- sim$agents$all
  for (i in 1:length(all)) {
    if (exists("step_data", where=sim$agents$all[[i]])){ 
      all_step_data <- rbind(all_step_data, sim$agents$all[[i]]$step_data)
    }  
  }
  return(all_step_data)
}

# CreateRunsList Function ------------------------------------------------------

###  Helper function for RunSimulation(). Creates a list of empty objects, each 
###    one named individually as "run_xx" where the number of leading zeroes 
###    matches the maximum digit width of the number of runs. 
###  Usage: CreateRunsList(run)
###  Arguments: runs = numeric, number of runs in RunSimulation() 
###  Returns: a list  
###  Notes: Where there are < 10 runs, there will have no leading zeros, 
###    when there are between 10 and 99 runs, there will be a leading zero for 
###    runs < 10.
###  Blake Massey
###  2015.03.29

CreateRunsList <- function(runs = runs) {
  format <- paste0("%0", nchar(runs), "d")
  runs_list <- list()
  for (i in 1:runs) {
    runs_list[[i]] <- NA
    names(runs_list)[[i]] <- paste0("run_", sprintf(format, i))
  }
  return(runs_list)
}

# CreateStepIntervals Function -------------------------------------------------

### Create step intervals within a summary interval
### Usage: CreateStepIntervals(sum_interval = sum_interval,
###   step_period = step_period)
### Arguments: sum_interval = an Interval object from lubridate that 
###              represents the current summary interval
###            step_period = a Period object specifying the length of the step 
###              period. The step_period should be of equal or greater length 
###              than the time_step_period. Note that if sum_interval is
###              not evenly divisible by step_period the last step_interval is 
###              removed so that all step_intervals returned are of equal 
###              duration.
### Returns: a list of intervals
### Notes:
### Javan Bauder
### 2015.03.29

CreateStepIntervals <- function(sum_interval = sum_interval,
                                step_period = sim$pars$global$step_period) {
  suppressPackageStartupMessages(require(lubridate))
  step_period <- step_period
  step_intervals <- list()
  interval_counter <- 1
  current_start <- int_start(sum_interval)
  current_end <- (current_start + step_period)
  stop_point<-int_end(sum_interval)
  while(current_start < (stop_point)) {
    current_end <- (current_start+step_period)
    step_intervals[[interval_counter]] <- new_interval(current_start,
      current_end)
    interval_counter <- interval_counter + 1
    current_start <- current_start + step_period
  }
  if(int_end(step_intervals[[length(step_intervals)]]) > stop_point){
    step_intervals[[length(step_intervals)]] <- NULL
  }
  return(step_intervals)
}

# CreateSumIntervals Function --------------------------------------------------

###  Creates summary intervals for the simulation input
###  Usage: CreateSumIntervals(sim)
###  Arguments: sim = a list that contains a "pars" list which contains a 
###               "global" list that contains the following objects:  
###               sim_start = a POSIXct object for simulation's start, required
###               sim_period = a Period object (from 'lubridate' package), e.g.
###                 period(10, "year"), period(10, "week"), period(10, "week")
###                 that sets the length of time the simulation runs
###               sim_end = a POSIXct object for simulation's end. Setting this 
###                 will override the 'sim_period' parameter. Default = NULL.
###               sum_period = a Period object that sets the summary interval  
###               sum_interval = vector or dateframe (with date in first column) 
###                 of the format %B%d, %b%d, %d%B, or %B%d (see ?strptime for 
###                 more details on formatting). Example: c("April01", "May15", 
###                 "Sep1", "Oct15"). Overrides 'sum_period' parameter.
###               sum_interval_custom = a vector of POSIXct objects that set the 
###                 start of the summary periods. If sim_start precedes the
###                 first value, the first value will go until from start_sim to 
###                 sum_interval_custom. Last interval will go from the last 
###                 date in sum_interval_custom to end of sim_period or end_sim. 
###                 Overrides 'sum_period' or 'sum_interval' parameters.
###  Returns: list of intervals  
###  Notes: Either 'sim_period' or 'sim_end' must be set, and either 
###    'sum_period', 'sum_interval', or 'sum_interval_custom' must be set. For
###    'sum_interval' the dates can be set as numeric for both day and month, 
###    but there is a possibility that the day and month will be inverted 
###    because dates with both values <12 can be confused (e.g. 02-10 can be 
###    either February 10th or October 2nd. Date format: "10Feb" or "02Oct" is 
###    preferred.  
###  Blake Massey
###  2015.03.10

CreateSumIntervals <- function(sim = sim) {  
  suppressPackageStartupMessages(require(lubridate))
  sim_start <- sim$pars$global$sim_start
  sim_period <- sim$pars$global$sim_period
  sim_end <- sim$pars$global$sim_end
  sum_period <- sim$pars$global$sum_period
  sum_interval <- sim$pars$global$sum_interval
  sum_interval_custom <- sim$pars$global$sum_interval_custom
  if (is.null(sim_period) && is.null(sim_end)) {
    stop("must provide either 'sim_period' or 'sim_end' parameter value")
  }
  if (!is.null(sim_period))  sim_end <- sim_start + sim_period   
  if (!is.null(sum_period) && !is.null(sum_interval)){
    print("'sum_period' used instead of 'sum_interval'")
  }
  if (!is.null(sum_period)) {
    sum_period_unit <- ExtractUnitFromPeriod(sum_period) 
    first_int <- as.interval(sim_start, ceiling_date(sim_start+period(1, 
      "second"), unit=sum_period_unit))
    sum_intervals <- sim_start
    end_int <- int_end(first_int)    
    sum_intervals <- with_tz(append(sum_intervals, end_int), tz(sim_start))      
    while (end_int < sim_end) {  
      end_int <- end_int + sum_period
      sum_intervals <- with_tz(append(sum_intervals, end_int), tz(sim_start))   
    }  
    if (tail(sum_intervals, 1) > sim_end) {
      sum_intervals[length(sum_intervals)] <- sim_end
    }
    sum_interval <- NULL
  } 
  if (!is.null(sum_interval)){
      if (is.data.frame(sum_interval)) sum_interval <- sum_interval[, 1]
      if (max(names(guess_formats(sum_interval, c("md", "dm")))) == "md") {
        md_format <- max(guess_formats(sum_interval, c("md", "dm")))
        sum_interval <- as.character(as.Date(sum_interval, format = md_format),
          "%d%b")  # date must be in day month order
      }                      
      dates <- dmy(paste0(sum_interval, 2000)) # creates POSIXct    
      sum_interval  <- sum_interval[order(dates)]  # orders sum_interval dates 
      first_sum_int <- FindFirstSumInterval(sim_start, sum_interval)
      end_int <- int_end(first_sum_int) 
      sum_intervals <- with_tz(append(sim_start, end_int), tz(sim_start))
      while (end_int < sim_end) {  
        end_int_jul <- yday(end_int)
        sum_int_jul <- yday(dmy(paste0(sum_interval, year(end_int))))
        sum_int_position <- match(end_int_jul, sum_int_jul)
        if (sum_int_position == length(sum_interval)) {
          next_int <- new_interval(end_int, dmy(paste0(sum_interval[1], 
            (year(end_int)))) + period(1, "year"))
        } else {
          next_int <- new_interval(end_int, dmy(paste0(sum_interval
            [sum_int_position+1], year(end_int))))
        }
        end_int <- int_end(next_int)
        sum_intervals <- with_tz(append(sum_intervals, end_int), tz(sim_start))   
      }      
      if (tail(sum_intervals, 1) > sim_end) {
        sum_intervals[length(sum_intervals)] <- sim_end
      }
      sum_interval_custom <- NULL
  }
  if(!is.null(sum_interval_custom)){
    if (is.data.frame(sum_interval)) sum_interval <- sum_interval_custom[,1]
    sum_interval <- sum_interval_custom # needed if sum_interval_custom != df
    sum_intervals <- as.POSIXct(sum_interval, tz=tz(sim_start))
    if (tail(sum_intervals, 1) > sim_end && !is.null(sum_period)) {
      sum_intervals[length(sum_intervals)] <- sim_end
      warning("'sum_inverval_custom' exceeds 'sim_period' time length")
    }
    if (tail(sum_intervals, 1) > sim_end && is.null(sum_period)) {
      sum_intervals[length(sum_intervals)] <- sim_end
      warning("'sum_inverval_custom' exceeds 'sim_end' time length")
    }
    if (tail(sum_intervals, 1) < sim_end) {
      sum_intervals[length(sum_intervals)+1] <- sim_end
    }
  }
  sum_intervals_list <- list()
  for (i in 1:(length(sum_intervals)-1)) {
    interval <- as.interval(sum_intervals[i], sum_intervals[i+1]-1)
    sum_intervals_list[[length(sum_intervals_list)+1]] <- interval
  }
  return(sum_intervals_list)
}

# CreateTimeSteps Function -----------------------------------------------------

### Create the time steps within a step interval
### Usage: CreateTimeSteps(step_interval)
### Arguments: step_interval = an Interval object from lubridate that 
###              represents the current step interval
###            time_step_period = a Period object specifying the length of the 
###              time step period. Note that if the length  of the 
###              current_step_interval is equal to the time_step_period the 
###              function will return a POSIXct object with the same date as the 
###              end of the current_step_interval.
### Returns: a list of intervals
### Notes: 
### Javan Bauder
### 2015.03.29

CreateTimeSteps <- function(step_interval = step_interval,
                            time_step_period = 
                              sim$pars$global$time_step_period) {
  suppressPackageStartupMessages(require(lubridate))
  options(lubridate.verbose=FALSE)
  step_interval_end <- int_end(step_interval)  
  steps <- int_start(step_interval) # creates output list object w/start time 
  end_int <- int_end(as.interval(time_step_period, int_start(step_interval))) 
  steps <- with_tz(append(steps, end_int), tz(int_end(step_interval))) # 2nd t     
  while (end_int < step_interval_end) {  
    end_int <- end_int + time_step_period
    steps <- with_tz(append(steps, end_int), tz(int_end(step_interval)))   
  }  
  if (tail(steps, 1) > step_interval_end) {
    steps[length(steps)] <- step_interval_end
  }
  time_step_list <- list()
  for (i in 1:(length(steps)-1)) {
    interval <- as.interval(steps[i], steps[i+1])
    time_step_list[[length(time_step_list)+1]] <- interval
  }
  return(time_step_list)
}

# ExtractUnitFromPeriod Function -----------------------------------------------

###  Helper function for CreateSumInterval(). Finds unit of Period object. 
###  Usage: ExtractUnitFromPeriod(period)
###  Arguments: sim_start = a POSIXct object for simulation's start, required
###  Returns: period object's unit as a chararcter 
###  Notes:
###  Blake Massey
###  2015.03.10

ExtractUnitFromPeriod <- function(period) {
  period <- period
  if (period$year > 0) unit <- "year"
  if (period$month > 0) unit <- "month"
  if (period$day == 7 ) unit <- "week"  
  if (period$day > 0 && period@day < 7) unit <- "day"
  if (period$day > 7) unit <- "day"
  if (period$hour > 0) unit <- "hour"  
  if (period$minute > 0) unit <- "minute"
  if (period$.Data > 0) unit <- "second"
  return(unit)
}

# FindFirstSumInterval Function ------------------------------------------------

###  Helper function for CreateSumInterval. Finds first interval based on  
###    sim_start datetime and sum_interval dates
###  Usage: FindFirstSumInterval(sim_start, sum_interval)
###  Arguments: sim_start = a POSIXct object for simulation's start, required
###             sum_interval = vector or dateframe (with date in first column) 
###               of the format %B%d, %b%d, %d%B, or %B%d (see ?strptime for 
###               more details on formatting). Example: c("April01", "May15", 
###               "Sep1", "Oct15"). Required.
###  Returns: an Interval object
###  Notes:
###  Blake Massey
###  2015.03.10

FindFirstSumInterval <- function(sim_start, 
                                 sum_interval){
  require(lubridate)
  start_year <- as.Date(0, origin=as.Date(floor_date(sim_start, "year")))
  if (length(sum_interval) == 1) {
    first_int <- as.interval(period(1, "year"), dmy(paste0(sum_interval[1], 
      year(start_year))))
    if (sim_start %within% first_int) { 
      first_sum_int <- new_interval(sim_start, dmy(paste0(sum_interval[1], 
        year(start_year + period(1, "year")))))
    } else {
      first_sum_int <- new_interval(sim_start, dmy(paste0(sum_interval[1], 
        year(start_year))))    
    }
  }
  if (length(sum_interval) == 2) {
    first_int <- as.interval(dmy(paste0("0101", year(start_year))),
      dmy(paste0(sum_interval[1], year(start_year))))
    if ((sim_start + period(1, "second")) %within% first_int) { 
      first_sum_int <- new_interval(sim_start, dmy(paste0(sum_interval[1], 
        year(start_year)))) 
    } 
    if (exists("first_sum_int") == FALSE) {
      second_int <- as.interval(dmy(paste0(sum_interval[1], year(start_year))),
        dmy(paste0(sum_interval[2], year(start_year))))
      if ((sim_start + period(1, "second")) %within% second_int) { 
        first_sum_int <- new_interval(sim_start, dmy(paste0(sum_interval[2], 
          year(start_year))))  
      } else { 
        first_sum_int <- new_interval(sim_start, dmy(paste0(sum_interval[1], 
          year(start_year+period(1, "year")))))   
      }
    }
  }
  if (length(sum_interval) > 2) { 
    first_int <- as.interval(dmy(paste0("0101", year(start_year))),
        dmy(paste0(sum_interval[1], year(start_year))))
    if ((sim_start + period(1, "second")) %within% first_int) { 
      first_sum_int <- new_interval(sim_start, dmy(paste0(sum_interval[1], 
        year(start_year))))
    }
    for (i in 2:(length(sum_interval)-1)) { 
      mid_int <- as.interval(dmy(paste0(sum_interval[i], year(start_year))),
        dmy(paste0(sum_interval[i+1], year(start_year))))
      if ((sim_start + period(1, "second")) %within% mid_int) { 
        first_sum_int <- new_interval(sim_start, dmy(paste0(sum_interval[i+1], 
          year(start_year))))  
      }
    }
    for (i in length(sum_interval)) { 
      last_int <- as.interval(dmy(paste0(sum_interval[i], year(start_year))),
        dmy(paste0("0101", (year(start_year)+1) )))
      if ((sim_start + period(1, "second")) %within% last_int) { 
        first_sum_int <- new_interval(sim_start, dmy(paste0(sum_interval[1], 
          (year(start_year)+1))))   
      }
    }
  }
  return(first_sum_int)
}

# FindSeasonFromDatetime Function ----------------------------------------------

###  Determines which season (from sim$pars$global$sim_seasons) a datetime is 
###    within
###  Usage: FindSeasonFromDatetime(datetime, season)
###  Arguments: datetime = a POSIXct object
###             seasons = a dataframe of starting dates ("start" column) and 
###               season names ("season" column), usually located at: 
###               sim$pars$global$sim_seasons
###  Returns: an atomic character object
###  Notes: 
###  Blake Massey
###  2015.03.28

FindSeasonFromDatetime <- function(datetime = datetime,
                                   seasons = sim$pars$global$sim_seasons) {
  require(lubridate)
  start_year <- as.Date(0, origin=as.Date(floor_date(datetime, "year")))
  if (max(names(guess_formats(seasons[, "start"], c("md", "dm")))) == "md") {
    md_format <- max(guess_formats(seasons[, "start"], c("md", "dm")))
    seasons[, "start"] <- as.character(as.Date(seasons[, "start"], format = 
      md_format), "%d%b")  # date must be in day month order
  } 
  dates <- dmy(paste0(seasons[,"start"], 2000)) # creates POSIXct 
  seasons  <- seasons[with(seasons, order(dates)), ]  # orders seasons by dates 
  if (nrow(seasons) == 2) {
    first_season <- as.interval(dmy(paste0(paste0(seasons[1, "start"], 
      year(start_year)))), dmy(paste0(seasons[2, "start"], year(start_year)))-1)
    if (datetime %within% first_season) { 
      season <- seasons[1, "season"] 
    } else {  
      season <- seasons[2, "season"]  
    }
  }
  if (nrow(seasons) > 2) { 
    last_season <- TRUE
    for (i in 1:(nrow(seasons)-1)) { 
      mid_season <- as.interval(dmy(paste0(seasons[i, "start"], year(
        start_year))), dmy(paste0(seasons[i+1, "start"], year(start_year)))-1)
      if (datetime %within% mid_season) { 
        season <- seasons[i, "season"]
        last_season = FALSE
      }
    }
    if(last_season == TRUE) season <- seasons[nrow(seasons), "season"]
  }
  return(season)
}

# NamedList Function -----------------------------------------------------------

###  Creates a list with named objects, and it tries not to replace any already-
###    named arguments.
###  Usage: NamedList(...)
###  Arguments: ... = objects to add to list
###  Returns: a list  
###  Notes: Original function came from Ben Bolker's answer on StackOverflow:
###    http://stackoverflow.com/questions/16951080
###  Blake Massey
###  2015.03.24

NamedList <- function(...) {
    list <- list(...)
    str_name <- sapply(substitute(list(...)),deparse)[-1]
    if (is.null(name <- names(list))) name <- str_name
    if (any(no_names <- name == "")) name[no_names] <- str_name[no_names]
    setNames(list, name)
}

# ReproductionSubModel Function ------------------------------------------------

###  Reproduction submodel 
###  Usage: MovementSubModel()
###  Arguments: agent_states, a 'agent_states' list object
###  Returns: 'agent_states' list object
###  Notes: Just a placeholder for now. 
###  Blake Massey
###  2015.03.29

ReproductionSubModel <- function(agent_states = agent_states,
                                 step_data = step_data) {
  agent_states <- agent_states
  return(agent_states)
}

# ReturnAliveIds Function --------------------------------------------------------

### Create vector of ids of alive agents
### Usage: ReturnAliveIds(all)
### Arguments: all = 'all' list
### Returns:
### Notes:
### Blake Massey
### 2015.03.28

ReturnAliveIds <- function(all = sim$agents$all){
  all <- all
  alive_ids <- vector()
  for (i in 1:length(all)) {
    agent <- all[[i]] 
    if (is.na(agent$states$died)) alive_ids <- append(alive_ids,agent$states$id)
  }
  return(alive_ids)
}
  

# RunSimulation Function -------------------------------------------------------

### The backbone for an agent-based model simulation.  
### Usage: RunSimulation()
### Arguments: sim = list of (agents, pars, and spatial)   
###            runs = number of runs
###            write = write  
###            output_dir = output files directory
### Returns: A list object 
### Notes:
### Blake Massey
### 2015.03.15

RunSimulation <- function(sim = sim, 
                          runs = 1,
                          write = FALSE,
                          output_dir = getwd()) {
  runs <- CreateRunsList(runs)
  for (i in 1:length(runs)) {
    sum_intervals <- CreateSumIntervals(sim)   
    sim$agents$all <- UpdateAgentStates(init=TRUE)
    sim$agents$all <- UpdateAgentStepData(init=TRUE)
    sim$agents$all <- UpdateAgentIntData(init=TRUE)
    sim$agents$all <- UpdateAgentParsData(init=TRUE)
    sim$agents$agents_sum_int_data <- UpdateAgentsSumIntData(init=TRUE)
    sim$agents$pop_sum_int_data <- UpdatePopSumIntData(init=TRUE)
    sim$spatial <- UpdateSpatial(init=TRUE)
    for (j in 1:length(sum_intervals)) {
      step_intervals <- CreateStepIntervals(sum_intervals[[j]])
      for (k in 1:length(step_intervals)) {
        time_steps <- CreateTimeSteps(step_intervals[[k]])
        for (m in 1:length(time_steps)) {
          step <- time_steps[[m]]
          alive_ids <- ReturnAliveIds(sim$agents$all)
          sim$agents$all <- UpdateAgentParsData(sim$agents$all)
          for (n in alive_ids) {
            agent_states <- sim$agents$all[[n]][["states"]] 
            step_data <- sim$agents$all[[n]][["step_data"]]
            int_data <- sim$agents$all[[n]][["int_data"]]            
            pars_data <- sim$agents$all[[n]][["pars_data"]]
            ##### Movement, Survival, Reproduction, etc. Submodels ########
            step_data <- MovementSubModel(agent_states, step_data, step)
            agent_states <- SurvivalSubModel(agent_states, step_data)
            agent_states <- ReproductionSubModel(agent_states, step_data) 
            ###############################################################
            sim$agents$all[[n]][["states"]] <- UpdateAgentStates(agent_states)
            sim$agents$all[[n]][["step_data"]] <- UpdateAgentStepData(step_data)
            sim$agents$all[[n]][["int_data"]] <- UpdateAgentIntData(int_data)
          } # end of alive_ids[[n]]
          sim$spatial <- UpdateSpatial(sim$spatial)
        } # end of time_steps[[m]]
      } # end of step_interval[[k]]
    sim$agents$agents_sum_int_data <- UpdateAgentsSumIntData()
    sim$agents$pop_sum_int_data <- UpdatePopSumIntData()
    } # end of sum_interval[[k]]
    runs[[i]] <- sim
    WriteSimList(write)
  }
  return(runs)
}


# SurvivalSubModel Function ----------------------------------------------------

###  Survival submodel 
###  Usage: SurvivalSubModel()
###  Arguments: agent_states, a 'agent_states' list object
###  Returns: 'agent_states' list object
###  Notes: Just a placeholder for now. 
###  Blake Massey
###  2015.03.29

SurvivalSubModel <- function(agent_states = agent_states,
                             step_data = step_data) {
  agent_states <- agent_states
  return(agent_states)
}

# UpdateAgentStates Function ---------------------------------------------------

###  Creates and updates agents list from the sim object
###  Usage: CreateAgentStates(sim)
###  Arguments: agent_states = a list from sim$agents$all$agent$states
###             init = logical, whether or not this is the initation step
###  Returns: an 'agents' list
###  Notes: New columns can be via the "add_columns" vector w  
###  Blake Massey
###  2015.03.24

UpdateAgentStates <- function(agent_states = NULL,
                              init = FALSE) {
  if (init == TRUE) {
    input <- sim$agents$input
    input_columns <- colnames(input)
    na_columns <- c("start_datetime", "died")
    all <- list()
    for (i in 1:nrow(input)) {
      states <- list()
      for (j in input_columns) states <- append(states, input[i, j])
      for (k in 1:length(na_columns)) states <- append(states, NA) 
      states <- setNames(states, c(input_columns, na_columns))    
      agent <- NamedList(states)
      all <- append(all, NamedList(agent))
    }
    return(all)
  } else {
    agent_states <- agent_states
    return(agent_states)
  }
}

# UpdateAgentStepData Function -------------------------------------------------

###  Creates and updates agents step data database
###  Usage: CreateAgentStepData(agents, init)
###  Arguments: agents = agents object
###             init = logical, whether or not this is the initation step
###  Returns: an 'agents' list
###  Notes: 
###  Blake Massey
###  2015.03.27

UpdateAgentStepData <- function(step_data = NULL,
                                init = FALSE) {
  source('C:/Work/R/Functions/gen.R')
  if (init == TRUE) {
    sim_start <- sim$pars$global$sim_start
    all <- sim$agents$all
    for (i in 1:length(all)) {
      agent <- all[[i]]
      step_data <- data.frame(id=agent$states$id, datetime=sim_start, 
        x=agent$states$start_x, y=agent$states$start_y, exp_angle=NA, 
        abs_angle=NA)
      agent  <-  append(agent, NamedList(step_data)) 
      all[[i]] <- agent  
    }
    return(all)
  } else {
    step_data <- step_data
    return(step_data)
  }
}

# UpdateAgentIntData Function --------------------------------------------------

###  Creates and updates agents interval data dataframe
###  Usage: CreateAgentInData(agents, init)
###  Arguments: int_data = agents object
###             init = logical, whether or not this is the initation step
###  Returns: an 'agents' list
###  Notes: 
###  Blake Massey
###  2015.03.27

UpdateAgentIntData <- function(agent_states = NULL,
                               step_data = NULL,
                               int_data = NULL,
                               init = FALSE,
                               j = j,
                               k = k,
                               m = m) {
  source('C:/Work/R/Functions/gen.R')
  all <- sim$agents$all
  if (init == TRUE) {
    input <- sim$agents$input
    for (i in 1:length(all)) {
       agent <- all[[i]]
       int_data <- data.frame(id=agent$states$id, sum_int = NA)
       all[[i]] <- append(agent, NamedList(int_data))   
    }
    return(all)
  } else {
#    if (k == length(sum_interval) && m == length(steps)){
#     int <- sum_intervals[[j]]
#     interval_step_data <- step_data[which(step_data$datetime %within% int), ] 
#     int_data <- int_data
#    }
#    return(int_data)
  }
}

# UpdateAgentParsData Function -------------------------------------------------

###  Creates and updates agents parameters dataframe
###  Usage: CreateAgentParsData(agents, init)
###  Arguments: pars_data = agents object
###             init = logical, whether or not this is the initation step
###  Returns: an 'agents' list
###  Notes: 
###  Blake Massey
###  2015.03.27

UpdateAgentParsData <- function(pars_data = NULL,
                                init = FALSE) {
  source('C:/Work/R/Functions/gen.R')
  all <- sim$agents$all
  if (init == TRUE) {
    input <- sim$agents$input
    for (i in 1:nrow(input)) {
      agent <- all[[i]]
      pars_data <- data.frame(id=agent$states$id, sum_int = NA)
      all[[i]] <- append(agent, NamedList(pars_data))   
    }
    return(all)
  } else {
    pars_data <- pars_data
    return(pars_data)
  }
}

# UpdateAgentsSumIntData Function ----------------------------------------------

###  Creates and updates summary interval dataframe
###  Usage: UpdateAgentsSumIntData(agents, init)
###  Arguments: agents_sum_int_data = agents summary interval data, usually from
###               sim$agents$agents_sum_int_data
###             init = logical, whether or not this is the initation step
###  Returns: an 'agents_sum_int_data' dataframe
###  Notes: 
###  Blake Massey
###  2015.03.27

UpdateAgentsSumIntData <- function(agents_sum_int_data = NULL,
                                   init = FALSE) {
  source('C:/Work/R/Functions/gen.R')
  if (init == TRUE) {
    agents_sum_int_data <- data.frame(sum_int = NA, id = NA, turbines_2km = NA)
  } else {
    agents_sum_int_data <- sim$agents$agents_sum_int_data
#   alive_ids <- ReturnAliveIds()
#   for (i in alive_ids) {
#   agent_sum_int_data <- sim$agents$all$ 
#   agents_sum_int_data <- rbind(agents_sum_int_data, agent_sum_int_data)
  }
  return(agents_sum_int_data)
}

# UpdatePopSumIntData Function -------------------------------------------------

###  Creates and updates population interval dataframe
###  Usage: UpdatePopSumData(init)
###  Arguments: init = logical, whether or not this is the initation step
###  Returns: a 'pop_sum_int_data' dataframe
###  Notes: 
###  Blake Massey
###  2015.03.27

UpdatePopSumIntData <- function(init = FALSE) {
  source('C:/Work/R/Functions/gen.R')
  if (init == TRUE) {
    pop_sum_int_data <- data.frame(sum_int = NA, total_n = NA, adults = NA, 
      mf_ratio = NA)    
  } else {
    pop_sum_int_data <- sim$agents$pop_sum_int_data
#    agents <- rbind(pop_sum_int_data)
#    sim$agents$pop_sum_int_data
  }
  return(pop_sum_int_data)
}

# UpdateSpatial Function ----------------------------------------------------

###  Creates and updates population interval dataframe
###  Usage: UpdateSpatial(spatial, init)
###  Arguments: spatial = spatial object
###             init = logical, whether or not this is the initation step
###  Returns: an 'agents' list
###  Notes: 
###  Blake Massey
###  2015.03.27

UpdateSpatial <- function(agents = agents, 
                          init = FALSE) {
  spatial <- sim$spatial
  if (init == TRUE) {
  
  } else {
#   spatial_timer <- UpdateSpatialTimer(spatial_timer)
#   if (spatial_timer == timer_number) UpdateSpatial(); rm(spatial_timer)    
#   spatial <- append()
  }
  return(spatial)
}

# WriteSimList Function --------------------------------------------------------

### Writes "sim" list to a directory, 
### Usage: WriteSimList(write, sim, output_dir, components)
### Arguments: write = logical, whether or not to write the sim list to a file.
###              Default is TRUE.
###            run = name of run. Default uses the name from the runs object to 
###              ensure that the proper number of leading zeros is used. If 
###              default is not use, the name will simply be: "sim_(run).RData"
###            sim = a 'sim' list object
###            output_dir = output directory, default is working directory
###            components = components of 'sim' to write. Options include: 
###              "agents", "pars", and "spatial". Default is all.
### Returns: Writes a file to the output_dir
### Notes:
### Blake Massey
### 2015.03.29

WriteSimList <- function(write = TRUE,
                         run = names(runs[j]),
                         sim = sim,
                         output_dir = getwd(),
                         components = "all") {
  if (write == TRUE) {
    if (components == "all") {
    file_path = file.path(output_dir, paste0("sim_",run,".RData"))
    print(paste0("Writing: ", file_path)) 
    save(sim, file = file_path)   
    }
  }
}
  