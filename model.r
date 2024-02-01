# Preliminary seabird colony SID model
# L.U.T. 2023-09-11
library(tidyverse)
library(patchwork)

# Biological parameters
SEASON_LENGTH          <- 121       # Number of days between May 1 and September 1
MIGRATION_INCOMING_DAY <- 0         # Season day that adults begin migrating into the breeding ground (May 1)
MIGRATION_OUTGOING_DAY <- 76        # Season day that adults begin migrating out of the breeding ground (July 15)
MIGRATION_IN_RATE      <- 0.05      # S adults arriving per day via migration (proportion of max population)
MIGRATION_OUT_RATE     <- 0.02      # S birds leaving per day via migration (proportion of current population)
INCOMING_I_RATE        <- 0.01     # I adults arriving on infection day (proportion of max population)

MAX_ADULTS             <- 1000      # Total population size of adults

HATCH_DAY              <- 42        # Season day that chicks start hatching (29 days after May 14)              

PRODUCTIVITY <- 0.5 * 2.33 * 0.81   # Productivity rate (== nests per parent * clutch size * hatch rate)

MAX_CHICKS             <- MAX_ADULTS * PRODUCTIVITY

# Infection rate function, scaling exponentially from 0 to 0.5 across the range from 0 to max individuals
ir <- function(n, maxInfectionRate) {
    b <- 0.00001  # Parameter b controls the rate of increase
    maxInfectionRate * (exp(b*n)-1) / (exp(b*(MAX_ADULTS+MAX_CHICKS))-1)
}

DEATH_RATE  <- 0.30       

# Function to simulate disease dynamics based on variable infection day
simulateSeason <- function(maxInfectionRate, infectionDay) {
    # Initialize model trackers
    S <- 0
    I <- 0
    D <- 0
    adultCounter <- 0

    for (day in 0:SEASON_LENGTH) {
        ## Migration and hatching

        # Incoming migration of susceptible adults
        if (day >= MIGRATION_INCOMING_DAY & adultCounter < MAX_ADULTS) {
            newAdults <- MAX_ADULTS * MIGRATION_IN_RATE
            adultCounter <- adultCounter + newAdults
            S <- S + newAdults
        }

        # Hatch chicks
        if (day == HATCH_DAY) {
            newChicks <- floor(S * PRODUCTIVITY)
            S <- S + newChicks
        }

        # Outgoing migration of susceptible (i.e., healthy) birds only
        if (day >= MIGRATION_OUTGOING_DAY) {
            S <- S - S * MIGRATION_OUT_RATE
        }

        # On infection day, add burst of infectious adults
        if (day == infectionDay) {
            I <- I + MAX_ADULTS * INCOMING_I_RATE
        } 

        ## Infection  updates

        infectionRate <- ir(S+I, maxInfectionRate)
        newInfections <- min(S, floor(infectionRate * (S*I)/(S+I)))
        newDeaths <- min(I, floor(DEATH_RATE * I))

        S <- S - newInfections
        I <- I + newInfections - newDeaths
        D <- D + newDeaths
    }

    monitor <- tibble(Max_Infection_Rate = maxInfectionRate,
                      Infection_Day = infectionDay,
                      Chicks_Born = newChicks,
                      D = D)
    return(monitor)
}

# Generate test conditions
maxInfectionRate <- seq(0, 1, by=0.01)
infectionDay  <- 0:99

variables <- expand.grid(maxInfectionRate, infectionDay)

# Test results from all possible infection days
results <- map2_dfr(variables[,1], variables[,2], simulateSeason)

summarized <- results |>
           group_by(Max_Infection_Rate, Infection_Day) |>
           summarize(Chicks_Born = max(Chicks_Born),
                     D = max(D))

plot_deaths <- ggplot(summarized) +
            geom_tile(aes(x=Infection_Day, y=Max_Infection_Rate, fill=D)) +
            scale_fill_gradient(low="#cfebf4", high="black", limits=c(0, (MAX_ADULTS+MAX_CHICKS))) +
            xlab("Date of infectious bird arrival") +
            ylab("Maximum density-dependent infection rate") +
            guides(fill=guide_legend(title="Total deaths")) +
            theme(legend.position = "top",
                  panel.grid=element_blank(),
                  panel.background=element_blank())

plot_chicks <- ggplot(summarized) +
            geom_tile(aes(x=Infection_Day, y=Max_Infection_Rate, fill=Chicks_Born)) +
            scale_fill_gradient(low="black", high="#d2f4cff4", limits=c(0, MAX_CHICKS)) +
            xlab("Date of infectious bird arrival") +
            ylab("Maximum density-dependent infection rate") +
            guides(fill=guide_legend(title="Chicks Born")) +
            theme(legend.position = "top",
                  panel.grid=element_blank(),
                  panel.background=element_blank())        

plot_combined <- plot_deaths + plot_chicks

ggsave(plot_combined, filename="C:/Users/ltayl/Desktop/temp_seabird_sid_plots.pdf",
       width=10, height=5)

# OLD
# plot <- ggplot(results) + 
#      geom_line(aes(x=Day, y=D/(MAX_ADULTS+MAX_CHICKS))) +
#      facet_wrap(facets=vars(Infection_Day), nrow=10, ncol=10) +
#      scale_x_continuous(limits=c(0, SEASON_LENGTH), breaks=seq(0, SEASON_LENGTH, by=50)) +
#      scale_y_continuous(limits=c(0, 0.5), breaks=seq(0, 1, by=0.25)) +
#      xlab("Day") +
#      ylab("Deaths (prop. of maximum pop.)") +
#      ggtitle("Day of infectious arrival") +
#      theme_bw()
# ggsave(filename="seabird_sid_plot.pdf", plot, width=10, height=10)
