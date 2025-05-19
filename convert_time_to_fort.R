
library(FortMyrmidon)
df <- read.csv("~/Downloads/trophallaxis_by_excluded.csv")

# create fort time to copy pasta in myrmidon
df$time_start_GMT <- as.POSIXct(df$start, format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
df$time_end_GMT <- as.POSIXct(df$end, format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
df$start_fort <- NULL
df$end_fort <- NULL
for (i in 1:nrow(df)) { #i = 1
  start_fort  <- capture.output(fmTimeCreate(df$time_start_GMT[i]))
  end_fort    <- capture.output(fmTimeCreate(df$time_end_GMT[i]))
  df$start_fort[i] <- start_fort
  df$end_fort[i]   <- end_fort
}

# bring new columns to the front of the df
new_order <- c(
  names(df)[1:3],              # Keep columns 1–3 as they are
  "start_fort", "end_fort",    # Move these to positions 4 and 5
  setdiff(names(df), c("start_fort", "end_fort", names(df)[1:3]))  # All others except the above
)
df <- df[, new_order]

#safe
write.csv(df, "~/Downloads/trophallaxis_by_excluded_DS.csv", row.names = FALSE)
