# Flowering-model utility functions.
# Original contributions: Flavian Tschurr (2024).

gapfill_and_smooth_function <- function(data_vect, timestamp){
  nr_knots = round(length(data_vect)/2)
  timestamp <- as.Date(timestamp)
  data_table <- data.frame(Date=timestamp, value=data_vect)
  data_table <- data_table[order(data_table$Date),]
  data_vect <- seq.Date(data_table$Date[1],data_table$Date[length(data_table$Date)],by="days")
  to_merge <- data.frame(Date =data_vect)
  data_table_long <- merge(to_merge,data_table,by.x="Date",all=TRUE)
  
  data_table_long$value_filled <- approx(data_table_long$value,n=length(data_table_long$value))$y
  data_table_long$value_smoothed  <- smooth.spline(data_table_long$value_filled, nknots=nr_knots)$y
  
  
  return(data_table_long)
  
}


DD <- function(expr, name, order = 1) {
  # Differentiate the expression recursively.
  if(order < 1) stop("'order' must be >= 1")
  if(order == 1) D(expr, name)
  else DD(D(expr, name), name, order - 1)
}
