require(data.table)
require(stAirPol)
require(ggplot2)
path = '~/data'
start_date = "2019-01-20"
end_date = "2019-01-20"

m.plz <- c(80539)
#Now you have to decide what the gridcellsize should be used, see
#?sf::st_make_grid for information’s about the cellsize parameter
m.grid_cellsize <- 0.005


#The specification which is needed is the aggregation interval and the
#time shift which is apply to the data. We choose an aggregation interval of 8
#hours, for more information about the aggregration_interval units see
#?lubridate::round_date. The time shift is applied to the data.
#Run the print method on the object for more information.

m.agg_info <- aggregation_information(timeshift = lubridate::hours(0),
                                      aggregation_interval = '24 hours')
print(m.agg_info)



# Data gathering ----------------------------------------------------------
m.date_pattern <- unique(substring(as.character(seq(as.Date(start_date),
                                             as.Date(end_date), 1)), 1,7))

#Collection the Information’s about the sensors and the collected data from
#the sensors in the chosen area.
sensors <- get_sensors(date_pattern = m.date_pattern, plz = m.plz, path = path)
sensor_age <- get_sensor_age(path = path)
# Traffic data ------------------------------------------------------------
sensor_data <- get_sensor_measured_values(sensors, m.date_pattern, path = path)


#Please note, if you want to gather data for a big area, I suggest to use a
#fixed value for lambda, e.g. lambda = 0.1, because the optimization of
#optim_lambda() will take a huge amount of computation costs.
lambda.p1 <- lambda.p2 <- 0.1
#If an optimization is applied, please check the validation plot.
#The local maximum should not be on the boarders, if it is please
#specify the lambda_range parameter.

estimate_grid_size(m.plz)
roads <- get_opentransportmap_data(m.plz, path = path, trafficvol_treshold = 1)
lambda.p1 <- optim_lambda(sensor_data[['P1']], sensors[['P1']], roads = roads,
                          validation_plot = TRUE)
lambda.p2 <- optim_lambda(sensor_data[['P2']], sensors[['P2']], roads = roads,
                          validation_plot = TRUE)
grid.traffic.p1 <- make_grid_traffic(lambda.p1, m.plz)
grid.traffic.p2 <- make_grid_traffic(lambda.p2, m.plz)
data.traffic.p1 <- make_data_traffic(sensors = sensors[['P1']], lambda = lambda.p1)
data.traffic.p2 <- make_data_traffic(sensors = sensors[['P2']], lambda = lambda.p2)

# Space Time data ---------------------------------------------------------
calculate_space_time_datasets()

# Time Variables ----------------------------------------------------------
# calculate_time_datasets()

setnames(data.rainhist.p1, "hist", 'rainhist')
setnames(grid.rainhist.p1, "hist", 'rainhist')
setnames(data.rainhist.p2, "hist", 'rainhist')
setnames(grid.rainhist.p2, "hist", 'rainhist')
setnames(data.windhist.p1, "hist", 'windhist')
setnames(grid.windhist.p1, "hist", 'windhist')
setnames(data.windhist.p2, "hist", 'windhist')
setnames(grid.windhist.p2, "hist", 'windhist')
setnames(data.humi.p2, "prediction", 'humi')
setnames(grid.humi.p2, "prediction", 'humi')
setnames(data.temp.p2, "prediction", 'temp')
setnames(grid.temp.p2, "prediction", 'temp')
setnames(data.humi.p1, "prediction", 'humi')
setnames(grid.humi.p1, "prediction", 'humi')
setnames(data.temp.p1, "prediction", 'temp')
setnames(grid.temp.p1, "prediction", 'temp')

# Combine all calculated datasets -----------------------------------------
# combine_datasets()
data.p1 <- get_model_frame(sensor_data[['P1']], sensors[['P1']], expand = TRUE)
data.final.p1 <- data.table(Reduce(function(x, y) dplyr::left_join(x, y),
                                   list(data.p1, data.humi.p1, data.temp.p1,
                                        #data.rainhist.p1, data.windhist.p1,
                                        data.traffic.p1, sensor_age)))
grid.final.p1 <- data.table(Reduce(function(x, y) dplyr::left_join(x, y),
                                   list(grid.humi.p1, grid.temp.p1,
                                        #grid.rainhist.p1, grid.windhist.p1,
                                        grid.traffic.p1)))

data.p2 <- get_model_frame(sensor_data[['P2']], sensors[['P2']], expand = TRUE)
data.final.p2 <- data.table(Reduce(function(x, y) dplyr::left_join(x, y),
                                   list(data.p2, data.humi.p2, data.temp.p2,
#                                        data.rainhist.p2, data.windhist.p2,
                                        data.traffic.p2, sensor_age)))
grid.final.p2 <- data.table(Reduce(function(x, y) dplyr::left_join(x, y),
                                   list(grid.humi.p2, grid.temp.p2,
#                                        grid.rainhist.p2, grid.windhist.p2,
                                        grid.traffic.p2)))
grid.final.p2[, sensor_id := .GRP, by = list(lon, lat)]
grid.final.p2[, sensor_age := 1.1]
grid.final.p1[, sensor_id := .GRP, by = list(lon, lat)]
grid.final.p1[, sensor_age := 1.1]
grid.final.p1 <<- grid.final.p1[substring(timestamp,1,7) %in% m.date_pattern]
grid.final.p2 <<- grid.final.p2[substring(timestamp,1,7) %in% m.date_pattern]
data.final.p1 <<- data.final.p1[substring(timestamp,1,7) %in% m.date_pattern]
data.final.p2 <<- data.final.p2[substring(timestamp,1,7) %in% m.date_pattern]


# Save RDS files ----------------------------------------------------------
saveRDS(data.final.p1, file = paste0(path, '/p1_model_data_',m.grid_cellsize, '.rds'))
saveRDS(grid.final.p1, file = paste0(path, '/p1_grid_data_',m.grid_cellsize, '.rds'))
saveRDS(data.final.p2, file = paste0(path, '/p2_model_data_',m.grid_cellsize, '.rds'))
saveRDS(grid.final.p2, file = paste0(path, '/p2_grid_data_',m.grid_cellsize, '.rds'))


