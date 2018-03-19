# Plotting Sam's NLP 1-week plots
# 2018-03-18

using DataFrames
using CSV
using Gadfly


NLP_data = CSV.read("NLP_one_week.csv", nullable = false)

test = DataFrame(test_x = 1:10, test_y = ones(10))

plot(test, x = "test_x", y = "test_y", Geom.line)

light_theme = Theme(
    background_color = "white",
    panel_fill = "white",
    default_color = "blue",
    key_position = :right
    )

CMH_plot =
    plot(NLP_data,
        x = "Timestep",
        y = "CMH",
        Geom.line,
        Guide.Title("Air Flow into Building"),
        Guide.XLabel("Time (hrs)"),
        Guide.YLabel("Air Flow (CMH)"),
        light_theme
    )

CO2_plot =
    plot(
        NLP_data,
        x = "Timestep",
        y = "CO2",
        Geom.line,
        Guide.Title("CO2 Concentration inside room (ppm)"),
        Guide.XLabel("Time (hrs)"),
        Guide.YLabel("CO2 Concentration (ppm)"),
        light_theme
    )

PM25_plot =
    plot(
        NLP_data,
        x = "Timestep",
        y = "PM25",
        Geom.line,
        Guide.Title("PM2.5 Concentration (µg)"),
        Guide.XLabel("Time (hrs)"),
        Guide.YLabel("PM2.5 Concentration (micrograms / m^3)"),
        light_theme
    )

PM25_absorbed_plot =
    plot(
        NLP_data,
        x = "Timestep",
        y = "PM25_removed",
        Geom.line,
        Guide.Title("PM2.5 Removed by Filters"),
        Guide.XLabel("Time (hrs)"),
        Guide.YLabel("PM2.5 Removed (micrograms)"),
        light_theme
    )

HUMID_plot =
    plot(
        NLP_data,
        x = "Timestep",
        y = "Humidity",
        Geom.line,
        Guide.Title("Humidity"),
        Guide.XLabel("Time (hrs)"),
        Guide.YLabel("Humidity Ratio (moisture content)"),
        light_theme
    )

HUMID_absorbed_plot =
    plot(
        NLP_data,
        x = "Timestep",
        y = "Humidity_removed",
        Geom.line,
        Guide.Title("Moisture Removed"),
        Guide.XLabel("Time (hrs)"),
        Guide.YLabel("Moisture Removed by Dehumidification (g H2O)"),
        light_theme
    )


final = vstack(
    hstack(CMH_plot, CO2_plot),
    hstack(PM25_plot, PM25_absorbed_plot),
    hstack(HUMID_plot, HUMID_absorbed_plot)
    )

img = SVG("NLP_one_week.svg", 12inch, 12inch)
draw(img, final)
