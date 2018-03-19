# Plotting Sam's NLP 1-week plots
# 2018-03-18

using DataFrames
using CSV

readcsv("NLP_one_week.csv")

df = CSV.read("NLP_one_week.csv")


CMH_plot =
    plot(
        x = 1:length(CMHresult),
        y = CMHresult,
        Geom.line,
        Guide.Title("CMH at timestep x")
    )

CO2_plot =
    plot(
        x = 0:length(CO2result) - 1,
        y = CO2result,
        Geom.line,
        Guide.Title("CO2 Concentration inside room (ppm)")
    )

PM25_plot =
    plot(
        x = 0:length(PM25result) - 1,
        y = PM25result,
        Geom.line,
        Guide.Title("PM2.5 Concentration (µg)")
    )

PM25_absorbed_plot =
    plot(
        x = 1:length(PM25Absorbedresult),
        y = PM25Absorbedresult,
        Geom.line,
        Guide.Title("PM2.5 Removed by Filters")
    )

HUMID_plot =
    plot(
        x = 0:length(HUMIDresult) - 1,
        y = HUMIDresult,
        Geom.line,
        Guide.Title("Humidity Ratio (kg Water / kg dry air)")
    )

HUMID_absorbed_plot =
    plot(
        x = 1:length(HUMIDAbsorbedresult),
        y = HUMIDAbsorbedresult,
        Geom.line,
        Guide.Title("Moisture Removed by Dehumidification")
    )


final = vstack(
    hstack(CMH_plot, CO2_plot),
    hstack(PM25_plot, PM25_absorbed_plot),
    hstack(HUMID_plot, HUMID_absorbed_plot)
    )

img = SVG("NLP_one_week.svg", 12inch, 12inch)
draw(img, final)
