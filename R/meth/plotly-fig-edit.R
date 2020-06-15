library(plotly)

Sys.setenv("plotly_username"="floris.barthel")
Sys.setenv("plotly_api_key"="xDVPnQwm21HVUZWA8iaD")

fig <- api_download_plot(143, "anderk")

fig@data
