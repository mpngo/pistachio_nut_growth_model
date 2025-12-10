
library(shiny)
library(caladaptr)
library(httr)
library(jsonlite)
library(dplyr)
library(lubridate)
library(plotly)
library(splines)
library(leaflet)
library(shinycssloaders)
library(ggplot2)
library(terra)
library(splines)
library(rsconnect)



# Load spline model
spline_model <- readRDS("kernel_dw_model.rds")
token <- "28260dffc4744631b9fddf06c34f6a10"

# Unified GDD threshold
GDD_MAX <- 5000

# Trait prediction function (single trait)
predict_trait <- function(gdd, trait) {
  if (is.na(gdd) || gdd > GDD_MAX) return(NA)
  log_GDD <- log(gdd)
  
  switch(trait,
         "kernel_dry_weight" = if (gdd >= min(spline_model$x) & gdd <= max(spline_model$x)) predict(spline_model, x = gdd)$y else NA,
         "nut_area" = { raw <- -10917.8286 + 4438.07059 * log_GDD - 587.23897 * log_GDD^2 + 25.97156 * log_GDD^3; if (raw < 0) NA else raw },
         "shell_texture" = exp(-2.570567 + 0.00449783 * gdd - 1.59005e-6 * gdd^2 + 2.61629e-10 * gdd^3 - 1.66428e-14 * gdd^4),
         "hull_texture" = exp(-0.356363 + 0.0013867 * gdd - 3.70287e-7 * gdd^2 + 2.005734e-11 * gdd^3),
         "nut_L" = exp(3.657422 + 0.001306558 * gdd - 1.105751e-6 * gdd^2 + 4.300462e-10 * gdd^3 - 7.604527e-14 * gdd^4 + 4.972184e-18 * gdd^5),
         "nut_A" = exp(2.108671 - 0.0027225 * gdd + 2.208453e-6 * gdd^2 - 5.732717e-10 * gdd^3 + 4.988507e-14 * gdd^4) - 14.8517,
         "nut_B" = ((2.33727 + 3.802293e-4 * gdd - 4.228952e-7 * gdd^2 + 1.750973e-10 * gdd^3 - 3.119361e-14 * gdd^4 + 1.993278e-18 * gdd^5) * -0.26263 + 1)^(1 / -0.26263),
         "kernel_area" = if (gdd >= 2000) -11001.714 + 2587.0618 * log_GDD - 149.9839 * log_GDD^2 else NA,
         "kernel_texture" = if (gdd >= 2000) ((-2.989603 + 0.0014403 * gdd - 1.35269e-7 * gdd^2) * -0.1414 + 1)^(1 / -0.1414) else NA,
         "kernel_L" = if (gdd >= 2000) ((3.739154 - 1.452817e-3 * gdd + 5.645595e-7 * gdd^2 - 9.746307e-11 * gdd^3 + 6.366107e-15 * gdd^4) * -0.3030303 + 1)^(1 / -0.3030303) else NA,
         "kernel_A" = if (gdd >= 2000) exp(8.945922 - 0.00507501 * gdd + 1.326560e-6 * gdd^2 - 1.067834e-10 * gdd^3) - 39.942 else NA,
         "kernel_B" = if (gdd >= 2000) exp(-0.862436 + 0.0038827 * gdd - 1.075414e-6 * gdd^2 + 9.733552e-11 * gdd^3) else NA,
         NA)
}

# Batch prediction for entire dataframe
apply_all_traits <- function(df) {
  df %>%
    rowwise() %>%
    mutate(
      kernel_dry_weight = predict_trait(GDD_cumulative, "kernel_dry_weight"),
      nut_area = predict_trait(GDD_cumulative, "nut_area"),
      shell_texture = predict_trait(GDD_cumulative, "shell_texture"),
      hull_texture = predict_trait(GDD_cumulative, "hull_texture"),
      nut_L = predict_trait(GDD_cumulative, "nut_L"),
      nut_A = predict_trait(GDD_cumulative, "nut_A"),
      nut_B = predict_trait(GDD_cumulative, "nut_B"),
      kernel_area = predict_trait(GDD_cumulative, "kernel_area"),
      kernel_texture = predict_trait(GDD_cumulative, "kernel_texture"),
      kernel_L = predict_trait(GDD_cumulative, "kernel_L"),
      kernel_A = predict_trait(GDD_cumulative, "kernel_A"),
      kernel_B = predict_trait(GDD_cumulative, "kernel_B")
    ) %>%
    mutate( # NA for negative values
      kernel_dry_weight = ifelse(!is.na(kernel_dry_weight) & kernel_dry_weight < 0, NA, kernel_dry_weight),
      nut_area = ifelse(!is.na(nut_area) & nut_area < 0, NA, nut_area),
      shell_texture = ifelse(!is.na(shell_texture) & shell_texture < 0, NA, shell_texture),
      hull_texture = ifelse(!is.na(hull_texture) & hull_texture < 0, NA, hull_texture),
      nut_L = ifelse(!is.na(nut_L) & nut_L < 0, NA, nut_L),
      nut_B = ifelse(!is.na(nut_B) & nut_B < 0, NA, nut_B),
      kernel_area = ifelse(!is.na(kernel_area) & kernel_area < 0, NA, kernel_area),
      kernel_texture = ifelse(!is.na(kernel_texture) & kernel_texture < 0, NA, kernel_texture),
      kernel_L = ifelse(!is.na(kernel_L) & kernel_L < 0, NA, kernel_L),
      kernel_B = ifelse(!is.na(kernel_B) & kernel_B < 0, NA, kernel_B)
    ) %>%
    ungroup()
}


# Unified plot generator function
generate_plotly <- function(df, y_var, x_axis_type = "date", custom_title = NULL) {
  y_unit <- switch(y_var,
                   "kernel_area" = "(mm²)", "nut_area" = "(mm²)", "shell_texture" = "(kg force)",
                   "hull_texture" = "(kg force)", "kernel_texture" = "(kg force)", "kernel_dry_weight" = "(g)", "")
  y_label <- tools::toTitleCase(gsub("_", " ", y_var))
  
  # Color traits block
  if (y_var %in% c("nut_color", "kernel_color")) {
    color_prefix <- if (y_var == "nut_color") "Nut" else "Kernel"
    p <- plot_ly()
    
    for (col in c("L", "A", "B")) {
      trait_col <- paste0(tolower(color_prefix), "_", col)
      dash_styles <- c("Live" = "solid", "Forecast" = "dash", "Historical" = "dot")
      color_values <- c("L" = "#3a3b3c", "A" = "#D55E00", "B" = "#0072B2")
      label_map <- c("L" = "Lightness (L*)", "A" = "Green-Red (a*)", "B" = "Blue-Yellow (b*)")
      
      for (src in unique(df$source)) {
        display_src <- if (src == "Historical") "Projected" else src
        temp_df <- df %>% filter(source == src)
        x_vals <- temp_df[[x_axis_type]]
        x_label <- if (x_axis_type == "GDD_cumulative") "GDD" else "Date"
        x_hover <- if (x_axis_type == "GDD_cumulative") round(x_vals, 1) else format(x_vals, "%b %d")
        
        hover_text <- paste0(
          "Date: ", format(temp_df$date, "%b %d"),
          "<br>GDD: ", round(temp_df$GDD_cumulative, 1),
          "<br>", label_map[[col]], ": ", round(temp_df[[trait_col]], 2)
        )
        
        p <- p %>%
          add_lines(
            data = temp_df,
            x = ~get(x_axis_type),
            y = temp_df[[trait_col]],
            name = paste0(if (color_prefix == "Nut") "Hull" else "Kernel", " ", label_map[[col]], " (", display_src, ")"),
            text = hover_text,
            line = list(width = 2, dash = dash_styles[[src]], color = color_values[[col]]),
            showlegend = TRUE,
            hoverinfo = "text"
          )
      }
    }
    
    p <- p %>%
      layout(
        title = list(
          text = if (!is.null(custom_title)) custom_title else paste("Predicted", y_label),
          y = 0.95, xanchor = "center", yanchor = "top"
        ),
        xaxis = list(
          title = if (x_axis_type == "GDD_cumulative") "Cumulative GDD" else "Date",
          titlefont = list(size = 16),
          tickfont = list(size = 14)
        ),
        yaxis = list(
          title = paste(y_label, y_unit),
          titlefont = list(size = 16),
          tickfont = list(size = 14)
        ),
        legend = list(font = list(size = 14)),
        margin = list(l = 60, r = 130, t = 50, b = 60),
        showlegend = TRUE
      )
    
    return(no_zoom(p))
  }
  
  # General traits block
  p <- plot_ly()
  
  for (src in unique(df$source)) {
    display_src <- if (src == "Historical") "Projected" else src
    temp_df <- df %>% filter(source == src)
    x_vals <- temp_df[[x_axis_type]]
    x_label <- if (x_axis_type == "GDD_cumulative") "GDD" else "Date"
    x_hover <- if (x_axis_type == "GDD_cumulative") round(x_vals, 1) else format(x_vals, "%b %d")
    
    hover_text <- paste0(
      "Date: ", format(temp_df$date, "%b %d"),
      "<br>GDD: ", round(temp_df$GDD_cumulative, 1),
      "<br>", y_label, ": ", round(temp_df[[y_var]], 2), " ", y_unit
    )
    
    
    
    line_style <- switch(src,
                         "Live" = list(width = 2, dash = "solid", color = "#FFBF00"),
                         "Forecast" = list(width = 2, dash = "dash", color = "#D55E00"),
                         "Historical" = list(width = 2, dash = "dot", color = "#002855")
    )
    
    p <- p %>%
      add_lines(
        data = temp_df,
        x = ~get(x_axis_type),
        y = temp_df[[y_var]],
        name = paste0(display_src, " - ", y_label),
        line = line_style,
        text = hover_text,
        showlegend = TRUE,
        hoverinfo = "text"
      )
  }
  
  p <- p %>%
    layout(
      title = list(
        text = if (!is.null(custom_title)) custom_title else paste("Predicted", y_label),
        y = 0.95, xanchor = "center", yanchor = "top"
      ),
      xaxis = list(
        title = if (x_axis_type == "GDD_cumulative") "Cumulative GDD" else "Date",
        titlefont = list(size = 16),
        tickfont = list(size = 14)
      ),
      yaxis = list(
        title = paste(y_label, y_unit),
        titlefont = list(size = 16),
        tickfont = list(size = 14)
      ),
      legend = list(font = list(size = 14)),
      margin = list(l = 60, r = 130, t = 50, b = 60),
      showlegend = TRUE
    )
  
  return(no_zoom(p))
}



process_uploaded_gdd <- function(df) {
  req("GDD_cumulative" %in% names(df))
  df %>% apply_all_traits()
}

no_zoom <- function(p) {
  p %>%
    layout(
      xaxis = list(fixedrange = TRUE),
      yaxis = list(fixedrange = TRUE)
    ) %>%
    config(
      scrollZoom = FALSE,
      displayModeBar = TRUE,
      modeBarButtonsToRemove = c(
        "zoom2d","zoomIn2d","zoomOut2d",
        "autoScale2d","resetScale2d",
        "pan2d","select2d","lasso2d"
      )
    )
}

fetch_openmeteo_hist <- function(lat, lon, year_start, year_end, base_temp_f = 44.6) {
  # Returns daily rows with date, high_temp, low_temp, GDD, year, jday
  resp <- httr::GET(
    url = "https://archive-api.open-meteo.com/v1/archive",
    query = list(
      latitude = lat,
      longitude = lon,
      start_date = sprintf("%d-01-01", year_start),
      end_date   = sprintf("%d-12-31", year_end),
      daily = "temperature_2m_max,temperature_2m_min",
      temperature_unit = "fahrenheit",
      timezone = "auto"
    )
  )
  if (httr::status_code(resp) != 200) {
    stop("Open-Meteo archive request failed: ", httr::content(resp, "text", encoding = "UTF-8"))
  }
  dat <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"))
  if (is.null(dat$daily$time)) stop("Open-Meteo returned no daily data")
  
  dplyr::tibble(
    date = as.Date(dat$daily$time),
    high_temp = dat$daily$temperature_2m_max,
    low_temp  = dat$daily$temperature_2m_min
  ) |>
    dplyr::mutate(
      GDD  = pmax((high_temp + low_temp) / 2 - base_temp_f, 0),
      year = lubridate::year(date),
      jday = lubridate::yday(date)
    )
}


desktop_layout <- function() {
  div(
    style = "min-height: 100vh; display: flex; flex-direction: column;",
    div(
      style = "flex-grow: 1;",
      navbarPage(
        title = div(
          "Nut Growth and Maturity Forecast Tool",
          style = "font-size: 22px; font-weight: bold; color: #002855; font-family: 'Source Sans Pro', sans-serif;"
        ),
        id = "main_navbar",
        
        tabPanel("Growth and Maturity Prediction Tool",
                 fluidPage(
                   # keep your styles and small scripts
                   tags$head(
                     tags$style(HTML("
                        body { background-color: #ffffff; color: #002855; font-family: 'Source Sans Pro', sans-serif; }
                        .well { background-color: #f9f9f9; border-radius: 12px; padding: 20px; }
                        h4 { font-weight: bold; font-size: 20px; }
                        .tabbable > .nav > li[class=active] > a { border: none; border-bottom: 3px solid #FFBF00 !important; background-color: transparent !important; }
                        .leaflet-container { border-radius: 12px; }
                        .form-control { border-radius: 8px; }
                        .shiny-download-link { color: #0066cc !important; text-decoration: underline !important; font-size: 14px; }
                        .btn-primary { background-color: #FFBF00; color: #ffffff !important; border: none; text-decoration: none !important; }
                        .btn-primary:hover { background-color: #e0aa00; color: #ffffff !important; text-decoration: none !important; }
                        .form-control, .selectize-input, .selectize-dropdown, .control-label { font-size: 16px !important; }
                        .nav-tabs > li > a, .navbar-nav > li > a { font-size: 16px !important; font-weight: 500; color: #002855; }
                        .note-under-plot { margin-top: 10px; font-size: 16px; color: #555; font-style: italic; text-align: center; line-height: 1.35; }
                      "))
                   ),
                   
                   fluidRow(
                     column(
                       width = 4,
                       div(
                         tags$h4("Select Orchard Location", style = "margin-bottom: 5px;"),
                         leafletOutput("map", height = 400)
                       ),
                       br(),
                       wellPanel(
                         tags$div(
                           tags$h4("Tool Options", style = "margin-bottom: 5px;"),
                           tags$div("Version 1.0 (Last Updated 7/18/2025)",
                                    style = "font-size: 15px; color: #666; margin-bottom: 10px;")
                         ),
                         
                         radioButtons(
                           "loc_mode", "Choose location input",
                           choices = c("Latitude & Longitude" = "coords",
                                       "Street Address (auto-geocode)" = "address"),
                           selected = "coords", inline = TRUE
                         ),
                         
                         # Address mode UI
                         conditionalPanel(
                           condition = "input.loc_mode === 'address'",
                           selectizeInput(
                             "address_pick", "Search address",
                             choices = NULL,
                             options = list(
                               placeholder = "Type an address (min 3 chars)",
                               create = FALSE,
                               maxOptions = 5,
                               # Send the typed query to server whenever it changes
                               onType = I("
                                function (str) {
                                  if (!str || str.length < 3) return;
                                  Shiny.setInputValue('address_type', str, {priority: 'event'});
                                }
                              ")
                             )
                           ),
                           # optional small hint
                           tags$div(style = "font-size: 13px; color: #666; margin-top: -6px;",
                                    "Pick a suggestion to auto-fill latitude and longitude.")
                         ),
                         
                         # Coordinates mode UI (unchanged)
                         conditionalPanel(
                           condition = "input.loc_mode === 'coords'",
                           fluidRow(
                             column(6, numericInput("latitude",  "Latitude",  value = 38.5449)),
                             column(6, numericInput("longitude", "Longitude", value = -121.7405))
                           )
                         ),
                         
                         
                         selectInput(
                           "crop_type", "Select Crop Type",
                           choices = c("Pistachio: Kerman", "Pistachio: Golden Hills", "Pistachio: Gumdrop","Walnut"),
                           selected = "Pistachio: Kerman"
                         ),
                         div(
                           style = "display: flex; gap: 15px; align-items: flex-start;",
                           div(
                             style = "flex: 1;",
                             tags$label("Bloom Date", style = "font-size: 15px; font-weight: 600;"),
                             dateInput("start_date", NULL,
                                       value = as.Date(paste0(format(Sys.Date(), "%Y"), "-04-04")),
                                       max = Sys.Date())
                           ),
                           div(
                             style = "flex: 1;",
                             tags$label("End Date", style = "font-size: 15px; font-weight: 600;"),
                             dateInput("end_date", NULL,
                                       value = as.Date(paste0(format(Sys.Date(), "%Y"), "-07-15")),
                                       max = Sys.Date())
                           )
                         ),
                         selectInput("model_var", "Select Trait to Plot:",
                                     choices = c("Areas: Kernel + Nut" = "areas_combined",
                                                 "Textures: Kernel + Shell + Hull" = "textures_combined",
                                                 "Kernel Dry Weight" = "kernel_dry_weight",
                                                 "Hull Color (L*, A*, B*)" = "nut_color",
                                                 "Kernel Color (L*, A*, B*)" = "kernel_color")
                         ),
                         conditionalPanel(
                           condition = "['nut_color', 'kernel_color'].includes(input.model_var)",
                           actionLink("show_lab_modal_link", "What is LAB Color Space?",
                                      style = "font-size: 15px; color: #0072B2; margin-top: 10px;")
                         ),
                         div(
                           style = "font-size: 16px; font-weight: 600;",
                           tags$label("View by:", style = "display: block; margin-bottom: 5px;"),
                           div(
                             style = "display: flex; align-items: center; gap: 10px;",
                             radioButtons("x_axis_type", NULL,
                                          choices = c("Date" = "date", "Growing Degree Days" = "GDD_cumulative"),
                                          selected = "date", inline = TRUE)
                           )
                         ),
                         actionButton("get_data", "Generate Model",
                                      class = "btn btn-primary btn-block",
                                      style = "font-size: 17px; padding: 6px 12px; color: #0066cc;")
                       ),
                       uiOutput("links_live")
                     ),
                     column(
                       width = 8,
                       tabsetPanel(
                         tabPanel("Real-time Estimate",
                                  div(
                                    plotlyOutput("growth_plot_live", height = "800px"),
                                    uiOutput("live_note")
                                  )
                         ),
                         tabPanel("Projected Estimate",
                                  div(
                                    plotlyOutput("growth_plot_hist", height = "800px"),
                                    uiOutput("hist_note")
                                  )
                         ),
                         tabPanel("Combined",
                                  plotlyOutput("growth_plot_combined", height = "800px"),
                                  uiOutput("combined_note"),
                                  uiOutput("reference_image_ui")
                         )
                       )
                     )
                   )
                 )
        ),
        
        tabPanel("Upload Custom GDD",
                 fluidPage(
                   h3("Upload GDD Dataset"),
                   tags$p(
                     "Upload a CSV with a column named ", tags$code("GDD_cumulative"), ". Optionally include ",
                     tags$code("date"), ". If you include a date it will be used for hover text and date axes.",
                     style = "font-size: 15px; color: #444; margin-bottom: 10px;"
                   ),
                   wellPanel(
                     tags$div(style = "font-size: 15px;",
                              tags$b("Requirements"),
                              tags$ul(
                                tags$li(tags$code("GDD_cumulative"), " must be numeric and non-decreasing."),
                                tags$li("Units assumed Fahrenheit base 44.6."),
                                tags$li("Optional ", tags$code("date"), " formats: ",
                                        tags$code("YYYY-MM-DD"), ", ", tags$code("MM/DD"), ", ", tags$code("MM/DD/YYYY"), "."),
                                tags$li("Header row required. Extra columns are ignored."),
                                tags$li("Missing or negative values are dropped.")
                              ),
                              tags$b("Example (first rows)"),
                              tags$pre("date,GDD_cumulative\n04/04,0\n04/05,12.4\n04/06,20.7\n04/07,27.1"),
                              div(style = "margin-top: 8px;", downloadButton("download_template", "Download CSV template", class = "btn btn-primary"))
                     )
                   ),
                   fileInput("gdd_file", "Upload CSV file:", accept = ".csv"),
                   div(style = "display: flex; gap: 10px; align-items: center;",
                       actionButton("process_file", "Generate Predictions", class = "btn btn-primary"),
                       uiOutput("download_button_ui")),
                   br(), br(),
                   DT::dataTableOutput("gdd_table"),
                   uiOutput("uploaded_var_selector"),
                   plotlyOutput("uploaded_growth_plot")
                 )
        ),
        
        tabPanel("Instructions", fluidPage(includeHTML("www/instructions.html"))),
        
        tabPanel("Resources", fluidPage(
          tags$h3("Resources", style = "color: #002855; font-weight: bold;"),
          tags$p("For more info, please see these additional resources:"),
          # paper
          fluidRow(style = "margin-bottom: 30px;",
                   column(2, tags$img(src = "paper_thumbnail.png",
                                      style = "max-width: 100%; border: 1px solid #000; display: block; margin-top: 5px;")),
                   column(10, style = "display: flex; align-items: center;",
                          tags$div(
                            tags$a("Nut Maturity Modeling Research Paper",
                                   href = "https://www.biorxiv.org/content/10.1101/2024.06.24.600444v1.full",
                                   target = "_blank",
                                   style = "font-size: 16px; font-weight: bold; color: #1a0dab;"),
                            tags$br(),
                            tags$i("Includes model validation, heat unit thresholds, and GDD framework for nut development"),
                            tags$br(), tags$br(),
                            tags$div(style = "font-size: 14px; color: #444;",
                                     tags$i("Adaskaveg et al. (2024)"),
                                     tags$br(),
                                     tags$i("Modeling nut maturity in pistachio with temperature-based growth predictors."),
                                     tags$br(),
                                     tags$i("bioRxiv. https://doi.org/10.1101/2024.06.24.600444"))
                          )
                   )
          ),
          # previous tool
          fluidRow(style = "margin-bottom: 30px;",
                   column(2, tags$img(src = "tool_thumbnail.png",
                                      style = "max-width: 100%; border: 1px solid #000; display: block; margin-top: 5px;")),
                   column(10, style = "display: flex; align-items: center;",
                          tags$div(
                            tags$a("Previous GDD-Based Forecast Tool (UC ANR IGIS)",
                                   href = "https://ucanr-igis.shinyapps.io/pist_gdd/",
                                   target = "_blank",
                                   style = "font-size: 16px; font-weight: bold; color: #1a0dab;"),
                            tags$br(),
                            tags$i("Legacy web app used to calculate GDD and visualize pistachio growth trends")
                          )
                   )
          ),
          # presentation
          fluidRow(style = "margin-bottom: 30px;",
                   column(2, tags$img(src = "presentation_thumbnail.png",
                                      style = "max-width: 100%; border: 1px solid #000; display: block; margin-top: 5px;")),
                   column(10, style = "display: flex; align-items: center;",
                          tags$div(
                            tags$a("Model Development Presentation",
                                   href = "https://drive.google.com/file/d/1sypxhM79EY5_AZPW68f9SDLbormHRLHF/view?usp=sharing",
                                   target = "_blank",
                                   style = "font-size: 16px; font-weight: bold; color: #1a0dab;"),
                            tags$br(),
                            tags$i("Overview of model design, GDD curve fitting, and app development")
                          )
                   )
          )
        )),
        
        tabPanel("Contact Us",
                 fluidPage(
                   tags$h3("Contact Us"),
                   tags$p("Please fill out the form below to get in touch with the project team."),
                   tags$iframe(
                     src = "https://forms.gle/C4xKwmtddoPLymE6A",
                     width = "100%", height = "800px", frameborder = "0",
                     marginheight = "0", marginwidth = "0", style = "border: none;"
                   )
                 )
        )
      )
    )
  )
}



### UI
ui <- tagList(
  desktop_layout(),
  # keep your exact same footer
  tags$footer(
    style = "background-color: #002855; color: white; padding: 20px 40px; font-family: 'Source Sans Pro', sans-serif;",
    div(
      style = "display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 30px;",
      div(
        style = "display: flex; gap: 20px; align-items: center;",
        tags$a(href = "https://blancolab.ucdavis.edu/", target = "_blank", tags$img(src = "Blancolab.png", height = "105px")),
        tags$a(href = "https://treesystemslab.ucdavis.edu/", target = "_blank", tags$img(src = "marino.png", height = "70px"))
      ),
      div(
        style = "text-align: left; font-size: 14px;",
        tags$b("Department of Plant Sciences:"),
        " General Inquiries - ",
        tags$a(href = "tel:5307520516", style = "color: #ffffff; text-decoration: underline;", "530-752-0516"),
        tags$br(),
        tags$a("University of California, Davis",
               href = "https://www.universityofcalifornia.edu/",
               style = "color: #ffffff; text-decoration: underline; margin-right: 10px;"),
        tags$br(),
        "Developed by ",
        tags$a("Tommy Ngo",
               href = "https://www.linkedin.com/in/tommyngo04/", target = "_blank",
               style = "color: #ffffff; text-decoration: underline;"),
        tags$br(),
        tags$span(tags$b("Funding for the development of this app was provided by the California Pistachio Research Board "))
      ),
      div(
        style = "display: flex; gap: 20px; align-items: center;",
        tags$a(href = "https://calpistachioresearch.org/", target = "_blank", tags$img(src = "cpb.png", height = "95px")),
        tags$a(href = "https://www.plantsciences.ucdavis.edu/", target = "_blank", tags$img(src = "ucd2.png", height = "60px"))
      )
    )
  )
)




### SERVER
server <- function(input, output, session) {
  
  # --- Reactive Values ---
  combined_data <- reactiveVal(NULL) 
  selected_station <- reactiveVal(NULL)
  uploaded_predictions <- reactiveVal(NULL)
  HIST_YEAR_START <- 2010
  HIST_YEAR_END   <- 2020
  coords <- reactiveVal(list(lat = 38.5449, lng = -121.7405))
  
  # Simple Nominatim geocoder
  geocode_address_nominatim <- function(address, limit = 5) {
    resp <- httr::GET(
      url = "https://nominatim.openstreetmap.org/search",
      query = list(q = address, format = "json", limit = limit, addressdetails = 0),
      httr::add_headers(`User-Agent` = "UC-Davis-Pistachio-App/1.0 (contact: your_email@ucdavis.edu)")
    )
    if (httr::status_code(resp) != 200) return(NULL)
    j <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"))
    if (length(j) == 0) return(NULL)
    # Return a compact data.frame we can store
    data.frame(
      label = j$display_name,
      lat   = as.numeric(j$lat),
      lon   = as.numeric(j$lon),
      stringsAsFactors = FALSE
    )
  }
  
  # Store the latest suggestion list
  address_suggestions <- reactiveVal(
    data.frame(label = character(), lat = numeric(), lon = numeric())
  )
  
  # Debounce the typed term to avoid hammering the API
  typed_term <- reactive(input$address_type)
  typed_term_deb <- debounce(typed_term, 400)
  
  
  
  # --- Leaflet Map Initialization and Interactivity ---
  DEFAULT_LAT <- 38.5449
  DEFAULT_LON <- -121.7405
  
  output$map <- renderLeaflet({
    leaflet() %>%
      addProviderTiles("CartoDB.Voyager") %>%
      setView(lng = DEFAULT_LON, lat = DEFAULT_LAT, zoom = 10) %>%
      # initial pointer so there’s always a marker even before inputs exist
      addMarkers(lng = DEFAULT_LON, lat = DEFAULT_LAT, layerId = "selected")
  })
  

  observe({
    # don’t try to add a marker until inputs are real numbers
    req(!is.null(input$longitude), !is.null(input$latitude))
    leafletProxy("map") %>%
      clearMarkers() %>%
      addMarkers(lng = input$longitude, lat = input$latitude, layerId = "selected")
  })
  
  observeEvent(input$map_click, {
    lat_rounded <- round(input$map_click$lat, 5)
    lon_rounded <- round(input$map_click$lng, 5)
    updateNumericInput(session, "latitude",  value = lat_rounded)
    updateNumericInput(session, "longitude", value = lon_rounded)
    
    leafletProxy("map") %>%
      clearMarkers() %>%
      addMarkers(lng = lon_rounded, lat = lat_rounded, layerId = "selected")
  })
  
  # 1) When the user types, fetch suggestions and update the dropdown
  observeEvent(typed_term_deb(), {
    req(typed_term_deb())
    q <- typed_term_deb()
    out <- tryCatch(geocode_address_nominatim(q, limit = 5), error = function(e) NULL)
    if (is.null(out) || nrow(out) == 0) {
      address_suggestions(data.frame(label = character(), lat = numeric(), lon = numeric()))
      updateSelectizeInput(session, "address_pick", choices = NULL, server = TRUE)
      return()
    }
    address_suggestions(out)
    
    # Show labels in the dropdown. Value is the label itself for easy lookup.
    updateSelectizeInput(
      session, "address_pick",
      choices = setNames(out$label, out$label),
      selected = NULL,
      server = TRUE
    )
  })
  
  # 2) When a suggestion is selected, paste lat/lon into inputs and update marker
  observeEvent(input$address_pick, {
    req(input$address_pick)
    df <- address_suggestions()
    if (nrow(df) == 0) return()
    row <- df[df$label == input$address_pick, , drop = FALSE]
    if (nrow(row) != 1 || any(is.na(c(row$lat, row$lon)))) return()
    
    updateNumericInput(session, "latitude",  value = round(row$lat, 6))
    updateNumericInput(session, "longitude", value = round(row$lon, 6))
    
    leafletProxy("map") %>%
      clearMarkers() %>%
      addMarkers(lng = row$lon, lat = row$lat, layerId = "selected")
  })
  
  
  # Legends at the bottom for mobile
  legend_bottom_if_mobile <- function(p) {
    p %>% layout(
      legend = list(
        orientation = "h",         # horizontal
        y = -0.3,                  # below plot area
        x = 0.5,                   # centered
        xanchor = "center",
        yanchor = "top"
      )
    )
  }
  
  
  
  # --- Station Info Rendering (used in multiple places) ---
  render_station_and_download_links <- function() {
    req(combined_data())
    req(selected_station())
    
    station_link <- paste0(
      "https://viewer.synopticdata.com/map/data/202506191823/air-temperature/",
      selected_station(),
      "/plots/temperature?layers=#map=7.44/38.32/-121.26&networks=66"
    )
    
    tags$div(
      style = "margin-top: 10px; display: flex; justify-content: center; gap: 50px; align-items: center;",
      
      tags$a(
        href = station_link,
        target = "_blank",
        style = "font-size: 15px; color: #0066cc; text-decoration: underline;",
        paste0("View CIMIS Station ", selected_station())
      ),
      
      downloadLink(
        outputId = "download_data",
        label = "Download Predictions (CSV)",
        class = "btn-link",
        style = "font-size: 15px; color: #0066cc; text-decoration: underline;"
      )
    )
  }
  
  
  # --- Download Links and Station Info ---
  output$links_live <- renderUI(render_station_and_download_links())
  output$links_hist <- renderUI(render_station_and_download_links())
  
  output$station_info_live <- renderUI({
    req(selected_station())
    station_link <- paste0(
      "https://viewer.synopticdata.com/map/data/202506191823/air-temperature/",
      selected_station(),
      "/plots/temperature?layers=#map=7.44/38.32/-121.26&networks=66"
    )
    tags$div(
      style = "margin-top: 10px; font-style: italic; color: #666666; font-size: 14px;",
      "Data source: ",
      tags$a(href = station_link, target = "_blank", paste0("CIMIS Station ", selected_station()))
    )
  })
  
  output$station_info_hist <- renderUI({
    req(selected_station())
    station_link <- paste0(
      "https://viewer.synopticdata.com/map/data/202506191823/air-temperature/",
      selected_station(),
      "/plots/temperature?layers=#map=7.44/38.32/-121.26&networks=66"
    )
    tags$div(
      style = "margin-top: 10px; font-style: italic; color = #666666; font-size: 14px;",
      "Data source: ",
      tags$a(href = station_link, target = "_blank", paste0("CIMIS Station ", selected_station()))
    )
  })
  
  # --- Main Data Fetch and Processing Workflow ---
  observeEvent(input$get_data, {
    
    progress <- Progress$new()
    on.exit(progress$close())
    progress$set(message = "Starting model generation...", value = 0)
    
    lat <- input$latitude
    lon <- input$longitude
    start <- as.Date(input$start_date)
    end <- as.Date(input$end_date)
    
    if (end < start) {
      showNotification("Error: End date cannot be earlier than Start date.", type = "error", duration = 5)
      return()
    }
    
    # Stage 1: Fetch Station Metadata/ Closest Station
    progress$inc(0.1, detail = "Fetching station metadata...")
    
    station_response <- GET("https://api.synopticdata.com/v2/stations/metadata", query = list(
      radius = paste(lat, lon, 100, sep = ","), network = 66, token = token))
    
    if (station_response$status_code != 200) {
      showNotification("Error: Failed to retrieve station metadata.", type = "error", duration = 5)
      return()
    }
    
    station_data <- content(station_response, as = "parsed", simplifyVector = TRUE)
    if (length(station_data$STATION) == 0) {
      showNotification("Error: No stations found nearby.", type = "error")
      return()
    }
    stid <- station_data$STATION$STID[1]
    selected_station(stid)
    print(stid)
    
    # Stage 2: Fetch Live Weather Data
    progress$inc(0.3, detail = "Fetching live station weather data...")
    
    response <- GET(
      "https://api.synopticdata.com/v2/stations/timeseries",
      query = list(
        stid = stid,
        start = format(start, "%Y%m%d0000"),
        end   = format(end,   "%Y%m%d2359"),
        vars  = "air_temp",
        units = "temp|F",
        obtimezone = "UTC",
        token = token
      )
    )
    
    if (response$status_code != 200) {
      showNotification("Error: Failed to retrieve station weather data.", type = "error", duration = 5)
      return()
    }
    
    raw <- content(response, as = "text", encoding = "UTF-8")
    parsed <- jsonlite::fromJSON(raw, flatten = TRUE)
    
    # Pull out times & temps safely
    times <- tryCatch(parsed$STATION$OBSERVATIONS.date_time[[1]], error = function(e) NULL)
    temps <- tryCatch(parsed$STATION$OBSERVATIONS.air_temp_set_1[[1]], error = function(e) NULL)
    
    # Guard: no observations returned
    if (is.null(times) || is.null(temps) || length(times) == 0 || length(temps) == 0) {
      showNotification(
        "No air temperature observations returned for this station/time window. Try a different end date or nearby station.",
        type = "error", duration = 7
      )
      return()
    }
    
    # Build live dataframe robustly
    df_live <- tibble::tibble(
      datetime    = as.POSIXct(times, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      temperature = suppressWarnings(as.numeric(temps))
    ) %>%
      dplyr::filter(!is.na(datetime), !is.na(temperature)) %>%
      dplyr::mutate(date = as.Date(datetime)) %>%
      dplyr::group_by(date) %>%
      dplyr::summarise(
        high_temp = max(temperature, na.rm = TRUE),
        low_temp  = min(temperature, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Guard: still empty after filtering
    if (nrow(df_live) == 0) {
      showNotification(
        "No valid temperature observations found after filtering. Try widening the date range.",
        type = "error", duration = 7
      )
      return()
    }
    
    df_live <- df_live %>%
      dplyr::mutate(
        GDD = pmax((high_temp + low_temp) / 2 - 44.6, 0),
        GDD_cumulative = cumsum(GDD),
        source = "Live"
      ) %>%
      apply_all_traits()
    
    print(df_live)
    
    # Stage 3: Forecast (if within 5 days from today)
    forecast_df <- NULL
    if (as.numeric(Sys.Date() - end) <= 5) {
      progress$inc(0.5, detail = "Fetching forecast data...")
      
      forecast_url <- "https://api.open-meteo.com/v1/forecast"
      forecast_params <- list(
        latitude = lat, longitude = lon,
        daily = "temperature_2m_max,temperature_2m_min",
        temperature_unit = "fahrenheit", timezone = "auto"
      )
      forecast_response <- GET(url = forecast_url, query = forecast_params)
      
      if (forecast_response$status_code == 200) {
        forecast_data <- fromJSON(content(forecast_response, "text", encoding = "UTF-8"))
        forecast_df <- data.frame(
          date = as.Date(forecast_data$daily$time),
          high_temp = forecast_data$daily$temperature_2m_max,
          low_temp = forecast_data$daily$temperature_2m_min
        ) %>%
          mutate(GDD = pmax((high_temp + low_temp) / 2 - 44.6, 0))
        
        last_cum_gdd <- tail(df_live$GDD_cumulative, 1)
        forecast_df <- forecast_df %>%
          mutate(GDD_cumulative = last_cum_gdd + cumsum(GDD),
                 source = "Forecast") %>%
          apply_all_traits()
      }
    }
    
    # Stage 4: Historical Data from Open-Meteo using Julian Date
    progress$inc(0.75, detail = "Fetching historical Open-Meteo data...")
    
    base_temp_f <- 44.6
    align_year <- lubridate::year(start)
    
    # Pull daily history for the chosen window, then average by day-of-year
    hist_daily <- fetch_openmeteo_hist(
      lat = lat, lon = lon,
      year_start = HIST_YEAR_START, year_end = HIST_YEAR_END,
      base_temp_f = base_temp_f
    )
    
    print(hist_daily)
    
    start_jday <- lubridate::yday(start)
    
    gdd_hist_avg <- hist_daily |>
      dplyr::group_by(jday) |>
      dplyr::summarise(GDD = mean(GDD, na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(jday) |>
      dplyr::filter(jday >= start_jday) |>
      dplyr::mutate(
        GDD_cumulative = cumsum(GDD),
        date = as.Date(jday - 1, origin = as.Date(sprintf("%d-01-01", align_year))),
        source = "Historical"
      ) |>
      apply_all_traits()
    
    
    # Stage 5: Combine All Sources
    progress$inc(0.95, detail = "Combining predictions and generating output...")
    combined <- bind_rows(df_live, forecast_df, gdd_hist_avg)
    
    # Store combined data
    combined_data(combined)
    
    
    # Stage 6: Plot Outputs
    ### Live
    output$growth_plot_live <- renderPlotly({
      req(combined_data())
      selected_data <- combined_data() %>% dplyr::filter(source %in% c("Live", "Forecast"))
      
      p <- if (input$model_var %in% c("areas_combined","textures_combined")) {
        generate_plotly_grouped(selected_data, input$model_var, x_axis_type = input$x_axis_type)
      } else {
        generate_plotly(selected_data, input$model_var, x_axis_type = input$x_axis_type)
      }
      
      # If the container is narrow (<768px), push legend to bottom and double the width
      w <- session$clientData$output_growth_plot_live_width
      if (!is.null(w) && w < 768) {
        p <- p %>% layout(
          width = w * 2,                          # double width on mobile
          legend = list(orientation = "h",
                        x = 0.5, xanchor = "center",
                        y = -0.25, yanchor = "top"),
          margin = list(b = 110)
        )
      }
      p
    })
    
    
    
    
    output$live_note <- renderUI({
      req(combined_data())
      used_forecast <- any(combined_data()$source == "Forecast")
      msg <- if (used_forecast) {
        tags$span(
          "Note: Real-time predictions use ",
          tags$a("CIMIS", href = "https://cimis.water.ca.gov/", target = "_blank"),
          " station observations for the selected station. Short-term forecasts are from ",
          tags$a("Open-Meteo", href = "https://open-meteo.com/", target = "_blank"), "."
        )
      } else {
        tags$span(
          "Note: Real-time predictions use ",
          tags$a("CIMIS", href = "https://cimis.water.ca.gov/", target = "_blank"),
          " station observations for the selected station. No forecast was added for this run."
        )
      }
      tags$div(class = "note-under-plot", msg)
    })
    
    
    ### History
    # --- HISTORICAL ---
    output$growth_plot_hist <- renderPlotly({
      req(combined_data())
      selected_data <- combined_data() %>% dplyr::filter(source == "Historical")
      
      custom_title <- paste0(
        "Predicted ",
        if (input$model_var == "areas_combined") "Areas"
        else if (input$model_var == "textures_combined") "Textures"
        else tools::toTitleCase(gsub("_", " ", input$model_var)),
        " (", HIST_YEAR_START, "–", HIST_YEAR_END, ")"
      )
      
      p <- if (input$model_var %in% c("areas_combined","textures_combined")) {
        generate_plotly_grouped(selected_data, input$model_var, x_axis_type = input$x_axis_type, custom_title = custom_title)
      } else {
        generate_plotly(selected_data, input$model_var, x_axis_type = input$x_axis_type, custom_title = custom_title)
      }
      
      p <- p %>% layout(
        xaxis = list(
          title = if (input$x_axis_type == "GDD_cumulative") "Cumulative GDD" else "Date",
          tickformat = if (input$x_axis_type == "date") "%b %d" else NULL,
          tickangle  = if (input$x_axis_type == "date") -45 else 0
        )
      )
      
      # If the container is narrow (<768px), push legend to bottom and double the width
      w <- session$clientData$output_growth_plot_hist_width
      if (!is.null(w) && w < 768) {
        p <- p %>% layout(
          width = w * 2,                          # double width on mobile
          legend = list(orientation = "h",
                        x = 0.5, xanchor = "center",
                        y = -0.25, yanchor = "top"),
          margin = list(b = 110)
        )
      }
      p
    })
    
    
    output$hist_note <- renderUI({
      tags$div(
        class = "note-under-plot",
        tags$span(
          "Note: These predictions use averaged historical weather data from ",
          HIST_YEAR_START, " to ", HIST_YEAR_END, " via ",
          tags$a("Open-Meteo", href = "https://open-meteo.com/", target = "_blank"), "."
        )
      )
    })
    
    
    
    ### Combined
    output$growth_plot_combined <- renderPlotly({
      req(combined_data())
      selected_data <- combined_data()
      
      p <- if (input$model_var %in% c("areas_combined","textures_combined")) {
        generate_plotly_grouped(
          selected_data,
          input$model_var,
          x_axis_type = input$x_axis_type
        )
      } else {
        generate_plotly(
          selected_data,
          input$model_var,
          x_axis_type = input$x_axis_type
        )
      }
      
      # Detect mobile width
      w <- session$clientData$output_growth_plot_combined_width
      if (!is.null(w) && w < 768) {
        p <- p %>% layout(
          width = w * 2,  # double the actual detected width
          legend = list(
            orientation = "h",
            x = 0.5, xanchor = "center",
            y = -0.25, yanchor = "top"
          ),
          margin = list(b = 110)
        )
      }
      
      p
    })
    
    
    progress$inc(1, detail = "Complete!")
  })
  
  # Grouped graphs
  generate_plotly_grouped <- function(df, group, x_axis_type = "date", custom_title = NULL) {
    if (group == "areas_combined") {
      vars <- c("kernel_area", "nut_area")
      labels <- c(kernel_area = "Kernel Area", nut_area = "Nut Area")
      unit <- "(mm²)"
      # use same color palette style
      var_colors <- c(kernel_area = "#3a3b3c",  # matches L*
                      nut_area    = "#D55E00")  # matches A*
      title_default <- "Predicted Areas"
    } else if (group == "textures_combined") {
      vars <- c("kernel_texture", "shell_texture", "hull_texture")
      labels <- c(kernel_texture = "Kernel Texture",
                  shell_texture  = "Shell Texture",
                  hull_texture   = "Hull Texture")
      unit <- "(kg force)"
      # use same palette style
      var_colors <- c(kernel_texture = "#3a3b3c",  # L*
                      shell_texture  = "#D55E00",  # A*
                      hull_texture   = "#0072B2")  # B*
      title_default <- "Predicted Textures"
    } else {
      stop("Unknown group for grouped plot")
    }
    
    line_style <- list(
      "Live"      = list(width = 2, dash = "solid"),
      "Forecast"  = list(width = 2, dash = "dash"),
      "Historical"= list(width = 2, dash = "dot")
    )
    
    p <- plotly::plot_ly()
    for (src in unique(df$source)) {
      display_src <- if (src == "Historical") "Projected" else src   # <— add this
      temp_df <- df %>% dplyr::filter(source == src)
      for (v in vars) {
        hover_text <- paste0(
          "Date: ", format(temp_df$date, "%b %d"),
          "<br>GDD: ", round(temp_df$GDD_cumulative, 1),
          "<br>", labels[[v]], ": ", round(temp_df[[v]], 2), " ", unit
        )
        p <- p %>%
          plotly::add_lines(
            data = temp_df,
            x = ~get(x_axis_type),
            y = temp_df[[v]],
            name = paste0(display_src, " - ", labels[[v]]),  # <— use display_src here
            text = hover_text,
            hoverinfo = "text",
            line = c(line_style[[src]], list(color = var_colors[[v]]))  # keep styling keyed to src
          )
      }
    }
    
    p %>%
      plotly::layout(
        title = list(text = if (!is.null(custom_title)) custom_title else title_default, y = 0.95),
        xaxis = list(title = if (x_axis_type == "GDD_cumulative") "Cumulative GDD" else "Date",
                     titlefont = list(size = 16), tickfont = list(size = 14)),
        yaxis = list(title = paste0(if (group == "areas_combined") "Area " else "Texture ", unit),
                     titlefont = list(size = 16), tickfont = list(size = 14)),
        legend = list(font = list(size = 14)),
        margin = list(l = 60, r = 130, t = 50, b = 60),
        showlegend = TRUE
      ) %>%
      no_zoom()
  }
  
  # Notes for combined legends
  output$combined_note <- renderUI({
    style_note <- tags$span(
      tags$b("Line style: "),
      "solid = Live (CIMIS), dashed = Forecast (Open-Meteo), dotted = Projected (historical average, Open-Meteo)."
    )
    color_note <- switch(
      input$model_var,
      "areas_combined" = tags$span(
        tags$b(" Colors: "),
        HTML("<span style='color:#3a3b3c;'>■</span> Kernel Area, "),
        HTML("<span style='color:#D55E00;'>■</span> Nut Area.")
      ),
      "textures_combined" = tags$span(
        tags$b(" Colors: "),
        HTML("<span style='color:#3a3b3c;'>■</span> Kernel Texture, "),
        HTML("<span style='color:#D55E00;'>■</span> Shell Texture, "),
        HTML("<span style='color:#0072B2;'>■</span> Hull Texture.")
      ),
      "kernel_color" = tags$span(
        tags$b(" Colors: "),
        HTML("<span style='color:#3a3b3c;'>■</span> L* (lightness), "),
        HTML("<span style='color:#D55E00;'>■</span> a* (green↔red), "),
        HTML("<span style='color:#0072B2;'>■</span> b* (blue↔yellow).")
      ),
      "nut_color" = tags$span(
        tags$b(" Colors: "),
        HTML("<span style='color:#3a3b3c;'>■</span> L* (lightness), "),
        HTML("<span style='color:#D55E00;'>■</span> a* (green↔red), "),
        HTML("<span style='color:#0072B2;'>■</span> b* (blue↔yellow).")
      ),
      tags$span(tags$b(" Colors: "), "Single variable shown (color indicates that variable).")
    )
    tags$div(class = "note-under-plot", style_note, tags$br(), color_note)
  })
  
  
  
  
  # Picture in combined colors
  output$reference_image_ui <- renderUI({
    req(input$model_var)
    
    if (input$model_var == "kernel_color") {
      tags$div(
        style = "padding-left: 90px; margin-top: 20px;",
        tags$img(
          src = "kernel_pic.png",
          style = "height: 120px; width: auto;"
        )
      )
    } else if (input$model_var == "nut_color") {
      tags$div(
        style = "padding-left: 105px; margin-top: 20px;",
        tags$img(
          src = "hull_pic.png",
          style = "height: 130px; width: auto;"
        )
      )
    } else {
      NULL
    }
  })
  
  # LAB color space
  observeEvent(input$show_lab_modal_link, {
    showModal(modalDialog(
      title = "Understanding the LAB Color Space",
      tags$img(src = "lab_colorspace.png", style = "width: 100%; height: auto; margin-bottom: 10px;"),
      tags$ul(
        tags$li(tags$b("L*:"), " Lightness — 0 (black) to 100 (white)"),
        tags$li(tags$b("a*:"), " Green (–) to Red (+)"),
        tags$li(tags$b("b*:"), " Blue (–) to Yellow (+)")
      ),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  
  # --- Upload & Render GDD Predictions ---
  observeEvent(input$process_file, {
    req(input$gdd_file)
    
    # Read and normalize headers
    df <- read.csv(input$gdd_file$datapath, check.names = FALSE, stringsAsFactors = FALSE)
    orig_names <- names(df)
    clean_names <- trimws(orig_names)
    names(df) <- clean_names
    
    # Allow a few common variants for the column name
    possible_gdd_names <- c("GDD_cumulative","gdd_cumulative","cumulative_gdd","GDD","gdd")
    gdd_col <- intersect(possible_gdd_names, names(df))
    if (length(gdd_col) == 0) {
      showNotification("CSV must include a 'GDD_cumulative' column.", type = "error", duration = 6)
      return()
    }
    gdd_col <- gdd_col[1]
    names(df)[names(df) == gdd_col] <- "GDD_cumulative"
    
    # Coerce to numeric, drop NA/negative, warn if fixed
    df$GDD_cumulative <- suppressWarnings(as.numeric(df$GDD_cumulative))
    if (any(is.na(df$GDD_cumulative))) {
      showNotification("Non-numeric values in GDD_cumulative were removed.", type = "warning", duration = 5)
    }
    df <- df[!is.na(df$GDD_cumulative) & df$GDD_cumulative >= 0, , drop = FALSE]
    
    # Optional date parsing, flexible formats
    if ("date" %in% names(df)) {
      dtry <- suppressWarnings(lubridate::ymd(df$date))
      dtry2 <- suppressWarnings(lubridate::mdy(df$date))
      # Allow MM/DD without year by assuming the current year
      dtry3 <- suppressWarnings({
        no_year <- grepl("^[0-1]?[0-9]/[0-3]?[0-9]$", df$date)
        out <- as.Date(NA)
        out[no_year] <- as.Date(paste0(format(Sys.Date(), "%Y"), "/", df$date[no_year]), format = "%Y/%m/%d")
        out
      })
      # Pick first successful parse
      parsed <- ifelse(!is.na(dtry), dtry,
                       ifelse(!is.na(dtry2), dtry2, dtry3))
      if (all(is.na(parsed))) {
        showNotification("Could not parse the 'date' column. It will be ignored.", type = "warning", duration = 5)
        df$date <- NULL
      } else {
        df$date <- as.Date(parsed, origin = "1970-01-01")
      }
    }
    
    if (nrow(df) < 2) {
      showNotification("Not enough rows to generate predictions.", type = "error", duration = 5)
      return()
    }
    
    # Monotonic check
    if (any(diff(df$GDD_cumulative) < 0, na.rm = TRUE)) {
      showNotification("GDD_cumulative must be non-decreasing by row. Please sort your file.", type = "error", duration = 7)
      return()
    }
    
    # Hard cap to your global max
    if (max(df$GDD_cumulative, na.rm = TRUE) > GDD_MAX) {
      showNotification(paste0("Values above ", GDD_MAX, " were trimmed."), type = "warning", duration = 5)
      df <- df[df$GDD_cumulative <= GDD_MAX, , drop = FALSE]
    }
    
    # Ensure an x-column for plotting by GDD
    if (!"GDD_cumulative" %in% names(df)) {
      showNotification("Internal error: GDD_cumulative missing after checks.", type = "error", duration = 5)
      return()
    }
    
    # Run predictions
    df_predicted <- df %>% apply_all_traits()
    uploaded_predictions(df_predicted)
  })
  
  
  output$gdd_table <- DT::renderDataTable({
    req(uploaded_predictions())
    DT::datatable(uploaded_predictions(), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  output$uploaded_var_selector <- renderUI({
    req(uploaded_predictions())
    selectInput("uploaded_model_var", "Select Variable to Plot:",
                choices = c("Kernel Area" = "kernel_area", 
                            "Nut Area" = "nut_area",
                            "Kernel Texture" = "kernel_texture", 
                            "Shell Texture" = "shell_texture",
                            "Hull Texture" = "hull_texture", 
                            "Kernel Dry Weight" = "kernel_dry_weight",
                            "Hull Color (L*, A*, B*)" = "nut_color", 
                            "Kernel Color (L*, A*, B*)" = "kernel_color")
    )
  })
  
  # Graphs for uploaded GDD
  output$uploaded_growth_plot <- renderPlotly({
    req(uploaded_predictions())
    df <- uploaded_predictions() %>% arrange(GDD_cumulative)
    selected_var <- if (!is.null(input$uploaded_model_var)) input$uploaded_model_var else "nut_area"
    
    if (selected_var %in% c("nut_color", "kernel_color")) {
      color_prefix <- if (selected_var == "nut_color") "nut" else "kernel"
      color_values <- c("L" = "#3a3b3c", "A" = "#D55E00", "B" = "#0072B2")
      dash_map <- c("L" = "solid", "A" = "dash", "B" = "dot")
      p <- plot_ly()
      for (col in c("L", "A", "B")) {
        trait_col <- paste0(color_prefix, "_", col)
        label_map <- c("L" = "Lightness (L*)", "A" = "Green-Red (a*)", "B" = "Blue-Yellow (b*)")
        p <- p %>% add_lines(
          data = df,
          x = ~GDD_cumulative,
          y = df[[trait_col]],
          name = label_map[[col]],
          text = ~paste0("GDD: ", round(GDD_cumulative, 1),
                         "<br>", label_map[[col]], ": ", round(get(trait_col), 2)),
          line = list(width = 2, dash = dash_map[[col]], color = color_values[[col]]),
          hoverinfo = "text"
        )
      }
      return(
        p %>%
          layout(
            title = paste("Predicted", if (color_prefix == "nut") "Hull" else "Kernel", "Color"),
            xaxis = list(title = "Cumulative GDD"),
            yaxis = list(title = "Color Value"),
            margin = list(t = 50),
            showlegend = TRUE
          ) %>%
          no_zoom()
      )
    }
    
    plot_ly(
      data = df,
      x = ~GDD_cumulative,
      y = as.formula(paste0("~", selected_var)),
      type = 'scatter',
      mode = 'lines+markers',
      text = ~paste0("GDD: ", round(GDD_cumulative, 1),
                     "<br>", selected_var, ": ", round(get(selected_var), 2)),
      hoverinfo = 'text'
    ) %>%
      layout(
        title = paste("Predicted", selected_var),
        xaxis = list(title = "Cumulative GDD"),
        yaxis = list(title = selected_var)
      ) %>%
      no_zoom()
  })
  
  # --- Download Handlers for Predictions ---
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("pistachio_predictions_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(combined_data())
      df <- combined_data()
      
      # Round numeric values
      num_cols <- sapply(df, is.numeric)
      df[num_cols] <- lapply(df[num_cols], function(x) round(x, 4))
      
      # Remove year for historical data
      df$date <- ifelse(
        df$source == "Historical",
        format(df$date, "%b %d"),
        as.character(df$date)
      )
      
      # Rename color columns for consistency
      names(df) <- names(df) %>%
        gsub("^nut_L$", "hull_L", .) %>%
        gsub("^nut_A$", "hull_A", .) %>%
        gsub("^nut_B$", "hull_B", .)
      
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  # --- Download Button for Uploaded Predictions ---
  output$download_button_ui <- renderUI({
    if (is.null(uploaded_predictions())) {
      actionButton("download_warning", "Download Predictions", class = "btn btn-primary", style = "color: white;")
    } else {
      downloadButton("download_gdd_predictions", "Download Predictions", class = "btn btn-primary", style = "color: white;")
    }
  })
  
  observeEvent(input$download_warning, {
    showNotification("Please generate predictions first before downloading.", type = "error", duration = 4)
  })
  
  # --- Download Handlers for Uploaded Predictions ---
  output$download_gdd_predictions <- downloadHandler(
    filename = function() {
      paste0("uploaded_gdd_predictions_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(uploaded_predictions())
      df <- uploaded_predictions()
      num_cols <- sapply(df, is.numeric)
      df[num_cols] <- lapply(df[num_cols], function(x) round(x, 4))
      names(df) <- names(df) %>%
        gsub("^nut_L$", "hull_L", .) %>%
        gsub("^nut_A$", "hull_A", .) %>%
        gsub("^nut_B$", "hull_B", .)
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  output$download_template <- downloadHandler(
    filename = function() "gdd_template.csv",
    content = function(file) {
      example <- data.frame(
        date = c("04/04","04/05","04/06","04/07","04/08"),
        GDD_cumulative = c(0, 12.4, 20.7, 27.1, 34.9)
      )
      write.csv(example, file, row.names = FALSE)
    }
  )
}






shinyApp(ui, server)

# %maxes, phone layout