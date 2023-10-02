
dt_output = function(title, id) {
  fluidRow(column(
    12, h1(paste0('Table ', sub('.*?([0-9]+)$', '\\1', id), ': ', title)),
    hr(), DTOutput(id)
  ))
}

bootstrapPage(
  plotOutput(outputId = "main_plot", height = "300px"),
  fluidRow(
    column(3,
           selectInput(inputId = "df",
                       label = "Knots df",
                       choices = selection_defaults$df,
                       selected = mspline_defaults$df),

           selectInput(inputId = "degree",
                       label = "Polynomial degree",
                       choices = 3,
                       selected = mspline_defaults$degree),

           checkboxInput(inputId = "bsmooth",
                         label = "Boundary smoothness",
                         value = mspline_defaults$bsmooth)
           ),
    column(4,
           DTOutput('knots_df',width="25%")
           ),
    column(4,
           DTOutput('coefs_df',width="25%")
           )
  )  
)

## TODO
## format
## why no x axis? 
## y axis scaling
## line colors, widths
## ggplot code to paste
## algebra explanations and language for model
