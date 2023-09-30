
dt_output = function(title, id) {
  fluidRow(column(
    12, h1(paste0('Table ', sub('.*?([0-9]+)$', '\\1', id), ': ', title)),
    hr(), DTOutput(id)
  ))
}

bootstrapPage(

  selectInput(inputId = "df",
      label = "Knots df",
      choices = 5:12,
      selected = mspline_defaults$df),

  selectInput(inputId = "degree",
      label = "Polynomial degree",
      choices = c(2,3,4),
      selected = mspline_defaults$degree),

  checkboxInput(inputId = "bsmooth",
      label = "Boundary smoothness",
      value = mspline_defaults$bsmooth),

  dt_output('Knots', 'knots_df'),
  dt_output('Coefs', 'coefs_df'),

  plotOutput(outputId = "main_plot", height = "300px"),

)

## TODO

## coefficients 

## knot positions 

## why no x axis? 
## y axis scaling
## line colors, widths
## ggplot code to paste
