library(shiny)
library(DT)

dt_output = function(title, id) {
  fluidRow(column(
    12, h1(paste0('Table ', sub('.*?([0-9]+)$', '\\1', id), ': ', title)),
    hr(), DTOutput(id)
  ))
}
render_dt = function(data, editable = 'cell', server = TRUE, ...) {
  renderDT(data, selection = 'none', server = server, editable = editable, ...)
}

shinyApp(
  ui = fluidPage(
    title = 'Double-click to edit table cells',
    dt_output('client-side processing (editable = "all")', 'x4'),
    dt_output('server-side processing (editable = "all")', 'x8'),
  ),

  server = function(input, output, session) {
    d1 = iris
    d1$Date = Sys.time() + seq_len(nrow(d1))
    d7 = d1

    options(DT.options = list(pageLength = 5))

    # server-side processing
    output$x7 = render_dt(d7, 'all')
    observeEvent(input$x7_cell_edit, {
      d8 <<- editData(d8, input$x8_cell_edit, 'x8')
    })

  }
)
