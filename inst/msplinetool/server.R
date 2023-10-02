render_dt = function(data, ...) {
  renderDT(data, selection = 'none', server = TRUE,
           editable = "column", options=list(info=FALSE, dom="t"), ...)
}

function(input, output, session) {

  ## default spline on startup 
  m <- reactiveValues()
  md <- msplinemodel_init()
  for (i in names(md)) m[[i]] <- md[[i]] # TODO nicer?
  m$knots_df <- data.frame(knots=md$knots)
  m$coefs_df <- data.frame(coefs=md$coefs)

  ## main plot 
  output$main_plot <- renderPlot({
    plot_mspline(knots = m$knots,
                 bknot = 10,
                 df = m$df,
                 degree = m$degree,
                 bsmooth = m$bsmooth,
                 coefs = m$coefs)
  })

  ## table of knots and coefficients (updated if spline changes)
  observeEvent(m$coefs_df, {
    output$coefs_df <- render_dt(m$coefs_df)
  })

  observeEvent(m$knots_df, {
    output$knots_df <- render_dt(m$knots_df)
  })
  
  ## make spline specification consistent if one component changes

  update_spline <- function(m, ...){
    an <- function(x){if(is.null(x)) NULL else as.numeric(x)}
    args <- list(...)
    for (i in names(args)) m[[i]] <- args[[i]]
    md <- msplinemodel_init(df = an(m$df), degree = an(m$degree),
                            bsmooth = m$bsmooth, knots = an(m$knots), coefs=an(m$coefs))
    for (i in names(md)) m[[i]] <- md[[i]]
    m$knots_df <- data.frame(knots=m$knots)
    m$coefs_df <- data.frame(coefs=m$coefs)
  }

  ## monitor if user changes spline components

  observeEvent(input$bsmooth, {
    degree_opts <- if (input$bsmooth) 3 else selection_defaults$degree_unsmooth
    df_opts <- selection_defaults$df
    update_spline(m, bsmooth = input$bsmooth, knots=NULL, coefs=NULL)
    updateSelectInput(session, "degree", choices=degree_opts, selected=3)
    updateSelectInput(session, "df", choices=df_opts, selected=m$df)
  })

  observeEvent(input$degree, {
    update_spline(m, degree = input$degree, knots=NULL, coefs=NULL)
  })

  observeEvent(input$df, {
    update_spline(m, df = input$df, knots=NULL, coefs=NULL)
  })

  observeEvent(input$knots_df_cell_edit, {
    m$knots_df <- editData(m$knots_df, input$knots_df_cell_edit)
    m$knots <- m$knots_df$knots
  })

  observeEvent(input$coefs_df_cell_edit, {
    m$coefs_df <- editData(m$coefs_df, input$coefs_df_cell_edit)
    m$coefs <- m$coefs_df$coefs
  })

}

