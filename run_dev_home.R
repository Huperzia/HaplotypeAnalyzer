library('shiny')
#library(reactlog)

#reactlog_enable()

#profvis::profvis(

  runApp(
  appDir = 'app',
  port = 7338,
  launch.browser = FALSE,
  host = getOption("shiny.host", "0.0.0.0"),
  workerId = "",
  quiet = FALSE,
  display.mode = "normal", 
  test.mode = getOption("shiny.testmode", FALSE)
)
#)

#shiny::reactlogShow()

  

  
