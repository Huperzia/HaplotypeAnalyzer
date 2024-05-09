#################
# Load packages #
#################

library(shiny)
library(bs4Dash)
library(waiter)
library(fresh)
library(stringi)


##############
# Aesthetics #
##############

#thematic_shiny()

drac <- c("#50fa7b", "#ffb86c", "#bd93f9", "#ff79c6", 
          "#ff5555", "#f1fa8c", "#6272a4", "#8be9fd", 
          "#282a36", "#44475a", "#44475a", "#f8f8f2")

cols <- drac

mytheme <- create_theme(
  bs4dash_vars(
    navbar_dark_color = drac[9],
    navbar_dark_active_color = drac[12],
    navbar_dark_hover_color = drac[12],
    navbar_light_color = drac[12],
    navbar_light_active_color = "#FFF",
    navbar_light_hover_color = "#FFF"
  ),
  bs4dash_layout(sidebar_width = "300px"),
  bs4dash_yiq(
    contrasted_threshold = 10,
    text_dark = "#FFF",
    text_light = "#000000"
  ),

  bs4dash_sidebar_dark(
    bg = drac[9], 
    color = drac[12],
    hover_color = "#FFF",
    hover_bg = drac[9],
    submenu_bg = drac[9], 
    submenu_hover_bg = drac[9],
    submenu_active_color = drac[10],
    submenu_active_bg = drac[10],
    header_color = drac[9],
    submenu_color = "#FFF",
    submenu_hover_color = "#FFF"
  ),
  bs4dash_status(
    primary = drac[10], 
    secondary = drac[9], 
    danger = "#BF616A", 
    light = drac[9], 
    dark = drac[9]
  ),
  bs4dash_color(
        white = drac[9],
        gray_900 = drac[9],
    fuchsia = drac[9],
    gray_800 = drac[9],
    gray_x_light = drac[9]
  )
)


#~~~~~~~~~~~



#############
# Tab Items #
#############

# Population Overview
tab_population_overview_samples <- tabItem(
  tabName = "population_overview_samples",    
  fluidRow(
    box(
      id = "card2",
      title = "Haplotype Overview", 
      width = 12,
      height = "1200",
      background = "fuchsia",
      maximizable = TRUE, 
      collapsible = TRUE,
      
      plotOutput("population_overview_samples", height = "100%")
    )
  )
)

tab_population_overview_samples2 <- tabItem(
  tabName = "population_overview_samples2",    
  fluidRow(
    box(
      id = "card2",
      title = "Population Overview", 
      width = 12,
      height = "9900",
      background = "fuchsia",
      maximizable = TRUE, 
      collapsible = TRUE,
      
      plotOutput("population_overview_samples2", height = "10000px")
    )
  )
)

# Population Overview
tab_population_overview_samples3 <- tabItem(
  tabName = "population_overview_samples3",    
  fluidRow(
    box(
      id = "card2",
      title = "Haplotype Overview", 
      width = 12,
      height = "1200",
      background = "fuchsia",
      maximizable = TRUE, 
      collapsible = TRUE,
      
      plotOutput("population_overview_samples3", height = "100%")
    )
  )
)


# Population Overview
tab_population_overview_samples4 <- tabItem(
  tabName = "population_overview_samples4",    
  fluidRow(
    box(
      id = "card4",
      title = "Haplotype Overview", 
      width = 12,
      height = "1200",
      background = "fuchsia",
      maximizable = TRUE, 
      collapsible = TRUE,
      
      plotOutput("population_overview_samples4", height = "100%")
    )
  )
)


tab_row_search <- tabItem(
  tabName = "rowSearch",
  h4("Top 12"),
  fluidRow(width = 12
           
  )
)

#~~~~~~~~~~~



##################
# User Interface #
##################

ui <- dashboardPage(
  preloader = list(html = tagList(spin_1(), "Loading ..."), color = "#343a40"),
  dark = TRUE,
  scrollToTop = TRUE,
  header = dashboardHeader(
    border = F,
    tags$style(HTML("
        .dark-mode .content-wrapper {
          background-color: #44475a;
        }
        .carousel-indicators {
          bottom:-50px;
          backcolor: #44475a;
        }
        ")
    ),
    tags$head(
      tags$style(type="text/css", "#inline label{ display: table-cell; 
                                                  padding-right: 10px;
                                                  padding-left: 30px; 
                                                  text-align: center; vertical-align: middle; }
                                                  
                #inline .form-group { display: table-row;}
                        .hovertext text {
                            font-size: 40px !important;
                        }
                        .legendtext {
                            font-size: 20px !important;
                        }")
    ),
    
    title = dashboardBrand(
      title = "Haplotype Browser",
      color = "olive",
      image = "icons/icon-144.png",
      opacity = 1
    ),
    fixed = TRUE,
    rightUi = tagList(
      userOutput("user")
    ),
    tags$div(id = "inline", textInput("range", "Locus:", value = "Chr2:236680000..236690000")),
    tags$div(id = "inline", textInput("NAME", "Trait:", value = "se1")),
    tags$div(id = "inline", textInput("VCF", "VCF:", value = "28.5M")),
    submitButton("Submit", icon("search")),
    verbatimTextOutput("value")
    
  ),
  sidebar = dashboardSidebar(
    fixed = TRUE,
    skin = "light",
    width = "600px",
    status = "primary",
    id = "sidebar",
    sidebarMenu(
      id = "current_tab",
      flat = FALSE,
      compact = FALSE,
      childIndent = TRUE,
      sidebarHeader("Citra Spring 2021"),
      menuItem(
        "Haplotype Analyzer",
        tabName = "population_overview_samples4",
        icon = icon("seedling")
      ),
      menuItem(
        "Haplotype Long View",
        tabName = "population_overview_samples2",
        icon = icon("seedling")
      )
    ),
    
downloadButton("downloadData", "Download CSV"),
downloadButton("downloadData2", "Download XLSX"),


    radioButtons("ID_TYPE", "Which ID are you using?",
                 c("Genotype ID" = "geno_id",
                   "vcf ID" = "vcf_id"),
                 selected = "vcf_id"),
    textAreaInput("trait_textinput", "Which traits to map?", "C19_Germ\n", height = "200px"),
    textAreaInput("populationA", "Population A", "P39\n", height = "200px"),
    textAreaInput("populationB", "Population B", "IL14H\n", height = "200px"),
    
    submitButton("Search", icon("search")),
    verbatimTextOutput("rowids"),
    verbatimTextOutput("rowids.notfound")
  ),
  body = dashboardBody(
    skin = "light",
    use_theme(mytheme),
    tabItems(
      tab_population_overview_samples,
      tab_population_overview_samples2,
      tab_population_overview_samples3,
      tab_population_overview_samples4
    )
  ),
  controlbar = dashboardControlbar(
    id = "controlbar",
    skin = "light",
    pinned = FALSE,
    overlay = FALSE
  ),
  title = "Haplotype Browser"
)




