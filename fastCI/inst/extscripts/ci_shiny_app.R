library(shiny)
library(gtools)


computeInversionsNumber <- function(x,N, unallowed_index){
  # order.y <- order(y)
  # y <- y[order.y]
  # x <- x[order.y]
  numInv <- 0
  unallowed_pairs <- as.data.frame(combn(N, 2)[,as.numeric(unlist(unallowed_index)), drop=FALSE])
  for(i in 1:length(x)){
    for(j in i:length(x)){
     if(x[i] > x[j] & !(list(c(x[j],x[i])) %in% unallowed_pairs)){
        numInv = numInv + 1
      }
    }
  }
  # browser()
  return(numInv)
}

computeDistN <- function(N, unallowed_index){
  
  myx <- permutations(N, N)
  x <- seq(1,N)
  # print(myx)
  inversions <- apply(myx, 1, function(myx){
    result <- computeInversionsNumber(x[myx], N, unallowed_index)
    return(result)
  })
  return(table(inversions))
  
}

ui <- fluidPage(
  titlePanel("Exploring the CI with arbitrary ties"),
  uiOutput("mainPage")
)

server <- function(input, output) {
  
 
  # N <- 3
  
  # toPlot <- 
  
  output$distribution <- renderPlot({
    
    count <- computeDistN(N(), input$invalidPairs)
    
    
    plot(count)})
  
  
  output$mainPage <- renderUI(fluidRow(
    fluidRow(column(10, sliderInput("numSamples", label = h3("Number of Samples"), min = 0,
                          max = 8, value = N()))),
    fluidRow(
      column(3,
             checkboxGroupInput("invalidPairs", label="Choose which Pairs Are Equivilent", choices = list_of_pairs())
             ),
      column(7, plotOutput("distribution"))
      )
    )
  )
  
  N <- reactive({
    N <- ifelse(is.null(input$numSamples), 4, input$numSamples)
    
  })
  
  list_of_pairs <- reactive({
    N <- N()
    list_of_pairs <- as.list(seq(1, choose(N,2)))
    names(list_of_pairs) <- apply(combn(N,2), 2, function(x) paste(x, collapse = " - "))
    list_of_pairs
  })
  
}

shinyApp(ui = ui, server = server)

