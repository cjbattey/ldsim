#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
source("global.R")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Linkage and Linkage Disequilibrium"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            actionButton("newpop","new population"),
            actionButton("meiosis","next generation"),
            br(),br(),
            sliderInput("r","recombination rate",0,0.5,value = 0.2),
            br(),br()
        ),

        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(type = "tabs",
                        tabPanel("chromosomes",
                                 helpText("Here's a diploid population where we're tracking two loci 
                                          (regions of the genome). At each locus there are two alleles -- 
                                          A/a and B/b. In the first generation all of the chromosomes are
                                          either AB or ab. Recombination generates new combinations of alleles
                                          by swapping pieces of chromosomes during meiosis -- shown here with
                                          crossing purple lines."),
                                  helpText("Try clicking 'next generation' a few times and 
                                          see what happens to the population. Does the number of alleles generally
                                          go up or down? What happens when we change the recombination rate (r)?"),
                                 plotOutput("pop_plot")),
                        tabPanel("frequency tables",
                                 helpText("Linkage disequilibrium measures how often we find two alleles on 
                                           the same chromosome, relative to random chance (that's the 'equilibrium'
                                           in the name). If two alleles are completely
                                           independent, then the probability of finding them together on
                                           one chromosome is just the product of finding each allele by itself: 
                                           p(a)p(b). To measure how far off from this expectation we are,
                                           we take our observed haplotype frequency and subtract this expectation: 
                                           D = p(ab) - p(a)p(b). Here you can find the haplotype and allele frequencies
                                           for your simulated population, as well as D calculated for the a and b alleles."),
                                 helpText("Try calculating D(a,b) yourself -- can you get the same answer? 
                                           What is D for the A and b alleles?"),
                                 tableOutput("hftable"),
                                 tableOutput("aftable"),
                                 textOutput("ld")),
                        tabPanel("Change Over Time",
                                 helpText("Linkage disequilbrium (D) decays over time as recombination 
                                           breaks up association among loci. In a random-mating 
                                           population with no selection, D should decrease by
                                           -rD each generation. Try a few values of r and n in the simulation 
                                           below. What happens when you increase r? What about n?"),
                                 sliderInput("n","population size",1,100,value = 50),
                                 sliderInput("gen",label = "generations",min=1,max=100,value=50),
                                 actionButton("timerun","run"),
                                 br(),br(),
                                 plotOutput("timeplot"))
                        )
        )
        
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    pop <- reactiveVal()
    
    observeEvent(input$newpop,{
        newpop <- data.frame(l1=c(rep("A",10),rep("a",10)),
                             l2=c(rep("B",10),rep("b",10)),
                             recomb=0)
        pop(newpop) #equivalent to pop <- newpop, but for reactive values. Note this must be in observeEvent(), not eventReactive()
    },ignoreNULL=F)
    
    observeEvent(input$meiosis,{
        newpop <- makegametes(pop(),input$r)
        pop(newpop)
    })
    
    p <- eventReactive(c(pop(),input$r),{
        plot_population(pop(),input$r)
    })
    output$pop_plot <- renderPlot(p())
    
    hf <- reactive({
        hf <- haplotype_freqs(pop())
    })
    output$hftable <- renderTable(hf())
    
    af <- reactive({
        allele_freqs(pop())[,1:2]
    })
    output$aftable <- renderTable(af())
    
    ld <- eventReactive(pop(),{
        a <- getLD("a","b",pop())
        paste("Linkage Disequilibrium (D) =",a)
    })
    output$ld <- renderText(ld())
    
    p2 <- eventReactive(input$timerun,{
        pop <- data.frame(l1=c(rep("A",input$n),rep("a",input$n)),
                          l2=c(rep("B",input$n),rep("b",input$n)),
                          recomb=0)
        d <- c(getLD("a","b",pop))
        for(i in 2:input$gen){
            pop <- makegametes(pop,input$r)
            d[i] <- getLD("a","b",pop)
        }
        d <- data.frame(gen=1:length(d),D=d)
        dpred <- LDpred(d$D[1],input$r,input$gen)
        d <- cbind(d,dpred)
        names(d) <- c("gen","simulation","theory")
        md <- melt(d,id.vars = "gen")
        md
    })
    
    output$test <- renderText(p2())
    
    p3 <- eventReactive(p2(),{
        a <- p2()
        ggplot(data=a,aes(x=gen,y=value))+
            theme_classic()+theme(legend.position=c(0.8,0.8),legend.text=element_text(size=12))+
            ylab("Linkage Disequilibrium (D)")+xlab("Generation")+
            scale_color_manual(values = c("black","red"),name="")+
            geom_line(aes(col=variable))+
            geom_point(data=subset(a,variable=="simulation"),color="black",size=1)
    })
    output$timeplot <- renderPlot(p3())
}

# Run the application 
shinyApp(ui = ui, server = server)
