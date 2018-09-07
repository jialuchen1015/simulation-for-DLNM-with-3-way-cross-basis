require(shiny)
require(plotrix)
plot.mypred <- function(ar,col="lightskyblue",
                        truth = NULL){

    if(!is.null(truth)){
        if(sum(dim(ar) != dim(truth)) > 0){
            stop("In `plot.mypred`: ar and truth should be of same D.")
        }
    }
    a <- ar$ar
    xname <- ar$xname
    uname <- ar$uname

    x <- ar$x
    u <- ar$u
    xdens <- length(x)
    udens <- length(u)

    lag <- ar$lag

    # CI
    lci <- ar$l
    hci <- ar$h

    # Plot range
    l <- min(a)
    h <- max(a)
    if(!is.null(truth)){
        l <- min(min(truth$ar),l)
        h <- max(max(truth$ar),h)
    }

    shinyApp(
        ui = fluidPage(
            titlePanel("3-way cross basis"),
            sidebarLayout(
              sidebarPanel(
                sliderInput('x', paste('i-th ',xname), min = 1,
                            max = xdens, step=1, value = 1),
                sliderInput('u', paste('i-th ',uname), min = 1,
                            max = udens, step=1, value = 1),
                sliderInput('lag', paste('i-th lag'), min = 0,
                            max = lag, step=1, value = 1)
              ),
              mainPanel("main panel",
                      fluidRow(
                             plotOutput('plot1',width = "850px", height = "650px")
                             )
                        ),
              position = "left"
            )


        ),

        server = function(input, output, session) {
            plot.arg <- list(ticktype = "detailed",
                           theta = 30,phi = 30,
                           zlab = "Outcome",
                           col = col,
                           zlim = c(l,h),
                           ltheta = 290, shade = 0.75, r = 1, d = 5)

            output$plot1 <- renderPlot({

              par(mfrow=c(2,3))
                i <- input$x
                # plot lag vs o3
                wantvar <- a[,,i]

                plot.arg.1 <- modifyList(plot.arg, list(
                                                      xlab = "Lag",
                                                      ylab = uname,
                                                      z = wantvar,
                                                      x = (0:lag),
                                                      y = u,
                                                      main=paste(xname," = ",round(x[i],2))))
                do.call("persp", plot.arg.1)

                j <- input$u

                wantvar <- a[,j,]
                plot.arg.2 <- modifyList(plot.arg, list(
                                                      xlab = "Lag",
                                                      ylab = xname,
                                                      z = wantvar,
                                                      x = (0:lag),
                                                      y = x,
                                                      main=paste(uname," = ",round(u[j],2))))
                do.call("persp", plot.arg.2)

                k <- input$lag + 1


                wantvar <- a[k,,]
                plot.arg.3 <- modifyList(plot.arg, list(
                                                      xlab = uname,
                                                      ylab = xname,
                                                      z = wantvar,
                                                      x = u,
                                                      y = x,
                                                      main=paste("Lag"," = ",k-1)))
                do.call("persp", plot.arg.3)


              k <- input$lag + 1
              j <- input$u

              wantvar <- a[k,j,]
              li <- lci[k,j,]
              hi <- hci[k,j,]
              plotCI(x=x,y=wantvar,ui=hi,li=li,xlab=xname,ylab="Outcome",
                     main=paste("Lag = ",k," ",
                                uname," = ",round(u[j],2)),
                     ylim = c(l,h))
              if(!is.null(truth)){
                  lines(x,truth$ar[k,j,])
              }

              k <- input$lag + 1
              i <- input$x

              # plot lag vs o3
              # plot lag vs o3
              wantvar <- a[k,,i]
              li <- lci[k,,i]
              hi <- hci[k,,i]
              plotCI(x=u,y=wantvar,ui=hi,li=li,xlab=uname,ylab="Outcome",
                     main=paste("Lag = ",k-1," ",
                                xname," = ",round(x[i],2)),
                     ylim = c(l,h))
              if(!is.null(truth)){
                  lines(u,truth$ar[k,,i])
              }

              j <- input$u
              i <- input$x

              wantvar <- a[,j,i]
              li <- lci[,j,i]
              hi <- hci[,j,i]
              plotCI(x=(0:lag),y=wantvar,ui=hi,li=li,xlab="Lag",ylab="Outcome",
                     main=paste(uname," = ",round(u[j],2)," ",
                                xname," = ",round(x[i],2)),
                     ylim = c(l,h))
              if(!is.null(truth)){
                  lines((0:lag),truth$ar[,j,i])
              }
            })
        }
    )
}
