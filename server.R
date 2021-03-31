
# install required packages not yet installed
pckgs <- c('shiny','reshape2','ggplot2','plyr','arm')
todo <- pckgs[!is.element(pckgs, names(installed.packages()[,"Package"]))]
# if(length(todo) > 0) install.packages(todo)  #,repos="http://cran.us.r-project.org")

# load in the packages
library(shiny)
library(ggplot2)  # ggplot2 -> plotting
library(dplyr)     # ddply -> data manipulation
library(pwr)


shifting <- function(nrs,newMu){
	scale(nrs, center = TRUE, scale = FALSE) + newMu
}
inflating <- function(nrs,newSd){
	scale(nrs, center = TRUE, scale = TRUE) * newSd + mean(nrs)
}
combined <- function(nrs,newMu,newSd){
	shifting(inflating(nrs,newSd),newMu)
}

# server file
shinyServer(function(input,output,session){

	set.seed(123)
	obsA <- rnorm(250)
	obsB <- rnorm(250)

	setDta <- reactive({

		sMuA <- as.numeric(as.character((input$sMuA)))
		sSdA <- as.numeric(as.character((input$sSdA)))
		sNrA <- as.numeric(as.character((input$sNrA)))
		sDAB <- as.numeric(as.character((input$sDAB)))
		sMuB <- as.numeric(as.character((input$sMuB)))
		sSdB <- as.numeric(as.character((input$sSdB)))
		sNrB <- as.numeric(as.character((input$sNrB)))
		nrsA <- sample(obsA,sNrA)
		nrsB <- sample(obsB,sNrB)
		ctrl <- combined(nrsA,sMuA,sSdA)
		treat <- combined(nrsB,sMuB,sSdB)	

		list(sMuA=sMuA,sMuB=sMuB,sSdA=sSdA,sSdB=sSdB,sNrA=sNrA,sNrB=sNrB,sDAB=sDAB,nrsA=nrsA,nrsB=nrsB,ctrl=ctrl,treat=treat)
	})
	setPwr <- reactive({
	
		sT1a <- as.numeric(as.character((input$sT1a)))

		list(sT1a=sT1a)
	})
	getPwr <- reactive({
	
		inc <- setDta()
		pwr <- setPwr()
		# total number
		n <- (inc$sNrA + inc$sNrB) #*(inc$sSdA^4+inc$sSdA^4)/(inc$sSdB^2+inc$sSdB^2)^2
		# effect size is difference on average standard deviation
		d <- (inc$sMuB-inc$sMuA)/sqrt((inc$sNrA/n)*inc$sSdA^2+(inc$sNrB/n)*inc$sSdB^2) 
		# degrees of freedom is Welch corrected combination of variances
		df <- ((inc$sSdA^2/inc$sNrA) + (inc$sSdB^2/inc$sNrB))^2 / ((inc$sSdA^4/(inc$sNrA^2*(inc$sNrA-1))) + (inc$sSdB^4/(inc$sNrB^2*(inc$sNrB-1)))) 
		# power for effect size, sample size (df+2), alpha
		pwrp <- pwr.t.test(d=d,n=round((df+2)/2),sig.level=pwr$sT1a,type="two.sample",alternative="greater")$power
		# non-centrality parameter as effect size times sample size (corrected for inequality)
		ncp <- d * sqrt((inc$sNrA * inc$sNrB)/(inc$sNrA + inc$sNrB))#(inc$sMuB-inc$sMuA) / ((inc$sSdA+inc$sSdA)/2 * sqrt(2)) * sqrt(n)
		
		list(d=d,n=n,pwr=pwrp,ncp=ncp,df=df)
	})
	plotRawDta <- function(){

		inc <- setDta()
		dtb <- rbind(data.frame(id="ctrl",score=inc$ctrl),data.frame(id="treat",score=inc$treat))
		dts <- dtb %>% group_by(id) %>% summarize(sd=sd(score),mean=mean(score))
		dtx <- dts %>% mutate(lower=mean-sd,upper=mean+sd)
		dtx$score <- dtx$mean
		cmp2grps <- ggplot(dtb,aes(y=score,x=id)) + geom_point(data=dtx,aes(y=mean,x=id,col=id,size=2)) 
		cmp2grps <- cmp2grps + geom_linerange(data=dtx,aes(ymin = lower, ymax = upper,col=id,size=1.25)) 
		cmp2grps <- cmp2grps + geom_jitter() + guides(col=F,size=F) + geom_boxplot(alpha=.1)

		print(cmp2grps)
	}
	plotNorm <- function(){

		inc <- setDta()
		xvalues <- data.frame(x=c(-3*max(inc$sSdA,inc$sSdB)-min(inc$sMuA,inc$sMuB),3*max(inc$sSdA,inc$sSdB)+max(inc$sMuA,inc$sMuB)))
		pNorm <- ggplot(xvalues,aes(x=x))#values))
		pn1 <- function(x){
			dnorm(x,mean=inc$sMuA,sd=inc$sSdA)
		}
		pn2 <- function(x){
			dnorm(x,mean=inc$sMuB,sd=inc$sSdB)
		}
		pNorm <- pNorm + 
			scale_x_continuous(labels=floor(xvalues[[1]][1]):ceiling(xvalues[[1]][2]),breaks=floor(xvalues[[1]][1]):ceiling(xvalues[[1]][2])) +
			stat_function(fun=pn1,colour="red") + 
			stat_function(fun=pn2, colour="darkgreen") + 
			stat_function(fun=pn1,geom="area",fill="red",alpha=.6) + 
			stat_function(fun=pn2,geom="area",fill="darkgreen",alpha=.6)

		print(pNorm)
	}
	plotDiff <- function(){

		inc <- setDta()
		pwr <- setPwr()
		xvalues <- data.frame(x=c(-3*2*max(inc$sSdA/sqrt(inc$sNrA),inc$sSdB/sqrt(inc$sNrB))-min(inc$sMuA,inc$sMuB),3*2*max(inc$sSdA/sqrt(inc$sNrA),inc$sSdB/sqrt(inc$sNrB))+max(inc$sMuA,inc$sMuB)))
		se <- sqrt(inc$sSdA^2/inc$sNrA + inc$sSdB^2/inc$sNrB)
		pn1 <- function(x){
			dnorm(x,mean=0,sd=se)
		}
		pn2 <- function(x){
			dnorm(x,mean=inc$sMuB-inc$sMuA,sd=se)
		}
		pa1 <- function(x){
			out <- dnorm(x,mean=0,sd=se)
			out[x < qnorm(1-pwr$sT1a,mean=0,sd=se)] <- NA
			return(out)
		}
		pb2 <- function(x){
			out <- dnorm(x,mean=inc$sMuB - inc$sMuA,sd=se)
			out[x >= qnorm(1-pwr$sT1a,mean=inc$sMuB - inc$sMuA,sd=se)] <- NA
			return(out)
		}
		pDiff <- ggplot(xvalues,aes(x=x))#values))
		pDiff <- pDiff + 
			scale_x_continuous(labels=floor(xvalues[[1]][1]):ceiling(xvalues[[1]][2]),breaks=floor(xvalues[[1]][1]):ceiling(xvalues[[1]][2])) +
			stat_function(fun=pn1,colour="red") + 
			stat_function(fun=pn2, colour="darkgreen") + 
			stat_function(fun=pa1,xlim=c(qnorm(1-pwr$sT1a,mean=0,sd=se),3*se),geom="area",fill="red",alpha=.2) + 
			stat_function(fun=pb2,xlim=c(-3*se,qnorm(1-pwr$sT1a,mean=0,sd=se)),geom="area",fill="darkgreen",alpha=.2) +
			geom_vline(xintercept=qnorm(1-pwr$sT1a,mean=0,sd=se),lwd=2,col="steelblue") + annotate("text", -Inf, Inf, label = as.character(as.expression("SE==sigma~sqrt(2/n)")),parse = TRUE, hjust = -.25, vjust = 3) + xlab("~ confidence intervals") + ylab("density")

		print(pDiff)
	}
	
	plotNcp <- function(){

		inc <- setDta()
		pwr <- setPwr()
		ncp <- getPwr()
		xvalues <- data.frame(x=c(-3,3+ncp$ncp))
		se <- 1
		pn1 <- function(x){
			dnorm(x,mean=0,sd=se)
		}
		pn2 <- function(x){
			dnorm(x,mean=ncp$ncp,sd=se)
		}
		pa1 <- function(x){
			out <- dnorm(x,mean=0,sd=se)
			out[x < qnorm(1-pwr$sT1a,mean=0,sd=se)] <- NA
			return(out)
		}
		pb2 <- function(x){
			out <- dnorm(x,mean=ncp$ncp,sd=se)
			out[x >= qnorm(1-pwr$sT1a,mean=ncp$ncp,sd=se)] <- NA
			return(out)
		}
		pNcp <- ggplot(xvalues,aes(x=x)) #values))
		pNcp <- pNcp + 
			scale_x_continuous(labels=floor(xvalues[[1]][1]):ceiling(xvalues[[1]][2]),breaks=floor(xvalues[[1]][1]):ceiling(xvalues[[1]][2])) +
			stat_function(fun=pn1,colour="red") + 
			stat_function(fun=pn2, colour="darkgreen") + 
			stat_function(fun=pa1,xlim=c(qnorm(1-pwr$sT1a,mean=0,sd=se),3*se),geom="area",fill="red",alpha=.2) + 
			stat_function(fun=pb2,xlim=c(-3*se,qnorm(1-pwr$sT1a,mean=0,sd=se)),geom="area",fill="darkgreen",alpha=.2) +
			geom_vline(xintercept=qnorm(1-pwr$sT1a,mean=0,sd=se),lwd=2,col="steelblue") + annotate("text", -Inf, Inf, label = "N(0,1)",parse = TRUE, hjust = -.25, vjust = 3) + xlab("Z-values") + ylab("density")
		print(pNcp)

	}
	output$sMuA <- renderUI({
		sliderInput("sMuA","average C",min=-10,max=10,value=0)
	})
	output$sDAB <- renderUI({
		sliderInput("sDAB","difference groups",min=0,max=10,value=2)
	})
	output$sMuB <- renderUI({
		sliderInput("sMuB","average T",min=-10,max=10,value=2)
	})
	observeEvent(input$sMuA,  {
		inc <- setDta()
		updateSliderInput(session = session, inputId = "sMuB", min = inc$sMuA, max = 10+inc$sMuA*10, value = inc$sMuA+inc$sDAB)
	})
	observeEvent(input$sDAB,  {
		inc <- setDta()
		updateSliderInput(session = session, inputId = "sMuB", min = inc$sMuA, max = 10+inc$sMuA*10, value = inc$sMuA+inc$sDAB)
	})
	observeEvent(input$sMuB,  {
		inc <- setDta()
		updateSliderInput(session = session, inputId = "sDAB", min = 0, max = 10+inc$sMuB*2, value = inc$sMuB-inc$sMuA)
	})
	output$sSdA <- renderUI({
		sliderInput("sSdA","standard deviation C",min=.01,max=16,value=4)
	})
	output$sSdB <- renderUI({
		sliderInput("sSdB","standard deviation T",min=.01,max=16,value=4)
	})
	output$sNrA <- renderUI({
		sliderInput("sNrA","number observations C",min=10,max=250,value=32)
	})
	output$sNrB <- renderUI({
		sliderInput("sNrB","number observations T",min=10,max=250,value=32)
	})
	output$sT1a <- renderUI({
		sliderInput("sT1a","Type I error (alpha)",min=0.01,max=.2,value=.05)
	})
	output$pRawDta <- renderPlot({
		req(input$sMuA,input$sMuB,input$sSdA,input$sSdB,input$sNrA,input$sNrB,input$sDAB)
		# plotOutput(plotRawDta())
		plotRawDta()
	})
	output$pNorm <- renderPlot({
		req(input$sMuA,input$sMuB,input$sSdA,input$sSdB,input$sNrA,input$sNrB,input$sDAB)
		# plotOutput(plotNorm())
		plotNorm()
	})
	output$pDiff <- renderPlot({
		req(input$sMuA,input$sMuB,input$sSdA,input$sSdB,input$sNrA,input$sNrB,input$sDAB)
		# plotOutput(plotDiff())
		plotDiff()
	})
	output$pNcp <- renderPlot({
		req(input$sMuA,input$sMuB,input$sSdA,input$sSdB,input$sNrA,input$sNrB,input$sDAB,input$sT1a)
		# plotOutput(plotNcp())
		plotNcp()
	})
	
	output$comments <- renderText({
		# if(any(is.null(input$sMuA),is.null(input$sMuB),is.null(input$sSdA),is.null(input$sSdB),is.null(input$sNrA),is.null(input$sNrB),is.null(input$sT1a),is.null(input$sDAB))){return("")}
		req(input$sMuA,input$sMuB,input$sSdA,input$sSdB,input$sNrA,input$sNrB,input$sDAB,input$sT1a)
		inc <- setDta()
		pwr <- setPwr()
		out <- getPwr()
		# list(d=d,n=n,pwr=pwrp,ncp=ncp)
		paste("<br>POWER = ",round(out$pwr,3),"<br><br>df = ",out$df,"<br><br>NCP = ",out$ncp,"<br>Cohen d = ",out$d,"<br>total n = ",out$n)
	})
})