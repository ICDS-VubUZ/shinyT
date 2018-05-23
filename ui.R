
shinyUI(fluidPage(

	titlePanel("Exploring relation data and power for t-test"),
		fluidRow(
			column(3,
				wellPanel(
					h3("control condition"),
					uiOutput("sMuA"),
					uiOutput("sSdA"),
					uiOutput("sNrA")
				)
			),
			column(3,
				wellPanel(
					h3("treatment condition"),
					uiOutput("sMuB"),
					uiOutput("sSdB"),
					uiOutput("sNrB")
				)
			),
			column(3,
				uiOutput("sDAB"),
				plotOutput(outputId="pRawDta",height="300")
			),
			column(3,
				plotOutput("pNorm")
			)
		),
		fluidRow(
			column(2,
				wellPanel(
					h3("one-sided t-test"),
					h4("independent means"),
					uiOutput("sT1a"),
					htmlOutput("comments")
				)
			),
			column(4,
				plotOutput("pDiff")
			),
			column(4,
				plotOutput("pNcp")
			),
			column(2,
				h5("")
			)
		)
	)
)