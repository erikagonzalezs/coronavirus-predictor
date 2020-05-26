#================================
# SIMULADOR COVID-19
# PROYECTO 
#================================
# Integrantes
# - Andres Camilo Giraldo Gil
# - Erika Alejandra Gonzalez
# - Leonel Steven Londono
#================================
# Analisis Numerico
#================================

#================================
# Instalar estas librerias
#================================
#install.packages("shinydashboard")
#install.packages("shinyjs")
#install.packages("deSolve")
#install.packages("phaseR")

#================================
# Librerias Usadas
#================================

library (shiny) 
library (shinydashboard) 
library (deSolve)
library (phaseR)

#================================
# Variables Globales
#================================

#================================
# Se carga base de datos
#================================

baseDeDatos <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv",
                        na.strings = "", fileEncoding = "UTF-8-BOM")


paisASimular <- "Colombia"
individuosMaximos <- 0
poblacionInfectadaActual <- 0
poblacionRecuperadaFallecidaActual <- 0
poblacionSusceptibleActual <- 0
casosPorDia <- c()

for(i in length(baseDeDatos$countriesAndTerritories):1)
{
  if(baseDeDatos$countriesAndTerritories[i]== paisASimular)
  {
    casosPorDia <- c(casosPorDia, baseDeDatos$cases[i])
    individuosMaximos <- baseDeDatos$popData2018[i]
    poblacionInfectadaActual <- poblacionInfectadaActual + baseDeDatos$cases[i]
    poblacionRecuperadaFallecidaActual <- poblacionRecuperadaFallecidaActual + baseDeDatos$deaths[i]
  }
  
}

poblacionSusceptibleActual <- individuosMaximos - poblacionInfectadaActual - poblacionRecuperadaFallecidaActual

cat("Maximos ", individuosMaximos, 
    "Infectados: ", poblacionInfectadaActual,
    "Fallecidos: ", poblacionRecuperadaFallecidaActual, 
    "susceptibles: ",poblacionSusceptibleActual)


#================================
# Funciones Globales
#================================
obtenerResultadosSIR <-function(input, metodo, nombreGrafica)
{
  
  #================================
  # Se obtienen los datos iniciales de la interfaz
  #================================
  entradaSusceptibles <- (input$suceptiblesInicialesSIR*1) / individuosMaximos
  entradaInfectados <- (input$infectadosInicialesSIR*1) / individuosMaximos
  entradaRecuperados <- (input$recuperadosFallecidosInicialesSIR*1) / individuosMaximos
  
  
  datosIniciales <- c( S = (entradaSusceptibles), I = (entradaInfectados), R = (entradaRecuperados))
  
  parameters <- c(beta = input$bethaSIR, gamma = input$gammaSIR)
  
  tiempoSimulacion <- seq(0, input$tiempoSimulacionSIR, by = 1)
  
  #================================
  # Se resuelve la ecuacion diferencial ordinaria
  #================================
  solucionEDO <- ode( y = datosIniciales,
                      times = tiempoSimulacion, 
                      func = funcionSIR, 
                      parms = parameters,
                      method = metodo)
  
  solucionEDO <- as.data.frame(solucionEDO)
  
  #================================
  # Se grafica la simulacion
  #================================
  plot(solucionEDO$time,solucionEDO$S*individuosMaximos,col="blue", type="l", xlab="Dias de simulacion", ylab="Numero de Individuos",ylim =c(0,input$poblacionInicialSIR), main = nombreGrafica)
  par(new=TRUE)
  plot(solucionEDO$time,solucionEDO$I*individuosMaximos,col="red", type="l", xlab="Dias de simulacion", ylab="Numero de Individuos", ylim = c(0,input$poblacionInicialSIR), main = nombreGrafica)
  par(new=T)
  plot(solucionEDO$time,solucionEDO$R*individuosMaximos, col="green", type="l", xlab="Dias de simulacion", ylab="Numero de Individuos" , ylim = c(0,input$poblacionInicialSIR), main = nombreGrafica)
  
  legend(x = "topright", legend=c("Susceptibles", "Infectados", "Muertos"), col=c("blue", "red","green"), lty=rep(1, 1))
  
  solucionEDO$time <- NULL 
  
  return (solucionEDO)
}


obtenerResultadosSIS <-function(input, metodo, nombreGrafica)
{
  
  #================================
  # Se obtienen los datos iniciales de la interfaz
  #================================
  entradaSusceptibles <- (input$suceptiblesInicialesSIS*1) / individuosMaximos
  entradaInfectados <- (input$infectadosInicialesSIS*1) / individuosMaximos
  
  
  datosIniciales <- c( S = (entradaSusceptibles), I = (entradaInfectados))
  
  parameters <- c(beta = input$bethaSIS, gamma = input$gammaSIS)
  
  tiempoSimulacion <- seq(0, input$tiempoSimulacionSIS, by = 1)
  
  #================================
  # Se resuelve la ecuacion diferencial ordinaria
  #================================
  solucionEDO <- ode( y = datosIniciales,
                      times = tiempoSimulacion, 
                      func = funcionSIS, 
                      parms = parameters,
                      method = metodo)
  
  solucionEDO <- as.data.frame(solucionEDO)
  
  #================================
  # Se grafica la simulacion
  #================================
  plot(solucionEDO$time,solucionEDO$S*individuosMaximos,col="blue", type="l", xlab="Dias de simulacion", ylab="Numero de Individuos",ylim =c(0,input$poblacionInicialSIS), main = nombreGrafica)
  par(new=TRUE)
  plot(solucionEDO$time,solucionEDO$I*individuosMaximos,col="red", type="l", xlab="Dias de simulacion", ylab="Numero de Individuos", ylim = c(0,input$poblacionInicialSIS), main = nombreGrafica)
  par(new=T)
  legend(x = "topright", legend=c("Susceptibles", "Infectados"), col=c("blue", "red"), lty=rep(1, 1))
  
  solucionEDO$time <- NULL 
  
  return (solucionEDO)
}


#================================
# Funcion SIR
#================================
funcionSIR <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- -beta * S * I
    dI <-  beta * S * I - gamma * I
    dR <-                 gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

#================================
# Funcion SIS
#================================
funcionSIS <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- gamma*S*I + beta*I
    dI <-  gamma*S*I - beta*I
    
    return(list(c(dS, dI)))
  })
}

#================================
# Contenido de informacion
#================================

contenidoInformacion <- fluidRow(
  box(
    title = "SIMULADOR DE CORONAVIRUS EN COLOMBIA",
    width = 8,
    h4("Descripcion: "),
    h5("Este es un aplicativo que realiza la simulacion del virus Covid-19 en colombia en un maximo de 100 dias, en este aplicativo,
       se utilizo una base de datos online en la cual se encuentran los datos actualizados del covid-19 en colombia dia a dia,se puede realizar la simulacion por dos modelos epidemologicos: "),
    h4("SIS: "),
    h5("Para este estudio se eligió un modelo discreto metapoblacional de transmisión de enfermedades
SIS. En este tipo de modelos, la población se divide en dos grupos de personas: las que
han sido infectadas por la enfermedad y son infecciosas y las que son susceptibles de ser
infectadas por la enfermedad.
Los modelos SIS se usan para enfermedades en las que no hay inmunidad, pues, una vez
que las personas infectadas se recuperan, pasan a ser de nuevo susceptibles. Por lo tanto, la
progresión de la enfermedad desde el punto de vista de un individuo es susceptible-infectado-susceptible"),
    h5("Ecuaciones utilizadas para el modelo: "),
    h5("dS = ???SI + ??I"),
    h5("dI = ???SI - ??I"),
    h4("SIR: "),
    h5("Este modelo se basa en ecuaciones diferenciales para describir la dinámica de los contagios en una población cerrada con N individuos que inicialmente son susceptibles (S) al patógeno y que, 
       a partir de un infectado inicial, van contagiándose a una determinada velocidad y pasando a ser infectados (I). Tras un período de enfermedad activa,
       los que no fallecen pasan al estado de inmunes: se han recuperado (R) y ya no contagiarán más. Por tanto, la población susceptible va disminuyendo hasta que ya no se produzcan más contagios."),
    h5("Ecuaciones utilizadas para el modelo: "),
    h5("1. dS = -??SI"),
    h5("2. dI = ??SI - ???I"),
    h5("3. dR = ???I")
  ),
  box(
    title = "Terminologia",
    width = 4,
    h5("S - Individuos susceptibles."),
    h5("I - Individuos infectados."),
    h5("R - Individuos recobrados o recuperados."),
    h5("?? - Tasa de contagios (probabilidad de que una persona enferme al estar en contacto con un infectado)."),
    h5("??? - Tasa de mortalidad (probabilidad de una persona infectada muera).")
  ),
  box(
    title = "Referencias",
    width = 4,
    h5("[1] Base de datos usada: 
       https://opendata.ecdc.europa.eu/
       covid19/casedistribution/csv"),
    h5("[2] Tesis de grado: Modelos epidemiológicos basados en ecuaciones diferenciales, 2016, universidad de la Rioja, Obtenido de: https://biblioteca.unirioja.es/tfe_e/TFE002211.pdf"),
    h5("[3] Reportes sobre la cantidad de casos a nivel mundial: https://www.who.int/emergencies/diseases/novel-coronavirus-2019/situation-reports "),
    h5("[4] Modelación SIR de la pandemia de COVID-19 en colombia 25 de marzo de 2020:
http://www.scielo.org.co/pdf/rsap/v22n1/0124-0064-rsap-22-01-e185977.pdf")
  )
)

#================================
# Contenido SIS
#================================
contenidoSIS <- fluidRow(
  
  box(
    title = "Datos iniciales S.I.S", 
    width = 4,
    shinyjs::useShinyjs(),
    shinyjs::disabled(
      numericInput("poblacionInicialSIS", "Poblacion Inicial:", individuosMaximos, min = 1, max = 100000)),
    sliderInput("infectadosInicialesSIS", "Infectados iniciales:", poblacionInfectadaActual,min = 1, max = individuosMaximos),
    sliderInput("suceptiblesInicialesSIS", "Susceptibles iniciales:", poblacionSusceptibleActual, min = 300, max = individuosMaximos),
    sliderInput("tiempoSimulacionSIS", "Tiempo limite:", 10, min = 10, max = 100),
    sliderInput("bethaSIS", "Tasa de recuperacion:", 0.2, min = 0.1, max = 1),
    sliderInput("gammaSIS", "Tasa de transmision:", 0.2, min = 0.1, max = 1),
    actionButton("botonCalcularSIS", "Calcular")
  ),
  tabBox(
    title = "Comportamiento del virus a traves del tiempo con Runge Kutta",
    width =8,
    id = "contenedorGraficaModeloSISRungeKutta", height = "300px",
    tabPanel(
      title = "Grafica modelo",
      width=200,
      plotOutput("graficaModeloSISRungeKutta", height = 200)
    ),
    tabPanel(
      title = "Grafica error",
      width=200,
      plotOutput("graficaErrorSISRungeKutta", height = 200)
    )
  ),
  tabBox(
    title = "Comportamiento del virus a traves del tiempo con Euler",
    width =8,
    id = "contenedorGraficaModeloSIS", height = "300px",
    tabPanel(
      title = "Grafica modelo",
      width=200,
      plotOutput("graficaModeloSISEuler", height = 200)
    ),
    tabPanel(
      title = "Grafica error",
      width=200,
      plotOutput("graficaErrorSISEuler", height = 200)
    )
  ),
  box(
    title = "Campo de pendientes del modelo SIS",
    width = 8,
    id = "contenedorPendientesModeloSIS", height = "300px",
    plotOutput("graficaPendientesSIS", height = 200)
  )
)

#================================
# Contenido SIR
#================================
contenidoSIR <- fluidRow(
  
  box(
    title = "Datos iniciales S.I.R", 
    width = 4,
    shinyjs::useShinyjs(),
    shinyjs::disabled(
      numericInput("poblacionInicialSIR", "Poblacion Inicial:", individuosMaximos, min = 1, max = individuosMaximos)),
    sliderInput("infectadosInicialesSIR", "Infectados iniciales:", poblacionInfectadaActual,min = 1, max = individuosMaximos),
    sliderInput("suceptiblesInicialesSIR", "Susceptibles iniciales:", poblacionSusceptibleActual, min = 300, max = individuosMaximos),
    sliderInput("recuperadosFallecidosInicialesSIR", "Recuperados/Fallecidos iniciales:", poblacionRecuperadaFallecidaActual, min = 0, max =  individuosMaximos),
    sliderInput("tiempoSimulacionSIR", "Tiempo limite:", 10, min = 10, max = 100),
    sliderInput("bethaSIR", "Tasa de transmision:", 0.2, min = 0.1, max = 1.0),
    sliderInput("gammaSIR", "Tasa de recuperacion:", 0.2, min = 0.1, max = 1.0),
    actionButton("botonCalcularSIR", "Calcular")
  ),
  tabBox(
    title = "Comportamiento del virus a traves del tiempo con Runge Kutta",
    width =8,
    id = "contenedorGraficaModeloSIRRungeKutta", height = "300px",
    tabPanel(
      title = "Grafica modelo",
      width=200,
      plotOutput("graficaModeloSIRRungeKutta", height = 200)
    ),
    tabPanel(
      title = "Grafica error",
      width=200,
      plotOutput("graficaErrorSIRRungeKutta", height = 200)
    )
  ),
  tabBox(
    title = "Comportamiento del virus a traves del tiempo con Euler",
    width =8,
    id = "contenedorGraficaModeloSIR", height = "300px",
    tabPanel(
      title = "Grafica modelo",
      width=200,
      plotOutput("graficaModeloSIREuler", height = 200)
    ),
    tabPanel(
      title = "Grafica error",
      width=200,
      plotOutput("graficaErrorSIREuler", height = 200)
    )
  ),
  box(
    title = "Campo de pendientes del modelo SIR",
    width = 8,
    id = "contenedorPendientesModeloSIR", height = "300px",
    plotOutput("graficaPendientesSIR", height = 200)
  )
  
)

#================================
# Elementos que tendra la aplicación
#================================
InterfazGrafica <- dashboardPage ( 
  dashboardSidebar(),
  dashboardHeader (),
  dashboardBody()
) 

#================================
# Barra de navegacion
#================================
barraNavegacion <- dashboardHeader(
  
  title = "Simulador del COVID-19 En Colombia",
  titleWidth = "500px"
  
)

#================================
# Menu lateral de la aplicacion
#================================
menuLateral <- dashboardSidebar(
  
  sidebarMenu(
    menuItem("Informacion", tabName = "Informacion", icon = icon("bell")),
    menuItem("S.I.S", tabName = "SIS", icon = icon("circle")    ),
    menuItem("S.I.R.", tabName = "SIR", icon = icon("circle") ),
    menuItem("Desarrolladores", icon = icon("gear"),
             menuSubItem("Erika Alejandra Gonzalez", tabName = "subitem1", icon = icon("user")),
             menuSubItem("Leonel Steven Londoño", tabName = "subitem2", icon = icon("user")),
             menuSubItem("Andres Camilo Giraldo", tabName = "subitem3", icon = icon("user"))
    )
  )
)

#================================
# Contenido de la apliacion
#================================
contenidoAplicacion <- dashboardBody(
  tabItems(
    tabItem("Informacion", contenidoInformacion),
    tabItem("SIS",contenidoSIS),
    tabItem("SIR", contenidoSIR)
  )
)

#================================
# Servidor de la aplicacion
#================================
Servidor <- function(input, output, session){
  
  #================================
  # Modelo del SIS
  #================================
  
  calcularSIS <- function(){
    output$graficaModeloSISRungeKutta <- renderPlot({
      
      solucionEDO = obtenerResultadosSIS(input,"rk4", "Metodo 1: RungeKutta")
      
    })
    
    output$graficaModeloSISEuler <- renderPlot({
      
      solucionEDO = obtenerResultadosSIS(input,"euler", "Metodo 1: Euler")
      
    })
    
    output$graficaPendientesSIS <- renderPlot({
      
      scopeField <- function(t, p, parameters){
        k <- parameters[1]
        n <- parameters[2]
        dp <- k*(p*(n-p))
        list(dp)
      }
      
      limiteY <- input$suceptiblesInicialesSIS / individuosMaximos
      
      scopeField.flowField <- flowField(scopeField, 
                                        xlim = c(0,
                                                 input$tiempoSimulacionSIS),
                                        ylim = c(0,limiteY), 
                                        parameters = c(0.0005,1000),
                                        system = "one.dim",
                                        add = FALSE, 
                                        xlab = "Dias de simulacion", 
                                        ylab = "Numero de individuos", 
                                        main = "Campo de pendientes")
      
    })
    
    
    output$graficaErrorSISRungeKutta <- renderPlot({
      
      #Experimental
      solucionEDORunge = obtenerResultadosSIS(input,"rk4", "Metodo 1: RungeKutta")
      
      #Teorico
      solucionEDOEuler = obtenerResultadosSIS(input,"euler", "Metodo 1: Euler")
      
      
      errores <- c()
      dias <- c()
      
      for(i in 1:length(solucionEDORunge$I-1)){
        errores[i] <- abs((solucionEDOEuler$I[i] - solucionEDORunge$I[i]) / (solucionEDOEuler$I[i]))*100
        dias <- c(dias, i)
      }
      
      plot(dias , errores, col="red", type="l", xlab="Dias Simulacion", ylab="Error", main = "Error de la simulacion", ylim = c(0,100))
    })
    
    output$graficaErrorSISEuler <- renderPlot({
      
      #Experimental
      solucionEDORunge = obtenerResultadosSIS(input,"rk4", "Metodo 1: RungeKutta")
      
      #Teorico
      solucionEDOEuler = obtenerResultadosSIS(input,"euler", "Metodo 1: Euler")
      
      
      errores <- c()
      dias <- c()
      
      for(i in 1:length(solucionEDORunge$I-1)){
        errores[i] <- abs((solucionEDORunge$I[i] - solucionEDOEuler$I[i]) / (solucionEDORunge$I[i]))*100
        dias <- c(dias, i)
      }
      
      plot(dias , errores, col="red", type="l", xlab="Dias Simulacion", ylab="Error", main = "Error de la simulacion", ylim = c(0,100))
    })
    
  }
  
  
  #================================
  # Modelo del SIR
  #================================
  calcularSIR <- function(){
    output$graficaModeloSIRRungeKutta <- renderPlot({
      
      solucionEDO = obtenerResultadosSIR(input,"rk4", "Metodo 1: RungeKutta")
      
    })
    
    output$graficaModeloSIREuler <- renderPlot({
      
      solucionEDO = obtenerResultadosSIR(input,"euler", "Metodo 1: Euler")
      
    })
    
    output$graficaPendientesSIR <- renderPlot({
      
      scopeField <- function(t, p, parameters){
        k <- parameters[1]
        n <- parameters[2]
        dp <- k*(p*(n-p))
        list(dp)
      }
      
      limiteY <- input$suceptiblesInicialesSIR / individuosMaximos
      
      scopeField.flowField <- flowField(scopeField, 
                                        xlim = c(0,input$tiempoSimulacionSIR),
                                        ylim = c(0,limiteY), 
                                        parameters = c(0.0005,1000),
                                        system = "one.dim",
                                        add = FALSE, 
                                        xlab = "Dias de simulacion", 
                                        ylab = "Numero de individuos", 
                                        main = "Campo de pendientes")
    })
    
    output$graficaErrorSIRRungeKutta <- renderPlot({
      
      #Experimental
      solucionEDORunge = obtenerResultadosSIR(input,"rk4", "Metodo 1: RungeKutta")
      
      #Teorico
      solucionEDOEuler = obtenerResultadosSIR(input,"euler", "Metodo 1: Euler")
      
      
      errores <- c()
      dias <- c()
      
      for(i in 1:length(solucionEDORunge$I-1)){
        errores[i] <- abs((solucionEDOEuler$I[i] - solucionEDORunge$I[i]) / (solucionEDOEuler$I[i]))*100
        dias <- c(dias, i)
      }
      
      plot(dias , errores, col="red", type="l", xlab="Dias Simulacion", ylab="Error", main = "Error de la simulacion", ylim = c(0,100))
    })
    
    output$graficaErrorSIREuler <- renderPlot({
      
      #Experimental
      solucionEDORunge = obtenerResultadosSIR(input,"rk4", "Metodo 1: RungeKutta")
      
      #Teorico
      solucionEDOEuler = obtenerResultadosSIR(input,"euler", "Metodo 1: Euler")
      
      
      errores <- c()
      dias <- c()
      
      for(i in 1:length(solucionEDORunge$I-1)){
        errores[i] <- abs((solucionEDORunge$I[i] - solucionEDOEuler$I[i]) / (solucionEDORunge$I[i]))*100
        dias <- c(dias, i)
      }
      
      plot(dias , errores, col="red", type="l", xlab="Dias Simulacion", ylab="Error", main = "Error de la simulacion", ylim = c(0,100))
    })
    
    
  }
  
  observeEvent(input$botonCalcularSIS, {
    
    calcularSIS()
    
  })
  
  observeEvent(input$botonCalcularSIR, {
    
    calcularSIR()
    
  })
  
}

#================================
# Se se pasan los elementos a la interfaz
#================================
InterfazGrafica <- dashboardPage(title = 'Analisis Numerico', header = barraNavegacion,sidebar = menuLateral, body = contenidoAplicacion, skin='purple')

#================================
# Se ejecuta la aplicación
#================================
shinyApp( ui = InterfazGrafica, server = Servidor)
