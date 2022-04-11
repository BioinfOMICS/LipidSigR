#' @title exp_profiling
#' @description This function three types of distribution plots for users to a simple view of sample variability, and compare the amount/expression difference of lipid between samples (i.e., patients vs. control).
#' \enumerate{
#' \item The first histogram depicts the numbers of lipids expressed in each sample.
#' \item The second histogram illustrates the total amount of lipid in each sample.
#' \item The last density plot visualizes the underlying probability distribution of the lipid expression in each sample (line).
#' }
#' @param exp_data A data frame of predictors, including features (molecule, lipid class, etc.) and their expression of each sample. NAs are not allowed. The name of the first column must be "feature" (lipid species).
#' @return Return 3 plots.
#' \enumerate{
#' \item histogram of number of expressed lipids.
#' \item histogram of the total amount of lipid in each sample.
#' \item density plot of the underlying probability distribution of the lipid expression in each sample.
#' }
#' @export
#' @examples
#' data("profiling_exp_data")
#' exp_data <- profiling_exp_data
#' exp_profiling(exp_data)
exp_profiling <- function(exp_data){
  ################################################
  ####                                        ####
  #### PLOT: total expressed lipid by samples ####
  ####                                        ####
  ################################################

  if(!is(exp_data[,1], 'character')){
    stop("The first column must contain a list of lipids names (features)")
  }
  if(ncol(exp_data)==2){
    if(!is(exp_data[,1], 'character') | sum(class(exp_data[,-1])%in%c("numeric","integer"))!=1){
      stop("First column type must be 'character',others must be 'numeric'")
    }
  }else{
    if(!is(exp_data[,1], 'character') | sum(sapply(exp_data[,-1], class)%in%c("numeric","integer"))!=ncol(exp_data[,-1])){
      stop("First column type must be 'character',others must be 'numeric'")
    }
  }
  if(nrow(exp_data)!=length(unique(exp_data[,1]))){
    stop("The lipids name (features) must be unique")
  }
  if(ncol(exp_data)<3){
    stop("At least 2 samples.")
  }
  if(sum(exp_data[,-1]<0,na.rm = TRUE)>0){
    stop("Variable must greater than zero")
  }
  if(sum(!is.na(exp_data[,-1]))==0 | sum(!is.null(exp_data[,-1]))==0){
    stop("Variables can not be all NULL/NA")
  }
  ## Count total expressed lipids ##
  num.lip <- exp_data %>%
    tidyr::gather(sample_name,value,-feature) %>%
    dplyr::mutate(is.expr = !is.na(value)) %>%
    dplyr::group_by(sample_name) %>%
    dplyr::summarise(expr.count = sum(is.expr)) %>%
    dplyr::arrange(sample_name)
  num.lip$sample_name <- factor(num.lip$sample_name,levels = unique(num.lip$sample_name)[order(num.lip$expr.count, decreasing = TRUE)])

  ## plot barplot ##
  i.expr.lip <- plotly::plot_ly(num.lip,x=~sample_name,y=~expr.count,
                        type="bar",color = ~sample_name,hoverinfo="text",
                        text=~paste("Sample name :",sample_name,"\nCount :",expr.count)) %>%
    plotly::layout(legend = list(title=list(text="Sample"),
                         font = list(family = 'arial')),
           xaxis = list(titlefont = list(family = 'arial'),
                        tickfont = list(family = 'arial'),
                        automargin = TRUE,
                        showgrid = FALSE),
           yaxis = list(title="Number of Expressed Lipids",
                        titlefont = list(family = 'arial'),
                        tickfont = list(family = 'arial'),
                        automargin = TRUE,
                        showgrid = FALSE))
  #############################################
  ####                                     ####
  #### PLOT: total lipid amount of samples ####
  ####                                     ####
  #############################################
  tot.lip <- exp_data %>%
    tidyr::gather(sample_name,value,-feature) %>%
    dplyr::group_by(sample_name) %>%
    dplyr::summarise(lipid_amount=sum(value, na.rm=TRUE))
  tot.lip$sample_name <- factor(tot.lip$sample_name,levels = unique(tot.lip$sample_name)[order(tot.lip$lipid_amount, decreasing = TRUE)])

  ## plot barplot ##
  i.p.amount <- plotly::plot_ly(tot.lip,x=~sample_name,y=~lipid_amount,
                        type="bar",color = ~sample_name,hoverinfo="text",
                        text=~paste("Sample name :",sample_name,"\nAmount :",round(lipid_amount,3))) %>%
    plotly::layout(legend = list(title=list(text="Sample"),
                         font = list(family = 'arial')),
           xaxis = list(titlefont = list(family = 'arial'),
                        tickfont = list(family = 'arial'),
                        automargin = TRUE,
                        showgrid = FALSE),
           yaxis = list(title="Lipid Amount",
                        titlefont = list(family = 'arial'),
                        tickfont = list(family = 'arial'),
                        automargin = TRUE,
                        showgrid = FALSE))

  #################################################
  ####                                         ####
  #### PLOT: hist of quantification per sample ####
  ####                                         ####
  #################################################
  p.hist.data <- exp_data %>%
    tidyr::gather(sample_name,value,-feature)
  p.hist.data <- stats::na.omit(p.hist.data)
  length_sample <- length(unique(p.hist.data$sample_name))
  for (i in seq_len(length_sample)) {
    if(i==1){
      density_sample <- p.hist.data[which(p.hist.data$sample_name==unique(p.hist.data$sample_name)[i]),]
      density_data_x <- data.frame(stats::density(log10(density_sample$value))$x)
      density_data_y <- data.frame(stats::density(log10(density_sample$value))$y)
    }else{
      density_sample <- p.hist.data[which(p.hist.data$sample_name==unique(p.hist.data$sample_name)[i]),]
      density_data_x <- data.frame(density_data_x,stats::density(log10(density_sample$value))$x)
      density_data_y <- data.frame(density_data_y,stats::density(log10(density_sample$value))$y)
    }
  }
  colnames(density_data_x) <- paste0("x",seq_len(length_sample))
  colnames(density_data_y) <- paste0("y",seq_len(length_sample))
  density_data <- cbind(density_data_x,density_data_y)
  ## plot densityplot ##
  for (i in seq_len(length_sample)) {
    if(i==1){
      p.hist.value <- plotly::plot_ly(data = density_data,
                              x=~x1,
                              y=~y1,
                              type = 'scatter',
                              mode = 'lines',
                              name = unique(p.hist.data$sample_name)[1],
                              hoverinfo="text",
                              text=~paste("Sample name :",
                                          unique(p.hist.data$sample_name)[1],
                                          "\nlog10(expression) :",
                                          round(x1,3),
                                          "\nDensity :",
                                          round(y1,3)))
    }else{
      text = paste0("p.hist.value %>% plotly::add_trace(x = ~x",
                    i,
                    ",",
                    "y = ~y",i,
                    ",name = unique(p.hist.data$sample_name)[",i,
                    "],text=~paste(\"Sample name :\",unique(p.hist.data$sample_name)[",i,
                    "],\"\nlog10(expression) :\",round(x",i,
                    ",3),\"\nDensity :\",round(y",i,",3)))")
      p.hist.value <- eval(parse(text=text))
    }
  }
  p.hist.value <- p.hist.value %>%
    plotly::layout(legend=list(title=list(text="Sample")),
           xaxis=list(title = "log10(expression)",zeroline= FALSE,showgrid = FALSE),
           yaxis=list(title = "Density",showgrid = FALSE))

  return(list(i.expr.lip = i.expr.lip,
              i.p.amount = i.p.amount,
              p.hist.value = p.hist.value))
}


