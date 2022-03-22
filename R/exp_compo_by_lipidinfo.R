#' @title exp_compo_by_lipidinfo
#' @description This function provides two plots depending on the "Lipid characteristics" table uploaded by users.
#' \enumerate{
#' \item A bar plot of lipid expression, which depicts the expression level of each sample within each group (e.g., PE, PC) of selected characteristics (e.g., class).
#' \item A stacked horizontal bar chart of the percentage of characteristics in each sample.
#' }
#' @param exp_data A data frame of predictors, including features (molecule, lipid class, etc.) and their expression of each sample. NAs are not allowed. The name of the first column must be "feature" (lipid species).
#' @param lipid_char_table A data frame of lipid characteristics such as name(feature) of lipid, class of lipid, the total length of lipid, and Fatty acid (FA_) related characteristics. NAs are allowed. The name of the first column must be "feature" (lipid species).
#' @param char_var A character string of the lipid characteristic selected by users from the column name of \bold{lipid_char_table}, such as 'class'.
#' @return Return 2 plots.
#' \enumerate{
#' \item A bar plot depicts the expression level of each sample within each group (e.g., PE, PC) of selected characteristics (e.g., class).
#' \item A stacked horizontal bar chart of lipid class composition by the user-inputed parameter, \bold{'char_var'} (such as class) .
#' }
#' @export
#' @examples
#' \dontrun{
#' data("profiling_exp_data")
#' exp_data <- profiling_exp_data
#' lipid_char_table <- data("profiling_lipid_char_table")
#' char_var <- colnames(lipid_char_table)[-1]
#' exp_compo_by_lipidinfo(exp_data, lipid_char_table, char_var[1])
#' }
exp_compo_by_lipidinfo <- function(exp_data, lipid_char_table, char_var){

  if(class(exp_data[,1])!="character"){
    stop("The first column must contain a list of lipids names (features)")
  }
  if(ncol(exp_data)==2){
    if(class(exp_data[,1])!="character" | sum(class(exp_data[,-1])%in%c("numeric","integer"))!=1){
      stop("First column type must be 'character',others must be 'numeric'")
    }
  }else{
    if(class(exp_data[,1])!="character" | sum(sapply(exp_data[,-1], class)%in%c("numeric","integer"))!=ncol(exp_data[,-1])){
      stop("First column type must be 'character',others must be 'numeric'")
    }
  }
  if(nrow(exp_data)!=length(unique(exp_data[,1]))){
    stop("The lipids name (features) must be unique")
  }
  if(ncol(exp_data)<3){
    stop("At least 2 samples.")
  }
  if(sum(exp_data[,-1]<0,na.rm = T)>0){
    stop("Variable must greater than zero")
  }
  if(sum(!is.na(exp_data[,-1]))==0 | sum(!is.null(exp_data[,-1]))==0){
    stop("Variables can not be all NULL/NA")
  }
  if(nrow(lipid_char_table)==nrow(exp_data)){
    if(sum(lipid_char_table[,1]%in%exp_data[,1])!=nrow(lipid_char_table)){
      stop("The lipids names (features) of lipid_char_table table must same as exp_data.")
    }
  }else{
    stop("The row number of lipid_char_table table must same as exp_data.")
  }
  if(class(lipid_char_table[,1])!="character"){
    stop("lipid_char_table first column must contain a list of lipids names (features).")
  }
  if(nrow(lipid_char_table)!=length(unique(lipid_char_table[,1]))){
    stop("lipid_char_table lipids names (features) must be unique.")
  }
  if("class" %in%colnames(lipid_char_table) & class(lipid_char_table[,'class'])!="character"){
    stop("lipid_char_table content of column 'class' must be characters")
  }
  if("totallength" %in%colnames(lipid_char_table) & !class(lipid_char_table[,'totallength'])%in%c("integer","numeric")){
    stop("lipid_char_table content of column 'totallength' must be numeric")
  }
  if("totaldb" %in%colnames(lipid_char_table) & !class(lipid_char_table[,'totaldb'])%in%c("integer","numeric")){
    stop("lipid_char_table content of column 'totaldb' must be numeric")
  }
  if("totaloh" %in%colnames(lipid_char_table) & !class(lipid_char_table[,'totaloh'])%in%c("integer","numeric")){
    stop("Thlipid_char_tablee content of column 'totaloh' must be numeric")
  }
  if(!char_var %in% colnames(lipid_char_table)){
      stop("char_var must be included in the lipid_char_table.")
  }
  if(ncol(dplyr::select(lipid_char_table,tidyselect::starts_with("FA_")))==0){
    warning("(OPTIONAL) lipid_char_table does not contain column names starting with 'FA_'")
  }else{
    FA_lipid_char_table <- lipid_char_table %>%dplyr::select(feature,tidyselect::starts_with("FA_"))
    FA_col <- grep("FA_",colnames(FA_lipid_char_table),value = TRUE)
    max_comma <- 0
    for(i in 1:length(FA_col)){
      col <- FA_col[i]
      comma_count <- max(stringr::str_count(FA_lipid_char_table[,col], ','),na.rm = T)
      if(comma_count>0){
        FA_lipid_char_table <- tidyr::separate(FA_lipid_char_table,col,c(col,paste0(col,"_",1:comma_count)),",", convert = TRUE)
      }
      if(comma_count>max_comma){max_comma <- comma_count}
    }
    FA_lipid_char_table <- FA_lipid_char_table %>% tidyr::gather(lipid.category, lipid.category.value,-feature)
    if(max_comma>0){
      for (i in 1:max_comma) {
        select_name <- paste0("_",i)
        FA_lipid_char_table <-FA_lipid_char_table[-intersect(grep(select_name,FA_lipid_char_table[,"lipid.category"]),which(is.na(FA_lipid_char_table$lipid.category.value))),]
      }
    }
    if(class(FA_lipid_char_table$lipid.category.value)=="character"){
      stop("In the 'FA_' related analyses, the values are positive integer or zero and separated by comma. i.e., 10,12,11")
    }else if(sum(stats::na.omit(as.numeric(FA_lipid_char_table$lipid.category.value))!=round(stats::na.omit(as.numeric(FA_lipid_char_table$lipid.category.value))))!=0 | min(stats::na.omit(as.numeric(FA_lipid_char_table$lipid.category.value)))<0){
      stop("In the 'FA_' related analyses, the values are positive integer or zero and separated by comma. i.e., 10,12,11")
    }
  }

  FA_col <- grep("FA_",colnames(lipid_char_table),value = TRUE)
  if(length(FA_col)>0){
    max_comma <- 0
    for(i in 1:length(FA_col)){
      col <- FA_col[i]
      comma_count <- max(stringr::str_count(lipid_char_table[,col], ','),na.rm = T)
      if(comma_count>0){
        lipid_char_table <- tidyr::separate(lipid_char_table,col,c(col,paste0(col,"_",1:comma_count)))
      }
      if(comma_count>max_comma){max_comma <- comma_count}
    }
    lipid_char_table <- lipid_char_table %>% tidyr::gather(lipid.category, lipid.category.value,-feature)
    if(max_comma>0){
      for (i in 1:max_comma) {
        select_name <- paste0("_",i)
        lipid_char_table <-lipid_char_table[-intersect(grep(select_name,lipid_char_table[,"lipid.category"]),which(is.na(lipid_char_table$lipid.category.value))),]
      }
      for(i in 1:length(FA_col)){
        col <- FA_col[i]
        lipid_char_table[grep(col,lipid_char_table[,"lipid.category"]),"lipid.category"]<-col
      }
    }

  }else{
    lipid_char_table <- lipid_char_table %>% tidyr::gather(lipid.category, lipid.category.value,-feature)
  }
  #######################################################
  ####                                                 ####
  #### PLOT: barplot of expression of each lipid class ####
  ####                                                 ####
  #########################################################


  p.barplot.p <-merge(exp_data %>% tidyr::gather(sample_name,value,-feature),
                      lipid_char_table,
                      by="feature") %>%
    dplyr::filter(lipid.category == char_var & !is.na(lipid.category.value)) %>%
    dplyr::group_by(sample_name, lipid.category, lipid.category.value) %>%
    dplyr::summarise(value=sum(value,na.rm=T)) %>%
    ungroup()%>%
    dplyr::group_by(sample_name) %>%
    dplyr::mutate(weight = 100/sum(value)) %>%
    dplyr::mutate(value=value*weight)%>%
    ungroup %>%
    plotly::plot_ly(x=~lipid.category.value,
            y=~value,
            color = ~sample_name,
            type="bar",
            hoverinfo="text",
            text=~paste(Hmisc::capitalize(char_var)," :",lipid.category.value,
                        "\nSample name :",sample_name,
                        "\nLipid Expression :",round(value,3))) %>%
    plotly::layout(legend = list(title=list(text="Sample name")),
           xaxis = list(title=Hmisc::capitalize(char_var)),
           yaxis = list(title='Lipid Expression',nticks=15),
           title = paste0('Lipid ', Hmisc::capitalize(char_var)))

  ######################################
  ####                              ####
  #### PLOT: lipid composition plot ####
  ####                              ####
  ######################################

  p.compos <- merge(exp_data %>% tidyr::gather(sample_name,value,-feature),
                    lipid_char_table,
                    by="feature") %>%
    dplyr::filter(lipid.category == char_var & !is.na(lipid.category.value)) %>%
    dplyr::group_by(sample_name, lipid.category, lipid.category.value) %>%
    dplyr::summarise(value=sum(value,na.rm=T)) %>%
    ungroup()%>%
    dplyr::group_by(sample_name) %>%
    dplyr::mutate(weight = 100/sum(value)) %>%
    dplyr::mutate(value=value*weight)%>%
    ungroup %>%
    plotly::plot_ly(x=~value,
            y=~sample_name,
            color = ~lipid.category.value,
            type="bar",
            hoverinfo="text",
            orientation = "h",
            text=~paste(Hmisc::capitalize(char_var),":",lipid.category.value,
                        "\nSample name :",sample_name,
                        "\nLipid Expression :",round(value,3))) %>%
    plotly::layout(barmode = 'stack',
           legend = list(title=list(text=Hmisc::capitalize(char_var))),
           xaxis = list(title="%",nticks=15),
           yaxis = list(title=''),
           title = paste0('Lipid ', Hmisc::capitalize(char_var), ' Composition'))
  return(list(p.barplot.p = p.barplot.p,
              p.compos = p.compos))
}






