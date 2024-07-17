

### WILL HAVE TO FIX THIS ###
get_monthly_filenames <- function(days,
                                  duration,
                                  direction,
                                  prefix = NULL,
                                  extension = NULL) {
  
  # Determine the minimum month (as a `Date`) for the model run
  if (direction == "backward") {
    min_month <-
      (lubridate::as_date(days) - (duration / 24) - lubridate::days(1)) %>%
      min() %>%
      lubridate::floor_date(unit = "month")
  } else if (direction == "forward") {
    min_month <-
      (lubridate::as_date(days)) %>%
      min() %>%
      lubridate::floor_date(unit = "month")
  }
  
  # Determine the maximum month (as a `Date`) for the model run
  if (direction == "backward") {
    max_month <-
      (lubridate::as_date(days)) %>%
      max() %>%
      lubridate::floor_date(unit = "month")
  } else if (direction == "forward") {
    max_month <-
      (lubridate::as_date(days) + (duration / 24) + lubridate::days(1)) %>%
      max() %>%
      lubridate::floor_date(unit = "month")
  }
  
  met_months <- seq(min_month, max_month, by = "1 month")
  
  months_short <- met_months %>% to_short_month()
  
  years_long <- lubridate::year(met_months)
  
  paste0(prefix, years_long, months_short, extension)
}


get_daily_filenames <- function(days,
                                duration,
                                direction,
                                prefix = NULL,
                                suffix = NULL) {
  
  lower_days <- lubridate::as_date(days) - 
    (duration/24) - 
    lubridate::days(1)
  
  upper_days <- lubridate::as_date(days) + 
    (duration/24) + 
    lubridate::days(1)
  
  if (direction == 'backward') {
    met_days <- union(days, lower_days)
    
  }
  else if (direction == 'forward') {
    met_days <- union(days, upper_days)
  }
 
  met_days <- met_days %>%
    as.Date() %>%
    sort() %>%
    as.character() %>%
    tidy_gsub("-", "")
  
  paste0(prefix, met_days, suffix)
}


get_met_files <- function(files, path_met_files, ftp_dir) {
  
  # Determine which met files are already locally available
  files_in_path <- list.files(path_met_files)
  
  # Download list of GFS0.25 met files by name
  if (!is.null(files)) {
    
    for (file in files) {
      
      if (!(file %in% files_in_path)) {
        
        downloader::download(
          url = file.path(ftp_dir, file),
          destfile = path.expand(file.path(path_met_files, file)),
          method = "auto",
          quiet = FALSE,
          mode = "wb",
          cacheOK = FALSE
        )
      }
    }
  }
  
  files
}

source('code/get_met_funs.R')

download_met_files <- function(days,
                               duration,
                               direction,
                               met_dir,
                               met_type = 'nams') {
  options(timeout = 800)
  
  if (met_type == "gdas1") {
    
    met_files <-
      get_met_gdas1(
        days = days,
        duration = duration,
        direction = direction,
        path_met_files = met_dir
      )
  }
  
  if (met_type == "gdas0.5") {
    
    met_files <-
      get_met_gdas0p5(
        days = days,
        duration = duration,
        direction = direction,
        path_met_files = met_dir
      )
  }
  
  if (met_type == "gfs0.25") {
    
    met_files <-
      get_met_gfs0p25(
        days = days,
        duration = duration,
        direction = direction,
        path_met_files = met_dir
      )
  }
  
  if (met_type == "reanalysis") {
    
    met_files <-
      get_met_reanalysis(
        days = days,
        duration = duration,
        direction = direction,
        path_met_files = met_dir
      )
  }
  
  if (met_type == "nam12") {
    
    met_files <-
      get_met_nam12(
        days = days,
        duration = duration,
        direction = direction,
        path_met_files = met_dir
      )
  }
  
  if (met_type == "nams") {
    
    met_files <-
      get_met_nams(
        days = days,
        duration = duration,
        direction = direction,
        path_met_files = met_dir
      )
  }
  
  met_files
}
