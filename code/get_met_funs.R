get_met_edas40 <- function(files = NULL,
                           years = NULL,
                           months = NULL,
                           path_met_files) {
  
  edas40_dir <- "ftp://arlftp.arlhq.noaa.gov/archives/edas40/"
  
  # Download the 'listing' file from NOAA server
  # It contains a list of EDAS40 files currently
  # available on the server
  download(
    url = paste0(edas40_dir, "listing"),
    destfile = paste0(getwd(), "/listing"),
    method = "auto",
    quiet = TRUE,
    mode = "wb",
    cacheOK = TRUE
  )
  
  edas40_listing <- readLines(paste0(getwd(), "/listing"))
  
  edas40_listing <- gsub(" ", "", edas40_listing)
  
  if (!is.null(years)) {
    if (length(years) > 1) {
      years <- seq.default(years[1], years[2])
    }
    
    years <- substr(years, 3, 4)
    
    for (i in 1:length(years)) {
      
      if (i == 1) {
        edas40_file_list <-vector(mode = "character")
      }
      
      edas40_file_list <-
        c(
          edas40_file_list,
          edas40_listing[
            which(
              grepl(
                paste0("[a-z][a-z][a-z]", years[i]),
                edas40_listing)
            )]
        )
    }
  }
  
  if (!is.null(months)) {
    
    months_3_letter <- 
      c("jan", "feb", "mar", "apr",
        "may", "jun", "jul", "aug",
        "sep", "oct", "nov", "dec")
    
    for (i in 1:length(months)) {
      if (i == 1) {
        edas40_file_list_month <- vector(mode = "character", length = 0)
      }
      
      edas40_file_list_month <-
        c(
          edas40_file_list_month,
          edas40_file_list[
            which(
              grepl(
                paste0("edas.", months_3_letter[months[i]]),
                edas40_file_list
              )
            )]
        )
    }
    
    edas40_file_list <- edas40_file_list_month
  }
  
  if (!is.null(files)) {
    edas40_file_list <- files
  }
  
  for (i in 1:length(edas40_file_list)) {
    download(
      url = paste0(edas40_dir, edas40_file_list[i]),
      destfile = paste0(path_met_files, edas40_file_list[i]),
      method = "auto",
      quiet = FALSE,
      mode = "wb",
      cacheOK = FALSE
    )
  }
}

get_met_era5 <- function(days,
                         duration,
                         direction,
                         path_met_files) {
  
  get_daily_filenames(
    days = days,
    duration = duration,
    direction = direction,
    suffix = ".ARL",
    prefix = "ERA5_"
  ) #%>%
  #get_met_files(
  #  path_met_files = path_met_files,
  #  ftp_dir = "ftp://arlftp.arlhq.noaa.gov/archives/gdas0p5"
  #)
}

get_met_forecast_nam <- function(path_met_files) {
  
  # Establish which forecast dirs are currently 
  # available on the server
  forecast_dirs <-
    grep(
      "^[0-9].*$",
      unlist(
        strsplit(
          RCurl::getURL(
            "ftp://arlftp.arlhq.noaa.gov/forecast/",
            dirlistonly = TRUE),
          "\n")
      ), 
      value = TRUE
    )
  
  # Get today's date and write in format equivalent
  # to FTP directories 
  today <- gsub("-", "", Sys.Date())
  
  # Download today's `namf` file
  # -- CONUS, 12 km, 3 hrly, pressure levels, 48 h forecast
  
  if (today %in% forecast_dirs) {
    
    downloader::download(
      url = paste0(
        "ftp://arlftp.arlhq.noaa.gov/forecast/",
        today, "/hysplit.t00z.namf"),
      destfile = paste0(path_met_files, paste0(today, ".t00z.namf")),
      method = "auto",
      quiet = FALSE,
      mode = "w",
      cacheOK = TRUE
    )
  } 
}

get_met_gdas0p5 <- function(days,
                            duration,
                            direction,
                            path_met_files) {
  
  get_daily_filenames(
    days = days,
    duration = duration,
    direction = direction,
    suffix = "_gdas0p5"
  ) %>%
    get_met_files(
      path_met_files = path_met_files,
      ftp_dir = "ftp://arlftp.arlhq.noaa.gov/archives/gdas0p5"
    )
}

get_met_gdas1 <- function(days,
                          duration,
                          direction,
                          path_met_files) {
  
  # Determine the minimum date (as a `Date`) for the model run
  if (direction == "backward") {
    min_date <- 
      (lubridate::as_date(days[1]) - (duration / 24)) %>%
      lubridate::floor_date(unit = "day")
  } else if (direction == "forward") {
    min_date <- 
      (lubridate::as_date(days[1])) %>%
      lubridate::floor_date(unit = "day")
  }
  
  # Determine the maximum date (as a `Date`) for the model run
  if (direction == "backward") {
    max_date <- 
      (lubridate::as_date(days[length(days)])) %>%
      lubridate::floor_date(unit = "day")
  } else if (direction == "forward") {
    max_date <- 
      (lubridate::as_date(days[length(days)]) + (duration / 24)) %>%
      lubridate::floor_date(unit = "day")
  }
  
  met_days <- 
    seq(min_date, max_date, by = "1 day") %>% 
    lubridate::day()
  
  # "en_US.UTF-8" is not a valid locale for windows.
  os_for_locale <- get_os()
  if(os_for_locale == "win"){
    month_names <- 
      seq(min_date, max_date, by = "1 day") %>%
      lubridate::month(label = TRUE, abbr = TRUE, locale = Sys.setlocale("LC_TIME", "English"))  %>%
      as.character() %>%
      tolower()
  } else {
    month_names <- 
      seq(min_date, max_date, by = "1 day") %>%
      lubridate::month(label = TRUE, abbr = TRUE, locale = "en_US.UTF-8")  %>%
      as.character() %>%
      tolower()
  }
  
  met_years <- 
    seq(min_date, max_date, by = "1 day") %>%
    lubridate::year() %>% 
    substr(3, 4)
  
  # Only consider the weeks of the month we need:
  #.w1 - days 1-7
  #.w2 - days 8-14
  #.w3 - days 15-21
  #.w4 - days 22-28
  #.w5 - days 29 - rest of the month 
  
  met_week <- ceiling(met_days / 7)
  
  files <- paste0("gdas1.", month_names, met_years, ".w", met_week) %>% unique()
  
  get_met_files(
    files = files,
    path_met_files = path_met_files,
    ftp_dir = "ftp://arlftp.arlhq.noaa.gov/archives/gdas1"
  )
  
}

get_met_gfs0p25 <- function(days,
                            duration,
                            direction,
                            path_met_files) {
  
  get_daily_filenames(
    days = days,
    duration = duration,
    direction = direction,
    suffix = "_gfs0p25"
  ) %>%
    get_met_files(
      path_met_files = path_met_files,
      ftp_dir = "ftp://arlftp.arlhq.noaa.gov/archives/gfs0p25"
    )
}

get_met_hrrr <- function(files = NULL,
                         path_met_files) {
  
  ftp_dir <- "ftp://arlftp.arlhq.noaa.gov/pub/archives/hrrr/"
  
  # Download list of reanalysis met files by name
  if (!is.null(files)) {
    
    for (i in 1:length(files)) {
      download(
        url = file.path(ftp_dir, files[i]),
        destfile = file.path(path_met_files, files[i]),
        method = "auto",
        quiet = FALSE,
        mode = "wb",
        cacheOK = FALSE) 
    }
  }
}

get_met_nam12 <- function(days,
                          duration,
                          direction,
                          path_met_files) {
  
  get_daily_filenames(
    days = days,
    duration = duration,
    direction = direction,
    suffix = "_nam12"
  ) %>%
    get_met_files(
      path_met_files = path_met_files,
      ftp_dir = "ftp://arlftp.arlhq.noaa.gov/archives/nam12"
    )
}

get_met_narr <- function(days,
                         duration,
                         direction,
                         path_met_files) {
  
  get_monthly_filenames(
    days = days,
    duration = duration,
    direction = direction,
    prefix = "NARR"
  ) %>%
    get_met_files(
      path_met_files = path_met_files,
      ftp_dir = "ftp://arlftp.arlhq.noaa.gov/archives/narr"
    )
}

get_met_reanalysis <- function(days,
                               duration,
                               direction,
                               path_met_files) {
  
  get_monthly_filenames(
    days = days,
    duration = duration,
    direction = direction,
    prefix = "RP",
    extension = ".gbl"
  ) %>%
    get_met_files(
      path_met_files = path_met_files,
      ftp_dir = "ftp://arlftp.arlhq.noaa.gov/archives/reanalysis"
    )
}

get_met_nams <- function(days,
                         duration,
                         direction,
                         path_met_files) {
  
  get_daily_filenames(
    days = days,
    duration = duration,
    direction = direction,
    suffix = "_hysplit.t00z.namsa"
  ) %>%
    get_met_files(
      path_met_files = path_met_files,
      ftp_dir = "ftp://arlftp.arlhq.noaa.gov/archives/nams"
    )
}