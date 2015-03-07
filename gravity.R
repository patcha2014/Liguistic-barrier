
if (Sys.info()["user"]=="Mint") {
  work.dir <- "/Users/Mint/Dropbox/LiguisticBarrierToTrade/data/"
}

library("foreign")
library("AER")

#-----------------------------------------
# Distance
#-----------------------------------------

# Record latitude and longitude for Indian states -----------------------

name = c("Andhra", "Assam", "Arunachal Pradesh", "Bihar", "Chattishgarh",
"Chandigarh", "Delhi", "Goa", "Gujarat", "Haryana", 
"Himachal Pradesh", "Jammu and Kasmir", "Jharkhand", "Karnataka", "Kerala", 
"Madhya Pradesh", "Maharashtra", "Manipur", "Meghalaya", "Mizoram",
"Nagaland", "Orissa", "Pondicheri and Karikal", "Punjab", "Rajasthan", 
"Tamil Nadu", "Tripura", "Uttaranchal", "Uttar Pradesh", "West Bengal")
# lat = north+ south- 
lat  = c( 16.5000, 26.1400, 27.0600, 25.3700, 21.2700, 
          30.7500, 28.6100, 15.4989, 23.2167, 30.7300, 
          31.1033, 33.4500, 23.3500, 12.9702, 8.5074, 
          23.2500, 18.9600, 24.8170, 25.5700, 23.3600,
          25.6700, 20.1500, 11.9310, 30.7900, 26.5727, 
          13.0900, 23.8400, 30.3300, 26.8500, 22.5667)
# lon = east+ west-
lon = c( 80.6400, 91.7700, 93.3700, 85.1300, 81.6000, 
          76.7800, 77.2300, 73.8278, 72.6833, 76.7800,
          77.1722, 76.2400, 85.3300, 77.5603, 76.9730, 
          77.4170, 72.8200, 93.9500, 91.8800, 92.0000, 
          94.1200, 85.5000, 79.7852, 76.7800, 73.8390, 
          80.2700, 91.2800, 78.0600, 80.9100, 88.3667) 
location.df <- data.frame(name,lat,lon)

# Function to calculate distance --------------------------------------------
#from http://www.statmethods.net/management/userfunctions.html

# Function to 
ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  
  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
  
  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints
  
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
  
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  
  return(mat.distances)
}

# Compute distance matrix --------------------------------------------
distance_m <- round(GeoDistanceInMetresMatrix(location.df) / 1000)

# Trade-pair vector (lower triangular pairs) -------------------------
# With help from Travis D. McArthur 

distance.df <- data.frame(t(combn(rownames(distance_m),2)), as.numeric(distance_m[lower.tri(distance_m)]))
names(distance.df) <- c("c1", "c2", "distance")

# Check if we get the right data 
for ( i in rownames(distance.df)) {
 for ( j in colnames(distance.df)) {
    if ( i == j) { next }
    # Report error if data cell is empty 
    stopifnot(length(distance.df[(distance.df$c1==i & distance.df$c2==j) |  
                                 (distance.df$c1==j & distance.df$c2==i), "distance"]) == 1) 
    # Report error if distance for each pair does NOT match data from distance matrix 
    stopifnot(distance_m[i, j] == distance.df[(distance.df$c1==i & distance.df$c2==j) | 
                                              (distance.df$c1==j & distance.df$c2==i), "distance"])
    cat(i, j, "\n")
  }
}

#-----------------------------------------
# Language commonality  
#-----------------------------------------

# 1. Common official language -------------------------------------------

# Record official languages for each state 
# From source[1] http://www.translationdirectory.com/articles/article2474.php
# Cross-checked with source[2] http://en.wikipedia.org/wiki/Languages_with_official_status_in_India
# accessed on (3/7/15)
o_lang.df <- read.table(text = "", 
                        colClasses = c("character", "character", "character"), 
                        col.names = c("name", "o_lang1", "o_lang2", "o_lang3") )

o_lang.df[1, ] <- c("Andhra", "Telugu", "Urdu", "")
o_lang.df[2, ] <- c("Assam", "Assamese", "Bengali", "Bodo")
o_lang.df[3, ] <- c("Arunachal Pradesh", "English", "", "")
o_lang.df[4, ] <- c("Bihar", "Hindi", "Urdu", "")
o_lang.df[5, ] <- c("Chattishgarh", "Chattisgarhi", "Hindi", "") 
o_lang.df[6, ] <- c("Chandigarh", "Punjabi", "Hindi", "English") 
#Chandigarh: from http://en.wikipedia.org/wiki/Chandigarh
o_lang.df[7, ] <- c("Delhi", "Punjabi", "Urdu", "Hindi")
#Delhi: from http://en.wikipedia.org/wiki/Delhi
o_lang.df[8, ] <- c("Goa", "Konkani", "Marathi", "English") 
#Goa: lang2-3 only in source[2]
o_lang.df[9, ] <- c("Gujarat", "Gujarati", "", "")
o_lang.df[10, ] <- c("Haryana", "Haryanvi", "Punjabi", "")
o_lang.df[11, ] <- c("Himachal Pradesh", "Hindi", "", "")
o_lang.df[12, ] <- c("Jammu and Kasmir", "Urdu", "English", "")
o_lang.df[13, ] <- c("Jharkhand", "Hindi", "Oriya", "Bengali") 
#Jharkhand: also have Santali from source[1]. Lang2 only from source[2]
o_lang.df[14, ] <- c("Karnataka", "Kannada", "", "")
o_lang.df[15, ] <- c("Kerala", "Malayalam", "English", "")
o_lang.df[16, ] <- c("Madhya Pradesh", "Hindi", "", "")
o_lang.df[17, ] <- c("Maharashtra", "Marathi", "", "")
o_lang.df[18, ] <- c("Manipur", "Meiteilon", "Manipuri", "English") 
# Manipur: Manipuri aka Meiteilon 
o_lang.df[19, ] <- c("Meghalaya", "English", "Khasi", "Garo")
o_lang.df[20, ] <- c("Mizoram", "Mizo", "", "")
o_lang.df[21, ] <- c("Nagaland", "English", "", "")
o_lang.df[22, ] <- c("Orissa", "Oriya", "", "") 
# Orissa: State name appeared as Odisha 
o_lang.df[23, ] <- c("Pondicheri and Karikal", "Tamil", "Telugu", "Malayalam") 
# Karaikal: http://en.wikipedia.org/wiki/Karaikal
# Pondicheri: http://en.wikipedia.org/wiki/Pondicherry
# Also have French 
o_lang.df[24, ] <- c("Punjab", "Punjabi", "", "")
o_lang.df[25, ] <- c("Rajasthan", "Hindi", "", "")
o_lang.df[26, ] <- c("Tamil Nadu", "Tamil", "English", "")
o_lang.df[27, ] <- c("Tripura", "Bengali", "Kokborok", "")
o_lang.df[28, ] <- c("Uttaranchal", "Hindi", "Sanskrit", "")
o_lang.df[29, ] <- c("Uttar Pradesh", "Hindi", "Urdu", "")
o_lang.df[30, ] <- c("West Bengal", "Bengali", "English", "Urdu") 
#Also have Nepali. lang3 only from source[2]

# Create common language(s) indicator 


# 2. Common spoken language ----------------------------------------------


