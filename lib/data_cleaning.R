loc_cleaning <- function(geocode_results){
    condition_a <- sapply(geocode_results, function(x) x["status"]=="OK")
    geocode_results <- geocode_results[condition_a]
    
    ### data cleaning
    # source("https://raw.githubusercontent.com/LucasPuente/geocoding/master/cleaning_geocoded_results.R")
    
    for(i in 1:length(geocode_results)){
        dynamic_j<-length(geocode_results[[i]]$results[[1]]$address_components)
        if(dynamic_j > 0) {
            for(j in 1:dynamic_j){
                if(length(geocode_results[[i]]$results[[1]]$address_components[[j]]$types)>2){
                    geocode_results[[i]]$results[[1]]$address_components[[j]]$types<-geocode_results[[i]]$results[[1]]$address_components[[j]]$types[(length(geocode_results[[i]]$results[[1]]$address_components[[j]]$types)-1):length(geocode_results[[i]]$results[[1]]$address_components[[j]]$types)]
                }
            }
            if(length(geocode_results[[i]]$results[[1]]$types)>2){
                geocode_results[[i]]$results[[1]]$types<-geocode_results[[i]]$results[[1]]$types[(length(geocode_results[[i]]$results[[1]]$types)-1):length(geocode_results[[i]]$results[[1]]$types)]
            }
            if(length(geocode_results[[i]]$results[[1]]$types)<1){
                geocode_results[[i]]$results[[1]]$types<-"Unknown"
            }
        }
        
        dynamic_k<-length(geocode_results[[i]]$results[[1]]$address_components)
        if(dynamic_k > 0) {
            for(k in 1:dynamic_k){
                if(length(geocode_results[[i]]$results[[1]]$address_components[[k]]$types)<1){
                    geocode_results[[i]]$results[[1]]$address_components[[k]]$types<-"Unknown"
                }
            }
            if(length(geocode_results[[i]]$results[[1]]$postcode_localities)>2){
                geocode_results[[i]]$results[[1]]$postcode_localities<-geocode_results[[i]]$results[[1]]$postcode_localities[(length(geocode_results[[i]]$results[[1]]$postcode_localities)-1):length(geocode_results[[i]]$results[[1]]$postcode_localities)]
            }
        }
        
    }
    
    #results_b <- lapply(geocode_results, as.data.frame)
    
    
    
    # Then, to simplify things, you should extract out only the columns you need to generate a map of this data.
    out <- data.frame()
    for(i in 1:length(geocode_results)){
        tem <- data.frame(address=geocode_results[[i]]$results[[1]]$formatted_address,
                          lat=geocode_results[[i]]$results[[1]]$geometry$location[['lat']],
                          lng=geocode_results[[i]]$results[[1]]$geometry$location[['lng']]
                          )
        out <- rbind(tem, out)
    }
    
    return(out)
    
}
