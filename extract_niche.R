setwd("C:/Users/Owner/Documents/data_analysis/niche_dynamics")

# Here are the main steps of the script:
# 1) spatial thinning of the occurrences (= spatial disaggregation)
# 2) delineation of the background area. Basically, the background area consists in 
# the biomes where the species is distributed. I splitted the invasive distribution 
# for each biogeographic realm
# 3) extraction of the environmental data
# 4) principal component analysis on the environmental data
# 5) niche comparison
# 6) figures and tables of the results



# load libraries and data -------------------------------------------------

library(geodata)
library(dismo)
library(ecospat)
library(terra)
library(sf)
library(dplyr)
library(mapview)
library(ggplot2)
library(RColorBrewer)
library(rnaturalearth)
library(ade4)



## set parameters ----
sp_list<-c('spermacoce', 
           'chromolaena',
           'prosopis',
           'mimosa',
           'sida')
th_ind = 0.01666667 # minimal distance between two points for spatial thinning
min_obs <- 5 # minimum observation per biome to be included in the analysis
n_bck <- 10000 # number of background points
sf_use_s2(FALSE) # import a map of the world
world <- rnaturalearth::ne_countries( returnclass = "sf") %>% 
  st_make_valid() %>% 
  st_union
sf_use_s2(TRUE)


for (spi in 1:length(sp_list)){ # loop accross all the species
  ## load species data ----
  sp_name<-sp_list[spi]
  
  # native
  sp_nat <- read.csv(paste0('data/sp_dist',spi,'/',sp_name,'_native.csv'))  #read data
  coordinates(sp_nat) <-  ~ long + lat # convert into spatial sp object
  sp_nat <- sp_nat %>%
    remove.duplicates(zero = th_ind) %>%  # remove close and aggregated observations
    st_as_sf () %>%   # convert into spatial sf object
    st_set_crs(4326) # set coordinate reference system
  
  #exotic
  sp_nov <- read.csv(paste0('data/sp_dist',spi,'/',sp_name,'_novel.csv'))  #read data
  coordinates(sp_nov) <-  ~ long + lat # convert into spatial sp object
  sp_nov <- sp_nov %>%
    sp::remove.duplicates(zero = th_ind) %>%  # remove close and aggregated observations
    st_as_sf () %>%   # convert into spatial sf object
    st_set_crs(4326) # set coordinate reference system
  
  # visualisation
  mapview::mapview (sp_nov, zcol = 'presence', col.region = 'red') +
    mapview::mapview (sp_nat, zcol = 'presence', col.region = 'darkgreen') # visualize distribution
  
  ggplot() +
    geom_sf(data = world)  +
    geom_sf(data = sp_nat,
            color = 'darkgreen')+
    geom_sf(data = sp_nov,
            color = 'darkred')# plot native distribution
  
  ggsave(file = paste0('data/sp_dist',spi,'/',sp_name,'.png'))
  
  
  
  ## load background data ----------------------------------------------------
  
  # sf_use_s2(TRUE)
  
  # # make the origninal shapefile valid
  # ecoreg <-
  #   sf::st_read('data/ecoregions/Ecoregions2017.shp') %>% # biomes and ecoregions from https://www.gislounge.com/terrestrial-ecoregions-gis-data/
  #   sf::st_make_valid() %>%
  #   sf::st_transform(st_crs(sp_nov)) %>%
  #   sf::st_write('data/ecoregions/Ecoregions2017_valid.shp')
  
  ecoreg <-
    sf::st_read('data/ecoregions/Ecoregions2017_valid.shp') # biomes and ecoregions clean and valid version
  
  ### generate extent and background data ------------------------------------------------
  
  #### for the native range ----
  my_dist <- sp_nat
  sp_range <- 'native'
  
  dist_cnt <-
    sf::st_intersects(ecoreg, my_dist) # intersects observations and ecoregions
  ecoreg_dist <- ecoreg
  ecoreg_dist$Nobs <- lengths(dist_cnt) # count points per ecoregion
  ecoreg_dist <-
    dplyr::mutate(ecoreg_dist, name = paste(REALM, BIOME_NUM)) # create a category name
  
  RealmxBiome <- sf::st_drop_geometry(ecoreg_dist) %>%  # count table
    dplyr::filter(Nobs > 0) %>% # select ecoregions with points
    dplyr::group_by(REALM, BIOME_NUM) %>% # group by biomes and realms
    dplyr::summarise(
      count = sum(Nobs),
      BIOME_NAME = unique(BIOME_NAME),
      name = unique (name)
    ) %>% # summarise by biomes and realms
    dplyr::arrange(dplyr::desc(count)) # sort descending
  
  RealmxBiome # Check result
  
  ggplot2::ggplot(data = RealmxBiome,
                  ggplot2::aes (
                    x = reorder(name, -count),
                    y = count,
                    fill = paste(REALM, BIOME_NUM)
                  )) +
    ggplot2::geom_bar(stat = 'identity') +
    ggplot2::labs(x = "Biomes", fill = "Biomes")# plot results
  
  if (!dir.exists(paste0('results/distributions/',sp_name))){# save plot
    dir.create(paste0('results/distributions/',sp_name),rec = T)
  }
  ggsave(file = paste0('results/distributions/',sp_name,'/native_biomes.png'))
  
  
  RealmxBiome <-
    filter(RealmxBiome, count >= min_obs) # keep biomes with 5 or more observations
  
  ecoreg_dist <- dplyr::filter(ecoreg_dist,
                               name %in% RealmxBiome$name) %>%  # select extent
    dplyr::group_by(REALM, BIOME_NUM) %>% # group by biomes and realms
    dplyr::summarise(
      count = sum(Nobs),
      BIOME_NAME = unique(BIOME_NAME),
      geometry = st_union (geometry)
    )  # summarise by biomes and realms
  
  my_dist <-
    my_dist[lengths(sf::st_intersects(my_dist, ecoreg_dist)) > 0, ]# keep only observation in the extent
  
  # my_dist_bck <-
  #   sf::st_sample(ecoreg_dist, size = n_bck)# generate background points in the extent
  
  # sampling strategy to avoid bug with st_sample
  grid.extent <- sf::st_bbox(ecoreg_dist) # set the grid extent for sampling
  env_grid <- terra::rast(
    res = 0.01666667,
    crs = "epsg:4326",
    xmin = grid.extent[1],
    ymin = grid.extent[2],
    xmax = grid.extent[3],
    ymax = grid.extent[4],
    vals = NA
  ) # create empty raster resolution is twice the env_grid to speed the process and lower the spatial autocorrelation
  env_grid<-terra::rasterize(x = ecoreg_dist,
                             y = env_grid,
                             fun = 'min',
                             field =1) # rasterize the extent for the background
  my_dist_bck<-which(values(env_grid)== 1) # select the cell with extent
  my_dist_bck<-terra::xyFromCell(env_grid,my_dist_bck) # convert raster to xy data frame
  if (nrow(my_dist_bck)>n_bck){ # sample the background (N max = 10000)
    my_dist_bck<-my_dist_bck[sample(1:nrow(my_dist_bck),n_bck),]
  }else{
    my_dist_bck<-my_dist_bck[1:nrow(my_dist_bck),]
  }
  colnames(my_dist_bck)<-c('X','Y')
  
  if (!dir.exists(paste0('data/sp_dist',spi,'/', sp_name))) {
    dir.create(paste0('data/sp_dist',spi,'/', sp_name), rec = TRUE)
  }
  
  write.table(
    sf::st_coordinates(my_dist),
    file = paste0('data/sp_dist',spi,'/',
                  sp_name, '/',
                  sp_range, '.txt'),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
  ) # write sp_dist
  
  write.table(
    my_dist_bck,
    file = paste0('data/sp_dist',spi,'/',
                  sp_name, '/', sp_range, '_bck_', '.txt'),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
  ) # write sp background
  
  sp_bbox<-st_bbox(st_bbox(ecoreg_dist))*1.05 # extent to plot
  ggplot() +
    geom_sf(data = world)  +
    geom_sf(data = st_union(ecoreg_dist),
            fill = 'lightyellow')+
    geom_sf(data = my_dist,
            color = 'darkgreen')+ 
    coord_sf(
      xlim = sp_bbox[c(1,3)],
      ylim = sp_bbox[c(2,4)],
      expand =FALSE
    ) # plot native distribution
  
  if (!dir.exists(paste0('results/distributions/',sp_name))){# save plot
    dir.create(paste0('results/distributions/',sp_name),rec = T)
  }
  ggsave(file = paste0('results/distributions/',sp_name,'/native_distributions.png'))
  
  
  #### for the novel range ----
  
  my_dist <- sp_nov
  sp_range <- 'novel'
  
  dist_cnt <-
    sf::st_intersects(ecoreg, my_dist) # intersects observations and ecoregions
  ecoreg_dist <- ecoreg
  ecoreg_dist$Nobs <- lengths(dist_cnt) # count points per ecoregion
  ecoreg_dist <-
    dplyr::mutate(ecoreg_dist, name = paste(REALM, BIOME_NUM)) # create a category name
  
  RealmxBiome <- sf::st_drop_geometry(ecoreg_dist) %>%  #count table
    dplyr::filter(Nobs > 0) %>% #select ecoregions with points
    dplyr::group_by(REALM, BIOME_NUM) %>% # groupb by biomes and realms
    dplyr::summarise(
      count = sum(Nobs),
      BIOME_NAME = unique(BIOME_NAME),
      name = unique (name)
    ) %>% # summarise by biomes and realms
    dplyr::arrange(dplyr::desc(count)) # sort descending
  
  RealmxBiome # Check result
  
  ggplot2::ggplot(data = RealmxBiome,
                  ggplot2::aes (
                    x = reorder(name, -count),
                    y = count,
                    fill = paste(REALM, BIOME_NUM)
                  )) +
    ggplot2::geom_bar(stat = 'identity') +
    ggplot2::labs(x = "Biomes", fill = "Biomes") # plot results
  
  if (!dir.exists(paste0('results/distributions/',sp_name))){# save plot
    dir.create(paste0('results/distributions/',sp_name),rec = T)
  }
  ggsave(file = paste0('results/distributions/',sp_name,'/novel_biomes.png'))
  
  RealmxBiome <-
    filter(RealmxBiome, count >= min_obs) # keep biomes with 5 or more observations
  
  for (i in unique(RealmxBiome$REALM)) {
    ecoreg_disti <- dplyr::filter(ecoreg_dist, REALM == i)
    ecoreg_disti <- dplyr::filter(ecoreg_disti,
                                  name %in% RealmxBiome$name) %>%  # select extent
      dplyr::group_by(REALM, BIOME_NUM) %>% # group by biomes and realms
      dplyr::summarise(
        count = sum(Nobs),
        BIOME_NAME = unique(BIOME_NAME),
        geometry = st_union (geometry)
      )  # summarise by biomes and realms
    
    my_disti <-
      my_dist[lengths(sf::st_intersects(my_dist, ecoreg_disti)) > 0, ] # keep only observation in the extent
    
    # my_dist_bck<-sf::st_sample(ecoreg_dist,size=n_bck)# generate background points in the extent
    # my_dist_bck <- sp::spsample(as_Spatial(ecoreg_disti),
    #                             n = n_bck,
    #                             type = 'random') %>%
    #   st_as_sf()
    
    # sampling strategy to avoid bug with st_sample
    grid.extent <- sf::st_bbox(ecoreg_disti) # set the grid extent for sampling
    env_grid <- terra::rast(
      res = 0.01666667,
      crs = "epsg:4326",
      xmin = grid.extent[1],
      ymin = grid.extent[2],
      xmax = grid.extent[3],
      ymax = grid.extent[4],
      vals = NA
    ) # create empty raster resolution is twice the env_grid to speed the process and lower the spatial autocorrelation
    env_grid<-terra::rasterize(ecoreg_disti,
                               env_grid,
                               fun = 'min',
                               field =1) # rasterize the extent for the background
    my_dist_bck<-which(values(env_grid)== 1) # select the cell with extent
    my_dist_bck<-terra::xyFromCell(env_grid,my_dist_bck) # convert raster to xy data frame
    if (nrow(my_dist_bck)>n_bck){ # sample the background (N max = 10000)
      my_dist_bck<-my_dist_bck[sample(1:nrow(my_dist_bck),n_bck),]
    }else{
      my_dist_bck<-my_dist_bck[1:nrow(my_dist_bck),]
    }
    colnames(my_dist_bck)<-c('X','Y')
    
    write.table(
      sf::st_coordinates(my_disti),
      file = paste0('data/sp_dist',spi,'/',
                    sp_name, '/',
                    sp_range, '_', i, '.txt'),
      sep = '\t',
      quote = FALSE,
      row.names = FALSE
    ) #write sp_dist
    
    write.table(
      my_dist_bck,
      file = paste0('data/sp_dist',spi,'/',
                    sp_name, '/', sp_range,
                    '_', i, '_bck', '.txt'),
      sep = '\t',
      quote = FALSE,
      row.names = FALSE
    ) # write sp background
    
    sp_bbox<-st_bbox(st_bbox(ecoreg_disti))# extent to plot
    ggplot() +
      geom_sf(data = world)  +
      geom_sf(data = st_union(ecoreg_disti),
              fill = 'lightyellow')+
      geom_sf(data = my_disti,
              color = 'darkred')+ 
      coord_sf(
        xlim = sp_bbox[c(1,3)],
        ylim = sp_bbox[c(2,4)],
        expand =FALSE
      ) # plot native distribution
    
    if (!dir.exists(paste0('results/distributions/',sp_name))){# save plot
      dir.create(paste0('results/distributions/',sp_name),rec = T)
    }
    ggsave(file = paste0('results/distributions/',sp_name,'/',i,'_novel_distributions.png'))
    
  }
  
  
  # niche dynamic analysis --------------------------------------------------
  
  ## load native and novel environmental data ----
  env<-terra::rast(c(list.files('data/clim_raster/', full.names = TRUE))) # load environmental variables
  names(env)<-list.files('data/clim_raster') # add the name of the layers
  sp_files<-list.files(paste0('data/sp_dist',spi,'/',sp_name),full.names = TRUE) # list of the distributions
  
  dist<-data.frame(X = numeric(),
                   Y = numeric(),
                   sp_range = character(),
                   presence = numeric()) # create empty data frame
  for (i in 1:length(sp_files)){ # create a data.frame with the different levels of distributions
    sp_file<-sp_files[i]
    presence = 1
    sp_range = 'native'
    if (length(grep('native',sp_file))==0){
      sp_range<-gsub('novel_','',sp_file) %>% 
        gsub('.txt','',x = .) %>% 
        strsplit(x = .,split  = '/')
      sp_range<-sp_range[[1]][4]
    }
    if (length(grep('bck',sp_file))==1){
      presence <- 0
      sp_range <- gsub('_bck','', sp_range)
    }
    
    disti<-read.table(sp_files[i],
                      header = TRUE,
                      sep = '\t')
    disti$sp_range<-sp_range
    disti$presence<-presence
    
    
    dist<-rbind(dist,disti)
  }
  
  dist<-st_as_sf(dist,coords = c('X','Y'),crs = 4326) # make it spatial
  dist<-cbind(dist,extract(env,vect(dist))) # extract environmental variable for each points
  dist<-na.exclude(dist) # remove NA values
  head(dist)
  
  ## PCA calibration ----
  env_col<-4:9 # index of the column to include in the Niche Analysis
  pca.env<-ade4::dudi.pca(st_drop_geometry(dist[,env_col]), scannf = FALSE, nf = 2) # make a PCA on the whole extent of the distribution
  
  scores.globclim <- pca.env$li # extract the score of the PCA
  
  # native niche
  scores.sp.nat<-pca.env$li[which(dist$sp_range == 'native' &
                                    dist$presence == 1),] # native species distribution
  scores.clim.nat<-pca.env$li[which(dist$sp_range == 'native' &
                                      dist$presence == 0),] # native background
  grid.clim.nat<-ecospat.grid.clim.dyn(glob=scores.globclim,
                                       glob1=scores.clim.nat,
                                       sp=scores.sp.nat, R=100,
                                       th.sp=0) #Grid the niche
  
  nov_ranges<-unique(dist$sp_range)[!(unique(dist$sp_range) %in% 'native')] # list of the different novel ranges
  
  ## Comparison native-novel (global) ----
  results<-data.frame(sp_range = c(nov_ranges,'all_novel'),
                      Nobs_native = nrow(scores.sp.nat),
                      Nobs = NA,
                      overlap_D = NA,
                      expansion = NA,
                      unfilling = NA,
                      stability = NA,
                      p_sim_cons = NA,
                      p_expansion = NA,
                      p_unfilling = NA,
                      p_stability = NA
  ) # empty dataframe
  
  
  # comparison with the whole novel ranges (all ranges pooled)
  scores.sp.inv<-pca.env$li[which(dist$sp_range != 'native' &
                                    dist$presence == 1),] # novel distribution
  scores.clim.inv<-pca.env$li[which(dist$sp_range != 'native' &
                                      dist$presence == 0),] # novel range background
  grid.clim.inv<-ecospat.grid.clim.dyn(glob=scores.globclim,
                                       glob1=scores.clim.inv,
                                       sp=scores.sp.inv, R=100,
                                       th.sp=0) # gridding the niche
  D.overlap <- ecospat.niche.overlap (grid.clim.nat, grid.clim.inv, cor = TRUE)$D # overlap D
  # To test for assumption of higher niche overlap
  sim.test <- ecospat.niche.similarity.test(grid.clim.nat, grid.clim.inv,rep=1000,
                                            overlap.alternative = "higher",
                                            expansion.alternative = "lower",
                                            stability.alternative = "higher",
                                            unfilling.alternative = "lower",
                                            intersection = 0.1,
                                            rand.type=2) # Similarity test
  
  ##### To test for assumption of lower niche overlap
  # sim.test <- ecospat.niche.similarity.test(grid.clim.nat, grid.clim.inv,rep=1000,
  #                                           overlap.alternative = "lower",
  #                                           expansion.alternative = "higher",
  #                                           stability.alternative = "lower",
  #                                           unfilling.alternative = "higher",
  #                                           intersection = 0.1,
  #                                           rand.type=2) # Similarity test
  
  niche.dyn <- ecospat.niche.dyn.index (grid.clim.nat, grid.clim.inv, intersection = 0.1) # niche dynamic indices
  
  # fill the results table
  results[results$sp_range== 'all_novel',]$Nobs<-nrow(scores.sp.inv)
  results[results$sp_range== 'all_novel',]$overlap_D<-D.overlap
  results[results$sp_range== 'all_novel',c('expansion','unfilling','stability')]<-niche.dyn$dynamic.index.w[c(1,3,2)]
  results[results$sp_range== 'all_novel',]$p_sim_cons<-sim.test$p.D
  results[results$sp_range== 'all_novel',]$p_expansion<-sim.test$p.expansion
  results[results$sp_range== 'all_novel',]$p_unfilling<-sim.test$p.unfilling
  results[results$sp_range== 'all_novel',]$p_stability<-sim.test$p.stability
  
  if (!dir.exists(paste0('results/niche/',sp_name))){ # save figures
    dir.create(paste0('results/niche/',sp_name),recursive = TRUE)
  }
  png (paste0('results/niche/',sp_name,'/all_novel.png'),
       width = 600, height = 400, units = "px")
  ecospat.plot.niche.dyn(grid.clim.nat, grid.clim.inv, quant=0.25, interest=2,
                         title= "Niche Overlap", name.axis1="PC1",
                         name.axis2="PC2")
  ecospat.shift.centroids(scores.sp.nat, scores.sp.inv, scores.clim.nat, scores.clim.inv)
  
  dev.off()
  
  png (paste0('results/niche/',sp_name,'/all_novel_tests.png'),
       width = 600, height = 600, units = "px" )
  par (mfrow = c(2,2)) 
  
  ecospat.plot.overlap.test(sim.test, "D", "Similarity tests")
  ecospat.plot.overlap.test(sim.test, "expansion", "")
  ecospat.plot.overlap.test(sim.test, "unfilling", "")
  ecospat.plot.overlap.test(sim.test, "stability", "")
  
  dev.off()
  
  ## Comparison native-novel (distinct ranges) ----
  sp_grid<-grid.clim.nat$z.uncor 
  extent_grid<-grid.clim.nat$Z
  for (i in nov_ranges){ #make niche analysis for each novel range
    scores.sp.inv<-pca.env$li[which(dist$sp_range == i &
                                      dist$presence == 1),]# novel distribution
    scores.clim.inv<-pca.env$li[which(dist$sp_range == i &
                                        dist$presence == 0),]#novel range background
    grid.clim.inv<-ecospat.grid.clim.dyn(glob=scores.globclim,
                                         glob1=scores.clim.inv,
                                         sp=scores.sp.inv, R=100,
                                         th.sp=0)#gridding the niche
    D.overlap <- ecospat.niche.overlap (grid.clim.nat, grid.clim.inv, cor = TRUE)$D #overlap D
    sim.test <- ecospat.niche.similarity.test(grid.clim.nat, grid.clim.inv,rep=100,
                                              overlap.alternative = "higher",
                                              expansion.alternative = "lower",
                                              stability.alternative = "higher",
                                              unfilling.alternative = "lower",
                                              intersection = 0.1,
                                              rand.type=2)# Similarity test
    niche.dyn <- ecospat.niche.dyn.index (grid.clim.nat, grid.clim.inv, intersection = 0.1) # niche dynamic indices
    
    # fill the results table
    results[results$sp_range== i,]$Nobs<-nrow(scores.sp.inv)
    results[results$sp_range== i,]$overlap_D<-D.overlap
    results[results$sp_range== i,c('expansion','unfilling','stability')]<-niche.dyn$dynamic.index.w[c(1,3,2)]
    results[results$sp_range== i,]$p_sim_cons<-sim.test$p.D
    results[results$sp_range== i,]$p_expansion<-sim.test$p.expansion
    results[results$sp_range== i,]$p_unfilling<-sim.test$p.unfilling
    results[results$sp_range== i,]$p_stability<-sim.test$p.stability
    
    if (!dir.exists(paste0('results/niche/',sp_name))){ #save figures
      dir.create(paste0('results/niche/',sp_name),recursive = TRUE)
    }
    png (paste0('results/niche/',sp_name,'/',i,'.png'),
         width = 600, height = 400, units = "px" )
    ecospat.plot.niche.dyn(grid.clim.nat, grid.clim.inv, quant=0.25, interest=2,
                           title= "Niche Overlap", name.axis1="PC1",
                           name.axis2="PC2")
    ecospat.shift.centroids(scores.sp.nat, scores.sp.inv, scores.clim.nat, scores.clim.inv)
    
    dev.off()
    
    png (paste0('results/niche/',sp_name,'/',i,'_tests.png'),
         width = 600, height = 600, units = "px" )
    par (mfrow = c(2,2)) 
    
    ecospat.plot.overlap.test(sim.test, "D", "Similarity tests")
    ecospat.plot.overlap.test(sim.test, "expansion", "")
    ecospat.plot.overlap.test(sim.test, "unfilling", "")
    ecospat.plot.overlap.test(sim.test, "stability", "")
    
    dev.off()
    
    sp_grid<-c(sp_grid,grid.clim.inv$z.uncor)
    extent_grid<-c(extent_grid,grid.clim.inv$Z)
  }
  
  ## Synthetize and save results ----
  write.table(results,file = paste0('results/niche/',sp_name,'/niche_dyn.txt'),
              sep = '\t', quote = F, row.names =F)  
  
  
  add.alpha <- function(col, alpha=1){ # function to add transprency to colors
    if(missing(col))
      stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
            rgb(x[1], x[2], x[3], alpha=alpha))  
  }
  
  # set color palette
  my.cols<-brewer.pal(8,'Set3')
  my.cols_trans <- add.alpha(my.cols[-1], alpha=0.3)
  my.cols_fill<-c(my.cols[1],my.cols_trans)
  
  #make figure
  png(paste0('results/niche/',sp_name,'/nich_dyn.png'),
      width = 800, height = 800, res = 150)
  
  for (i in 1:length(sp_grid)){
    add = FALSE
    main = paste0(sp_name)
    if (i >1){add = TRUE}
    my.rast<-rast(sp_grid[[i]])>0
    plot(my.rast,
         col = c(NA,my.cols_fill[i]),
         add = add,
         type = 'classes',
         legend = FALSE,
         mar = c(3.1, 3.1, 2.1, 7.1),
         main = main
    )
    plot(as.contour(my.rast, levels = TRUE), col = my.cols[i],add = T) 
  }
  
  # add legend
  xmax<-ext(sp_grid[[1]])[2]
  ymax<-ext(sp_grid[[1]])[4]
  y.length<- ext(sp_grid[[1]])[4]-ext(sp_grid[[1]])[3]
  x.length<-ext(sp_grid[[1]])[2]-ext(sp_grid[[1]])[1]
  y.leg<-seq(ext(sp_grid[[1]])[4]- (0.05*y.length),
             ext(sp_grid[[1]])[3]+(0.05*y.length),
             length.out = length(sp_grid))
  leg.txt<-c('Native',nov_ranges)
  par(xpd=TRUE)
  legend(xmax+0.05*x.length, 
         y.leg[1], 
         leg.txt, 
         col = my.cols, 
         cex = 1,
         pch = 15, 
         bty = 'n')
  par(xpd=FALSE)
  dev.off()
  
  print(results)
}


