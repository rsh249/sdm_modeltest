# Distribution modeling for model comparison
library(ENMeval)
library(spocc)
library(spThin)
library(stringr)
library(maxnet)
library(usdm)
library(parallel)
library(rgdal)
library(ggsci)
library(ggplot2)
library(reshape)


#parameters
thin.dist = 45
parallel = FALSE
nclus = 60
method = 'block'
kfolds = 10
buffer = 2000
ext = extent(c(-180,-40, 0, 80))

#setup 

#devtools::install_github("rsh249/cRacle")
library(cRacle)


#setup RSpatial
#devtools::install_github('oshea-patrick/RSpatial')
library(RSpatial)

#devtools::install_github('rsh249/vegdistmod')
library(vegdistmod)


#get Git repo with Little shapefiles:
system('git clone https://github.com/wpetry/USTreeAtlas')

shapefiles = list.files('USTreeAtlas/SHP', full.names = TRUE)

#data table
dat = read.csv('USTreeAtlas/Little_datatable.csv', stringsAsFactors = FALSE)



##TEST tax.list for GBIF hits:
test = vector()
for (i in 1:nrow(dat)) {
  test[i] = occ(dat[i, 'Latin.Name'], from = 'gbif', limit = 1)$gbif$meta$found
}
dat = cbind(dat, test)
targets = dat[which(dat$test >= 500 & dat$test <= 10000), ]

tax.list = targets$Latin.Name

#CLIMATE DATA
if (!file.exists('envirem_clim.gri')) {
  envir = get_envirem_clim()
  writeRaster(envir, filename = 'envirem_clim')
  envir.topo = get_envirem_elev()
  writeRaster(envir.topo, filename = 'envirem_topo')
} else {
  envir = stack('envirem_clim.gri')
  envir.topo = stack('envirem_topo.gri')
}

r = stack(envir, envir.topo)
r2 = crop(r, ext)
r2 = r2[[c(1, 2, 3, 4, 6, 8, 9, 10, 13, 16)]]




evalfun = function(tax) {
  out = tryCatch({
    #DOWNLOAD DATA
    
    occ = occ2df(occ(tax, from = 'gbif', limit = 20000))
    occ.ex = extract(r2, occ[, c('longitude', 'latitude')])
    occ = cbind(occ[c('longitude', 'latitude')], occ.ex)
    occ = na.omit(occ)
    
    #ext = extent(c(min(occ$longitude)-5, max(occ$longitude) + 5, min(occ$latitude)-5, max(occ$latitude)+5))
    
    
    #background sampling
    bg.pts = cRacle::rad_bg(occ[, c('longitude', 'latitude')],
                            r2,
                            radius = buffer,
                            n = 100000 / nrow(occ))
    #bg.pts = poThin(bg.pts,spacing = thin.dist, dimension = nrow(bg.pts), lat = 'lat', lon = 'lon')
    
    # ENMeval with unthinned data
    
    
    #leave this for adding more feature classes back in later
    #in.fc = c("L", "Q", "H", "P", "T")
    #fc = vector();
    #for (i in 1:length(in.fc)){
    #  comb = combn(in.fc, i)
    #  fc1 = do.call(paste, as.data.frame(t(comb), stringsAsFactors=FALSE));
    #  fc = c(fc, fc1)
    #}
    
    #fc = gsub(" ", "", fc)
    #fc = fc[grep("[LQ]", fc)]
    fc = c("L", "LQ")
    set.eval = ENMevaluate(
      occ[, c('longitude', 'latitude')],
      r2,
      rasterPreds = TRUE,
      parallel = parallel,
      fc = fc,
      #fc = fc[1:3],
      numCores = nclus,
      method = method,
      kfolds = kfolds,
      bg.coords = bg.pts[, c('lon', 'lat')],
      clamp = TRUE,
      RMvalues = c(0.5, 1, 1.5, 2, 2.5, 3)
      #RMvalues = c(0.5, 2)
    )
    
    
    
    best = which(set.eval@results[, 'AICc'] == min(na.omit(set.eval@results[, 'AICc'])))
    ev.set <-
      evaluate(occ[, c('longitude', 'latitude')], set.eval@bg.pts, set.eval@models[[best]], r2)
    thr.set <- threshold(ev.set)
    
    # For picking model parameters on the complete set
    best_param = set.eval@results[best, 1]
    best_arr = strsplit(as.character(best_param), "_")
    
    rm = best_arr[[1]][length(best_arr[[1]])]
    
    fc1 = best_arr[[1]][1:(length(best_arr[[1]]) - 1)]
    
    maxmatr = rbind(set.eval@occ.pts, set.eval@bg.pts)
    pres = c(rep(1, nrow(set.eval@occ.pts)), rep(0, nrow(set.eval@bg.pts)))
    maxmatr = cbind(maxmatr, pres)
    
    maxextr = extract(r2, maxmatr[, c('LON', 'LAT')])
    best_mod = maxnet(
      p = maxmatr[, 'pres'],
      data = as.data.frame(maxextr),
      maxnet.formula(
        p = maxmatr[, 'pres'],
        data = as.data.frame(maxextr),
        classes = stringr::str_to_lower(fc1)
      ),
      regmult = as.numeric(rm)
    )
    m = predict(r2, best_mod, type = 'cloglog')
    p.test = extract(m, set.eval@occ.pts)
    ab.test = extract(m, set.eval@bg.pts)
    e.sub = evaluate(p.test, ab.test,
                     best_mod)
    th.sub = threshold(e.sub)
    thr = m > th.sub$equal_sens_spec
  #  plot(thr, col = c('black', 'blue'), main = 'no thin')
    bin.occ = extract(thr, set.eval@occ.pts)
    bin.abs = extract(thr, set.eval@bg.pts)
    ev.set.bin <- evaluate(bin.occ, bin.abs)
    
    # thin with poThin
    p = proc.time()
    occ2thin = poThin(
      df = occ[, c('longitude', 'latitude')],
      spacing = thin.dist,
      dimension = nrow(occ),
      lon = 'longitude',
      lat = 'latitude'
    )
    if (length(occ2thin) < 1) {
      occ.thin = occ[, c('longitude', 'latitude')]
    } else {
      occ.thin = occ[-occ2thin, c('longitude', 'latitude')]
    }
    poThin.time = proc.time() - p
    
    
    colnames(occ.thin) = c('LON', 'LAT')
    maxmatr.thin = rbind(occ.thin, set.eval@bg.pts)
    pres.thin = c(rep(1, nrow(occ.thin)), rep(0, nrow(set.eval@bg.pts)))
    maxmatr.thin = cbind(maxmatr.thin, pres.thin)
    maxextr.thin = extract(r2, maxmatr.thin[, c('LON', 'LAT')])
    best_mod.thin = maxnet(
      p = maxmatr.thin[, 'pres.thin'],
      data = as.data.frame(maxextr.thin),
      maxnet.formula(
        p = maxmatr.thin[, 'pres.thin'],
        data = as.data.frame(maxextr.thin),
        classes = stringr::str_to_lower(fc1)
      ),
      regmult = as.numeric(rm)
    )
    m.thin = predict(r2, best_mod.thin, type = 'cloglog')
    
    
    
    p.test.thin = extract(m.thin, occ.thin)
    ab.test.thin = extract(m.thin, set.eval@bg.pts)
    e.sub.thin = evaluate(p.test.thin,
                          ab.test.thin,
                          best_mod.thin)
    
    th.sub.thin = threshold(e.sub.thin)
    thr.thin = m.thin > th.sub.thin$equal_sens_spec
   # plot(thr.thin,
     #    col = c('black', 'blue'),
        # main = 'poThin')
    bin.occ.thin = extract(thr.thin, occ.thin)
    bin.abs.thin = extract(thr.thin, set.eval@bg.pts)
    ev.set.thin.bin <- evaluate(bin.occ.thin, bin.abs.thin)
    
    #thin with spThin
    occ = cbind(rep(tax, nrow(occ)), occ)
    colnames(occ)[1] = 'tax'
    p = proc.time()
    spthinned = thin(
      loc.data = occ,
      lat.col = "latitude",
      long.col = "longitude",
      spec.col = "tax",
      thin.par = thin.dist,
      reps = 1,
      locs.thinned.list.return = TRUE,
      write.files = FALSE,
      max.files = 5,
      out.dir = "test/",
      out.base = "thinned",
      write.log.file = FALSE,
      log.file = "log_file.txt"
    )
    spthin.time = proc.time() - p
    
    
    #maxent with spThinned data
    sp.thin = spthinned[[1]]
    colnames(sp.thin) = c('LON', 'LAT')
    maxmatr.spthin = rbind(sp.thin, set.eval@bg.pts)
    pres.spthin = c(rep(1, nrow(sp.thin)), rep(0, nrow(set.eval@bg.pts)))
    maxmatr.spthin = cbind(maxmatr.spthin, pres.spthin)
    maxextr.spthin = extract(r2, maxmatr.spthin[, c('LON', 'LAT')])
    best_mod.spthin = maxnet(
      p = maxmatr.spthin[, 'pres.spthin'],
      data = as.data.frame(maxextr.spthin),
      maxnet.formula(
        p = maxmatr.spthin[, 'pres.spthin'],
        data = as.data.frame(maxextr.spthin),
        classes = stringr::str_to_lower(fc1)
      ),
      regmult = as.numeric(rm)
    )
    m.spthin = predict(r2, best_mod.spthin, type = 'cloglog')
    p.test.spthin = extract(m, sp.thin)
    ab.test.spthin = extract(m, set.eval@bg.pts)
    e.sub.spthin = evaluate(p.test.spthin,
                            ab.test.spthin,
                            best_mod.spthin)
    th.sub.spthin = threshold(e.sub.spthin)
    thr.spthin = m.spthin > th.sub.spthin$equal_sens_spec
    #plot(thr.spthin,
       #  col = c('black', 'blue'),
         #main = 'spThin')
    
    bin.occ.spthin = extract(thr.spthin, sp.thin)
    bin.abs.spthin = extract(thr.spthin, set.eval@bg.pts)
    ev.set.spthin.bin <- evaluate(bin.occ.spthin, bin.abs.spthin)
    
    # report timings nad model evaluation
    print(poThin.time)
    print(spthin.time)
    
    print(e.sub)
    print(e.sub.thin)
    print(e.sub.spthin)
    
    un2po = nicheOverlap(
      m,
      m.thin,
      stat = 'D',
      mask = FALSE,
      checkNegatives = FALSE
    )
    un2sp = nicheOverlap(
      m,
      m.spthin,
      stat = 'D',
      mask = FALSE,
      checkNegatives = FALSE
    )
    sp2po = nicheOverlap(
      m.spthin,
      m.thin,
      stat = 'D',
      mask = FALSE,
      checkNegatives = FALSE
    )
    
    un2po.thr = nicheOverlap(
      thr,
      thr.thin,
      stat = 'D',
      mask = FALSE,
      checkNegatives = FALSE
    )
    un2sp.thr = nicheOverlap(
      thr,
      thr.spthin,
      stat = 'D',
      mask = FALSE,
      checkNegatives = FALSE
    )
    sp2po.thr = nicheOverlap(
      thr.spthin,
      thr.thin,
      stat = 'D',
      mask = FALSE,
      checkNegatives = FALSE
    )
    filebase = paste('rasters/', gsub(" ", "_", tax), 'base', sep = '_')
    writeRaster(m, filename = filebase, overwrite = TRUE)
    filethin = paste('rasters/', gsub(" ", "_", tax), 'poThin', sep = '_')
    writeRaster(m.thin, filename = filethin, overwrite = TRUE)
    filespThin = paste('rasters/', gsub(" ", "_", tax), 'spThin', sep = '_')
    writeRaster(m.spthin, filename = filespThin, overwrite = TRUE)
    
    distsbase = paste('dists/', gsub(" ", "_", tax), 'base', sep = '_')
    write.csv(occ, file = distsbase)
    diststhin = paste('dists/', gsub(" ", "_", tax), 'poThin', sep = '_')
    write.csv(occ.thin, file = diststhin)
    distsspThin = paste('dists/', gsub(" ", "_", tax), 'spThin', sep = '_')
    write.csv(sp.thin, file = distsspThin)
    
    return(
      c(
        tax,
        poThin.time[3],
        spthin.time[3],
        nrow(occ),
        nrow(occ.thin),
        nrow(sp.thin),
        e.sub@auc,
        e.sub.thin@auc,
        e.sub.spthin@auc,
        un2po,
        un2sp,
        sp2po,
        un2sp.thr,
        un2po.thr,
        sp2po.thr
      )
    )
  },
  error = function(cond) {
    # Choose a return value in case of error
    return(NA)
  })
  
}





cl = makeCluster(nclus, type = "FORK")
p = proc.time()

coll = parLapply(cl, tax.list, evalfun)

proc.time() - p

stopCluster(cl)

#Collect results
colldf = (coll[[1]])


for (i in 2:length(coll)) {
  if (!is.na(coll[i])) {
    colldf = rbind(colldf, coll[[i]])
    
  }
}

colnames(colldf) = c(
  'tax',
  'poThin.time',
  'spthin.time',
  'occ.count',
  'occ.thin.count',
  'occ.sp.thin.count',
  'e.sub.auc',
  'e.sub.thin.auc',
  'e.sub.spthin.auc',
  'un2po',
  'un2sp',
  'sp2po',
  'un2sp.thr',
  'un2po.thr',
  'sp2po.thr'
)
colldf = as.data.frame(colldf)
for (co in 2:ncol(colldf)) {
  colldf[, co] = as.numeric(as.character(colldf[, co]))
}


#plotting

ggplot(data = colldf) +
  geom_point(aes(x = poThin.time, y = spthin.time)) +
  scale_y_log10() +
  theme_linedraw() +
  xlab('poThin time (s)') +
  ylab('spThin time (s)')

niche.overlap = melt(colldf,
                     id.vars = "tax",
                     measure.vars = c('un2po', 'un2sp', 'sp2po'))
ggplot(data = niche.overlap, aes(x = variable, y = value)) +
  geom_violin() +
  geom_jitter(position = position_jitter(0.02)) +
  theme_linedraw()
ggsave('niche_overlap.png', dpi=500, width = 7.25, height = 4)

##Test niche overlap differences
aov.overlap = aov(formula = value ~ variable, data = niche.overlap)
niche.tukeyHSD = TukeyHSD(aov.overlap)


niche.overlap.thr = melt(
  colldf,
  id.vars = "tax",
  measure.vars = c('un2po.thr', 'un2sp.thr', 'sp2po.thr')
)
ggplot(data = niche.overlap.thr, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(0.05))



### REMEMBER TO ADD LEGEND
ggplot(data = colldf) +
  geom_point(aes(x = occ.count, y = poThin.time/60)) +
  geom_point(aes(x = occ.count, y = spthin.time/60), col = 'darkblue')  +
  geom_smooth(aes(x = occ.count, y = poThin.time/60), method='loess', col = 'black') +
  geom_smooth(aes(x = occ.count, y = spthin.time/60), col = 'darkblue', method = 'loess')  +
  scale_y_sqrt() +
  theme_linedraw() +
  xlab('Starting Occurrence Count') +
  ylab('Time (min)') +
  theme(legend.position = "right")
ggsave('thin_timings.png', dpi = 500, width = 7.25, height = 4)



## Compare to Little shapefile "Expert" range maps
little.map.eval <- function(n) {
  out = tryCatch({
    taxstub = targets[n, 'SHP..']
    tax = targets[n, 'Latin.Name']
    filebase = paste('rasters/', gsub(" ", "_", tax), 'base', sep = '_')
    m = raster::raster(filebase)
    filethin = paste('rasters/', gsub(" ", "_", tax), 'poThin', sep = '_')
    m.thin = raster::raster(filethin)
    filespThin = paste('rasters/', gsub(" ", "_", tax), 'spThin', sep = '_')
    m.spThin = raster::raster(filespThin)
    
    shp = readOGR(paste('USTreeAtlas/SHP/', taxstub, sep = ''), taxstub)
    shpr = rasterize(shp, m) >= 1
    shpr[is.na(shpr[])] <- 0
    shpr = mask(shpr, r2[[1]])
    
    
    #shpr.df = as.data.frame(shpr, xy = TRUE)
    #mdf = as.data.frame(m, xy = TRUE)
    #m.thindf = as.data.frame(m.thin, xy = TRUE)
    #m.spthindf = as.data.frame(m.spThin, xy = TRUE)
    # plot(density(mdf$layer[which(shpr.df$layer > 0)], na.rm = T), ylim = c(0, 3.5))
    # points(density(m.thindf$layer[which(shpr.df$layer > 0)], na.rm = T), type =
    #          'l')
    # points(density(m.spthindf$layer[which(shpr.df$layer > 0)], na.rm = T), type =
    #          'l')
    distsbase = paste('dists/', gsub(" ", "_", tax), 'base', sep = '_')
    
    dist = read.csv(distsbase, stringsAsFactors = F)
    dist$longitude = as.numeric(dist$longitude)
    dist$latitude = as.numeric(dist$latitude)
    dist.bg = vegdistmod::rad_bg(dist[,3:4], r2, buffer/2,
                                 n = 100000 / nrow(dist))
    bg.t = poThin(dist.bg, spacing = thin.dist, 
                  dimension = nrow(dist.bg), 
                  lat = 'lat', 
                  lon='lon')
    if (length(bg.t) < 1) {
      dist.bg = dist.bg
    } else {
      dist.bg = dist.bg[-bg.t,]
    }
    
    m.ex = raster::extract(stack(m, shpr), dist.bg[,c('lon', 'lat')])
    m.ex = na.omit(m.ex)
    m.ev = evaluate(m.ex[m.ex[, 2] == 1, 1], m.ex[m.ex[, 2] == 0, 1])
    
    m.spthin.ex = raster::extract(stack(m.spThin, shpr), dist.bg[,c('lon', 'lat')])
    m.spthin.ex = na.omit(m.spthin.ex)
    m.spthin.ev = evaluate(m.spthin.ex[m.spthin.ex[, 2] == 1, 1], m.spthin.ex[m.spthin.ex[, 2] == 0, 1])
    
    m.thin.ex = raster::extract(stack(m.thin, shpr), dist.bg[,c('lon', 'lat')])
    m.thin.ex = na.omit(m.thin.ex)
    m.thin.ev = evaluate(m.thin.ex[m.thin.ex[, 2] == 1, 1], m.thin.ex[m.thin.ex[, 2] == 0, 1])
    return(c(m.ev@auc, m.thin.ev@auc, m.spthin.ev@auc))
  },
  error = function(cond) {
    # Choose a return value in case of error
    return(NA)
  })
}
minclus = nclus/8 ; #Smaller cluster for memory intensive steps
cl2 = makeCluster(minclus, type = 'FORK')
little.auc.l = parLapplyLB(cl2,  1:nrow(targets), little.map.eval, chunk.size = 1)

#little.auc.l = parLapply(cl2,  1:minclus, little.map.eval)
stopCluster(cl2)

little.auc = little.auc.l[[1]]
for (co in 2:length(little.auc.l)) {
  little.auc= rbind(little.auc, little.auc.l[[co]])
}

colnames(little.auc) = c('base.auc', 'poThin.auc', 'spThin.auc')
little.auc.m = melt(little.auc)

ggplot(data=little.auc.m) + 
  geom_boxplot(aes(x=X2, y=value))

##Test niche overlap differences
aov.little.auc = aov(formula = value ~ X2, data = little.auc.m)
little.auc.tukeyHSD = TukeyHSD(aov.little.auc)
