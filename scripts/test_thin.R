# Distribution modeling for model comparison
library(ENMeval)
library(spocc)
library(spThin)
library(stringr)
library(maxnet)
library(usdm)

#parameters
thin.dist = 40
tax = 'Tsuga canadensis'
parallel = TRUE
nclus = 16
method = 'checkerboard2'

#setup ENMTools
#devtools::install_github("rsh249/cRacle")
library(cRacle)


#setup RSpatial
#devtools::install_github('oshea-patrick/RSpatial')
library(RSpatial)

#devtools::install_github('rsh249/vegdistmod')
library(vegdistmod)

ext = extent(c(-125,-50, 10, 60))

# load climate parameters
# r = stack(list.files(
#   '/usr/share/data/wc1.4/wc2-5/',
#   pattern = '.bil',
#   full.names = TRUE
# ))

if (!file.exists('envirem_clim.gri')){
  envir = get_envirem_clim()
  writeRaster(envir,filename = 'envirem_clim')
  envir.topo = get_envirem_elev()
  writeRaster(envir.topo, filename = 'envirem_topo')
} else {
  envir = stack('envirem_clim.gri')
  envir.topo = stack('envirem_topo.gri')
}

r= stack(envir, envir.topo)
r2 = crop(r, ext)

#DOWNLOAD DATA

occ = occ2df(occ(tax, from = 'gbif', limit = 20000))
occ.ex = extract(r2, occ[, c('longitude', 'latitude')])
occ = cbind(occ[c('longitude', 'latitude')], occ.ex)
occ = na.omit(occ)


#ext = extent(c(min(occ$longitude)-5, max(occ$longitude) + 5, min(occ$latitude)-5, max(occ$latitude)+5))


#background sampling
bg.pts = cRacle::rad_bg(occ[, c('longitude', 'latitude')],
                        r2,
                        radius = 2500,
                        n = 40000 / nrow(occ))


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
  bg.coords = bg.pts[, c('lon', 'lat')],
  clamp = TRUE,
  RMvalues = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)
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
fc = sapply(strsplit(fc1, ""), tolower)

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
plot(thr, col = c('black', 'blue'), main = 'no thin')
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
occ.thin = occ[-occ2thin, c('longitude', 'latitude')]
poThin.time = proc.time() - p


###TEMPORARY ADDITION####
occ.thin.dists = (rasterExtras:::distance(as.matrix(occ.thin), as.matrix(occ.thin)))
############


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
plot(thr.thin, col = c('black', 'blue'), main = 'poThin')
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
plot(thr.spthin, col = c('black', 'blue'), main = 'spThin')

bin.occ.spthin = extract(thr.spthin, sp.thin)
bin.abs.spthin = extract(thr.spthin, set.eval@bg.pts)
ev.set.spthin.bin <- evaluate(bin.occ.spthin, bin.abs.spthin)

# report timings nad model evaluation
print(poThin.time)
print(spthin.time)

print(e.sub)
print(e.sub.thin)
print(e.sub.spthin)
