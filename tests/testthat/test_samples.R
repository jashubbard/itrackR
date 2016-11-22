context('samples')

z <- itrackr(path=itrackR.data('path'))

test_that('Sample Functions', {



  z <- set_index(z,varnames=c('Block','Trial'), patterns=c('^BLOCK [0-9]*','^TRIAL [0-9]*'), numeric.only=T)
  z <- find_messages(z,varnames = c('STIMONSET'),patterns = c('STIMONSET'),timestamp = T)

  z <- add_behdata(z,itrackR.data('beh'),append=F)

  z <- load_samples(z,force=T, parallel=F)
  expect_equal(length(z$samples),length(z$subs))
  expect_true(file.exists(z$samples[[1]]))


  z <- remove_blinks(z, interpolate = T)
  samp <- get_samples(z,z$subs[1])

  expect_is(samp,'data.frame')
  expect_false(any(is.na(samp$pa)))



  expect_true('STIMONSET' %in% names(z$header))

  epochs <- get_sample_epochs(z, factors = c('Task','Conflict'), parallel=T,
                              event='STIMONSET', epoch = c(-500,500), baseline= c(-100,0),aggregate=T)

  expect_is(epochs,'data.frame')
  expect_true(all(c('Task','Conflict') %in% names(epochs)))






})
