library(itrackR)

context("samples")

z <- itrackr(edfs=itrackR.data('edfs')[1])
z <- load_samples(z,force=T, parallel=F)

test_that("Checking sample loading (serial)",
          {
            expect_equal(length(z$samples),1)
            expect_true(file.exists(z$samples[[1]]))

          })


test_that("sample blink removal", {

  z <- remove_blinks(z, interpolate = T)
  samp <- readRDS(z$samples[[1]])

  expect_is(samp,'data.table')
  expect_false(any(is.na(samp$pa)))
  expect_false(any(samp$pa<500))


})


test_that("epoching samples", {

  z <- find_messages(z,varnames = c('STIMONSET'),patterns = c('STIMONSET'),timestamp = T)

  expect_true('STIMONSET' %in% names(z$header))

  z <- epoch_samples(z,'STIMONSET',epoch = c(-500,500))

  expect_match(names(z$epochs$samples)[[1]],'STIMONSET')
  expect_equal(length(z$epochs$samples$STIMONSET),1)
  expect_is(z$epochs$samples$STIMONSET[[1]]$epochs,'matrix')
  expect_equal(ncol(z$epochs$samples$STIMONSET[[1]]$epochs),1000)
  expect_gt(var(z$epochs$samples$STIMONSET[[1]]$epochs[100,],na.rm=T),0)


  expect_warning(allepochs <- get_all_epochs(z,'STIMONSET',shape='wide',baseline=NULL,beh=T))

  expect_is(allepochs,'data.frame')
  expect_equal(ncol(allepochs),1002)

  check1 <- z$epochs$samples$STIMONSET[[1]]$epochs
  check2 <- allepochs[1:nrow(check1), 3:ncol(allepochs)]

  expect_true(all(check1==check2))


})
